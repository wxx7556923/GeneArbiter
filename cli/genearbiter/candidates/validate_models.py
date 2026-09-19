#!/usr/bin/env python3
# script_id_md5: f3b92224c2c9ce7179b7de0212b5103d
# created: 2026-06-23
# modified: 2026-06-26
# owner: project
# status: project_code
# purpose: 在证据卡前对当前注释/从头预测候选编码模型做致命质量控制准入筛查。
# inputs: 命令行参数指定的 TSV/JSONL/GFF3/FASTA 输入。
# outputs: 命令行参数指定的 TSV/JSONL/GFF3/Markdown 输出。
# notes: 该步骤仅做候选模型结构和编码检查。

"""在卡片构建前对当前注释/从头预测候选编码模型运行致命质量控制。

输出是确定性的逐模型致命质量控制摘要。它会通过 --validation-summary
传给 build_arbitration_cards.py，使存在硬性编码结构失败的模型能在
证据仲裁前被排除。本脚本有意不做证据排序，也不评估 RNA/protein 支持。
"""

from __future__ import annotations

import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_COMMON_DIR = _SCRIPT_DIR.parent / "common"
if str(_COMMON_DIR) not in sys.path:
    sys.path.insert(0, str(_COMMON_DIR))

import argparse
import csv
import gzip
from typing import Dict, List, Sequence, Set, Tuple

from gff_utils import Transcript, parse_annotation_gff, write_tsv


STOP_CODONS = {"TAA", "TAG", "TGA"}

FIELDS = [
    "model_id",
    "source",
    "gene_id",
    "transcript_id",
    "seqid",
    "start",
    "end",
    "strand",
    "exon_count",
    "cds_count",
    "cds_length",
    "fatal_qc_status",
    "validation_status",
    "overall_validation_status",
    "cds_status",
    "cds_length_mod3",
    "internal_stop_count",
    "internal_stop_status",
    "phase_status",
    "assembly_gap_or_N_in_CDS",
    "sequence_status",
    "fatal_fail_reasons",
    "validation_fail_reasons",
    "blocking_reasons",
]


def parse_tool_arg(text: str) -> Tuple[str, str]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("--tool-gff must be name=path")
    name, path = text.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("--tool-gff must be name=path")
    return name, path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--current-gff", required=True)
    parser.add_argument("--tool-gff", action="append", type=parse_tool_arg, default=[], help="Repeatable: name=path")
    parser.add_argument("--gene-list", default="", help="Optional newline/TSV gene or transcript ID list for pilot subset selection")
    parser.add_argument("--region-table", default="", help="Optional TSV with seqid/start/end columns for pilot subset selection")
    parser.add_argument("--region-flank", type=int, default=0, help="Flank bp added when selecting by --region-table")
    parser.add_argument("--output-validation", required=True)
    return parser.parse_args()


def safe_id(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text)


def split_ids(text: str) -> List[str]:
    ids: List[str] = []
    for item in text.replace(";", ",").split(","):
        item = item.strip()
        if item and item != ".":
            ids.append(item)
    return ids


def load_gene_filter(path: str) -> Set[str]:
    if not path:
        return set()
    ids: Set[str] = set()
    with open(path, "r") as probe:
        first = probe.readline().rstrip("\n")
    looks_like_tsv = "\t" in first and any(field in first.split("\t") for field in ["gene_id", "current_gene_id", "transcript_id", "locus_id"])
    if looks_like_tsv:
        with open(path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for field in ["gene_id", "current_gene_id", "transcript_id", "locus_id"]:
                    if field in row:
                        ids.update(split_ids(row.get(field, "")))
        return ids
    with open(path, "r") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                ids.update(split_ids(line.strip()))
    return ids


def load_regions(path: str, flank: int) -> List[Tuple[str, int, int]]:
    if not path:
        return []
    regions: List[Tuple[str, int, int]] = []
    with open(path, "r", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            seqid = row.get("seqid") or row.get("chrom") or row.get("chromosome") or ""
            try:
                start = int(float(row.get("start") or row.get("region_start") or 0))
                end = int(float(row.get("end") or row.get("region_end") or 0))
            except ValueError:
                continue
            if not seqid or not start or not end:
                continue
            if start > end:
                start, end = end, start
            regions.append((seqid, max(1, start - flank), end + flank))
    return regions


def transcript_overlaps_regions(tx: Transcript, regions: Sequence[Tuple[str, int, int]]) -> bool:
    for seqid, start, end in regions:
        if tx.seqid == seqid and tx.end >= start and tx.start <= end:
            return True
    return False


def transcript_selected(tx: Transcript, model_id: str, gene_id: str, gene_ids: Set[str], regions: Sequence[Tuple[str, int, int]]) -> bool:
    if not gene_ids and not regions:
        return True
    if gene_id in gene_ids or tx.transcript_id in gene_ids or model_id in gene_ids:
        return True
    return transcript_overlaps_regions(tx, regions) if regions else False


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def load_fasta(path: str) -> Dict[str, str]:
    genome: Dict[str, List[str]] = {}
    current_id = ""
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                genome.setdefault(current_id, [])
            elif current_id:
                genome[current_id].append(line.upper())
    return {seqid: "".join(parts) for seqid, parts in genome.items()}


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")
    return seq.translate(table)[::-1].upper()


def get_genomic_seq(genome: Dict[str, str], seqid: str, start: int, end: int) -> str:
    seq = genome.get(seqid, "")
    if not seq or start < 1 or end > len(seq) or start > end:
        return ""
    return seq[start - 1:end].upper()


def cds_blocks_in_transcript_order(tx: Transcript) -> List[Tuple[int, int, str]]:
    blocks = tx.sorted_cds()
    if tx.strand == "-":
        return sorted(blocks, key=lambda item: (item[0], item[1]), reverse=True)
    return sorted(blocks, key=lambda item: (item[0], item[1]))


def extract_cds_sequence(genome: Dict[str, str], tx: Transcript) -> Tuple[str, bool]:
    parts = []
    missing = False
    for start, end, _phase in cds_blocks_in_transcript_order(tx):
        part = get_genomic_seq(genome, tx.seqid, start, end)
        if not part:
            missing = True
            continue
        if tx.strand == "-":
            part = revcomp(part)
        parts.append(part)
    return "".join(parts).upper(), missing


def check_phase(tx: Transcript) -> Tuple[str, List[str]]:
    blocks = cds_blocks_in_transcript_order(tx)
    if not blocks:
        return "not_assessable", []
    blocking_reasons = []
    metadata_reasons = []
    cumulative = 0
    for index, (start, end, phase) in enumerate(blocks):
        if phase in {"", "."}:
            metadata_reasons.append("phase_missing_metadata:{0}-{1}:{2}".format(start, end, phase or "."))
        elif phase not in {"0", "1", "2"}:
            blocking_reasons.append("invalid_phase_value:{0}-{1}:{2}".format(start, end, phase))
        else:
            expected = 0 if index == 0 else (3 - (cumulative % 3)) % 3
            if int(phase) != expected:
                blocking_reasons.append("phase_mismatch:{0}-{1}:observed_{2}:expected_{3}".format(start, end, phase, expected))
        cumulative += end - start + 1
    if blocking_reasons:
        return "phase_inconsistent", blocking_reasons + metadata_reasons
    if metadata_reasons:
        return "phase_missing_metadata", metadata_reasons
    return "phase_consistent", []



def validate_transcript(model_id: str, source: str, gene_id: str, tx: Transcript, genome: Dict[str, str]) -> Dict[str, object]:
    cds_seq, missing_cds_sequence = extract_cds_sequence(genome, tx)
    cds_length = len(cds_seq)
    cds_blocks = tx.sorted_cds()
    exons = tx.sorted_exons()
    blocking: List[str] = []

    if not cds_blocks:
        blocking.append("no_CDS_features")
    if missing_cds_sequence:
        blocking.append("CDS_sequence_unavailable_or_out_of_bounds")

    if cds_length:
        cds_length_mod3 = "valid" if cds_length % 3 == 0 else "invalid"
        if cds_length % 3 != 0:
            blocking.append("CDS_length_mod3_invalid")
        codons = [cds_seq[i:i + 3] for i in range(0, max(0, len(cds_seq) - 3), 3)]
        internal_stop_count = sum(1 for codon in codons if codon in STOP_CODONS)
        internal_stop_status = "absent" if internal_stop_count == 0 else "present"
        if internal_stop_count:
            blocking.append("internal_stop_present")
    else:
        cds_length_mod3 = "not_assessable"
        internal_stop_count = 0
        internal_stop_status = "not_assessable"

    phase_status, phase_reasons = check_phase(tx)
    if phase_status == "phase_inconsistent":
        blocking.append("phase_inconsistent")
    blocking.extend(reason for reason in phase_reasons if reason.startswith(("invalid_phase_value", "phase_mismatch")))

    fatal_qc_status = "fail" if blocking else "pass"
    fail_reasons = ",".join(sorted(set(blocking)))

    return {
        "model_id": model_id,
        "source": source,
        "gene_id": gene_id,
        "transcript_id": tx.transcript_id,
        "seqid": tx.seqid,
        "start": tx.start,
        "end": tx.end,
        "strand": tx.strand,
        "exon_count": len(exons),
        "cds_count": len(cds_blocks),
        "cds_length": cds_length,
        "fatal_qc_status": fatal_qc_status,
        "validation_status": fatal_qc_status,
        "overall_validation_status": fatal_qc_status,
        "cds_status": "CDS_present" if cds_blocks else "CDS_absent",
        "cds_length_mod3": cds_length_mod3,
        "internal_stop_count": internal_stop_count,
        "internal_stop_status": internal_stop_status,
        "phase_status": phase_status,
        "assembly_gap_or_N_in_CDS": "yes" if "N" in cds_seq else "no",
        "sequence_status": "unavailable_or_out_of_bounds" if missing_cds_sequence else "available",
        "fatal_fail_reasons": fail_reasons,
        "validation_fail_reasons": fail_reasons,
        "blocking_reasons": fail_reasons,
    }

def validate_gff(
    path: str,
    source: str,
    genome: Dict[str, str],
    gene_ids: Set[str],
    regions: Sequence[Tuple[str, int, int]],
) -> List[Dict[str, object]]:
    _genes, transcripts, tx_to_gene = parse_annotation_gff(path, source)
    rows = []
    for tx_id, tx in sorted(transcripts.items()):
        gene_id = tx_to_gene.get(tx_id, tx.gene_id) or tx_id
        if not tx.seqid or not tx.start or not tx.end:
            continue
        model_id = "{0}:{1}:{2}".format(safe_id(source), safe_id(gene_id), safe_id(tx_id))
        if not transcript_selected(tx, model_id, gene_id, gene_ids, regions):
            continue
        rows.append(validate_transcript(model_id, source, gene_id, tx, genome))
    return rows


def main() -> None:
    args = parse_args()
    genome = load_fasta(args.genome_fasta)
    gene_ids = load_gene_filter(args.gene_list)
    regions = load_regions(args.region_table, args.region_flank)
    rows = []
    for source, path in [("current", args.current_gff)] + list(args.tool_gff):
        rows.extend(validate_gff(path, source, genome, gene_ids, regions))
    write_tsv(args.output_validation, FIELDS, rows)


if __name__ == "__main__":
    main()
