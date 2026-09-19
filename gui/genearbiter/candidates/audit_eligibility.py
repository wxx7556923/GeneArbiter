#!/usr/bin/env python3
# script_id_md5: 0e98830471edd1a8865ed92f1cc7c373
# created: 2026-06-25
# modified: 2026-06-26
# owner: project
# status: project_code
# purpose: 审查候选模型是否可进入 GeneArbiter 可选候选池。
# inputs: 基因组 FASTA；当前注释 GFF3；可选工具 GFF3；可选 RNA junction TSV。
# outputs: 逐模型准入资格 TSV 和可选的来源/原因汇总 TSV。
# notes: 不读取真值注释；硬检查限于致命 CDS/phase/坐标问题；start/stop 和非典型剪接只作警告。

"""审查候选编码模型是否具备进入 AI 前候选池的资格。

本脚本把模型准入资格和证据强度分开处理。未通过确定性致命编码/坐标
检查的模型会被标记为不能进入 AI 可选候选池。RNA/protein
证据强度、start/stop codon 缺失以及非典型剪接 motif 在这里
有意不作为硬过滤条件。
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
import os
from collections import defaultdict
from typing import Dict, List, Sequence, Set, Tuple

from gff_utils import Transcript, iter_gff_rows, parse_annotation_gff, write_tsv


STOP_CODONS = {"TAA", "TAG", "TGA"}
CANONICAL_SPLICE_PAIRS = {"GT-AG", "GC-AG", "AT-AC"}

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
    "eligible_for_ai",
    "eligibility_status",
    "hard_fail_reasons",
    "warning_reasons",
    "cds_length_mod3",
    "start_codon_status",
    "stop_codon_status",
    "internal_stop_count",
    "phase_status",
    "splice_motif_status",
    "noncanonical_splice_count",
    "noncanonical_splice_junction_supported_count",
    "coordinate_status",
    "parent_status",
    "neighbor_conflict_status",
    "source_file",
]

SUMMARY_FIELDS = [
    "source",
    "eligibility_status",
    "reason",
    "count",
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
    parser.add_argument("--junctions", action="append", default=[], help="RNA junction TSV; repeatable")
    parser.add_argument("--min-junction-support-count", type=int, default=3)
    parser.add_argument("--min-junction-sample-count", type=int, default=1)
    parser.add_argument("--allow-partial-models", action="store_true", help="Compatibility flag; missing start/stop codons are always warnings, not hard eligibility failures.")
    parser.add_argument("--same-source-conflict-fraction", type=float, default=0.90, help="Hard fail same-source different-gene high-overlap conflicts; set 0 to disable.")
    parser.add_argument("--gene-list", default="", help="Optional newline/TSV gene or transcript ID list for pilot subset selection")
    parser.add_argument("--region-table", default="", help="Optional TSV with seqid/start/end columns for pilot subset selection")
    parser.add_argument("--region-flank", type=int, default=0)
    parser.add_argument("--output-eligibility", required=True)
    parser.add_argument("--output-summary", default="")
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
    return any(tx.seqid == seqid and tx.end >= start and tx.start <= end for seqid, start, end in regions)


def transcript_selected(tx: Transcript, model_id: str, gene_id: str, gene_ids: Set[str], regions: Sequence[Tuple[str, int, int]]) -> bool:
    if not gene_ids and not regions:
        return True
    if gene_id in gene_ids or tx.transcript_id in gene_ids or model_id in gene_ids:
        return True
    return transcript_overlaps_regions(tx, regions) if regions else False


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


def check_phase(tx: Transcript) -> Tuple[str, List[str], List[str]]:
    blocks = cds_blocks_in_transcript_order(tx)
    if not blocks:
        return "not_assessable", [], []
    hard: List[str] = []
    warnings: List[str] = []
    cumulative = 0
    for index, (start, end, phase) in enumerate(blocks):
        if phase in {"", "."}:
            warnings.append("phase_missing_metadata:{0}-{1}".format(start, end))
        elif phase not in {"0", "1", "2"}:
            hard.append("invalid_phase_value:{0}-{1}:{2}".format(start, end, phase))
        else:
            expected = 0 if index == 0 else (3 - (cumulative % 3)) % 3
            if int(phase) != expected:
                hard.append("phase_mismatch:{0}-{1}:observed_{2}:expected_{3}".format(start, end, phase, expected))
        cumulative += end - start + 1
    if hard:
        return "phase_inconsistent", hard, warnings
    if warnings:
        return "phase_missing_metadata", [], warnings
    return "phase_consistent", [], []


def load_junctions(paths: Sequence[str], min_support_count: int, min_sample_count: int) -> Set[Tuple[str, int, int, str]]:
    junctions: Set[Tuple[str, int, int, str]] = set()
    merged: Dict[Tuple[str, int, int, str], Dict[str, object]] = {}
    for path in paths:
        with open(path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    start = int(row.get("intron_start_1based") or row.get("intron_start") or 0)
                    end = int(row.get("intron_end_1based") or row.get("intron_end") or 0)
                    support = int(float(row.get("support_count") or row.get("junction_read_count") or row.get("unique_junction_read_count") or 0))
                except ValueError:
                    continue
                if support < min_support_count:
                    continue
                seqid = row.get("seqid") or row.get("chrom") or row.get("chromosome") or ""
                strand = row.get("strand") or "."
                key = (seqid, start, end, strand)
                entry = merged.setdefault(key, {"support": 0, "samples": set()})
                entry["support"] += support
                sample_ids = row.get("sample_ids") or row.get("supporting_samples") or ""
                if sample_ids:
                    for sample_id in sample_ids.split(","):
                        if sample_id:
                            entry["samples"].add(sample_id)
                else:
                    try:
                        sample_count = int(float(row.get("sample_count") or row.get("sample_support_count") or 1))
                    except ValueError:
                        sample_count = 1
                    for index in range(max(1, sample_count)):
                        entry["samples"].add("{0}:sample_{1}".format(os.path.basename(path), index + 1))
    for key, entry in merged.items():
        if len(entry["samples"]) >= min_sample_count:
            junctions.add(key)
    return junctions


def exact_junction_supported(junctions: Set[Tuple[str, int, int, str]], tx: Transcript, start: int, end: int) -> bool:
    return (tx.seqid, start, end, tx.strand) in junctions or (tx.seqid, start, end, ".") in junctions


def check_splice_motifs(
    genome: Dict[str, str],
    tx: Transcript,
    junctions: Set[Tuple[str, int, int, str]],
) -> Tuple[str, int, int, List[str], List[str]]:
    introns = tx.introns()
    if not introns:
        return "single_exon_no_splice_to_assess", 0, 0, [], []
    hard: List[str] = []
    warnings: List[str] = []
    noncanonical = 0
    supported_noncanonical = 0
    for intron_start, intron_end in introns:
        intron_seq = get_genomic_seq(genome, tx.seqid, intron_start, intron_end)
        if len(intron_seq) < 2:
            noncanonical += 1
            hard.append("intron_sequence_unavailable:{0}-{1}".format(intron_start, intron_end))
            continue
        if tx.strand == "-":
            intron_seq = revcomp(intron_seq)
        motif = "{0}-{1}".format(intron_seq[:2], intron_seq[-2:])
        if motif in CANONICAL_SPLICE_PAIRS:
            continue
        noncanonical += 1
        detail = "noncanonical_splice:{0}-{1}:{2}".format(intron_start, intron_end, motif)
        if junctions and exact_junction_supported(junctions, tx, intron_start, intron_end):
            supported_noncanonical += 1
            warnings.append(detail + ":exact_junction_supported")
        else:
            warnings.append(detail + ":no_exact_junction_support")
    if noncanonical == 0:
        return "canonical_splice_site", 0, 0, hard, warnings
    if supported_noncanonical == noncanonical:
        return "noncanonical_splice_site_junction_supported", noncanonical, supported_noncanonical, hard, warnings
    return "noncanonical_splice_site_warning", noncanonical, supported_noncanonical, hard, warnings


def audit_raw_parent_status(path: str) -> Dict[str, str]:
    transcript_ids: Set[str] = set()
    child_parent_ids: Set[str] = set()
    for _seqid, _source, feature, _start, _end, _score, _strand, _phase, attrs in iter_gff_rows(path):
        feature_l = feature.lower()
        if feature_l in {"mrna", "transcript", "rna"}:
            tx_id = attrs.get("ID") or attrs.get("Name")
            if tx_id:
                transcript_ids.add(tx_id)
        elif feature_l in {"exon", "cds"}:
            child_parent_ids.update(parent for parent in attrs.get("Parent", "").split(",") if parent)
    status: Dict[str, str] = {}
    for parent in child_parent_ids:
        status[parent] = "parent_transcript_present" if parent in transcript_ids else "child_parent_transcript_feature_missing"
    return status


def blocks_inside_exons(tx: Transcript) -> bool:
    exons = tx.sorted_exons()
    if not exons or not tx.sorted_cds():
        return True
    for cds_start, cds_end, _phase in tx.sorted_cds():
        if not any(exon_start <= cds_start and cds_end <= exon_end for exon_start, exon_end in exons):
            return False
    return True


def same_source_neighbor_conflicts(records: Sequence[Dict[str, object]], threshold: float) -> Dict[str, List[str]]:
    conflicts: Dict[str, List[str]] = defaultdict(list)
    if threshold <= 0:
        return conflicts
    by_source_seqid_strand: Dict[Tuple[str, str, str], List[Dict[str, object]]] = defaultdict(list)
    for record in records:
        by_source_seqid_strand[(str(record["source"]), str(record["seqid"]), str(record["strand"]))].append(record)
    for rows in by_source_seqid_strand.values():
        rows = sorted(rows, key=lambda row: (int(row["start"]), int(row["end"]), str(row["model_id"])))
        for i, left in enumerate(rows):
            left_start = int(left["start"])
            left_end = int(left["end"])
            left_len = max(1, left_end - left_start + 1)
            for right in rows[i + 1:]:
                right_start = int(right["start"])
                right_end = int(right["end"])
                if right_start > left_end:
                    break
                if str(left["gene_id"]) == str(right["gene_id"]):
                    continue
                overlap = min(left_end, right_end) - max(left_start, right_start) + 1
                if overlap <= 0:
                    continue
                right_len = max(1, right_end - right_start + 1)
                if overlap / left_len >= threshold or overlap / right_len >= threshold:
                    reason = "same_source_different_gene_high_overlap:{0}:overlap_bp_{1}".format(right["model_id"], overlap)
                    conflicts[str(left["model_id"])].append(reason)
                    reason = "same_source_different_gene_high_overlap:{0}:overlap_bp_{1}".format(left["model_id"], overlap)
                    conflicts[str(right["model_id"])].append(reason)
    return conflicts


def validate_transcript(
    model_id: str,
    source: str,
    gene_id: str,
    tx: Transcript,
    source_file: str,
    genome: Dict[str, str],
    junctions: Set[Tuple[str, int, int, str]],
    parent_status_by_tx: Dict[str, str],
    allow_partial_models: bool,
) -> Dict[str, object]:
    hard: List[str] = []
    warnings: List[str] = []
    exons = tx.sorted_exons()
    cds_blocks = tx.sorted_cds()
    cds_seq, missing_cds_sequence = extract_cds_sequence(genome, tx)
    cds_length = len(cds_seq)

    coordinate_status = "ok"
    if tx.seqid not in genome:
        hard.append("seqid_not_found_in_genome")
        coordinate_status = "fail"
    elif tx.start < 1 or tx.end > len(genome[tx.seqid]) or tx.start > tx.end:
        hard.append("transcript_coordinates_out_of_bounds")
        coordinate_status = "fail"
    if tx.strand not in {"+", "-"}:
        hard.append("invalid_or_missing_strand")
        coordinate_status = "fail"
    if not exons:
        hard.append("no_exon_features")
    if not cds_blocks:
        hard.append("no_CDS_features")
    if not blocks_inside_exons(tx):
        hard.append("CDS_not_nested_within_exon")
        coordinate_status = "fail"
    if missing_cds_sequence:
        hard.append("CDS_sequence_unavailable_or_out_of_bounds")
        coordinate_status = "fail"

    parent_status = parent_status_by_tx.get(tx.transcript_id, "parent_status_not_observed")
    if parent_status == "child_parent_transcript_feature_missing":
        hard.append(parent_status)
    elif parent_status == "parent_status_not_observed":
        warnings.append(parent_status)

    if cds_length:
        cds_length_mod3 = "valid" if cds_length % 3 == 0 else "invalid"
        if cds_length % 3 != 0:
            hard.append("CDS_length_mod3_invalid")
        start_codon = cds_seq[:3] if len(cds_seq) >= 3 else ""
        stop_codon = cds_seq[-3:] if len(cds_seq) >= 3 else ""
        start_codon_status = "present" if start_codon == "ATG" else "absent_or_noncanonical"
        stop_codon_status = "present" if stop_codon in STOP_CODONS else "absent_or_noncanonical"
        if start_codon_status != "present":
            warnings.append("start_codon_absent_or_noncanonical")
        if stop_codon_status != "present":
            warnings.append("stop_codon_absent_or_noncanonical")
        codons = [cds_seq[i:i + 3] for i in range(0, max(0, len(cds_seq) - 3), 3)]
        internal_stop_count = sum(1 for codon in codons if codon in STOP_CODONS)
        if internal_stop_count:
            hard.append("internal_stop_present")
    else:
        cds_length_mod3 = "not_assessable"
        start_codon_status = "not_assessable"
        stop_codon_status = "not_assessable"
        internal_stop_count = 0

    phase_status, phase_hard, phase_warnings = check_phase(tx)
    hard.extend(phase_hard)
    warnings.extend(phase_warnings)
    if phase_status == "phase_inconsistent":
        hard.append("phase_inconsistent")

    splice_status, noncanonical_count, supported_noncanonical_count, splice_hard, splice_warnings = check_splice_motifs(genome, tx, junctions)
    hard.extend(splice_hard)
    warnings.extend(splice_warnings)

    status = "eligible" if not hard else "ineligible"
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
        "eligible_for_ai": "true" if status == "eligible" else "false",
        "eligibility_status": status,
        "hard_fail_reasons": ",".join(sorted(set(hard))),
        "warning_reasons": ",".join(sorted(set(warnings))),
        "cds_length_mod3": cds_length_mod3,
        "start_codon_status": start_codon_status,
        "stop_codon_status": stop_codon_status,
        "internal_stop_count": internal_stop_count,
        "phase_status": phase_status,
        "splice_motif_status": splice_status,
        "noncanonical_splice_count": noncanonical_count,
        "noncanonical_splice_junction_supported_count": supported_noncanonical_count,
        "coordinate_status": coordinate_status,
        "parent_status": parent_status,
        "neighbor_conflict_status": "not_assessed",
        "source_file": source_file,
    }


def audit_gff(
    path: str,
    source: str,
    genome: Dict[str, str],
    junctions: Set[Tuple[str, int, int, str]],
    gene_ids: Set[str],
    regions: Sequence[Tuple[str, int, int]],
    allow_partial_models: bool,
) -> List[Dict[str, object]]:
    _genes, transcripts, tx_to_gene = parse_annotation_gff(path, source)
    parent_status_by_tx = audit_raw_parent_status(path)
    rows: List[Dict[str, object]] = []
    for tx_id, tx in sorted(transcripts.items()):
        gene_id = tx_to_gene.get(tx_id, tx.gene_id) or tx_id
        if not tx.seqid or not tx.start or not tx.end:
            continue
        model_id = "{0}:{1}:{2}".format(safe_id(source), safe_id(gene_id), safe_id(tx_id))
        if not transcript_selected(tx, model_id, gene_id, gene_ids, regions):
            continue
        rows.append(validate_transcript(model_id, source, gene_id, tx, path, genome, junctions, parent_status_by_tx, allow_partial_models))
    return rows


def apply_neighbor_conflicts(rows: List[Dict[str, object]], threshold: float) -> None:
    conflicts = same_source_neighbor_conflicts(rows, threshold)
    for row in rows:
        model_id = str(row["model_id"])
        if model_id not in conflicts:
            row["neighbor_conflict_status"] = "none"
            continue
        existing = [item for item in str(row.get("hard_fail_reasons", "")).split(",") if item]
        all_reasons = sorted(set(existing + conflicts[model_id]))
        row["hard_fail_reasons"] = ",".join(all_reasons)
        row["eligible_for_ai"] = "false"
        row["eligibility_status"] = "ineligible"
        row["neighbor_conflict_status"] = "same_source_high_overlap_conflict"


def build_summary(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    counts: Dict[Tuple[str, str, str], int] = defaultdict(int)
    for row in rows:
        source = str(row.get("source", ""))
        status = str(row.get("eligibility_status", ""))
        reasons = [item for item in str(row.get("hard_fail_reasons", "")).split(",") if item] or ["none"]
        for reason in reasons:
            counts[(source, status, reason)] += 1
    return [
        {"source": source, "eligibility_status": status, "reason": reason, "count": count}
        for (source, status, reason), count in sorted(counts.items())
    ]


def main() -> None:
    args = parse_args()
    genome = load_fasta(args.genome_fasta)
    junctions = load_junctions(args.junctions, args.min_junction_support_count, args.min_junction_sample_count)
    gene_ids = load_gene_filter(args.gene_list)
    regions = load_regions(args.region_table, args.region_flank)

    rows: List[Dict[str, object]] = []
    for source, path in [("current", args.current_gff)] + list(args.tool_gff):
        rows.extend(audit_gff(path, source, genome, junctions, gene_ids, regions, args.allow_partial_models))
    apply_neighbor_conflicts(rows, args.same_source_conflict_fraction)
    write_tsv(args.output_eligibility, FIELDS, rows)
    if args.output_summary:
        write_tsv(args.output_summary, SUMMARY_FIELDS, build_summary(rows))


if __name__ == "__main__":
    main()
