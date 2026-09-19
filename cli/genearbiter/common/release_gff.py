#!/usr/bin/env python3
# script_id_md5: 5d35a7ed6f3f4eb7a3db9e865b722a3c
# created: 2026-07-14
# modified: 2026-07-14
# owner: project
# status: project_code
# purpose: 对 GeneArbiter 命名层 GFF 执行 release-like 结构检查和保守修复，输出可发布兼容 GFF。
# inputs: policy_named_annotation.gff3 或 final_named_annotation.gff3。
# outputs: 精简溯源版 release_annotation.gff3；结构版 release_annotation.clean.gff3；release_gff_qc_summary.tsv；release_gff_id_mapping.tsv；release_gff_report.md。
# notes: 不读取 truth annotation；不重跑 AI；不发明新候选坐标。新增 gene/exon 只由已有 transcript/CDS/UTR 坐标推导。

"""Validate and repair GeneArbiter GFF3 for release-like downstream use.

The upstream GeneArbiter proposal and policy exporters preserve source GFF dialects
where possible. That is useful internally, but some current/backbone files
do not use a release-like ``gene -> transcript -> exon/CDS`` hierarchy. This
post-processing step creates a stricter compatibility layer:

* top-level transcript records are promoted to a gene plus renamed transcript;
* transcript spans are expanded to cover their children;
* missing exon rows are synthesized from existing exon-like/CDS intervals;
* invalid phase values are made syntactically legal and reported;
* Parent references are rewritten after transcript/gene ID repair;
* public clean and compact-trace profiles share identical structures and IDs.

The script does not inspect a genome FASTA and does not infer coding frames. CDS
phase repair is therefore a syntax repair, not biological frame validation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, MutableMapping, Sequence, Tuple
from urllib.parse import quote, unquote


TRANSCRIPT_FEATURES = {
    "mrna",
    "transcript",
    "rna",
    "lnc_rna",
    "ncrna",
    "rrna",
    "trna",
    "mirna",
    "snrna",
    "snorna",
    "pre_mirna",
    "primary_transcript",
}
EXON_LIKE_FEATURES = {"exon", "five_prime_utr", "three_prime_utr", "utr"}
CHILD_FEATURES = EXON_LIKE_FEATURES | {"cds", "start_codon", "stop_codon"}
PUBLIC_SOURCE = "GeneArbiter"
COMMON_GENE_ATTRS = ("ID", "Name", "biotype")
COMMON_TRANSCRIPT_ATTRS = ("ID", "Parent", "Name", "biotype")
COMMON_CDS_ATTRS = ("ID", "Parent", "protein_id")
COMMON_CHILD_ATTRS = ("ID", "Parent")
RELATION_TO_CHANGE_TYPE = {
    "one_to_one_replacement": "replacement",
    "split_replacement": "split",
    "merge_replacement": "merge",
    "complex_replacement": "complex",
    "novel_gene": "novel",
}
FEATURE_RANK = {
    "gene": 0,
    "mrna": 1,
    "transcript": 1,
    "rna": 1,
    "lnc_rna": 1,
    "ncrna": 1,
    "rrna": 1,
    "trna": 1,
    "mirna": 1,
    "snrna": 1,
    "snorna": 1,
    "exon": 2,
    "five_prime_utr": 3,
    "three_prime_utr": 3,
    "utr": 3,
    "cds": 4,
    "start_codon": 5,
    "stop_codon": 5,
}
TRUTH_LIKE_RE = re.compile(r"truth|manual", re.IGNORECASE)


@dataclass
class Record:
    seqid: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: Dict[str, str]
    idx: int
    synthetic: bool = False
    raw_feature: str = ""

    @property
    def feature_l(self) -> str:
        return self.feature.lower()

    @property
    def row_id(self) -> str:
        return self.attrs.get("ID", "")

    @property
    def parents(self) -> List[str]:
        return split_list(self.attrs.get("Parent", ""))


@dataclass
class TranscriptBlock:
    tx_id: str
    record: Record
    gene_id: str = ""
    child_records: List[Record] = field(default_factory=list)

    def child_span(self) -> Tuple[int, int]:
        starts = [self.record.start]
        ends = [self.record.end]
        starts.extend(child.start for child in self.child_records)
        ends.extend(child.end for child in self.child_records)
        return min(starts), max(ends)

    def exon_intervals(self) -> List[Tuple[int, int]]:
        intervals = [(child.start, child.end) for child in self.child_records if child.feature_l == "exon"]
        return merge_intervals(intervals)

    def exon_seed_intervals(self) -> List[Tuple[int, int]]:
        intervals = [
            (child.start, child.end)
            for child in self.child_records
            if child.feature_l in EXON_LIKE_FEATURES or child.feature_l == "cds"
        ]
        return merge_intervals(intervals)


@dataclass
class GeneBlock:
    gene_id: str
    record: Record
    transcript_ids: List[str] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-gff", required=True, help="Named GeneArbiter GFF3 to repair.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--output-gff-name", default="release_annotation.gff3")
    parser.add_argument("--clean-gff-name", default="release_annotation.clean.gff3")
    parser.add_argument("--summary-name", default="release_gff_qc_summary.tsv")
    parser.add_argument("--mapping-name", default="release_gff_id_mapping.tsv")
    parser.add_argument("--report-name", default="release_gff_report.md")
    parser.add_argument("--source-label", default="", help="Optional source label for synthesized gene/exon rows.")
    parser.add_argument("--transcript-template", default="{gene_id}.t{index:02d}")
    parser.add_argument("--exon-template", default="{transcript_id}.exon{index:03d}")
    parser.add_argument("--clean", action="store_true", default=True, help="Write a clean GFF without GeneArbiter internal attributes.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    parser.add_argument("--fail-on-truth-path", action="store_true", help="Refuse truth/manual-like paths.")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open(encoding="utf-8")


def split_list(text: object) -> List[str]:
    return [item.strip() for item in str(text or "").replace(";", ",").split(",") if item.strip()]


def attr_value(text: object) -> str:
    value = str(text or "")
    for _ in range(3):
        decoded = unquote(value)
        if decoded == value:
            break
        value = decoded
    return quote(value, safe="._:-|,")


def attrs_to_text(items: Iterable[Tuple[str, object]]) -> str:
    return ";".join("{0}={1}".format(key, attr_value(value)) for key, value in items if value not in {"", None})


def clean_attrs(record: Record) -> Dict[str, str]:
    """Return the shared, feature-aware attributes for both public GFF profiles."""
    if record.feature_l == "gene":
        allowed = COMMON_GENE_ATTRS
    elif record.feature_l in TRANSCRIPT_FEATURES:
        allowed = COMMON_TRANSCRIPT_ATTRS
    elif record.feature_l == "cds":
        allowed = COMMON_CDS_ATTRS
    else:
        allowed = COMMON_CHILD_ATTRS

    attrs = dict(record.attrs)
    if record.feature_l == "gene" and not attrs.get("biotype"):
        attrs["biotype"] = attrs.get("gene_biotype", "")
    elif record.feature_l in TRANSCRIPT_FEATURES and not attrs.get("biotype"):
        attrs["biotype"] = attrs.get("transcript_biotype", "")

    result = {key: attrs[key] for key in allowed if attrs.get(key)}
    if result.get("Name") and unquote(result["Name"]) == unquote(result.get("ID", "")):
        del result["Name"]
    return result


def trace_review_reason(record: Record) -> str:
    """Collapse internal review/risk text into a short controlled public vocabulary."""
    values = [
        record.attrs.get("review_flags", ""),
        record.attrs.get("risk_tags", ""),
        record.attrs.get("blocking_reasons", ""),
        record.attrs.get("final_call", ""),
    ]
    blob = " ".join(unquote(value).lower() for value in values if value)
    if not blob:
        return ""

    reasons: List[str] = []
    checks = [
        ("novel_candidate", ("manual_review_novel", "novel_gene_candidate")),
        ("unsupported_locus", ("unsupported_locus", "unsupported_junction")),
        ("junction_not_assessed", ("not_assessed", "not_junction_assessable")),
        ("partial_junction_support", ("partial_junction",)),
        ("complex_structure", ("over_split", "fusion", "fragment", "complex_locus")),
        ("model_conflict", ("conflict", "disagree", "mutually_exclusive")),
        ("limited_evidence", ("uncertain", "weak", "low_confidence", "insufficient")),
    ]
    for reason, markers in checks:
        if any(marker in blob for marker in markers):
            reasons.append(reason)
    if not reasons:
        reasons.append("workflow_flag")
    return ",".join(reasons)


def trace_attrs(record: Record) -> Dict[str, str]:
    """Return compact provenance/review attributes for the public trace GFF."""
    result = clean_attrs(record)
    if record.feature_l == "gene":
        relation = record.attrs.get("final_id_relation", "")
        result["change_type"] = RELATION_TO_CHANGE_TYPE.get(relation, "unchanged")
        origin_source = record.attrs.get("selected_source", "")
        if origin_source:
            result["origin_source"] = origin_source
        if record.attrs.get("source_gene_id"):
            result["origin_gene_id"] = record.attrs["source_gene_id"]
        current_ids = record.attrs.get("current_gene_ids") or record.attrs.get("replaces_current_genes", "")
        if current_ids:
            result["current_gene_ids"] = current_ids
        if record.attrs.get("retired_current_gene_ids"):
            result["retired_current_gene_ids"] = record.attrs["retired_current_gene_ids"]
        review_reason = trace_review_reason(record)
        if review_reason:
            result["review_recommended"] = "true"
            result["review_reason"] = review_reason
    elif record.feature_l in TRANSCRIPT_FEATURES and record.attrs.get("source_transcript_id"):
        result["origin_transcript_id"] = record.attrs["source_transcript_id"]
    return result


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    merged: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def natural_seq_key(seqid: str) -> Tuple[object, ...]:
    return tuple(int(token) if token.isdigit() else token for token in re.split(r"(\d+)", seqid))


def fail_if_truth_like(paths: Sequence[Path]) -> None:
    for path in paths:
        if TRUTH_LIKE_RE.search(str(path)):
            raise SystemExit("Refusing truth/manual-like path for release GFF QC: {0}".format(path))


def read_gff(path: Path) -> Tuple[List[str], List[Record], Counter]:
    headers: List[str] = []
    records: List[Record] = []
    counts: Counter = Counter()
    with open_text(path) as handle:
        for idx, line in enumerate(handle):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                counts["malformed_rows"] += 1
                continue
            seqid, source, feature, start_s, end_s, score, strand, phase, attr_text = parts
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                counts["bad_coordinate_rows"] += 1
                continue
            if start < 1 or end < start:
                counts["bad_coordinate_rows"] += 1
                continue
            records.append(Record(seqid, source, feature, start, end, score or ".", strand or ".", phase or ".", parse_attrs(attr_text), idx, raw_feature=feature))
            counts["input_rows"] += 1
            counts["input_feature:" + feature] += 1
    return headers, records, counts


def unique_id(base: str, used: set[str]) -> str:
    candidate = base
    idx = 1
    while not candidate or candidate in used:
        idx += 1
        candidate = "{0}.dup{1:02d}".format(base or "genearbiter_feature", idx)
    used.add(candidate)
    return candidate


def feature_type(feature_l: str) -> str:
    if feature_l == "gene":
        return "gene"
    if feature_l in TRANSCRIPT_FEATURES:
        return "transcript"
    if feature_l in CHILD_FEATURES:
        return "child"
    return "other"


def allocate_transcript_id(gene_id: str, index: int, template: str, used: set[str]) -> str:
    preferred = template.format(gene_id=gene_id, index=index)
    return unique_id(preferred, used)


def build_release_records(records: List[Record], transcript_template: str, exon_template: str, source_label: str) -> Tuple[List[Record], List[Dict[str, str]], Counter]:
    counts: Counter = Counter()
    mapping_rows: List[Dict[str, str]] = []
    original_gene_ids = {rec.row_id for rec in records if rec.feature_l == "gene" and rec.row_id}
    used_gene_ids: set[str] = set()
    used_tx_ids: set[str] = set()
    used_child_ids: set[str] = set()

    gene_records: Dict[str, GeneBlock] = {}
    transcript_records: Dict[str, TranscriptBlock] = {}
    child_records: List[Record] = []
    other_records: List[Record] = []

    for rec in records:
        ftype = feature_type(rec.feature_l)
        if ftype == "gene":
            old_id = rec.row_id or "gene_{0}".format(rec.idx)
            new_id = unique_id(old_id, used_gene_ids)
            if new_id != rec.row_id:
                rec.attrs["original_release_gene_id"] = rec.row_id
                rec.attrs["ID"] = new_id
                counts["renamed_duplicate_or_missing_gene_ids"] += 1
                mapping_rows.append(mapping_row("gene", rec.row_id, new_id, "", "rename_duplicate_or_missing_gene_id"))
            else:
                rec.attrs["ID"] = new_id
            gene_records[new_id] = GeneBlock(new_id, rec, [])
        elif ftype == "transcript":
            old_id = rec.row_id or "transcript_{0}".format(rec.idx)
            parent_ids = rec.parents
            gene_id = parent_ids[0] if parent_ids else ""
            transcript_records[old_id] = TranscriptBlock(old_id, rec, gene_id, [])
        elif ftype == "child":
            child_records.append(rec)
        else:
            other_records.append(rec)

    tx_old_to_new: Dict[str, str] = {}
    tx_old_to_gene: Dict[str, str] = {}
    tx_index_by_gene: Counter = Counter()

    for old_tx_id, block in sorted(transcript_records.items(), key=lambda item: item[1].record.idx):
        rec = block.record
        gene_id = block.gene_id
        if not gene_id:
            proposed_gene_id = old_tx_id if old_tx_id not in original_gene_ids and old_tx_id not in gene_records else old_tx_id + ".gene"
            gene_id = unique_id(proposed_gene_id, used_gene_ids)
            gene_rec = Record(
                rec.seqid,
                source_label or rec.source,
                "gene",
                rec.start,
                rec.end,
                ".",
                rec.strand,
                ".",
                {"ID": gene_id, "Name": gene_id, "genearbiter_release_repair": "synthetic_gene_from_top_level_transcript", "source_transcript_id": old_tx_id},
                rec.idx - 0.1,  # type: ignore[arg-type]
                synthetic=True,
                raw_feature="gene",
            )
            gene_records[gene_id] = GeneBlock(gene_id, gene_rec, [])
            counts["synthetic_gene_records"] += 1
            mapping_rows.append(mapping_row("gene", "", gene_id, old_tx_id, "synthetic_gene_from_top_level_transcript"))
        elif gene_id not in gene_records:
            gene_id = unique_id(gene_id, used_gene_ids)
            gene_rec = Record(
                rec.seqid,
                source_label or rec.source,
                "gene",
                rec.start,
                rec.end,
                ".",
                rec.strand,
                ".",
                {"ID": gene_id, "Name": gene_id, "genearbiter_release_repair": "synthetic_gene_from_transcript_parent"},
                rec.idx - 0.1,  # type: ignore[arg-type]
                synthetic=True,
                raw_feature="gene",
            )
            gene_records[gene_id] = GeneBlock(gene_id, gene_rec, [])
            counts["synthetic_gene_records"] += 1
            mapping_rows.append(mapping_row("gene", "", gene_id, old_tx_id, "synthetic_gene_from_missing_parent_gene"))

        tx_index_by_gene[gene_id] += 1
        if old_tx_id == gene_id or old_tx_id in used_tx_ids or old_tx_id in used_gene_ids:
            new_tx_id = allocate_transcript_id(gene_id, tx_index_by_gene[gene_id], transcript_template, used_tx_ids)
            rec.attrs["original_release_transcript_id"] = old_tx_id
            rec.attrs["ID"] = new_tx_id
            counts["renamed_transcript_ids"] += 1
            mapping_rows.append(mapping_row("transcript", old_tx_id, new_tx_id, gene_id, "rename_transcript_for_release_gene_hierarchy"))
        else:
            new_tx_id = unique_id(old_tx_id, used_tx_ids)
            if new_tx_id != old_tx_id:
                rec.attrs["original_release_transcript_id"] = old_tx_id
                rec.attrs["ID"] = new_tx_id
                counts["renamed_transcript_ids"] += 1
                mapping_rows.append(mapping_row("transcript", old_tx_id, new_tx_id, gene_id, "rename_duplicate_transcript_id"))
        rec.attrs["Parent"] = gene_id
        tx_old_to_new[old_tx_id] = new_tx_id
        tx_old_to_gene[old_tx_id] = gene_id
        block.gene_id = gene_id
        gene_records[gene_id].transcript_ids.append(old_tx_id)

    children_by_tx: Dict[str, List[Record]] = defaultdict(list)
    for rec in child_records:
        parents = rec.parents
        new_parents: List[str] = []
        for parent in parents:
            if parent in tx_old_to_new:
                new_parents.append(tx_old_to_new[parent])
                children_by_tx[parent].append(rec)
            elif parent in gene_records:
                new_parents.append(parent)
                counts["child_parent_is_gene"] += 1
            else:
                counts["missing_parent_refs"] += 1
        if new_parents:
            rec.attrs["Parent"] = ",".join(dict.fromkeys(new_parents))
        if rec.feature_l == "cds" and rec.phase not in {"0", "1", "2"}:
            rec.attrs["original_release_phase"] = rec.phase
            rec.phase = "0"
            counts["fixed_invalid_cds_phase_to_zero"] += 1
        elif rec.feature_l != "cds" and rec.phase != ".":
            rec.attrs["original_release_phase"] = rec.phase
            rec.phase = "."
            counts["fixed_non_cds_phase_to_dot"] += 1
        if rec.row_id:
            if rec.feature_l in {"exon", "cds"}:
                new_id = unique_id(rec.row_id, used_child_ids)
                if new_id != rec.row_id:
                    rec.attrs["original_release_feature_id"] = rec.row_id
                    rec.attrs["ID"] = new_id
                    counts["renamed_duplicate_child_ids"] += 1
                    mapping_rows.append(mapping_row(rec.feature, rec.row_id, new_id, rec.attrs.get("Parent", ""), "rename_duplicate_child_id"))

    parent_records: Dict[str, Record] = {block.record.row_id: block.record for block in transcript_records.values()}
    parent_records.update({gene_id: block.record for gene_id, block in gene_records.items()})
    for rec in child_records:
        parent_strands = {parent_records[parent].strand for parent in rec.parents if parent in parent_records}
        if len(parent_strands) == 1 and rec.strand not in parent_strands:
            rec.attrs["original_release_strand"] = rec.strand
            rec.strand = next(iter(parent_strands))
            counts["fixed_child_strand_to_parent"] += 1

    for old_tx_id, block in transcript_records.items():
        block.child_records = children_by_tx.get(old_tx_id, [])
        start, end = block.child_span()
        if start != block.record.start or end != block.record.end:
            block.record.attrs["original_release_span"] = "{0}-{1}".format(block.record.start, block.record.end)
            block.record.start = start
            block.record.end = end
            counts["expanded_transcript_spans"] += 1
        if not block.exon_intervals():
            seeds = block.exon_seed_intervals()
            for exon_idx, (start, end) in enumerate(seeds, start=1):
                new_tx_id = tx_old_to_new[old_tx_id]
                exon_id = unique_id(exon_template.format(transcript_id=new_tx_id, index=exon_idx), used_child_ids)
                exon = Record(
                    block.record.seqid,
                    source_label or block.record.source,
                    "exon",
                    start,
                    end,
                    ".",
                    block.record.strand,
                    ".",
                    {"ID": exon_id, "Parent": new_tx_id, "genearbiter_release_repair": "synthetic_exon_from_cds_or_utr"},
                    block.record.idx + exon_idx / 100000.0,  # type: ignore[arg-type]
                    synthetic=True,
                    raw_feature="exon",
                )
                block.child_records.append(exon)
                counts["synthetic_exon_records"] += 1
                mapping_rows.append(mapping_row("exon", "", exon_id, new_tx_id, "synthetic_exon_from_cds_or_utr"))

    for gene_id, block in gene_records.items():
        starts = [block.record.start]
        ends = [block.record.end]
        for old_tx_id in block.transcript_ids:
            tx = transcript_records.get(old_tx_id)
            if tx:
                starts.append(tx.record.start)
                ends.append(tx.record.end)
                starts.extend(child.start for child in tx.child_records)
                ends.extend(child.end for child in tx.child_records)
        new_start, new_end = min(starts), max(ends)
        if new_start != block.record.start or new_end != block.record.end:
            block.record.attrs["original_release_span"] = "{0}-{1}".format(block.record.start, block.record.end)
            block.record.start = new_start
            block.record.end = new_end
            counts["expanded_gene_spans"] += 1

    release_records: List[Record] = []
    emitted_records: set[int] = set()

    def append_once(record: Record) -> None:
        """A multi-parent child belongs to several transcript blocks but is one GFF row."""
        object_key = id(record)
        if object_key in emitted_records:
            return
        emitted_records.add(object_key)
        release_records.append(record)

    for gene_id, gene in gene_records.items():
        append_once(gene.record)
        for old_tx_id in gene.transcript_ids:
            tx = transcript_records.get(old_tx_id)
            if not tx:
                continue
            append_once(tx.record)
            for child in tx.child_records:
                append_once(child)
    for record in other_records:
        append_once(record)

    release_records.sort(key=record_sort_key)
    counts.update(validate_release_records(release_records, source_label))
    counts["output_rows"] = len(release_records)
    return release_records, mapping_rows, counts


def mapping_row(feature: str, old_id: str, new_id: str, parent_id: str, reason: str) -> Dict[str, str]:
    return {"feature": feature, "old_id": old_id, "new_id": new_id, "parent_id": parent_id, "reason": reason}


def record_sort_key(rec: Record) -> Tuple[object, int, int, int, object]:
    return (
        natural_seq_key(rec.seqid),
        rec.start,
        FEATURE_RANK.get(rec.feature_l, 9),
        rec.end,
        rec.idx,
    )


def validate_release_records(records: Sequence[Record], generated_source_label: str = "") -> Counter:
    counts: Counter = Counter()
    ids: Dict[str, Record] = {}
    for rec in records:
        counts["output_feature:" + rec.feature] += 1
        rid = rec.row_id
        if rid:
            if rid in ids and rec.feature_l in {"gene"} | TRANSCRIPT_FEATURES | {"exon"}:
                counts["postcheck_duplicate_gene_tx_exon_ids"] += 1
            ids.setdefault(rid, rec)
        if rec.feature_l == "cds" and rec.phase not in {"0", "1", "2"}:
            counts["postcheck_invalid_cds_phase"] += 1
        elif rec.feature_l != "cds" and rec.phase != ".":
            counts["postcheck_invalid_non_cds_phase"] += 1
    for rec in records:
        for parent in rec.parents:
            if parent not in ids:
                counts["postcheck_missing_parent_refs"] += 1
                continue
            par = ids[parent]
            known_strand_conflict = rec.strand in {"+", "-"} and par.strand in {"+", "-"} and rec.strand != par.strand
            hierarchy_conflict = rec.seqid != par.seqid or known_strand_conflict or rec.start < par.start or rec.end > par.end
            if hierarchy_conflict:
                generated = bool(generated_source_label) and (
                    rec.source == generated_source_label or par.source == generated_source_label
                )
                if generated or rec.synthetic or par.synthetic:
                    counts["postcheck_child_outside_parent"] += 1
                else:
                    counts["postcheck_inherited_source_child_outside_parent"] += 1
    return counts


def write_gff(path: Path, records: Sequence[Record], input_gff: Path, clean: bool) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for rec in records:
            row_attrs: MutableMapping[str, str] = clean_attrs(rec) if clean else trace_attrs(rec)
            handle.write(
                "\t".join(
                    [
                        rec.seqid,
                        PUBLIC_SOURCE,
                        rec.feature,
                        str(rec.start),
                        str(rec.end),
                        rec.score or ".",
                        rec.strand or ".",
                        rec.phase or ".",
                        attrs_to_text(row_attrs.items()) or ".",
                    ]
                )
                + "\n"
            )


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_report(path: Path, args: argparse.Namespace, summary_rows: Sequence[Dict[str, object]]) -> None:
    lines = [
        "# GeneArbiter Release GFF QC Report",
        "",
        "Generated: {0}".format(datetime.now(timezone.utc).isoformat()),
        "",
        "- input_gff: `{0}`".format(args.input_gff),
        "- output_gff: `{0}`".format(Path(args.out_dir) / args.output_gff_name),
        "- clean_gff: `{0}`".format(Path(args.out_dir) / args.clean_gff_name),
        "- truth_usage: none",
        "",
        "## Summary",
        "",
    ]
    for row in summary_rows:
        lines.append("- {0}: {1}".format(row["metric"], row["value"]))
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- This step repairs GFF3 structure for release-like compatibility.",
            "- New gene records are inferred only from existing transcript spans.",
            "- New exon records are inferred only from existing exon-like/CDS intervals.",
            "- CDS phase syntax repair does not validate biological reading frame.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_gff = Path(args.input_gff)
    out_dir = Path(args.out_dir)
    if args.fail_on_truth_path:
        fail_if_truth_like([input_gff, out_dir])
    if not input_gff.exists():
        raise SystemExit("Missing input GFF: {0}".format(input_gff))
    out_dir.mkdir(parents=True, exist_ok=True)

    output_gff = out_dir / args.output_gff_name
    clean_gff = out_dir / args.clean_gff_name if args.clean_gff_name else None
    summary_path = out_dir / args.summary_name
    mapping_path = out_dir / args.mapping_name
    report_path = out_dir / args.report_name
    outputs = [output_gff, summary_path, mapping_path, report_path]
    if clean_gff:
        outputs.append(clean_gff)
    if not args.force:
        existing = [str(path) for path in outputs if path.exists()]
        if existing:
            raise SystemExit("Output exists; pass --force to overwrite: {0}".format(",".join(existing)))

    _headers, records, input_counts = read_gff(input_gff)
    release_records, mapping_rows, repair_counts = build_release_records(
        records,
        args.transcript_template,
        args.exon_template,
        args.source_label,
    )
    counts = Counter()
    counts.update(input_counts)
    counts.update(repair_counts)

    summary_rows = [{"metric": key, "value": value} for key, value in sorted(counts.items())]
    write_gff(output_gff, release_records, input_gff, clean=False)
    if clean_gff:
        write_gff(clean_gff, release_records, input_gff, clean=True)
    write_tsv(summary_path, ["metric", "value"], summary_rows)
    write_tsv(mapping_path, ["feature", "old_id", "new_id", "parent_id", "reason"], mapping_rows)
    write_report(report_path, args, summary_rows)

    print("release_gff\t{0}".format(output_gff))
    if clean_gff:
        print("release_clean_gff\t{0}".format(clean_gff))
    print("summary\t{0}".format(summary_path))
    print("mapping\t{0}".format(mapping_path))
    print("synthetic_gene_records\t{0}".format(counts.get("synthetic_gene_records", 0)))
    print("synthetic_exon_records\t{0}".format(counts.get("synthetic_exon_records", 0)))
    print("postcheck_missing_parent_refs\t{0}".format(counts.get("postcheck_missing_parent_refs", 0)))
    print("postcheck_invalid_cds_phase\t{0}".format(counts.get("postcheck_invalid_cds_phase", 0)))


if __name__ == "__main__":
    main()
