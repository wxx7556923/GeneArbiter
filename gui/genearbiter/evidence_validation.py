#!/usr/bin/env python3
# script_id_md5: 4cb3e6f4df64517bf45c74d879b4a063
# created: 2026-07-08
# modified: 2026-07-09
# owner: project
# status: project_code
# purpose: 校验 GeneArbiter evidence TSV/GFF 的最小 schema 和关键坐标字段。
# inputs: short_read_evidence directory or individual evidence TSV/protein GFF paths。
# outputs: validation messages on stdout and process exit status。
# notes: 只做输入结构和基础坐标检查；多份 protein GFF 只报警，因底层会合并且不去重。

"""Evidence schema validation for GeneArbiter."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


REQUIRED_COLUMNS = {
    "sample_manifest": ["sample_id"],
    "merged_splice_junctions": [
        "chrom",
        "intron_start_1based",
        "intron_end_1based",
        "strand",
        "junction_read_count",
        "sample_support_count",
        "supporting_samples",
    ],
    "junction_support_by_model": [
        "transcript_id",
        "gene_id",
        "n_model_junctions",
        "n_supported_junctions",
        "junction_support_fraction",
        "unsupported_junctions",
        "novel_supported_junctions_nearby",
        "donor_supported",
        "acceptor_supported",
        "full_intron_chain_supported",
    ],
}


def read_header(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        first = handle.readline().rstrip("\n")
    return first.split("\t") if first else []


def check_columns(path: Path, table_kind: str) -> List[str]:
    messages: List[str] = []
    if not path.exists():
        return [f"missing:{table_kind}:{path}"]
    header = read_header(path)
    missing = [col for col in REQUIRED_COLUMNS[table_kind] if col not in header]
    if missing:
        messages.append(f"missing_columns:{table_kind}:{','.join(missing)}:{path}")
    else:
        messages.append(f"ok_columns:{table_kind}:{path}")
    return messages


def iter_rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def check_junction_coordinates(path: Path, limit: int = 20) -> List[str]:
    messages: List[str] = []
    if not path.exists():
        return messages
    bad = 0
    total = 0
    for row in iter_rows(path):
        total += 1
        try:
            start = int(row.get("intron_start_1based") or 0)
            end = int(row.get("intron_end_1based") or 0)
            bed_start = row.get("bed_start_0based")
            bed_end = row.get("bed_end_0based")
            strand = row.get("strand") or ""
            read_count = int(row.get("junction_read_count") or 0)
            sample_count = int(row.get("sample_support_count") or 0)
            sample_list = [item.strip() for item in (row.get("supporting_samples") or "").split(",") if item.strip()]
            sample_ids = sorted(set(sample_list))
            if start <= 0 or end <= 0 or start >= end:
                bad += 1
            elif bed_start not in {None, ""} and int(bed_start) != start - 1:
                bad += 1
            elif bed_end not in {None, ""} and int(bed_end) != end:
                bad += 1
            elif strand not in {"+", "-", "."} or read_count < 0 or sample_count < 1:
                bad += 1
            elif sample_list != sample_ids or len(sample_ids) != sample_count:
                bad += 1
        except ValueError:
            bad += 1
        if bad >= limit:
            break
    if bad:
        messages.append(f"bad_junction_coordinates:{bad}_examples_or_more:{path}")
    else:
        messages.append(f"ok_junction_coordinates:rows={total}:{path}")
    return messages


def check_numeric_bounds(path: Path, table_kind: str) -> List[str]:
    messages: List[str] = []
    if not path.exists():
        return messages
    bad = 0
    total = 0
    fraction_columns = ["junction_support_fraction"]
    for row in iter_rows(path):
        total += 1
        for col in fraction_columns:
            value = row.get(col, "")
            if value in {"", "not_assessable"}:
                continue
            try:
                number = float(value)
            except ValueError:
                bad += 1
                continue
            if number < 0 or number > 1:
                bad += 1
    if bad:
        messages.append(f"bad_fraction_values:{table_kind}:{bad}:{path}")
    else:
        messages.append(f"ok_fraction_values:{table_kind}:rows={total}:{path}")
    return messages


PROTEIN_HIT_FEATURES = {"mrna", "match", "protein_match", "cdna_match"}
PROTEIN_PART_FEATURES = {"cds", "match_part", "hsp"}


def split_named_path(value: str) -> tuple[str, Path]:
    if "=" in value and not Path(value).exists():
        name, path = value.split("=", 1)
        return name.strip() or Path(path).stem, Path(path).expanduser().resolve()
    path = Path(value).expanduser().resolve()
    return path.stem, path


def parse_gff_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def check_protein_gff(label: str, path: Path, limit: int = 20) -> List[str]:
    messages: List[str] = []
    if not path.exists():
        return [f"missing:protein_gff:{label}:{path}"]
    hit_count = 0
    part_count = 0
    target_count = 0
    identity_count = 0
    bad = 0
    rows = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                bad += 1
                if bad >= limit:
                    break
                continue
            feature = fields[2].lower()
            if feature not in PROTEIN_HIT_FEATURES and feature not in PROTEIN_PART_FEATURES:
                continue
            try:
                start = int(float(fields[3]))
                end = int(float(fields[4]))
                if start <= 0 or end <= 0 or start > end:
                    bad += 1
            except ValueError:
                bad += 1
            attrs = parse_gff_attrs(fields[8])
            if feature in PROTEIN_HIT_FEATURES:
                hit_count += 1
                if attrs.get("Target"):
                    target_count += 1
                if attrs.get("Identity") or attrs.get("identity"):
                    identity_count += 1
            else:
                part_count += 1
            if bad >= limit:
                break
    if bad:
        messages.append(f"bad_protein_gff_coordinates_or_columns:{label}:{bad}_examples_or_more:{path}")
    elif hit_count == 0:
        messages.append(f"bad_protein_gff_no_hit_features:{label}:expected_one_of={','.join(sorted(PROTEIN_HIT_FEATURES))}:{path}")
    else:
        messages.append(f"ok_protein_gff:{label}:hits={hit_count}:parts={part_count}:rows={rows}:{path}")
        if target_count == 0:
            messages.append(f"warning_protein_gff_no_target:{label}:Target attribute is recommended:{path}")
        if identity_count == 0:
            messages.append(f"warning_protein_gff_no_identity:{label}:Identity attribute is recommended:{path}")
    return messages


def paths_from_args(args: argparse.Namespace) -> Dict[str, Path]:
    base = Path(args.evidence_dir).expanduser().resolve() if args.evidence_dir else None
    paths: Dict[str, Path] = {}
    for kind, filename in [
        ("sample_manifest", "sample_manifest.tsv"),
        ("merged_splice_junctions", "merged_splice_junctions.tsv"),
        ("junction_support_by_model", "junction_support_by_model.tsv"),
    ]:
        value = getattr(args, kind)
        if value:
            paths[kind] = Path(value).expanduser().resolve()
        elif base:
            paths[kind] = base / filename
    return paths


def validate_evidence(args: argparse.Namespace) -> int:
    paths = paths_from_args(args)
    protein_gffs = [split_named_path(value) for value in args.protein_gff]
    required_kinds: List[str] = list(paths)
    if args.evidence_dir:
        required_kinds = ["sample_manifest", "merged_splice_junctions", "junction_support_by_model"]
    elif not paths and not protein_gffs:
        required_kinds = ["sample_manifest", "merged_splice_junctions", "junction_support_by_model"]
    messages: List[str] = []
    if len(protein_gffs) > 1:
        labels = ",".join(label for label, _path in protein_gffs)
        messages.append(f"warning_multiple_protein_gff:{labels}:protein hits are merged without deduplication")
    for kind in required_kinds:
        path = paths.get(kind)
        if not path:
            messages.append(f"missing_argument:{kind}")
            continue
        messages.extend(check_columns(path, kind))
        if kind == "merged_splice_junctions":
            messages.extend(check_junction_coordinates(path))
        if kind == "junction_support_by_model":
            messages.extend(check_numeric_bounds(path, kind))
    for label, path in protein_gffs:
        messages.extend(check_protein_gff(label, path))
    failed = any(msg.startswith(("missing", "bad_")) for msg in messages)
    for msg in messages:
        print(msg)
    if failed and args.strict:
        return 1
    return 0


def add_validate_evidence_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("validate-evidence", help="Validate GeneArbiter evidence TSV schemas.")
    parser.add_argument("--evidence-dir", default="", help="Directory containing sample_manifest.tsv, merged_splice_junctions.tsv and junction_support_by_model.tsv.")
    parser.add_argument("--sample-manifest", dest="sample_manifest", default="")
    parser.add_argument("--merged-splice-junctions", dest="merged_splice_junctions", default="")
    parser.add_argument("--junction-support-by-model", dest="junction_support_by_model", default="")
    parser.add_argument("--protein-gff", action="append", default=[], help="Optional protein alignment GFF3/GTF evidence; repeatable; accepts path or name=path.")
    parser.add_argument("--strict", action="store_true", help="Return non-zero if required checks fail.")
    parser.set_defaults(func=validate_evidence)
