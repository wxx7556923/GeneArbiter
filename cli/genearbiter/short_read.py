#!/usr/bin/env python3
# script_id_md5: d71e1885c43e6fb8d94d3d81e597c8cc
# created: 2026-07-08
# modified: 2026-07-08
# owner: project
# status: project_code
# purpose: 从短读长 RNA-seq junction 表和 GFF 提取 GeneArbiter 可用的 short-read evidence TSV。
# inputs: sample_manifest.tsv; STAR SJ.out.tab 或已标准化 junction TSV; one or more GFF files。
# outputs: sample_manifest.tsv; merged_splice_junctions.tsv; junction_support_by_model.tsv。
# notes: 构建 splice-junction evidence，并汇总到转录本模型层级。

"""Short-read splice-junction evidence extraction utilities for GeneArbiter."""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


JUNCTION_COLUMNS = [
    "chrom",
    "intron_start_1based",
    "intron_end_1based",
    "bed_start_0based",
    "bed_end_0based",
    "strand",
    "donor_site",
    "acceptor_site",
    "splice_motif",
    "junction_read_count",
    "unique_junction_read_count",
    "max_overhang",
    "sample_support_count",
    "supporting_samples",
    "coordinate_note",
]

SUPPORT_COLUMNS = [
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
]


def parse_key_value(value: str) -> Tuple[str, str]:
    if "=" not in value:
        raise SystemExit(f"Expected sample=path or name=path, got: {value}")
    left, right = value.split("=", 1)
    if not left or not right:
        raise SystemExit(f"Expected non-empty sample=path, got: {value}")
    return left, right


def open_tsv(path: Path):
    return path.open("r", encoding="utf-8", newline="")


def parse_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in text.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def split_parents(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def load_gff_transcripts(gff_specs: Sequence[str]) -> Dict[str, Dict[str, object]]:
    transcripts: Dict[str, Dict[str, object]] = {}
    gene_by_tx: Dict[str, str] = {}
    for spec in gff_specs:
        name, raw_path = parse_key_value(spec) if "=" in spec else (Path(spec).stem, spec)
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise SystemExit(f"GFF does not exist: {path}")
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue
                seqid, _source, feature, start, end, _score, strand, _phase, attr_text = parts
                feature_l = feature.lower()
                attrs = parse_attrs(attr_text)
                if feature_l in {"mrna", "transcript", "rna"}:
                    tx_id = attrs.get("ID") or attrs.get("transcript_id")
                    if not tx_id:
                        continue
                    parent = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("gene") or ""
                    gene_id = split_parents(parent)[0] if parent else tx_id
                    key = f"{name}|{tx_id}"
                    gene_by_tx[tx_id] = gene_id
                    transcripts.setdefault(
                        key,
                        {"source": name, "transcript_id": tx_id, "gene_id": gene_id, "seqid": seqid, "strand": strand, "exons": []},
                    )
                elif feature_l == "exon":
                    parents = split_parents(attrs.get("Parent") or attrs.get("transcript_id") or "")
                    for tx_id in parents:
                        key = f"{name}|{tx_id}"
                        row = transcripts.setdefault(
                            key,
                            {
                                "source": name,
                                "transcript_id": tx_id,
                                "gene_id": gene_by_tx.get(tx_id, attrs.get("gene_id") or tx_id),
                                "seqid": seqid,
                                "strand": strand,
                                "exons": [],
                            },
                        )
                        row["exons"].append((int(start), int(end)))  # type: ignore[index]
    return transcripts


def star_strand(value: str) -> str:
    return {"0": ".", "1": "+", "2": "-"}.get(value, ".")


def star_motif(value: str) -> str:
    return {
        "0": "non_canonical",
        "1": "GT/AG",
        "2": "CT/AC",
        "3": "GC/AG",
        "4": "CT/GC",
        "5": "AT/AC",
        "6": "GT/AT",
    }.get(value, value)


def load_sample_ids(sample_manifest: Path) -> List[str]:
    if not sample_manifest.exists():
        raise SystemExit(f"sample manifest does not exist: {sample_manifest}")
    with open_tsv(sample_manifest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise SystemExit("sample_manifest.tsv must contain sample_id column")
        return [row["sample_id"] for row in reader if row.get("sample_id")]


def read_junction_file(sample_id: str, path: Path) -> Iterable[Dict[str, object]]:
    with open_tsv(path) as handle:
        first = handle.readline()
        if not first:
            return
        handle.seek(0)
        if first.startswith("chrom\t") or "intron_start_1based" in first:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                yield {
                    "chrom": row.get("chrom") or row.get("seqid") or row.get("chr") or "",
                    "start": int(float(row.get("intron_start_1based") or row.get("intron_start") or 0)),
                    "end": int(float(row.get("intron_end_1based") or row.get("intron_end") or 0)),
                    "strand": row.get("strand") or ".",
                    "motif": row.get("splice_motif") or "",
                    "read_count": int(float(row.get("junction_read_count") or row.get("support_count") or 0)),
                    "unique_count": int(float(row.get("unique_junction_read_count") or row.get("junction_read_count") or 0)),
                    "max_overhang": int(float(row.get("max_overhang") or 0)),
                    "sample_id": sample_id,
                    "note": row.get("coordinate_note") or "standard_junction_tsv",
                }
        else:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) < 9:
                    continue
                unique = int(float(row[6] or 0))
                multi = int(float(row[7] or 0))
                yield {
                    "chrom": row[0],
                    "start": int(float(row[1])),
                    "end": int(float(row[2])),
                    "strand": star_strand(row[3]),
                    "motif": star_motif(row[4]),
                    "read_count": unique + multi,
                    "unique_count": unique,
                    "max_overhang": int(float(row[8] or 0)),
                    "sample_id": sample_id,
                    "note": "STAR_SJ.out.tab_1based_intron_coordinates",
                }


def merge_junctions(junction_specs: Sequence[str]) -> Dict[Tuple[str, int, int, str], Dict[str, object]]:
    merged: Dict[Tuple[str, int, int, str], Dict[str, object]] = {}
    for spec in junction_specs:
        sample_id, raw_path = parse_key_value(spec)
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise SystemExit(f"junction file does not exist: {path}")
        for row in read_junction_file(sample_id, path):
            chrom = str(row["chrom"])
            start = int(row["start"])
            end = int(row["end"])
            strand = str(row["strand"] or ".")
            if not chrom or start <= 0 or end <= 0 or start >= end:
                continue
            key = (chrom, start, end, strand)
            item = merged.setdefault(
                key,
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "motifs": set(),
                    "read_count": 0,
                    "unique_count": 0,
                    "max_overhang": 0,
                    "samples": set(),
                    "notes": set(),
                },
            )
            item["motifs"].add(str(row.get("motif") or ""))  # type: ignore[index]
            item["read_count"] = int(item["read_count"]) + int(row["read_count"])
            item["unique_count"] = int(item["unique_count"]) + int(row["unique_count"])
            item["max_overhang"] = max(int(item["max_overhang"]), int(row["max_overhang"]))
            item["samples"].add(str(row["sample_id"]))  # type: ignore[index]
            item["notes"].add(str(row.get("note") or ""))  # type: ignore[index]
    return merged


def write_merged_junctions(path: Path, merged: Mapping[Tuple[str, int, int, str], Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=JUNCTION_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(merged):
            item = merged[key]
            samples = sorted(item["samples"])  # type: ignore[index]
            motifs = sorted(m for m in item["motifs"] if m)  # type: ignore[index]
            notes = sorted(n for n in item["notes"] if n)  # type: ignore[index]
            writer.writerow(
                {
                    "chrom": item["chrom"],
                    "intron_start_1based": item["start"],
                    "intron_end_1based": item["end"],
                    "bed_start_0based": int(item["start"]) - 1,
                    "bed_end_0based": item["end"],
                    "strand": item["strand"],
                    "donor_site": "",
                    "acceptor_site": "",
                    "splice_motif": ";".join(motifs),
                    "junction_read_count": item["read_count"],
                    "unique_junction_read_count": item["unique_count"],
                    "max_overhang": item["max_overhang"],
                    "sample_support_count": len(samples),
                    "supporting_samples": ",".join(samples),
                    "coordinate_note": ";".join(notes),
                }
            )


def supported_junctions(
    merged: Mapping[Tuple[str, int, int, str], Mapping[str, object]],
    min_read_count: int,
    min_sample_count: int,
) -> Dict[Tuple[str, int, int, str], Mapping[str, object]]:
    out = {}
    for key, item in merged.items():
        if int(item["read_count"]) >= min_read_count and len(item["samples"]) >= min_sample_count:  # type: ignore[arg-type]
            out[key] = item
    return out


def introns_from_exons(exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    ordered = sorted(exons)
    introns = []
    for left, right in zip(ordered, ordered[1:]):
        start = left[1] + 1
        end = right[0] - 1
        if start < end:
            introns.append((start, end))
    return introns


def junction_match(
    supported: Mapping[Tuple[str, int, int, str], Mapping[str, object]],
    seqid: str,
    start: int,
    end: int,
    strand: str,
) -> bool:
    strands = [strand] if strand and strand != "." else ["+", "-", "."]
    return any((seqid, start, end, s) in supported for s in strands) or (seqid, start, end, ".") in supported


def endpoint_supported(
    supported: Mapping[Tuple[str, int, int, str], Mapping[str, object]],
    seqid: str,
    pos: int,
    strand: str,
    endpoint: str,
) -> bool:
    for j_seqid, j_start, j_end, j_strand in supported:
        if j_seqid != seqid:
            continue
        if strand not in {"", "."} and j_strand not in {strand, "."}:
            continue
        if endpoint == "donor" and j_start == pos:
            return True
        if endpoint == "acceptor" and j_end == pos:
            return True
    return False


def nearby_supported(
    supported: Mapping[Tuple[str, int, int, str], Mapping[str, object]],
    seqid: str,
    start: int,
    end: int,
    strand: str,
    window: int,
) -> List[str]:
    hits = []
    for j_seqid, j_start, j_end, j_strand in supported:
        if j_seqid != seqid:
            continue
        if strand not in {"", "."} and j_strand not in {strand, "."}:
            continue
        if j_start == start and j_end == end:
            continue
        if abs(j_start - start) <= window or abs(j_end - end) <= window:
            hits.append(f"{j_seqid}:{j_start}-{j_end}:{j_strand}")
    return sorted(set(hits))[:20]


def write_model_junction_support(
    path: Path,
    transcripts: Mapping[str, Mapping[str, object]],
    merged: Mapping[Tuple[str, int, int, str], Mapping[str, object]],
    min_read_count: int,
    min_sample_count: int,
    nearby_window: int,
) -> None:
    supported = supported_junctions(merged, min_read_count, min_sample_count)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUPPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(transcripts):
            tx = transcripts[key]
            exons = tx.get("exons") or []
            introns = introns_from_exons(exons)  # type: ignore[arg-type]
            seqid = str(tx.get("seqid") or "")
            strand = str(tx.get("strand") or ".")
            unsupported = []
            nearby = []
            supported_count = 0
            donor_flags = []
            acceptor_flags = []
            for start, end in introns:
                exact = junction_match(supported, seqid, start, end, strand)
                if exact:
                    supported_count += 1
                else:
                    unsupported.append(f"{seqid}:{start}-{end}:{strand}")
                    nearby.extend(nearby_supported(supported, seqid, start, end, strand, nearby_window))
                donor_flags.append(endpoint_supported(supported, seqid, start, strand, "donor"))
                acceptor_flags.append(endpoint_supported(supported, seqid, end, strand, "acceptor"))
            n_introns = len(introns)
            if n_introns == 0:
                fraction = ""
                donor_supported = "not_assessable"
                acceptor_supported = "not_assessable"
                full_chain = "not_assessable"
            else:
                fraction = f"{supported_count / n_introns:.6f}"
                donor_supported = "true" if all(donor_flags) else "false"
                acceptor_supported = "true" if all(acceptor_flags) else "false"
                full_chain = "true" if supported_count == n_introns else "false"
            writer.writerow(
                {
                    "transcript_id": tx.get("transcript_id") or key,
                    "gene_id": tx.get("gene_id") or "",
                    "n_model_junctions": n_introns,
                    "n_supported_junctions": supported_count,
                    "junction_support_fraction": fraction,
                    "unsupported_junctions": ";".join(unsupported),
                    "novel_supported_junctions_nearby": ";".join(sorted(set(nearby))),
                    "donor_supported": donor_supported,
                    "acceptor_supported": acceptor_supported,
                    "full_intron_chain_supported": full_chain,
                }
            )


def extract_short_read(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_manifest = Path(args.sample_manifest).expanduser().resolve()
    sample_ids = set(load_sample_ids(sample_manifest))
    junction_specs = list(args.sj_out or []) + list(args.junction_tsv or [])
    if not junction_specs:
        raise SystemExit("Provide at least one --sj-out sample=SJ.out.tab or --junction-tsv sample=path.")
    for spec in junction_specs:
        sample_id, _path = parse_key_value(spec)
        if sample_ids and sample_id not in sample_ids:
            raise SystemExit(f"junction sample not found in sample_manifest.tsv: {sample_id}")
    if not args.gff:
        raise SystemExit("Provide at least one --gff name=path to compute junction_support_by_model.tsv.")

    shutil.copyfile(sample_manifest, out_dir / "sample_manifest.tsv")
    merged = merge_junctions(junction_specs)
    write_merged_junctions(out_dir / "merged_splice_junctions.tsv", merged)
    transcripts = load_gff_transcripts(args.gff)
    write_model_junction_support(
        out_dir / "junction_support_by_model.tsv",
        transcripts,
        merged,
        args.min_junction_read_count,
        args.min_sample_support_count,
        args.nearby_window,
    )
    print(f"out_dir\t{out_dir}")
    print(f"sample_manifest\t{out_dir / 'sample_manifest.tsv'}")
    print(f"merged_splice_junctions\t{out_dir / 'merged_splice_junctions.tsv'}")
    print(f"junction_support_by_model\t{out_dir / 'junction_support_by_model.tsv'}")
    return 0


def add_extract_short_read_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("extract-short-read", help="Build short-read junction evidence tables from STAR SJ.out.tab or junction TSV files.")
    parser.add_argument("--sample-manifest", required=True, help="TSV with at least sample_id column.")
    parser.add_argument("--gff", action="append", default=[], help="Repeatable name=path GFF used to compute junction_support_by_model.tsv.")
    parser.add_argument("--sj-out", action="append", default=[], help="Repeatable sample_id=STAR_SJ.out.tab.")
    parser.add_argument("--junction-tsv", action="append", default=[], help="Repeatable sample_id=existing_junction.tsv with GeneArbiter-style columns.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-junction-read-count", type=int, default=3)
    parser.add_argument("--min-sample-support-count", type=int, default=1)
    parser.add_argument("--nearby-window", type=int, default=500)
    parser.set_defaults(func=extract_short_read)
