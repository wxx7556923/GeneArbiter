#!/usr/bin/env python3
# script_id_md5: 4d6a08cfbce1c7cf271155e1c6687b33
# created: 2026-07-10
# modified: 2026-07-10
# owner: project
# status: project_code
# purpose: 将不同工具输出的候选 GFF 清洗为统一 gene/mRNA/exon/CDS 层级。
# inputs: Repeatable name=path candidate GFF entries.
# outputs: normalized GFF3 files plus manifest/report TSV.
# notes: 只清洗候选来源 GFF；current/backbone 默认不改，避免影响 current ID 替换映射。

"""Normalize candidate GFF files before GeneArbiter card construction."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")


@dataclass
class Feature:
    seqid: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: Dict[str, str]


@dataclass
class TxRecord:
    tx_id: str
    gene_id: str
    seqid: str = ""
    source: str = ""
    start: int = 0
    end: int = 0
    strand: str = "."
    score: str = "."
    attrs: Dict[str, str] = field(default_factory=dict)
    exons: List[Feature] = field(default_factory=list)
    cds: List[Feature] = field(default_factory=list)


@dataclass
class GeneRecord:
    gene_id: str
    seqid: str = ""
    source: str = ""
    start: int = 0
    end: int = 0
    strand: str = "."
    score: str = "."
    attrs: Dict[str, str] = field(default_factory=dict)
    tx_ids: List[str] = field(default_factory=list)
    inferred: bool = False


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def parse_attrs(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for match in ATTR_RE.finditer(text):
        attrs[match.group(1).strip()] = match.group(2).strip().strip('"')
    if not attrs:
        for item in text.strip().split(";"):
            if not item:
                continue
            parts = item.strip().split(None, 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
    return attrs


def format_attrs(items: Sequence[Tuple[str, object]]) -> str:
    parts = []
    for key, value in items:
        value_s = str(value)
        value_s = value_s.replace(";", ",").replace("\t", " ").replace("\n", " ")
        parts.append(f"{key}={value_s}")
    return ";".join(parts)


def iter_features(path: Path) -> Iterable[Feature]:
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start_s, end_s, score, strand, phase, attr_text = fields[:9]
            try:
                start = int(float(start_s))
                end = int(float(end_s))
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            yield Feature(seqid, source, feature, start, end, score, strand, phase, parse_attrs(attr_text))


def safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text) or "source"


def output_name(index: int, source_name: str) -> str:
    return f"{index:02d}_{safe_name(source_name)}.normalized.gff3"


def split_named_path(value: str) -> Tuple[str, Path]:
    if "=" not in value:
        path = Path(value).expanduser().resolve()
        return path.stem, path
    name, path = value.split("=", 1)
    return name.strip(), Path(path).expanduser().resolve()


def infer_gene_id(tx_id: str) -> str:
    match = re.match(r"^(.+)\.t\d+(?:\..*)?$", tx_id)
    if match:
        return match.group(1)
    match = re.match(r"^(.+?)[._-]mRNA\d*$", tx_id, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return tx_id


def update_span(obj, feature: Feature) -> None:
    obj.seqid = obj.seqid or feature.seqid
    obj.source = obj.source or feature.source
    obj.strand = obj.strand if obj.strand != "." else feature.strand
    obj.score = obj.score if obj.score != "." else feature.score
    obj.start = min(obj.start, feature.start) if obj.start else feature.start
    obj.end = max(obj.end, feature.end)


def normalize_one(source_name: str, in_path: Path, out_path: Path) -> Dict[str, object]:
    genes: Dict[str, GeneRecord] = {}
    txs: Dict[str, TxRecord] = {}
    orphan_tx = 0
    inferred_genes = set()
    original_features = defaultdict(int)

    for feature in iter_features(in_path):
        feature_l = feature.feature.lower()
        original_features[feature.feature] += 1
        row_source = source_name or feature.source
        feature.source = row_source
        if feature_l == "gene":
            gene_id = feature.attrs.get("ID") or feature.attrs.get("Name")
            if not gene_id:
                gene_id = f"anonymous_gene_{len(genes) + 1}"
            gene = genes.setdefault(gene_id, GeneRecord(gene_id))
            gene.attrs.update(feature.attrs)
            update_span(gene, feature)
        elif feature_l in {"mrna", "transcript", "rna"}:
            tx_id = feature.attrs.get("ID") or feature.attrs.get("Name")
            if not tx_id:
                continue
            parent = feature.attrs.get("Parent", "").split(",")[0].strip()
            if not parent:
                parent = infer_gene_id(tx_id)
                orphan_tx += 1
                inferred_genes.add(parent)
            tx = txs.setdefault(tx_id, TxRecord(tx_id=tx_id, gene_id=parent))
            tx.gene_id = parent
            tx.attrs.update(feature.attrs)
            update_span(tx, feature)
            gene = genes.setdefault(parent, GeneRecord(parent, inferred=True))
            if parent in inferred_genes:
                gene.inferred = True
            update_span(gene, feature)
            if tx_id not in gene.tx_ids:
                gene.tx_ids.append(tx_id)
        elif feature_l in {"exon", "cds"}:
            parents = [parent for parent in feature.attrs.get("Parent", "").split(",") if parent]
            if not parents:
                continue
            for tx_id in parents:
                tx = txs.setdefault(tx_id, TxRecord(tx_id=tx_id, gene_id=infer_gene_id(tx_id)))
                update_span(tx, feature)
                if feature_l == "exon":
                    tx.exons.append(feature)
                else:
                    tx.cds.append(feature)
                gene = genes.setdefault(tx.gene_id, GeneRecord(tx.gene_id, inferred=True))
                update_span(gene, feature)
                if tx_id not in gene.tx_ids:
                    gene.tx_ids.append(tx_id)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        out.write(f"# normalized_by=genearbiter_gff_normalize;source={source_name};input={in_path}\n")
        for gene in sorted(genes.values(), key=lambda g: (g.seqid, g.start, g.end, g.gene_id)):
            attrs = [("ID", gene.gene_id), ("Name", gene.attrs.get("Name", gene.gene_id)), ("source_name", source_name)]
            if gene.inferred:
                attrs.append(("genearbiter_inferred_gene", "true"))
            out.write("\t".join([gene.seqid, gene.source or source_name, "gene", str(gene.start), str(gene.end), gene.score or ".", gene.strand or ".", ".", format_attrs(attrs)]) + "\n")
            for tx_id in sorted(gene.tx_ids, key=lambda t: (txs.get(t, TxRecord(t, '')).start, t)):
                tx = txs.get(tx_id)
                if not tx:
                    continue
                tx_attrs = [("ID", tx.tx_id), ("Parent", gene.gene_id), ("Name", tx.attrs.get("Name", tx.tx_id)), ("source_name", source_name)]
                out.write("\t".join([tx.seqid, tx.source or source_name, "mRNA", str(tx.start), str(tx.end), tx.score or ".", tx.strand or ".", ".", format_attrs(tx_attrs)]) + "\n")
                exons = tx.exons or [Feature(tx.seqid, tx.source or source_name, "exon", start, end, ".", tx.strand, ".", {}) for start, end, _phase in tx.cds]
                for i, exon in enumerate(sorted(exons, key=lambda f: (f.start, f.end)), start=1):
                    out.write("\t".join([exon.seqid, exon.source or source_name, "exon", str(exon.start), str(exon.end), exon.score or ".", exon.strand or tx.strand or ".", ".", format_attrs([("ID", f"{tx.tx_id}.exon.{i}"), ("Parent", tx.tx_id)])]) + "\n")
                for i, cds in enumerate(sorted(tx.cds, key=lambda f: (f.start, f.end)), start=1):
                    phase = cds.phase if cds.phase in {"0", "1", "2"} else "."
                    out.write("\t".join([cds.seqid, cds.source or source_name, "CDS", str(cds.start), str(cds.end), cds.score or ".", cds.strand or tx.strand or ".", phase, format_attrs([("ID", f"{tx.tx_id}.CDS.{i}"), ("Parent", tx.tx_id)])]) + "\n")

    return {
        "source": source_name,
        "input_gff": str(in_path),
        "output_gff": str(out_path),
        "gene_count": len(genes),
        "transcript_count": len(txs),
        "orphan_transcripts_fixed": orphan_tx,
        "inferred_gene_count": len([gene for gene in genes.values() if gene.inferred]),
        "original_gene_features": original_features.get("gene", 0),
        "original_mrna_features": sum(original_features.get(item, 0) for item in ["mRNA", "transcript", "RNA", "rna"]),
        "original_exon_features": original_features.get("exon", 0),
        "original_cds_features": original_features.get("CDS", 0),
    }


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fields = ["source", "input_gff", "output_gff", "gene_count", "transcript_count", "orphan_transcripts_fixed", "inferred_gene_count", "original_gene_features", "original_mrna_features", "original_exon_features", "original_cds_features"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", action="append", default=[], help="Repeatable source=path candidate GFF input.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--report", required=True)
    args = parser.parse_args(argv)
    if not args.gff:
        raise SystemExit("At least one --gff source=path is required.")
    out_dir = Path(args.out_dir).expanduser().resolve()
    rows = []
    for index, value in enumerate(args.gff, start=1):
        source_name, in_path = split_named_path(value)
        out_path = out_dir / output_name(index, source_name)
        rows.append(normalize_one(source_name, in_path, out_path))
    write_tsv(Path(args.manifest).expanduser().resolve(), rows)
    write_tsv(Path(args.report).expanduser().resolve(), rows)
    for row in rows:
        print(f"normalized_gff\t{row['source']}\t{row['output_gff']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
