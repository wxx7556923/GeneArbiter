#!/usr/bin/env python3
# script_id_md5: 9bd83a83c82b773611ded72ee332f10a
# created: 2026-06-22
# modified: 2026-06-25
# owner: project
# status: project_code
# purpose: 提供 GeneArbiter 候选模型复用的 GFF/interval 工具函数。
# inputs: 命令行参数指定的 TSV/JSONL/GFF3/FASTA 输入。
# outputs: 命令行参数指定的 TSV/JSONL/GFF3/Markdown 输出。
# notes: 只提供坐标和 GFF 基础操作，不做裁决。

"""Small GFF/interval helpers for GeneArbiter candidate workflows.

The helpers keep only traceable coordinates from input GFF-like files.  They do
not use benchmark truth labels and do not infer new coordinates by themselves.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
from bisect import bisect_left
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str = ""
    seqid: str = ""
    source: str = ""
    strand: str = "."
    start: int = 0
    end: int = 0
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int, str]] = field(default_factory=list)

    def span(self) -> Tuple[int, int]:
        starts = [self.start] if self.start else []
        ends = [self.end] if self.end else []
        starts.extend(start for start, _end in self.exons)
        ends.extend(end for _start, end in self.exons)
        starts.extend(start for start, _end, _phase in self.cds)
        ends.extend(end for _start, end, _phase in self.cds)
        return (min(starts), max(ends)) if starts and ends else (self.start, self.end)

    def sorted_exons(self) -> List[Tuple[int, int]]:
        blocks = self.exons or [(start, end) for start, end, _phase in self.cds]
        return sorted(set(blocks))

    def sorted_cds(self) -> List[Tuple[int, int, str]]:
        return sorted(set(self.cds), key=lambda item: (item[0], item[1], item[2]))

    def introns(self) -> List[Tuple[int, int]]:
        blocks = self.sorted_exons()
        return [(left[1] + 1, right[0] - 1) for left, right in zip(blocks, blocks[1:])]

    def intron_chain_key(self) -> str:
        introns = self.introns()
        if not introns:
            return "single_exon"
        return ",".join("{0}-{1}".format(start, end) for start, end in introns)

    def cds_key(self) -> str:
        cds = self.sorted_cds()
        if not cds:
            return "no_cds"
        return ",".join("{0}-{1}:{2}".format(start, end, phase) for start, end, phase in cds)


@dataclass
class Gene:
    gene_id: str
    seqid: str = ""
    source: str = ""
    strand: str = "."
    start: int = 0
    end: int = 0
    transcript_ids: List[str] = field(default_factory=list)


def ensure_dir(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for match in ATTR_RE.finditer(attr_text):
        attrs[match.group(1)] = match.group(2).strip('"')
    if not attrs:
        for item in attr_text.strip().split(";"):
            if not item:
                continue
            parts = item.strip().split(None, 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
    return attrs


def iter_gff_rows(path: str) -> Iterator[Tuple[str, str, str, int, int, str, str, str, Dict[str, str]]]:
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start_s, end_s, score, strand, phase, attr_text = fields[:9]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            yield seqid, source, feature, start, end, score, strand, phase, parse_attrs(attr_text)


def parse_annotation_gff(path: str, source_label: str = "") -> Tuple[Dict[str, Gene], Dict[str, Transcript], Dict[str, str]]:
    genes: Dict[str, Gene] = {}
    transcripts: Dict[str, Transcript] = {}
    tx_to_gene: Dict[str, str] = {}
    anonymous_gene_index = 0

    for seqid, source, feature, start, end, _score, strand, phase, attrs in iter_gff_rows(path):
        feature_l = feature.lower()
        row_source = source_label or source
        if feature_l == "gene":
            gene_id = attrs.get("ID") or attrs.get("Name")
            if not gene_id:
                anonymous_gene_index += 1
                gene_id = "anonymous_gene_{0}".format(anonymous_gene_index)
            genes[gene_id] = Gene(gene_id, seqid, row_source, strand, start, end, [])
        elif feature_l in {"mrna", "transcript", "rna"}:
            tx_id = attrs.get("ID") or attrs.get("Name")
            if not tx_id:
                continue
            parent = attrs.get("Parent", "").split(",")[0]
            if not parent:
                parent = tx_id
            tx = transcripts.setdefault(tx_id, Transcript(tx_id))
            tx.gene_id = parent
            tx.seqid = seqid
            tx.source = row_source
            tx.strand = strand
            tx.start = start
            tx.end = end
            tx_to_gene[tx_id] = parent
            genes.setdefault(parent, Gene(parent, seqid, row_source, strand, start, end, []))
            if tx_id not in genes[parent].transcript_ids:
                genes[parent].transcript_ids.append(tx_id)
        elif feature_l in {"exon", "cds"}:
            parents = [parent for parent in attrs.get("Parent", "").split(",") if parent]
            for tx_id in parents:
                tx = transcripts.setdefault(tx_id, Transcript(tx_id))
                if not tx.seqid:
                    tx.seqid = seqid
                    tx.source = row_source
                    tx.strand = strand
                tx.start = min(tx.start, start) if tx.start else start
                tx.end = max(tx.end, end)
                if feature_l == "exon":
                    tx.exons.append((start, end))
                else:
                    tx.cds.append((start, end, phase))

    for tx in transcripts.values():
        start, end = tx.span()
        tx.start, tx.end = start, end
        if tx.gene_id:
            gene = genes.setdefault(tx.gene_id, Gene(tx.gene_id, tx.seqid, tx.source, tx.strand, start, end, []))
            gene.seqid = gene.seqid or tx.seqid
            gene.source = gene.source or tx.source
            gene.strand = gene.strand if gene.strand != "." else tx.strand
            gene.start = min(gene.start, start) if gene.start else start
            gene.end = max(gene.end, end)
            if tx.transcript_id not in gene.transcript_ids:
                gene.transcript_ids.append(tx.transcript_id)
            tx_to_gene[tx.transcript_id] = tx.gene_id

    return genes, transcripts, tx_to_gene


def transcript_signature(tx: Transcript, mode: str = "intron_chain") -> str:
    if mode == "cds":
        return "{0}|{1}|{2}|{3}".format(tx.seqid, tx.strand, tx.cds_key(), tx.intron_chain_key())
    if mode == "exon":
        exons = tx.sorted_exons()
        key = ",".join("{0}-{1}".format(start, end) for start, end in exons) if exons else "no_exon"
        return "{0}|{1}|{2}".format(tx.seqid, tx.strand, key)
    introns = tx.introns()
    if not introns:
        exons = tx.sorted_exons()
        exon_key = ",".join("{0}-{1}".format(start, end) for start, end in exons) if exons else "{0}-{1}".format(tx.start, tx.end)
        return "{0}|{1}|single_exon:{2}".format(tx.seqid, tx.strand, exon_key)
    intron_key = ",".join("{0}-{1}".format(start, end) for start, end in introns)
    return "{0}|{1}|{2}".format(tx.seqid, tx.strand, intron_key)


def gene_signatures(gene: Gene, transcripts: Dict[str, Transcript], mode: str = "intron_chain") -> List[str]:
    signatures = []
    for tx_id in gene.transcript_ids:
        tx = transcripts.get(tx_id)
        if tx:
            signatures.append(transcript_signature(tx, mode))
    return sorted(set(signatures))


def build_interval_index(items: Iterable[Tuple[str, int, int, str]]) -> Tuple[Dict[str, List[Tuple[int, int, str]]], Dict[str, List[int]]]:
    by_seqid: Dict[str, List[Tuple[int, int, str]]] = {}
    for seqid, start, end, item_id in items:
        by_seqid.setdefault(seqid, []).append((start, end, item_id))
    for seqid in by_seqid:
        by_seqid[seqid].sort(key=lambda item: (item[0], item[1], item[2]))
    starts = {seqid: [item[0] for item in rows] for seqid, rows in by_seqid.items()}
    return by_seqid, starts


def query_overlaps(
    seqid: str,
    start: int,
    end: int,
    index: Dict[str, List[Tuple[int, int, str]]],
    starts_by_seqid: Dict[str, List[int]],
) -> List[Tuple[str, int, int, int]]:
    rows = index.get(seqid, [])
    starts = starts_by_seqid.get(seqid, [])
    if not rows:
        return []
    pos = max(0, bisect_left(starts, start) - 1)
    hits: List[Tuple[str, int, int, int]] = []
    while pos < len(rows):
        hit_start, hit_end, hit_id = rows[pos]
        if hit_start > end:
            break
        ov_start = max(start, hit_start)
        ov_end = min(end, hit_end)
        if ov_start <= ov_end:
            hits.append((hit_id, ov_start, ov_end, ov_end - ov_start + 1))
        pos += 1
    hits.sort(key=lambda item: (-item[3], item[0]))
    return hits


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_jsonl(path: str, rows: Iterable[Dict[str, object]]) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def pct(numerator: int, denominator: int) -> str:
    if not denominator:
        return "0.0000"
    return "{0:.4f}".format(float(numerator) / float(denominator))
