#!/usr/bin/env python3
# script_id_md5: 816dc384c2e14e2717b808f011e462f3
# created: 2026-06-24
# modified: 2026-07-09
# owner: project
# status: project_code
# purpose: 构建当前注释/从头预测模型仲裁的位点级完整卡。
# inputs: 命令行参数指定的 TSV/JSONL/GFF3/FASTA 输入。
# outputs: 命令行参数指定的 TSV/JSONL/GFF3/Markdown 输出。
# notes: 新增 --source-priority；protein evidence 支持纯路径或 name=path，避免 final config 传参失败。

"""为当前注释/从头预测模型仲裁构建位点级卡。

本脚本把当前注释和从头预测结果作为可比较的候选模型。它不选择最终模型，
也不编辑 GFF 文件；只为受约束 AI 判断层准备可追溯的可选模型 ID、
RNA/protein 证据摘要和复杂度标记。
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Set, Tuple

from gff_utils import (
    Transcript,
    build_interval_index,
    iter_gff_rows,
    parse_annotation_gff,
    pct,
    query_overlaps,
    read_tsv,
    transcript_signature,
    write_jsonl,
    write_tsv,
)


MODEL_FIELDS = [
    "locus_id",
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
    "intron_count",
    "signature",
    "rna_junction_status",
    "exact_junction_supported_introns",
    "not_observed_introns",
    "junction_support_fraction",
    "conflicting_junction_count",
    "nearby_alternative_junction_count",
    "protein_overlap_count",
    "protein_best_overlap_bp",
    "protein_best_target",
    "protein_best_identity",
    "long_read_status",
    "long_read_support_class",
    "long_read_exact_chain_support",
    "long_read_partial_junction_support",
    "long_read_conflicting_chain_count",
    "long_read_risk_tags",
    "model_level_junction_status",
    "model_level_junction_support_fraction",
    "model_level_full_intron_chain_supported",
    "validation_status",
    "validation_fail_reasons",
    "hard_flags",
    "risk_tags",
    "recommended_model_role",
]

COMPLEX_FIELDS = [
    "locus_id",
    "seqid",
    "start",
    "end",
    "candidate_model_count",
    "source_count",
    "current_model_count",
    "tool_model_count",
    "unique_structure_count",
    "max_gene_count_per_source",
    "complexity_class",
    "complexity_reasons",
    "risk_flagged",
    "risk_triggers",
    "decision_required",
    "recommended_action_space",
    "conflict_class",
    "recommended_action",
]

SUMMARY_FIELDS = [
    "locus_id",
    "seqid",
    "start",
    "end",
    "candidate_model_count",
    "source_count",
    "complexity_class",
    "risk_flagged",
    "risk_triggers",
    "decision_required",
    "decision_requirement",
    "recommended_action_space",
    "conflict_class",
    "selectable_model_ids",
    "recommended_next_gate",
    "risk_tags",
]

UNIT_FIELDS = [
    "locus_id",
    "seqid",
    "start",
    "end",
    "strand",
    "candidate_model_count",
    "source_count",
    "current_gene_count",
    "model_ids",
    "unit_edge_count",
    "edge_types",
    "construction_mode",
]

EDGE_FIELDS = [
    "locus_id",
    "model_id_a",
    "model_id_b",
    "edge_type",
    "edge_strength",
    "seqid",
    "strand",
    "evidence_detail",
]

SET_FIELDS = [
    "locus_id",
    "set_id",
    "source",
    "set_role",
    "gene_count",
    "transcript_count",
    "model_ids",
    "gene_ids",
    "seqid",
    "start",
    "end",
    "strand",
    "validation_pass_count",
    "validation_warning_count",
    "validation_fail_count",
    "validation_not_assessed_count",
    "rna_full_junction_support_count",
    "rna_partial_junction_support_count",
    "rna_conflicting_junction_count",
    "rna_none_observed_count",
    "protein_supported_count",
    "long_read_exact_supported_transcripts",
    "long_read_partial_supported_transcripts",
    "long_read_unsupported_transcripts",
    "long_read_mono_exon_not_assessable_count",
    "long_read_conflicting_chain_count",
    "long_read_support_level_counts",
    "long_read_best_supported_model_ids",
    "risk_tags",
    "hard_flags",
    "candidate_set_mode",
    "supporting_sources",
    "representative_source",
    "equivalent_source_set_ids",
    "structure_signature_key",
]


@dataclass
class CandidateModel:
    model_id: str
    source: str
    gene_id: str
    transcript_id: str
    transcript: Transcript
    source_file: str
    signature: str


@dataclass
class JunctionEvidence:
    junction_id: str
    seqid: str
    intron_start: int
    intron_end: int
    strand: str
    support_count: int
    sample_count: int
    sample_ids: str
    source_file: str


@dataclass
class ProteinHit:
    hit_id: str
    seqid: str
    start: int
    end: int
    strand: str
    target: str = ""
    identity: str = ""
    source_tool: str = ""



_JUNCTION_INDEX_CACHE = {}
DEFAULT_SOURCE_PRIORITY = ["current"]
SOURCE_PRIORITY = list(DEFAULT_SOURCE_PRIORITY)


def parse_source_priority(values: Sequence[str]) -> List[str]:
    order: List[str] = []
    for value in values or []:
        for item in str(value).split(","):
            source = item.strip()
            if source and source not in order:
                order.append(source)
    if "current" not in order:
        order.insert(0, "current")
    return order or list(DEFAULT_SOURCE_PRIORITY)


def source_priority_rank(source: object) -> int:
    source_name = str(source)
    try:
        return SOURCE_PRIORITY.index(source_name)
    except ValueError:
        return len(SOURCE_PRIORITY)


def parse_tool_arg(text: str) -> Tuple[str, str]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("--tool-gff must be name=path")
    name, path = text.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("--tool-gff must be name=path")
    return name, path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-mode", choices=["correction", "de_novo_annotation"], default="correction")
    parser.add_argument("--current-gff", default="", help="Required in correction mode; not allowed in de_novo_annotation mode.")
    parser.add_argument("--tool-gff", action="append", type=parse_tool_arg, default=[], help="Repeatable: name=path")
    parser.add_argument("--source-priority", action="append", default=[], help="Optional comma-separated source priority order used to break evidence ties.")
    parser.add_argument("--junctions", action="append", default=[], help="RNA junction TSV; repeatable")
    parser.add_argument("--junction-support-summary", default="", help="Optional model-level RNA junction support summary TSV keyed by transcript_id")
    parser.add_argument("--long-read-support-dir", default="", help="Optional directory with per-source long_read_model_support_by_model.tsv and long_read_set_support_by_locus.tsv")
    parser.add_argument("--protein-gff", action="append", default=[], help="Optional miniprot/spaln-like protein alignment GFF; repeatable")
    parser.add_argument("--signature-mode", choices=["intron_chain", "exon", "cds"], default="intron_chain")
    parser.add_argument("--min-overlap-fraction", type=float, default=0.20)
    parser.add_argument("--min-junction-support-count", type=int, default=3)
    parser.add_argument("--min-junction-sample-count", type=int, default=1)
    parser.add_argument("--boundary-window", type=int, default=30)
    parser.add_argument("--nearby-window", type=int, default=500)
    parser.add_argument("--protein-flank", type=int, default=0)
    parser.add_argument("--validation-summary", default="", help="Optional per-model fatal QC TSV keyed by model_id or object_id")
    parser.add_argument("--eligibility-summary", default="", help="Optional pre-AI eligibility TSV keyed by model_id; ineligible models are excluded from AI-selectable sets.")
    parser.add_argument("--candidate-set-mode", choices=["structure_consensus", "source"], default="structure_consensus", help="structure_consensus collapses source-level sets with identical model structure signatures; source keeps one set per input source.")
    parser.add_argument("--max-candidate-models", type=int, default=12)
    parser.add_argument("--max-unique-structures", type=int, default=8)
    parser.add_argument("--max-transcripts-per-source", type=int, default=6)
    parser.add_argument("--max-locus-bp", type=int, default=500000)
    parser.add_argument("--ai-soft-locus-bp", type=int, default=50000, help="Window span above this is checked with gene/model/weak-edge complexity gates; 0 disables.")
    parser.add_argument("--ai-hard-locus-bp", type=int, default=100000, help="Window span above this is always human-review-only before AI decision; 0 disables.")
    parser.add_argument("--max-ai-current-genes-in-long-window", type=int, default=3, help="For windows above --ai-soft-locus-bp, current gene count above this is human-review-only; 0 disables.")
    parser.add_argument("--max-ai-models-per-source-in-long-window", type=int, default=3, help="For windows above --ai-soft-locus-bp, candidate model count above source_count * this value is human-review-only; 0 disables.")
    parser.add_argument("--gene-list", default="", help="Optional newline/TSV gene or transcript ID list for pilot subset selection")
    parser.add_argument("--region-table", default="", help="Optional TSV with seqid/start/end columns for pilot subset selection")
    parser.add_argument("--region-flank", type=int, default=0, help="Flank bp added when selecting by --region-table")
    parser.add_argument("--output-units", default="", help="Optional evidence-aware 仲裁单元 summary TSV")
    parser.add_argument("--output-unit-edges", default="", help="Optional evidence-aware unit graph edge TSV")
    parser.add_argument("--output-model-sets", default="", help="Optional gene-set arbitration candidates TSV")
    parser.add_argument("--output-models", required=True)
    parser.add_argument("--output-complex-loci", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-jsonl", required=True)
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
        rows = read_tsv(path)
        preferred = ["gene_id", "current_gene_id", "transcript_id", "locus_id"]
        for row in rows:
            for field in preferred:
                if field in row:
                    ids.update(split_ids(row.get(field, "")))
            if not any(field in row for field in preferred):
                for value in row.values():
                    ids.update(split_ids(value))
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
    for row in read_tsv(path):
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


def model_overlaps_regions(model: CandidateModel, regions: Sequence[Tuple[str, int, int]]) -> bool:
    tx = model.transcript
    for seqid, start, end in regions:
        if tx.seqid == seqid and tx.end >= start and tx.start <= end:
            return True
    return False


def filter_models_for_subset(
    models: Sequence[CandidateModel],
    gene_ids: Set[str],
    regions: Sequence[Tuple[str, int, int]],
) -> List[CandidateModel]:
    if not gene_ids and not regions:
        return list(models)
    selected = []
    for model in models:
        by_id = model.gene_id in gene_ids or model.transcript_id in gene_ids or model.model_id in gene_ids
        by_region = model_overlaps_regions(model, regions) if regions else False
        if by_id or by_region:
            selected.append(model)
    return selected


def load_candidate_models(current_gff: str, tool_gffs: Sequence[Tuple[str, str]], signature_mode: str, task_mode: str) -> List[CandidateModel]:
    all_sources: List[Tuple[str, str]] = []
    if task_mode == "correction" and current_gff:
        all_sources.append(("current", current_gff))
    all_sources.extend(list(tool_gffs))
    models: List[CandidateModel] = []
    for source, path in all_sources:
        _genes, transcripts, tx_to_gene = parse_annotation_gff(path, source)
        for tx_id, tx in sorted(transcripts.items()):
            gene_id = tx_to_gene.get(tx_id, tx.gene_id) or tx_id
            if not tx.seqid or not tx.start or not tx.end:
                continue
            model_id = "{0}:{1}:{2}".format(safe_id(source), safe_id(gene_id), safe_id(tx_id))
            models.append(
                CandidateModel(
                    model_id=model_id,
                    source=source,
                    gene_id=gene_id,
                    transcript_id=tx_id,
                    transcript=tx,
                    source_file=path,
                    signature=transcript_signature(tx, signature_mode),
                )
            )
    return models


def interval_overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    start = max(a_start, b_start)
    end = min(a_end, b_end)
    return max(0, end - start + 1)


def interval_list_overlap_bp(a_items: Sequence[Tuple[int, int]], b_items: Sequence[Tuple[int, int]]) -> int:
    total = 0
    for a_start, a_end in a_items:
        for b_start, b_end in b_items:
            total += interval_overlap_bp(a_start, a_end, b_start, b_end)
    return total


def compatible_strand(a: str, b: str) -> bool:
    return a == b or "." in {a, b} or not a or not b


def splice_sites(tx: Transcript) -> Set[int]:
    sites: Set[int] = set()
    for intron_start, intron_end in tx.introns():
        sites.add(intron_start)
        sites.add(intron_end)
    return sites


def edge_key(a: str, b: str, edge_type: str, detail: str) -> Tuple[str, str, str, str]:
    left, right = sorted([a, b])
    return left, right, edge_type, detail


def add_edge(
    graph: Dict[str, Set[str]],
    edges: Dict[Tuple[str, str, str, str], Dict[str, object]],
    a: CandidateModel,
    b: CandidateModel,
    edge_type: str,
    edge_strength: str,
    detail: str,
) -> None:
    if a.model_id == b.model_id:
        return
    left, right, e_type, e_detail = edge_key(a.model_id, b.model_id, edge_type, detail)
    if edge_type in {"protein_bridge", "RNA_junction_bridge"}:
        e_detail = edge_type
    graph[left].add(right)
    graph[right].add(left)
    edges.setdefault((left, right, e_type, e_detail), {
        "model_id_a": left,
        "model_id_b": right,
        "edge_type": edge_type,
        "edge_strength": edge_strength,
        "seqid": a.transcript.seqid,
        "strand": a.transcript.strand if a.transcript.strand == b.transcript.strand else "compatible_or_mixed",
        "evidence_detail": detail,
    })


def structural_edges_for_pair(
    a: CandidateModel,
    b: CandidateModel,
    min_overlap_fraction: float,
) -> List[Tuple[str, str, str]]:
    if a.transcript.seqid != b.transcript.seqid:
        return []
    if not compatible_strand(a.transcript.strand, b.transcript.strand):
        return []

    tx_a = a.transcript
    tx_b = b.transcript
    out: List[Tuple[str, str, str]] = []

    span_bp = interval_overlap_bp(tx_a.start, tx_a.end, tx_b.start, tx_b.end)
    if span_bp:
        len_a = max(1, tx_a.end - tx_a.start + 1)
        len_b = max(1, tx_b.end - tx_b.start + 1)
        frac_a = span_bp / float(len_a)
        frac_b = span_bp / float(len_b)
        if frac_a >= min_overlap_fraction or frac_b >= min_overlap_fraction:
            out.append(("span_overlap", "moderate", "overlap_bp={0};frac_a={1:.4f};frac_b={2:.4f}".format(span_bp, frac_a, frac_b)))

    exon_a = tx_a.sorted_exons()
    exon_b = tx_b.sorted_exons()
    exon_bp = interval_list_overlap_bp(exon_a, exon_b)
    if exon_bp:
        out.append(("exon_overlap", "strong", "overlap_bp={0}".format(exon_bp)))

    cds_a = [(start, end) for start, end, _phase in tx_a.sorted_cds()]
    cds_b = [(start, end) for start, end, _phase in tx_b.sorted_cds()]
    cds_bp = interval_list_overlap_bp(cds_a, cds_b)
    if cds_bp:
        out.append(("CDS_overlap", "strong", "overlap_bp={0}".format(cds_bp)))

    introns_a = tx_a.introns()
    introns_b = tx_b.introns()
    if introns_a and introns_a == introns_b:
        out.append(("identical_intron_chain", "strong", "intron_count={0}".format(len(introns_a))))
    elif introns_a and introns_b:
        shared_introns = sorted(set(introns_a) & set(introns_b))
        if shared_introns:
            out.append(("shared_introns", "strong", "shared_intron_count={0}".format(len(shared_introns))))
        shared_sites = sorted(splice_sites(tx_a) & splice_sites(tx_b))
        if shared_sites:
            out.append(("shared_splice_site", "moderate", "shared_site_count={0}".format(len(shared_sites))))

    if introns_a and introns_a == introns_b and (tx_a.start != tx_b.start or tx_a.end != tx_b.end):
        out.append(("terminal_boundary_difference_only", "moderate", "same_intron_chain_with_different_span"))

    return out


def build_same_source_gene_edges(
    models: Sequence[CandidateModel],
    graph: Dict[str, Set[str]],
    edges: Dict[Tuple[str, str, str, str], Dict[str, object]],
) -> None:
    """保证同一来源同一基因的所有模型完整保留在同一个 仲裁单元 中。"""
    by_source_gene: Dict[Tuple[str, str], List[CandidateModel]] = defaultdict(list)
    for model in models:
        if model.source and model.gene_id:
            by_source_gene[(model.source, model.gene_id)].append(model)
    for (source, gene_id), gene_models in by_source_gene.items():
        if len(gene_models) < 2:
            continue
        gene_models = sorted(gene_models, key=lambda model: model.model_id)
        detail = "source={0};gene_id={1}".format(source, gene_id)
        for index, left in enumerate(gene_models):
            for right in gene_models[index + 1:]:
                add_edge(graph, edges, left, right, "same_source_gene", "strong", detail)


def build_protein_bridge_edges(
    models: Sequence[CandidateModel],
    protein_hits: Sequence[ProteinHit],
    protein_flank: int,
    graph: Dict[str, Set[str]],
    edges: Dict[Tuple[str, str, str, str], Dict[str, object]],
) -> None:
    if not protein_hits:
        return
    protein_items = []
    protein_by_id: Dict[str, ProteinHit] = {}
    for hit in protein_hits:
        start = max(1, hit.start - protein_flank)
        end = hit.end + protein_flank
        protein_items.append((hit.seqid, start, end, hit.hit_id))
        protein_by_id[hit.hit_id] = hit
    protein_idx, protein_starts = build_interval_index(protein_items)
    models_by_protein: Dict[str, List[CandidateModel]] = defaultdict(list)
    for model in models:
        tx = model.transcript
        hits = query_overlaps(tx.seqid, tx.start, tx.end, protein_idx, protein_starts)
        for hit_id, _ov_start, _ov_end, _ov_bp in hits:
            models_by_protein[hit_id].append(model)
    for hit_id, hit_models in models_by_protein.items():
        if len(hit_models) < 2:
            continue
        hit = protein_by_id[hit_id]
        hit_models = sorted(hit_models, key=lambda model: model.model_id)
        for index, left in enumerate(hit_models):
            for right in hit_models[index + 1:]:
                if not compatible_strand(left.transcript.strand, right.transcript.strand):
                    continue
                detail = "protein_hit={0};target={1}".format(hit.hit_id, hit.target)
                add_edge(graph, edges, left, right, "protein_bridge", "moderate", detail)

def junction_bridges_models(junction: JunctionEvidence, left: CandidateModel, right: CandidateModel, boundary_window: int, nearby_window: int) -> bool:
    if left.transcript.seqid != right.transcript.seqid or junction.seqid != left.transcript.seqid:
        return False
    if not compatible_strand(left.transcript.strand, right.transcript.strand):
        return False
    if junction.strand not in {left.transcript.strand, right.transcript.strand, "."} and left.transcript.strand != "." and right.transcript.strand != ".":
        return False
    a, b = sorted([left, right], key=lambda model: (model.transcript.start, model.transcript.end))
    if a.transcript.end >= b.transcript.start:
        return False
    gap = b.transcript.start - a.transcript.end - 1
    if gap > nearby_window:
        return False
    left_near = abs(junction.intron_start - a.transcript.end) <= boundary_window or a.transcript.start <= junction.intron_start <= a.transcript.end
    right_near = abs(junction.intron_end - b.transcript.start) <= boundary_window or b.transcript.start <= junction.intron_end <= b.transcript.end
    return left_near and right_near and junction.intron_start < junction.intron_end


def build_junction_bridge_edges(
    models: Sequence[CandidateModel],
    junctions: Dict[Tuple[str, int, int, str], JunctionEvidence],
    boundary_window: int,
    nearby_window: int,
    graph: Dict[str, Set[str]],
    edges: Dict[Tuple[str, str, str, str], Dict[str, object]],
) -> None:
    if not junctions:
        return
    by_seqid: Dict[str, List[CandidateModel]] = defaultdict(list)
    for model in models:
        by_seqid[model.transcript.seqid].append(model)
    junction_items = [(j.seqid, j.intron_start, j.intron_end, j.junction_id) for j in junctions.values()]
    junction_index, junction_starts = build_interval_index(junction_items)
    junction_by_id = {j.junction_id: j for j in junctions.values()}

    for seq_models in by_seqid.values():
        seq_models.sort(key=lambda model: (model.transcript.start, model.transcript.end, model.model_id))
        for index, left in enumerate(seq_models):
            for right in seq_models[index + 1:]:
                if right.transcript.start - left.transcript.end - 1 > nearby_window:
                    break
                if left.transcript.end >= right.transcript.start:
                    continue
                query_start = max(1, left.transcript.start - boundary_window)
                query_end = right.transcript.end + boundary_window
                hits = query_overlaps(left.transcript.seqid, query_start, query_end, junction_index, junction_starts)
                for junction_id, _hit_start, _hit_end, _ov_bp in hits:
                    junction = junction_by_id[junction_id]
                    if junction_bridges_models(junction, left, right, boundary_window, nearby_window):
                        detail = "{0};support={1};samples={2}".format(junction.junction_id, junction.support_count, junction.sample_count)
                        add_edge(graph, edges, left, right, "RNA_junction_bridge", "strong", detail)
                        break



def connected_components(
    models: Sequence[CandidateModel],
    graph: Dict[str, Set[str]],
) -> List[List[CandidateModel]]:
    by_id = {m.model_id: m for m in models}
    loci: List[List[CandidateModel]] = []
    seen: Set[str] = set()
    for model in sorted(models, key=lambda item: (item.transcript.seqid, item.transcript.start, item.transcript.end, item.model_id)):
        if model.model_id in seen:
            continue
        component_ids = []
        queue = deque([model.model_id])
        seen.add(model.model_id)
        while queue:
            model_id = queue.popleft()
            component_ids.append(model_id)
            for next_id in sorted(graph[model_id]):
                if next_id not in seen:
                    seen.add(next_id)
                    queue.append(next_id)
        loci.append([by_id[model_id] for model_id in sorted(component_ids)])
    return loci


def build_loci(
    models: Sequence[CandidateModel],
    min_overlap_fraction: float,
    junctions: Dict[Tuple[str, int, int, str], JunctionEvidence],
    boundary_window: int,
    nearby_window: int,
    protein_hits: Sequence[ProteinHit],
    protein_flank: int,
    regions: Sequence[Tuple[str, int, int]],
) -> Tuple[List[List[CandidateModel]], Dict[Tuple[str, str], List[Dict[str, object]]]]:
    """基于可追溯、非 truth 证据构建 evidence-aware 仲裁单元s。"""
    items = [(m.transcript.seqid, m.transcript.start, m.transcript.end, m.model_id) for m in models]
    index, starts = build_interval_index(items)
    by_id = {m.model_id: m for m in models}
    graph: Dict[str, Set[str]] = {m.model_id: set() for m in models}
    edge_records: Dict[Tuple[str, str, str, str], Dict[str, object]] = {}

    for model in models:
        tx = model.transcript
        for hit_id, _ov_start, _ov_end, _ov_bp in query_overlaps(tx.seqid, tx.start, tx.end, index, starts):
            if hit_id == model.model_id:
                continue
            other = by_id[hit_id]
            for edge_type, strength, detail in structural_edges_for_pair(model, other, min_overlap_fraction):
                add_edge(graph, edge_records, model, other, edge_type, strength, detail)

    # region 只作为子集筛选上下文，不能用来合并生物学意义上的 仲裁单元s。
    build_same_source_gene_edges(models, graph, edge_records)
    build_protein_bridge_edges(models, protein_hits, protein_flank, graph, edge_records)
    build_junction_bridge_edges(models, junctions, boundary_window, nearby_window, graph, edge_records)

    loci = connected_components(models, graph)
    component_by_model: Dict[str, str] = {}
    for unit_index, locus_models in enumerate(loci, 1):
        unit_id = "locus_{0:08d}".format(unit_index)
        for model in locus_models:
            component_by_model[model.model_id] = unit_id

    edges_by_pair: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for record in edge_records.values():
        unit_a = component_by_model.get(str(record["model_id_a"]), "")
        unit_b = component_by_model.get(str(record["model_id_b"]), "")
        if unit_a and unit_a == unit_b:
            record = dict(record)
            record["locus_id"] = unit_a
            edges_by_pair[(str(record["model_id_a"]), str(record["model_id_b"]))].append(record)
    return loci, edges_by_pair


def unit_rows_for_loci(loci: Sequence[Sequence[CandidateModel]], edges_by_pair: Dict[Tuple[str, str], List[Dict[str, object]]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    edge_rows = [edge for edges in edges_by_pair.values() for edge in edges]
    edges_by_locus: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for edge in edge_rows:
        edges_by_locus[str(edge["locus_id"])].append(edge)
    for locus_index, locus_models in enumerate(loci, 1):
        locus_id = "locus_{0:08d}".format(locus_index)
        seqid = locus_models[0].transcript.seqid
        start = min(model.transcript.start for model in locus_models)
        end = max(model.transcript.end for model in locus_models)
        strands = sorted(set(model.transcript.strand for model in locus_models))
        sources = sorted(set(model.source for model in locus_models))
        current_genes = sorted(set(model.gene_id for model in locus_models if model.source == "current"))
        unit_edges = edges_by_locus.get(locus_id, [])
        rows.append(
            {
                "locus_id": locus_id,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strands[0] if len(strands) == 1 else ",".join(strands),
                "candidate_model_count": len(locus_models),
                "source_count": len(sources),
                "current_gene_count": len(current_genes),
                "model_ids": ",".join(sorted(model.model_id for model in locus_models)),
                "unit_edge_count": len(unit_edges),
                "edge_types": ",".join(sorted(set(str(edge["edge_type"]) for edge in unit_edges))),
                "construction_mode": "evidence_aware_arbitration_unit",
            }
        )
    return rows


def load_junctions(paths: Sequence[str], min_support_count: int, min_sample_count: int) -> Dict[Tuple[str, int, int, str], JunctionEvidence]:
    merged: Dict[Tuple[str, int, int, str], Dict[str, object]] = {}
    for path in paths:
        with open(path, "r", newline="") as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            for row in rows:
                try:
                    start = int(row.get("intron_start_1based") or row.get("intron_start") or 0)
                    end = int(row.get("intron_end_1based") or row.get("intron_end") or 0)
                    support = int(float(row.get("support_count") or row.get("junction_read_count") or row.get("unique_junction_read_count") or 0))
                except ValueError:
                    continue
                if support < min_support_count:
                    continue
                seqid = row.get("seqid") or row.get("chrom") or row.get("chromosome") or ""
                strand = row.get("strand", ".")
                key = (seqid, start, end, strand)
                entry = merged.setdefault(key, {"support": 0, "samples": set(), "source_files": set()})
                entry["support"] += support
                sample_ids = row.get("sample_ids") or row.get("supporting_samples") or ""
                if sample_ids:
                    for sample_id in sample_ids.split(","):
                        if sample_id:
                            entry["samples"].add(sample_id)
                else:
                    sample_count = int(float(row.get("sample_count") or row.get("sample_support_count") or 1))
                    for index in range(max(1, sample_count)):
                        entry["samples"].add("{0}:sample_{1}".format(os.path.basename(path), index + 1))
                entry["source_files"].add(path)

    junctions = {}
    for (seqid, start, end, strand), entry in merged.items():
        sample_count = len(entry["samples"])
        if sample_count < min_sample_count:
            continue
        junction_id = "rna_junction_{0}:{1}-{2}:{3}".format(seqid, start, end, strand)
        junctions[(seqid, start, end, strand)] = JunctionEvidence(
            junction_id=junction_id,
            seqid=seqid,
            intron_start=start,
            intron_end=end,
            strand=strand,
            support_count=int(entry["support"]),
            sample_count=sample_count,
            sample_ids=",".join(sorted(entry["samples"])),
            source_file=",".join(sorted(entry["source_files"])),
        )
    return junctions


def nearby_junctions(
    junctions: Dict[Tuple[str, int, int, str], JunctionEvidence],
    seqid: str,
    strand: str,
    start: int,
    end: int,
    window: int,
) -> List[JunctionEvidence]:
    cache_key = id(junctions)
    cached = _JUNCTION_INDEX_CACHE.get(cache_key)
    if cached is None:
        items = [(j.seqid, j.intron_start, j.intron_end, j.junction_id) for j in junctions.values()]
        index, starts = build_interval_index(items)
        by_id = {j.junction_id: j for j in junctions.values()}
        cached = (index, starts, by_id)
        _JUNCTION_INDEX_CACHE[cache_key] = cached
    index, starts_by_seqid, by_id = cached
    hits = []
    for junction_id, _ov_start, _ov_end, _ov_bp in query_overlaps(seqid, max(1, start - window), end + window, index, starts_by_seqid):
        junction = by_id[junction_id]
        if junction.strand not in {strand, "."} and strand != ".":
            continue
        hits.append(junction)
    hits.sort(key=lambda item: (-item.sample_count, -item.support_count, item.intron_start, item.intron_end))
    return hits


def split_named_path(value: str) -> Tuple[str, str]:
    if "=" in value and not os.path.exists(value):
        name, path = value.split("=", 1)
        return safe_id(name.strip()), path.strip()
    return "", value


def parse_protein_gff(path: str, source_label: str = "") -> List[ProteinHit]:
    if not path:
        return []
    if not source_label:
        source_label = os.path.basename(path).replace(".gff3", "").replace(".gff", "")
    hits: Dict[str, ProteinHit] = {}
    current_hit_id = ""
    for line_number, (seqid, source, feature, start, end, _score, strand, _phase, attrs) in enumerate(iter_gff_rows(path), 1):
        feature_l = feature.lower()
        if feature_l in {"mrna", "match", "protein_match", "cdna_match"}:
            raw_hit_id = attrs.get("ID") or attrs.get("Name") or "protein_hit_line_{0}".format(line_number)
            hit_id = "{0}|{1}".format(source_label, raw_hit_id)
            hits[hit_id] = ProteinHit(
                hit_id=hit_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                target=attrs.get("Target", ""),
                identity=attrs.get("Identity", attrs.get("identity", "")),
                source_tool=source,
            )
            current_hit_id = hit_id
        elif feature_l in {"cds", "match_part", "hsp"}:
            parent = attrs.get("Parent", "").split(",")[0] or current_hit_id
            if parent in hits:
                hits[parent].start = min(hits[parent].start, start)
                hits[parent].end = max(hits[parent].end, end)
    return list(hits.values())


def parse_protein_gffs(paths: Sequence[str]) -> List[ProteinHit]:
    hits: List[ProteinHit] = []
    for value in paths:
        source_label, path = split_named_path(str(value))
        hits.extend(parse_protein_gff(path, source_label))
    return hits


def protein_index(protein_hits: Sequence[ProteinHit], flank: int):
    items = []
    for hit in protein_hits:
        start = max(1, hit.start - flank)
        end = hit.end + flank
        items.append((hit.seqid, start, end, hit.hit_id))
    return build_interval_index(items), {hit.hit_id: hit for hit in protein_hits}


def load_validation_summary(path: str) -> Dict[str, Dict[str, str]]:
    if not path:
        return {}
    rows = {}
    for row in read_tsv(path):
        model_id = row.get("model_id") or row.get("object_id") or row.get("candidate_id") or ""
        if model_id:
            rows[model_id] = row
    return rows


def validation_for_model(model_id: str, validation_rows: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    row = validation_rows.get(model_id, {})
    status = row.get("fatal_qc_status") or row.get("validation_status") or row.get("overall_validation_status") or "not_assessed"
    fail_reasons = row.get("fatal_fail_reasons") or row.get("validation_fail_reasons") or row.get("blocking_reasons") or ""
    return {
        "fatal_qc_status": status,
        "validation_status": status,
        "validation_fail_reasons": fail_reasons,
        "cds_status": row.get("cds_status") or "",
        "cds_length_mod3": row.get("cds_length_mod3") or "",
        "internal_stop_status": row.get("internal_stop_status") or "",
        "phase_status": row.get("phase_status") or row.get("phase_consistency") or "",
        "sequence_status": row.get("sequence_status") or "",
    }


def load_eligibility_summary(path: str) -> Dict[str, Dict[str, str]]:
    if not path:
        return {}
    rows = {}
    for row in read_tsv(path):
        model_id = row.get("model_id") or row.get("object_id") or row.get("candidate_id") or ""
        if model_id:
            rows[model_id] = row
    return rows


def eligibility_for_model(model_id: str, eligibility_rows: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    row = eligibility_rows.get(model_id, {})
    if not row:
        return {
            "eligible_for_ai": "true",
            "eligibility_status": "not_assessed_allowed_for_legacy_compatibility",
            "eligibility_fail_reasons": "",
        }
    eligible = str(row.get("eligible_for_ai", "")).strip().lower() in {"true", "1", "yes", "y", "eligible"}
    return {
        "eligible_for_ai": "true" if eligible else "false",
        "eligibility_status": row.get("eligibility_status") or ("eligible" if eligible else "ineligible"),
        "eligibility_fail_reasons": row.get("hard_fail_reasons") or row.get("eligibility_fail_reasons") or row.get("blocking_reasons") or "",
    }


def count_junction_list(text: str) -> int:
    if not text or text == ".":
        return 0
    return len([item for item in text.split(";") if item and item != "."])


def preview_junction_list(text: str, limit: int = 5) -> List[str]:
    if not text or text == ".":
        return []
    return [item for item in text.split(";") if item and item != "."][:limit]


def parse_bool(text: str) -> bool:
    return str(text).strip().lower() in {"true", "1", "yes", "y"}


def empty_junction_support_summary(status: str) -> Dict[str, object]:
    return {
        "status": status,
        "n_model_junctions": "",
        "n_supported_junctions": "",
        "junction_support_fraction": "",
        "full_intron_chain_supported": "",
        "donor_supported": "",
        "acceptor_supported": "",
        "unsupported_junction_count": "",
        "unsupported_junction_preview": [],
        "nearby_novel_supported_junction_count": "",
        "nearby_novel_supported_junction_preview": [],
    }


def load_junction_support_summary(path: str, transcript_ids: Set[str]) -> Dict[str, Dict[str, object]]:
    if not path:
        return {}
    summaries: Dict[str, Dict[str, object]] = {}
    with open(path, "r", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx_id = row.get("transcript_id") or ""
            if transcript_ids and tx_id not in transcript_ids:
                continue
            unsupported = row.get("unsupported_junctions") or ""
            nearby = row.get("novel_supported_junctions_nearby") or ""
            summaries[tx_id] = {
                "status": "direct_current_summary_available",
                "n_model_junctions": row.get("n_model_junctions") or "",
                "n_supported_junctions": row.get("n_supported_junctions") or "",
                "junction_support_fraction": row.get("junction_support_fraction") or "",
                "full_intron_chain_supported": parse_bool(row.get("full_intron_chain_supported") or ""),
                "donor_supported": parse_bool(row.get("donor_supported") or ""),
                "acceptor_supported": parse_bool(row.get("acceptor_supported") or ""),
                "unsupported_junction_count": count_junction_list(unsupported),
                "unsupported_junction_preview": preview_junction_list(unsupported),
                "nearby_novel_supported_junction_count": count_junction_list(nearby),
                "nearby_novel_supported_junction_preview": preview_junction_list(nearby),
            }
    return summaries




def _safe_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return default


def _split_csv(text: object) -> List[str]:
    out: List[str] = []
    for item in str(text or "").replace(";", ",").split(","):
        item = item.strip()
        if item and item != ".":
            out.append(item)
    return out


def load_long_read_support_dir(path: str) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], Dict[Tuple[str, str], Dict[str, str]], Set[str]]:
    """Load optional long-read support tables keyed by original source IDs."""
    model_rows: Dict[Tuple[str, str], Dict[str, str]] = {}
    set_rows: Dict[Tuple[str, str], Dict[str, str]] = {}
    loaded_sources: Set[str] = set()
    if not path:
        return model_rows, set_rows, loaded_sources
    if not os.path.isdir(path):
        raise SystemExit("--long-read-support-dir does not exist: {0}".format(path))
    for source in sorted(os.listdir(path)):
        source_dir = os.path.join(path, source)
        if not os.path.isdir(source_dir):
            continue
        model_path = os.path.join(source_dir, "long_read_model_support_by_model.tsv")
        set_path = os.path.join(source_dir, "long_read_set_support_by_locus.tsv")
        if not os.path.exists(model_path) and not os.path.exists(set_path):
            continue
        loaded_sources.add(source)
        if os.path.exists(model_path):
            for row in read_tsv(model_path):
                model_set = row.get("model_set") or source
                model_id = row.get("model_id") or ""
                if model_id:
                    model_rows[(model_set, model_id)] = row
        if os.path.exists(set_path):
            for row in read_tsv(set_path):
                model_set = row.get("model_set") or source
                source_locus_id = row.get("locus_id") or ""
                if source_locus_id:
                    set_rows[(model_set, source_locus_id)] = row
    return model_rows, set_rows, loaded_sources


def empty_long_read_model_support(status: str) -> Dict[str, object]:
    return {
        "status": status,
        "support_class": "",
        "exact_chain_support": 0,
        "partial_junction_support": "",
        "conflicting_chain_count": 0,
        "best_chain_id": "",
        "risk_tags": [],
    }


def long_read_for_model(
    model: CandidateModel,
    long_read_model_rows: Dict[Tuple[str, str], Dict[str, str]],
    long_read_sources: Set[str],
) -> Dict[str, object]:
    if not long_read_model_rows and not long_read_sources:
        return empty_long_read_model_support("not_configured")
    if model.source not in long_read_sources:
        return empty_long_read_model_support("not_available_for_source")
    row = long_read_model_rows.get((model.source, model.transcript_id))
    if not row:
        return empty_long_read_model_support("missing_for_model")
    return {
        "status": "available",
        "support_class": row.get("support_class", ""),
        "exact_chain_support": _safe_int(row.get("exact_chain_support")),
        "partial_junction_support": row.get("partial_junction_support", ""),
        "conflicting_chain_count": _safe_int(row.get("conflicting_chain_count")),
        "best_chain_id": row.get("best_chain_id", ""),
        "risk_tags": _split_csv(row.get("risk_tags")),
    }


def add_long_read_evidence_flags(
    hard_flags: List[str],
    risk_tags: List[str],
    long_read_summary: Dict[str, object],
) -> Tuple[List[str], List[str]]:
    status = str(long_read_summary.get("status", ""))
    support_class = str(long_read_summary.get("support_class", ""))
    if status == "not_configured":
        return sorted(set(hard_flags)), sorted(set(risk_tags))
    if status == "not_available_for_source":
        risk_tags.append("long_read_not_available_for_source")
        return sorted(set(hard_flags)), sorted(set(risk_tags))
    if status == "missing_for_model":
        risk_tags.append("long_read_missing_for_model")
        return sorted(set(hard_flags)), sorted(set(risk_tags))

    if support_class == "exact_chain_supported":
        hard_flags.append("long_read_exact_chain_supported")
    elif support_class == "single_read_exact_chain":
        hard_flags.append("long_read_single_read_exact_chain")
        risk_tags.append("long_read_single_read_exact_chain")
    elif support_class == "all_junctions_supported_no_exact_chain":
        hard_flags.append("long_read_all_junctions_supported_no_exact_chain")
        risk_tags.append("long_read_no_exact_chain")
    elif support_class == "partial_junction_supported":
        hard_flags.append("long_read_partial_junction_supported")
        risk_tags.append("long_read_partial_junction_supported")
    elif support_class == "unsupported_by_long_reads":
        risk_tags.append("long_read_unsupported_by_long_reads")
    elif support_class == "mono_exon_not_assessable":
        risk_tags.append("long_read_mono_exon_not_assessable")

    for tag in long_read_summary.get("risk_tags") or []:
        value = "long_read_" + str(tag)
        risk_tags.append(value)
    return sorted(set(hard_flags)), sorted(set(risk_tags))


def aggregate_long_read_set_summary(
    source: str,
    models: Sequence[Dict[str, object]],
    long_read_set_rows: Dict[Tuple[str, str], Dict[str, str]],
) -> Dict[str, object]:
    status_values = [str((m.get("long_read_model_support") or {}).get("status", "not_configured")) for m in models]
    if not status_values or all(status == "not_configured" for status in status_values):
        return {"status": "not_configured", "support_level_counts": {}, "model_support_class_counts": {}, "risk_tags": []}
    if all(status == "not_available_for_source" for status in status_values):
        return {"status": "not_available_for_source", "support_level_counts": {}, "model_support_class_counts": {}, "risk_tags": ["long_read_not_available_for_source"]}

    gene_ids = sorted({str(m.get("source_gene_id", "")) for m in models if m.get("source_gene_id")})
    gene_rows = [long_read_set_rows[(source, gene_id)] for gene_id in gene_ids if (source, gene_id) in long_read_set_rows]
    support_levels = [row.get("support_level", "") for row in gene_rows if row.get("support_level")]
    model_support_classes = [str((m.get("long_read_model_support") or {}).get("support_class", "")) for m in models]
    support_level_counts = value_counts(support_levels)
    model_support_class_counts = value_counts([item for item in model_support_classes if item])
    best_supported: List[str] = []
    risk_tags: List[str] = []
    exact_supported = 0
    partial_supported = 0
    unsupported = 0
    conflicting = 0
    for row in gene_rows:
        exact_supported += _safe_int(row.get("exact_supported_transcripts"))
        partial_supported += _safe_int(row.get("partial_supported_transcripts"))
        unsupported += _safe_int(row.get("unsupported_transcripts"))
        conflicting += _safe_int(row.get("conflicting_chain_count"))
        best_supported.extend(_split_csv(row.get("best_supported_model_ids")))
        for tag in _split_csv(row.get("risk_tags")):
            risk_tags.append("long_read_" + tag)

    if not gene_rows:
        for model in models:
            summary = model.get("long_read_model_support") or {}
            support_class = str(summary.get("support_class", ""))
            if support_class in {"exact_chain_supported", "single_read_exact_chain"}:
                exact_supported += 1
            elif support_class in {"all_junctions_supported_no_exact_chain", "partial_junction_supported"}:
                partial_supported += 1
            elif support_class == "unsupported_by_long_reads":
                unsupported += 1
            conflicting += _safe_int(summary.get("conflicting_chain_count"))
            for tag in summary.get("risk_tags") or []:
                risk_tags.append("long_read_" + str(tag))

    mono_count = support_level_counts.get("mono_exon_not_assessable", 0) + model_support_class_counts.get("mono_exon_not_assessable", 0)
    return {
        "status": "available" if gene_rows or any(status == "available" for status in status_values) else "missing_for_set",
        "source_gene_count_with_long_read_rows": len(gene_rows),
        "exact_supported_transcripts": exact_supported,
        "partial_supported_transcripts": partial_supported,
        "unsupported_transcripts": unsupported,
        "mono_exon_not_assessable_count": mono_count,
        "best_supported_model_ids": sorted(set(best_supported))[:20],
        "conflicting_chain_count": conflicting,
        "support_level_counts": support_level_counts,
        "model_support_class_counts": model_support_class_counts,
        "risk_tags": sorted(set(risk_tags)),
    }


def empty_optional_evidence_slots() -> Dict[str, object]:
    return {
        "te_overlap_summary": {"status": "not_configured"},
        "synteny_summary": {"status": "not_configured"},
        "domain_summary": {"status": "not_configured"},
    }


def junction_support_for_model(model: CandidateModel, support_rows: Dict[str, Dict[str, object]]) -> Dict[str, object]:
    if not support_rows:
        return empty_junction_support_summary("not_configured")
    if model.source != "current":
        return empty_junction_support_summary("not_available_for_de_novo_model")
    return support_rows.get(model.transcript_id, empty_junction_support_summary("missing_for_current_model"))


def add_auxiliary_evidence_flags(
    hard_flags: List[str],
    risk_tags: List[str],
    junction_support_summary: Dict[str, object],
) -> Tuple[List[str], List[str]]:
    if junction_support_summary.get("status") == "direct_current_summary_available":
        fraction = junction_support_summary.get("junction_support_fraction")
        try:
            fraction_f = float(fraction)
        except (TypeError, ValueError):
            fraction_f = -1.0
        if junction_support_summary.get("full_intron_chain_supported"):
            hard_flags.append("model_level_full_intron_chain_supported")
        elif fraction_f == 0:
            risk_tags.append("model_level_no_junction_support")
        elif 0 < fraction_f < 1:
            hard_flags.append("model_level_partial_junction_support")
            risk_tags.append("model_level_partial_junction_support")
    elif junction_support_summary.get("status") in {"missing_for_current_model", "not_available_for_de_novo_model"}:
        risk_tags.append(str(junction_support_summary.get("status")))
    return sorted(set(hard_flags)), sorted(set(risk_tags))


def model_junction_summary(
    model: CandidateModel,
    junctions: Dict[Tuple[str, int, int, str], JunctionEvidence],
    boundary_window: int,
    nearby_window: int,
) -> Dict[str, object]:
    """汇总 RNA 剪接连接证据，但不创建坐标补丁。"""
    tx = model.transcript
    introns = tx.introns()
    if not introns or not junctions:
        return {
            "status": "not_assessable",
            "exact_supported": 0,
            "not_observed": 0,
            "fraction": "0.0000",
            "conflicting_count": 0,
            "nearby_alternative_count": 0,
        }

    exact_supported = 0
    not_observed = 0
    conflicting_count = 0
    nearby_alt_count = 0
    for intron_start, intron_end in introns:
        exact_key = (tx.seqid, intron_start, intron_end, tx.strand)
        exact_key_unstranded = (tx.seqid, intron_start, intron_end, ".")
        if exact_key in junctions or exact_key_unstranded in junctions:
            exact_supported += 1
            continue
        nearby = nearby_junctions(junctions, tx.seqid, tx.strand, intron_start, intron_end, nearby_window)
        observed_conflict = False
        for junction in nearby[:5]:
            left_delta = junction.intron_start - intron_start
            right_delta = junction.intron_end - intron_end
            if junction.intron_start == intron_start and junction.intron_end == intron_end:
                continue
            observed_conflict = True
            if not (abs(left_delta) <= boundary_window and abs(right_delta) <= boundary_window):
                nearby_alt_count += 1
        if observed_conflict:
            conflicting_count += 1
        else:
            not_observed += 1

    if conflicting_count:
        status = "conflicting_junction"
    elif exact_supported == len(introns):
        status = "full_junction_support"
    elif exact_supported:
        status = "partial_junction_support"
    else:
        status = "none_observed"

    return {
        "status": status,
        "exact_supported": exact_supported,
        "not_observed": not_observed,
        "fraction": pct(exact_supported, len(introns)),
        "conflicting_count": conflicting_count,
        "nearby_alternative_count": nearby_alt_count,
    }


def model_protein_summary(model: CandidateModel, protein_indexes, proteins_by_id: Dict[str, ProteinHit]) -> Dict[str, object]:
    if not proteins_by_id:
        return {"count": 0, "best_overlap_bp": 0, "best_target": "", "best_identity": ""}
    (index, starts) = protein_indexes
    tx = model.transcript
    hits = query_overlaps(tx.seqid, tx.start, tx.end, index, starts)
    if not hits:
        return {"count": 0, "best_overlap_bp": 0, "best_target": "", "best_identity": ""}
    best_id, _start, _end, best_bp = hits[0]
    best = proteins_by_id[best_id]
    return {
        "count": len(hits),
        "best_overlap_bp": best_bp,
        "best_target": best.target,
        "best_identity": best.identity,
    }


def classify_model_evidence(
    junction_summary: Dict[str, object],
    protein_summary: Dict[str, object],
    model: CandidateModel,
    validation_summary: Dict[str, str],
) -> Tuple[List[str], List[str], str]:
    intron_count = len(model.transcript.introns())
    protein_count = int(protein_summary["count"])
    cds_count = len(model.transcript.sorted_cds())
    exon_count = len(model.transcript.sorted_exons())
    junction_status = str(junction_summary["status"])
    risk_tags: List[str] = []
    hard_flags: List[str] = []
    if junction_status == "full_junction_support":
        hard_flags.append("all_introns_junction_supported")
    elif junction_status == "partial_junction_support":
        hard_flags.append("partial_introns_junction_supported")
    elif junction_status == "conflicting_junction":
        risk_tags.append("conflicting_RNA_junction")
    elif junction_status == "none_observed":
        risk_tags.append("RNA_junction_none_observed_not_negative_evidence")
    else:
        if intron_count:
            risk_tags.append("RNA_junction_not_assessable")
        else:
            risk_tags.append("single_exon_model_no_splice_junction_to_assess")
    if int(junction_summary["nearby_alternative_count"]):
        risk_tags.append("nearby_alternative_junctions")
    if protein_count:
        hard_flags.append("protein_overlap_present")
    else:
        risk_tags.append("protein_evidence_none_observed")
    if cds_count:
        hard_flags.append("CDS_present")
    else:
        risk_tags.append("no_CDS_features")
    if exon_count == 0:
        risk_tags.append("no_exon_features")
    if model.source == "current":
        recommended = "base_model_candidate"
    elif junction_status in {"full_junction_support", "partial_junction_support"} or protein_count:
        recommended = "supported_alternative_candidate"
    else:
        recommended = "available_candidate_needs_validation"
    return sorted(set(hard_flags)), sorted(set(risk_tags)), recommended


def model_sort_key(model_record: Dict[str, object]) -> Tuple[int, int, int, int, str]:
    junction_rank = {
        "full_junction_support": 0,
        "partial_junction_support": 1,
        "conflicting_junction": 2,
        "none_observed": 3,
        "not_assessable": 4,
    }.get(str(model_record["rna_junction_summary"]["status"]), 5)
    validation_rank = {
        "pass": 0,
        "warning": 1,
        "not_assessed": 2,
        "missing": 3,
        "fail": 4,
    }.get(str(model_record["validation_summary"]["validation_status"]), 3)
    protein_rank = 0 if int(model_record["protein_summary"]["count"]) else 1
    source_rank = source_priority_rank(model_record["source"])
    return (validation_rank, junction_rank, protein_rank, source_rank, str(model_record["model_id"]))


def value_counts(values: Sequence[object]) -> Dict[str, int]:
    counter: Dict[str, int] = defaultdict(int)
    for value in values:
        counter[str(value)] += 1
    return dict(counter)


def build_source_candidate_model_sets(
    locus_id: str,
    selectable_models: Sequence[Dict[str, object]],
    long_read_set_rows: Dict[Tuple[str, str], Dict[str, str]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    by_source: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for model in selectable_models:
        source = str(model.get("source", "")) or "unknown"
        by_source[source].append(model)

    sets: List[Dict[str, object]] = []
    rows: List[Dict[str, object]] = []
    for source in sorted(by_source, key=lambda s: (s != "current", s)):
        models = sorted(by_source[source], key=lambda m: str(m.get("model_id", "")))
        model_ids = [str(m.get("model_id", "")) for m in models if m.get("model_id")]
        gene_ids = sorted({str(m.get("source_gene_id", "")) for m in models if m.get("source_gene_id")})
        starts = [int(m.get("start") or 0) for m in models if m.get("start")]
        ends = [int(m.get("end") or 0) for m in models if m.get("end")]
        seqids = sorted({str(m.get("seqid", "")) for m in models if m.get("seqid")})
        strands = sorted({str(m.get("strand", "")) for m in models if m.get("strand")})
        validation_statuses = [str((m.get("validation_summary") or {}).get("validation_status", "not_assessed")) for m in models]
        rna_statuses = [str((m.get("rna_junction_summary") or {}).get("status", "not_assessable")) for m in models]
        protein_supported = [m for m in models if int((m.get("protein_summary") or {}).get("count") or 0) > 0]
        risk_tags = sorted({str(tag) for m in models for tag in (m.get("risk_tags") or []) if str(tag)})
        hard_flags = sorted({str(tag) for m in models for tag in (m.get("hard_flags") or []) if str(tag)})
        validation_counts = value_counts(validation_statuses)
        rna_counts = value_counts(rna_statuses)
        long_read_summary = aggregate_long_read_set_summary(source, models, long_read_set_rows)
        risk_tags = sorted(set(risk_tags + [str(tag) for tag in (long_read_summary.get("risk_tags") or []) if str(tag)]))
        set_id = "{0}_set".format(safe_id(source))
        set_role = "current_gene_set" if source == "current" else "candidate_gene_set"
        record = {
            "set_id": set_id,
            "source": source,
            "set_role": set_role,
            "model_ids": model_ids,
            "gene_ids": gene_ids,
            "seqid": seqids[0] if len(seqids) == 1 else ",".join(seqids),
            "start": min(starts) if starts else "",
            "end": max(ends) if ends else "",
            "strand": strands[0] if len(strands) == 1 else ",".join(strands),
            "gene_count": len(gene_ids),
            "transcript_count": len(model_ids),
            "validation_summary": {
                "pass": validation_counts.get("pass", 0),
                "warning": validation_counts.get("warning", 0),
                "fail": validation_counts.get("fail", 0),
                "not_assessed": validation_counts.get("not_assessed", 0),
            },
            "rna_junction_summary": {
                "full_junction_support": rna_counts.get("full_junction_support", 0),
                "partial_junction_support": rna_counts.get("partial_junction_support", 0),
                "conflicting_junction": rna_counts.get("conflicting_junction", 0),
                "none_observed": rna_counts.get("none_observed", 0),
                "not_assessable": rna_counts.get("not_assessable", 0),
            },
            "protein_supported_model_count": len(protein_supported),
            "short_read_junction_summary": {
                "full_junction_support": rna_counts.get("full_junction_support", 0),
                "partial_junction_support": rna_counts.get("partial_junction_support", 0),
                "conflicting_junction": rna_counts.get("conflicting_junction", 0),
                "none_observed": rna_counts.get("none_observed", 0),
                "not_assessable": rna_counts.get("not_assessable", 0),
            },
            "homolog_protein_summary": {
                "supported_model_count": len(protein_supported),
                "status": "available" if protein_supported else "none_observed",
            },
            "long_read_set_support_summary": long_read_summary,
            "optional_evidence_slots": empty_optional_evidence_slots(),
            "evidence_summary": {
                "short_read_transcriptome": {
                    "junction_summary": {
                        "full_junction_support": rna_counts.get("full_junction_support", 0),
                        "partial_junction_support": rna_counts.get("partial_junction_support", 0),
                        "conflicting_junction": rna_counts.get("conflicting_junction", 0),
                        "none_observed": rna_counts.get("none_observed", 0),
                        "not_assessable": rna_counts.get("not_assessable", 0),
                    },
                },
                "long_read_transcriptome": long_read_summary,
                "homolog_protein": {
                    "supported_model_count": len(protein_supported),
                    "status": "available" if protein_supported else "none_observed",
                },
                "future_optional": empty_optional_evidence_slots(),
            },
            "risk_tags": risk_tags,
            "hard_flags": hard_flags,
        }
        sets.append(record)
        rows.append(
            {
                "locus_id": locus_id,
                "set_id": set_id,
                "source": source,
                "set_role": set_role,
                "gene_count": len(gene_ids),
                "transcript_count": len(model_ids),
                "model_ids": ",".join(model_ids),
                "gene_ids": ",".join(gene_ids),
                "seqid": record["seqid"],
                "start": record["start"],
                "end": record["end"],
                "strand": record["strand"],
                "validation_pass_count": validation_counts.get("pass", 0),
                "validation_warning_count": validation_counts.get("warning", 0),
                "validation_fail_count": validation_counts.get("fail", 0),
                "validation_not_assessed_count": validation_counts.get("not_assessed", 0),
                "rna_full_junction_support_count": rna_counts.get("full_junction_support", 0),
                "rna_partial_junction_support_count": rna_counts.get("partial_junction_support", 0),
                "rna_conflicting_junction_count": rna_counts.get("conflicting_junction", 0),
                "rna_none_observed_count": rna_counts.get("none_observed", 0),
                "protein_supported_count": len(protein_supported),
                "long_read_exact_supported_transcripts": long_read_summary.get("exact_supported_transcripts", ""),
                "long_read_partial_supported_transcripts": long_read_summary.get("partial_supported_transcripts", ""),
                "long_read_unsupported_transcripts": long_read_summary.get("unsupported_transcripts", ""),
                "long_read_mono_exon_not_assessable_count": long_read_summary.get("mono_exon_not_assessable_count", ""),
                "long_read_conflicting_chain_count": long_read_summary.get("conflicting_chain_count", ""),
                "long_read_support_level_counts": ";".join("{0}:{1}".format(k, v) for k, v in sorted((long_read_summary.get("support_level_counts") or {}).items())),
                "long_read_best_supported_model_ids": ",".join(long_read_summary.get("best_supported_model_ids") or []),
                "risk_tags": ",".join(risk_tags),
                "hard_flags": ",".join(hard_flags),
                "candidate_set_mode": "source",
                "supporting_sources": source,
                "representative_source": source,
                "equivalent_source_set_ids": set_id,
                "structure_signature_key": "|".join(sorted(str(m.get("signature", "")) for m in models if m.get("signature"))),
            }
        )
    return sets, rows


def source_set_priority(model_set: Dict[str, object]) -> Tuple[int, str]:
    source = str(model_set.get("source", ""))
    return (source_priority_rank(source), source)


def source_set_structure_key(model_set: Dict[str, object], models_by_id: Dict[str, Dict[str, object]]) -> str:
    signatures = []
    for model_id in model_set.get("model_ids") or []:
        model = models_by_id.get(str(model_id), {})
        signatures.append(str(model.get("signature", "")))
    return "|".join(sorted(sig for sig in signatures if sig))


def collapse_structure_consensus_sets(
    locus_id: str,
    source_sets: Sequence[Dict[str, object]],
    source_rows: Sequence[Dict[str, object]],
    selectable_models: Sequence[Dict[str, object]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    models_by_id = {str(model.get("model_id", "")): model for model in selectable_models if model.get("model_id")}
    row_by_set = {str(row.get("set_id", "")): row for row in source_rows}
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for model_set in source_sets:
        key = source_set_structure_key(model_set, models_by_id)
        groups[key or str(model_set.get("set_id", ""))].append(model_set)

    collapsed_sets: List[Dict[str, object]] = []
    collapsed_rows: List[Dict[str, object]] = []
    ordered_groups = sorted(groups.items(), key=lambda item: min(source_set_priority(row) for row in item[1]))
    for index, (structure_key, members) in enumerate(ordered_groups, start=1):
        members = sorted(members, key=source_set_priority)
        canonical = dict(members[0])
        canonical_source = str(canonical.get("source", ""))
        supporting_sources = sorted({str(member.get("source", "")) for member in members if member.get("source")})
        equivalent_set_ids = sorted(str(member.get("set_id", "")) for member in members if member.get("set_id"))
        if len(members) > 1:
            set_id = "structure_set_{0:02d}".format(index)
            if canonical.get("set_role") == "current_gene_set":
                set_id = "current_structure_set_{0:02d}".format(index)
        else:
            set_id = str(canonical.get("set_id", ""))

        canonical.update(
            {
                "set_id": set_id,
                "candidate_set_mode": "structure_consensus",
                "supporting_sources": supporting_sources,
                "representative_source": canonical_source,
                "equivalent_source_set_ids": equivalent_set_ids,
                "structure_signature_key": structure_key,
            }
        )
        risk_tags = sorted(set(str(tag) for tag in (canonical.get("risk_tags") or [])))
        hard_flags = sorted(set(str(flag) for flag in (canonical.get("hard_flags") or [])))
        if len(members) > 1:
            hard_flags.append("structure_consensus_multi_source_support")
            risk_tags = [tag for tag in risk_tags if tag != "similar_structure_seen_in_multiple_candidate_sets"]
        canonical["risk_tags"] = sorted(set(risk_tags))
        canonical["hard_flags"] = sorted(set(hard_flags))
        collapsed_sets.append(canonical)

        base_row = dict(row_by_set.get(str(members[0].get("set_id", "")), {}))
        base_row.update(
            {
                "locus_id": locus_id,
                "set_id": set_id,
                "source": canonical.get("source", ""),
                "set_role": canonical.get("set_role", ""),
                "risk_tags": ",".join(canonical.get("risk_tags") or []),
                "hard_flags": ",".join(canonical.get("hard_flags") or []),
                "candidate_set_mode": "structure_consensus",
                "supporting_sources": ",".join(supporting_sources),
                "representative_source": canonical_source,
                "equivalent_source_set_ids": ",".join(equivalent_set_ids),
                "structure_signature_key": structure_key,
            }
        )
        collapsed_rows.append(base_row)
    return collapsed_sets, collapsed_rows


def build_candidate_model_sets(
    locus_id: str,
    selectable_models: Sequence[Dict[str, object]],
    long_read_set_rows: Dict[Tuple[str, str], Dict[str, str]],
    candidate_set_mode: str,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    source_sets, source_rows = build_source_candidate_model_sets(locus_id, selectable_models, long_read_set_rows)
    if candidate_set_mode == "source":
        return source_sets, source_rows
    if candidate_set_mode != "structure_consensus":
        raise ValueError("Unsupported candidate_set_mode: {0}".format(candidate_set_mode))
    return collapse_structure_consensus_sets(locus_id, source_sets, source_rows, selectable_models)


def classify_complexity(
    locus_id: str,
    models: Sequence[CandidateModel],
    max_candidate_models: int,
    max_unique_structures: int,
    max_transcripts_per_source: int,
    max_locus_bp: int,
    task_mode: str,
) -> Tuple[Dict[str, object], bool]:
    seqid = models[0].transcript.seqid
    start = min(m.transcript.start for m in models)
    end = max(m.transcript.end for m in models)
    sources = sorted(set(m.source for m in models))
    signatures = sorted(set(m.signature for m in models))
    by_source = defaultdict(set)
    for model in models:
        by_source[model.source].add(model.gene_id)
    max_gene_count_per_source = max((len(genes) for genes in by_source.values()), default=0)
    current_count = sum(1 for m in models if m.source == "current")
    current_gene_count = len({m.gene_id for m in models if m.source == "current" and m.gene_id})
    tool_count = len(models) - current_count

    hard_risk_reasons = []
    ai_arbitrable_reasons = []
    if len(models) > max_candidate_models:
        ai_arbitrable_reasons.append("many_candidate_models_AI_arbitrable")
    if len(signatures) > max_unique_structures:
        ai_arbitrable_reasons.append("many_unique_structures_AI_arbitrable")
    if max_gene_count_per_source > max_transcripts_per_source:
        ai_arbitrable_reasons.append("multi_gene_source_set_arbitration")
    if end - start + 1 > max_locus_bp:
        hard_risk_reasons.append("locus_span_too_large")
    if current_count > 1 and tool_count > 1 and max_gene_count_per_source > 1:
        ai_arbitrable_reasons.append("multi_current_multi_tool_gene_set_arbitration")

    reasons = sorted(set(hard_risk_reasons + ai_arbitrable_reasons))
    if hard_risk_reasons:
        complexity_class = "complex_risk_flagged"
    elif ai_arbitrable_reasons:
        complexity_class = "ai_arbitrable_complex"
    else:
        complexity_class = "auto_arbitration_eligible"
    action_space = "choose_candidate_set" if task_mode == "de_novo_annotation" else "keep_current_set,choose_candidate_set"
    row = {
        "locus_id": locus_id,
        "seqid": seqid,
        "start": start,
        "end": end,
        "candidate_model_count": len(models),
        "source_count": len(sources),
        "current_model_count": current_count,
        "current_gene_count": current_gene_count,
        "tool_model_count": tool_count,
        "unique_structure_count": len(signatures),
        "max_gene_count_per_source": max_gene_count_per_source,
        "complexity_class": complexity_class,
        "complexity_reasons": ",".join(reasons),
        "risk_flagged": "true" if hard_risk_reasons else "false",
        "risk_triggers": ",".join(sorted(set(hard_risk_reasons))),
        "decision_required": "false" if hard_risk_reasons else "true",
        "recommended_action_space": "risk_flagged" if hard_risk_reasons else action_space,
        "conflict_class": "complex_risk_flagged" if hard_risk_reasons else ("ai_arbitrable_complex" if ai_arbitrable_reasons else "none"),
        "recommended_action": "risk_flagged" if hard_risk_reasons else "AI_must_select_best_gene_set",
    }
    return row, bool(hard_risk_reasons)


WEAK_LONG_WINDOW_EDGE_TYPES = {"span_overlap", "protein_bridge", "RNA_junction_bridge"}


def apply_window_risk_gate(
    complexity_row: Dict[str, object],
    edge_records: Sequence[Dict[str, object]],
    ai_soft_locus_bp: int,
    ai_hard_locus_bp: int,
    max_ai_current_genes_in_long_window: int,
    max_ai_models_per_source_in_long_window: int,
) -> Tuple[Dict[str, object], List[str]]:
    """标记不应进入普通 AI 仲裁的窗口。

    V1 保留类似 IGV 的多轨道上下文。只要窗口仍有清楚的生物学边界，
    长窗口是允许的；但如果长窗口包含过多当前注释基因、相对于来源数
    过多的模型，或只由弱桥接边连接，则在进入 AI 前标记为 risk。
    """
    row = dict(complexity_row)
    risks: List[str] = []
    try:
        span_bp = int(row.get("end") or 0) - int(row.get("start") or 0) + 1
    except (TypeError, ValueError):
        span_bp = 0
    current_gene_count = int(row.get("current_gene_count") or 0)
    candidate_model_count = int(row.get("candidate_model_count") or 0)
    source_count = max(1, int(row.get("source_count") or 1))
    edge_types = {str(edge.get("edge_type", "")) for edge in edge_records if edge.get("edge_type")}

    if ai_hard_locus_bp and span_bp > ai_hard_locus_bp:
        risks.append("window_span_above_hard_limit_for_AI")

    long_window = bool(ai_soft_locus_bp and span_bp > ai_soft_locus_bp)
    if long_window:
        if max_ai_current_genes_in_long_window and current_gene_count > max_ai_current_genes_in_long_window:
            risks.append("long_window_too_many_current_genes_for_AI")
        dynamic_model_limit = source_count * max_ai_models_per_source_in_long_window if max_ai_models_per_source_in_long_window else 0
        if dynamic_model_limit and candidate_model_count > dynamic_model_limit:
            risks.append("long_window_too_many_candidate_models_for_AI")
        if edge_types and edge_types.issubset(WEAK_LONG_WINDOW_EDGE_TYPES):
            risks.append("long_window_weak_bridge_only_for_AI")

    if not risks:
        return row, []

    existing_reasons = [item for item in str(row.get("complexity_reasons", "")).split(",") if item]
    existing_risks = [item for item in str(row.get("risk_triggers", "")).split(",") if item]
    all_reasons = sorted(set(existing_reasons + risks))
    risk_triggers = sorted(set(existing_risks + risks))
    row.update(
        {
            "complexity_class": "complex_risk_flagged",
            "complexity_reasons": ",".join(all_reasons),
            "risk_flagged": "true",
            "risk_triggers": ",".join(risk_triggers),
            "decision_required": "false",
            "recommended_action_space": "risk_flagged",
            "conflict_class": "window_boundary_risk",
            "recommended_action": "risk_flagged",
        }
    )
    return row, sorted(set(risks))


def is_supported_model(model_record: Dict[str, object]) -> bool:
    junction_status = str(model_record["rna_junction_summary"]["status"])
    validation_status = str(model_record["validation_summary"]["validation_status"])
    if validation_status == "fail":
        return False
    if junction_status in {"full_junction_support", "partial_junction_support"}:
        return True
    return int(model_record["protein_summary"]["count"]) > 0


def apply_risk_router(
    complexity_row: Dict[str, object],
    selectable_models: Sequence[Dict[str, object]],
    task_mode: str,
) -> Dict[str, object]:
    triggers = [item for item in str(complexity_row.get("complexity_reasons", "")).split(",") if item]
    conflict_class = "complex_locus" if triggers else "none"
    hard_risk_triggers = [item for item in str(complexity_row.get("risk_triggers", "")).split(",") if item]
    hard_risk = bool(hard_risk_triggers)

    supported = [model for model in selectable_models if is_supported_model(model)]
    supported_signatures = set(str(model["signature"]) for model in supported)
    if len(supported_signatures) > 1:
        triggers.append("multiple_supported_distinct_structures")
        conflict_class = "supported_model_conflict"

    rna_supported = [
        model for model in selectable_models
        if str(model["rna_junction_summary"]["status"]) in {"full_junction_support", "partial_junction_support"}
    ]
    protein_supported = [
        model for model in selectable_models
        if int(model["protein_summary"]["count"]) > 0
    ]
    if rna_supported and protein_supported:
        rna_signatures = set(str(model["signature"]) for model in rna_supported)
        protein_signatures = set(str(model["signature"]) for model in protein_supported)
        if rna_signatures.isdisjoint(protein_signatures):
            triggers.append("RNA_protein_support_different_structures")
            conflict_class = "RNA_protein_conflict"

    if any(str(model["validation_summary"]["validation_status"]) == "fail" for model in selectable_models):
        triggers.append("one_or_more_candidate_validation_fail")
        if conflict_class == "none":
            conflict_class = "validation_conflict"

    triggers = sorted(set(triggers))
    hard_risk_triggers = sorted(set(hard_risk_triggers))
    risk_flagged = hard_risk
    if hard_risk:
        decision_required = "false"
        decision_requirement = "risk_flagged_no_ai_decision"
        action_space = "risk_flagged"
        next_gate = "risk_review_after_pipeline"
    elif triggers:
        decision_required = "true"
        if task_mode == "de_novo_annotation":
            decision_requirement = "AI_must_select_best_candidate_gene_set_with_complexity_or_conflict_flags"
            action_space = "choose_candidate_set"
        else:
            decision_requirement = "AI_must_select_best_gene_set_with_complexity_or_conflict_flags"
            action_space = "keep_current_set,choose_candidate_set"
        next_gate = "AI_gene_set_arbitration_with_complexity_flags_then_set_validation"
    else:
        decision_required = "true"
        if task_mode == "de_novo_annotation":
            decision_requirement = "AI_must_select_one_best_candidate_gene_set_id"
            action_space = "choose_candidate_set"
        else:
            decision_requirement = "AI_must_select_one_best_gene_set_id"
            action_space = "keep_current_set,choose_candidate_set"
        next_gate = "AI_gene_set_arbitration_then_set_validation"

    row = dict(complexity_row)
    row.update(
        {
            "risk_flagged": "true" if risk_flagged else "false",
            "risk_triggers": ",".join(hard_risk_triggers if hard_risk else triggers),
            "decision_required": decision_required,
            "decision_requirement": decision_requirement,
            "recommended_action_space": action_space,
            "conflict_class": conflict_class,
            "recommended_next_gate": next_gate,
        }
    )
    return row


def main() -> None:
    args = parse_args()
    global SOURCE_PRIORITY
    SOURCE_PRIORITY = parse_source_priority(args.source_priority)
    if args.task_mode == "correction" and not args.current_gff:
        raise SystemExit("--current-gff is required when --task-mode correction")
    if args.task_mode == "de_novo_annotation" and args.current_gff:
        raise SystemExit("--current-gff is not allowed when --task-mode de_novo_annotation; provide all inputs with --tool-gff.")
    if not args.tool_gff:
        raise SystemExit("At least one --tool-gff candidate source is required.")
    models = load_candidate_models(args.current_gff, args.tool_gff, args.signature_mode, args.task_mode)
    gene_filter = load_gene_filter(args.gene_list)
    region_filter = load_regions(args.region_table, args.region_flank)
    models = filter_models_for_subset(models, gene_filter, region_filter)
    if not models:
        raise SystemExit("No candidate models remained after subset filtering.")
    junctions = load_junctions(args.junctions, args.min_junction_support_count, args.min_junction_sample_count)
    protein_hits = parse_protein_gffs(args.protein_gff)
    loci, edges_by_pair = build_loci(
        models,
        args.min_overlap_fraction,
        junctions,
        args.boundary_window,
        args.nearby_window,
        protein_hits,
        args.protein_flank,
        region_filter,
    )
    protein_indexes, proteins_by_id = protein_index(protein_hits, args.protein_flank)
    validation_rows = load_validation_summary(args.validation_summary)
    eligibility_rows = load_eligibility_summary(args.eligibility_summary)
    current_transcript_ids = {model.transcript_id for model in models if model.source == "current"}
    junction_support_rows = load_junction_support_summary(args.junction_support_summary, current_transcript_ids)
    long_read_model_rows, long_read_set_rows, long_read_sources = load_long_read_support_dir(args.long_read_support_dir)

    model_rows = []
    set_rows = []
    complex_rows = []
    summary_rows = []
    cards = []
    all_edge_rows = [edge for edges in edges_by_pair.values() for edge in edges]
    edges_by_locus: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for edge in all_edge_rows:
        edges_by_locus[str(edge.get("locus_id", ""))].append(edge)

    for locus_index, locus_models in enumerate(loci, 1):
        locus_id = "locus_{0:08d}".format(locus_index)
        locus_models = sorted(locus_models, key=lambda item: (item.source != "current", item.source, item.gene_id, item.transcript_id))
        complexity_row, is_complex = classify_complexity(
            locus_id,
            locus_models,
            args.max_candidate_models,
            args.max_unique_structures,
            args.max_transcripts_per_source,
            args.max_locus_bp,
            args.task_mode,
        )

        selectable_models = []
        excluded_candidate_models = []
        locus_risks: List[str] = []
        for model in locus_models:
            junction_summary = model_junction_summary(model, junctions, args.boundary_window, args.nearby_window)
            protein_summary = model_protein_summary(model, protein_indexes, proteins_by_id)
            validation_summary = validation_for_model(model.model_id, validation_rows)
            junction_support_summary = junction_support_for_model(model, junction_support_rows)
            long_read_summary = long_read_for_model(model, long_read_model_rows, long_read_sources)
            eligibility_summary = eligibility_for_model(model.model_id, eligibility_rows)
            hard_flags, risk_tags, recommended_role = classify_model_evidence(
                junction_summary,
                protein_summary,
                model,
                validation_summary,
            )
            hard_flags, risk_tags = add_auxiliary_evidence_flags(
                hard_flags,
                risk_tags,
                junction_support_summary,
            )
            hard_flags, risk_tags = add_long_read_evidence_flags(
                hard_flags,
                risk_tags,
                long_read_summary,
            )
            locus_risks.extend(risk_tags)
            model_record = {
                "model_id": model.model_id,
                "source": model.source,
                "source_file": model.source_file,
                "source_gene_id": model.gene_id,
                "source_transcript_id": model.transcript_id,
                "seqid": model.transcript.seqid,
                "start": model.transcript.start,
                "end": model.transcript.end,
                "strand": model.transcript.strand,
                "exons": model.transcript.sorted_exons(),
                "cds": model.transcript.sorted_cds(),
                "introns": model.transcript.introns(),
                "signature": model.signature,
                "hard_flags": hard_flags,
                "risk_tags": risk_tags,
                "recommended_model_role": recommended_role,
                "rna_junction_summary": junction_summary,
                "model_level_junction_support_summary": junction_support_summary,
                "short_read_transcriptome_summary": {
                    "junction_summary": junction_summary,
                    "model_level_junction_support_summary": junction_support_summary,
                },
                "long_read_model_support": long_read_summary,
                "protein_summary": protein_summary,
                "homolog_protein_summary": protein_summary,
                "optional_evidence_slots": empty_optional_evidence_slots(),
                "validation_summary": validation_summary,
                "eligibility_summary": eligibility_summary,
            }
            fatal_qc_status = str(validation_summary.get("fatal_qc_status") or validation_summary.get("validation_status") or "not_assessed")
            if fatal_qc_status == "fail":
                excluded = dict(model_record)
                excluded["exclusion_reason"] = "fatal_qc_fail"
                excluded_candidate_models.append(excluded)
            elif eligibility_summary["eligible_for_ai"] == "true":
                selectable_models.append(model_record)
            else:
                excluded = dict(model_record)
                excluded["exclusion_reason"] = "pre_ai_eligibility_fail"
                excluded_candidate_models.append(excluded)
                locus_risks.append("pre_ai_eligibility_excluded_model")
            model_rows.append(
                {
                    "locus_id": locus_id,
                    "model_id": model.model_id,
                    "source": model.source,
                    "gene_id": model.gene_id,
                    "transcript_id": model.transcript_id,
                    "seqid": model.transcript.seqid,
                    "start": model.transcript.start,
                    "end": model.transcript.end,
                    "strand": model.transcript.strand,
                    "exon_count": len(model.transcript.sorted_exons()),
                    "cds_count": len(model.transcript.sorted_cds()),
                    "intron_count": len(model.transcript.introns()),
                    "signature": model.signature,
                    "rna_junction_status": junction_summary["status"],
                    "exact_junction_supported_introns": junction_summary["exact_supported"],
                    "not_observed_introns": junction_summary["not_observed"],
                    "junction_support_fraction": junction_summary["fraction"],
                    "conflicting_junction_count": junction_summary["conflicting_count"],
                    "nearby_alternative_junction_count": junction_summary["nearby_alternative_count"],
                    "protein_overlap_count": protein_summary["count"],
                    "protein_best_overlap_bp": protein_summary["best_overlap_bp"],
                    "protein_best_target": protein_summary["best_target"],
                    "protein_best_identity": protein_summary["best_identity"],
                    "long_read_status": long_read_summary["status"],
                    "long_read_support_class": long_read_summary["support_class"],
                    "long_read_exact_chain_support": long_read_summary["exact_chain_support"],
                    "long_read_partial_junction_support": long_read_summary["partial_junction_support"],
                    "long_read_conflicting_chain_count": long_read_summary["conflicting_chain_count"],
                    "long_read_risk_tags": ",".join(long_read_summary.get("risk_tags") or []),
                    "model_level_junction_status": junction_support_summary["status"],
                    "model_level_junction_support_fraction": junction_support_summary["junction_support_fraction"],
                    "model_level_full_intron_chain_supported": junction_support_summary["full_intron_chain_supported"],
                    "validation_status": validation_summary["validation_status"],
                    "validation_fail_reasons": validation_summary["validation_fail_reasons"],
                    "hard_flags": ",".join(hard_flags),
                    "risk_tags": ",".join(risk_tags),
                    "recommended_model_role": recommended_role,
                }
            )

        selectable_models.sort(key=model_sort_key)
        selectable_model_ids = [item["model_id"] for item in selectable_models]
        candidate_model_sets, candidate_set_rows = build_candidate_model_sets(locus_id, selectable_models, long_read_set_rows, args.candidate_set_mode)
        set_rows.extend(candidate_set_rows)
        complexity_row, window_risks = apply_window_risk_gate(
            complexity_row,
            edges_by_locus.get(locus_id, []),
            args.ai_soft_locus_bp,
            args.ai_hard_locus_bp,
            args.max_ai_current_genes_in_long_window,
            args.max_ai_models_per_source_in_long_window,
        )
        locus_risks.extend(window_risks)
        router_row = apply_risk_router(complexity_row, selectable_models, args.task_mode)
        if not selectable_models:
            router_row.update(
                {
                    "risk_flagged": "true",
                    "risk_triggers": "no_ai_eligible_candidate_models",
                    "decision_required": "false",
                    "decision_requirement": "risk_flagged_no_ai_eligible_candidate_models",
                    "recommended_action_space": "risk_flagged",
                    "conflict_class": "no_ai_eligible_candidate_models",
                    "recommended_next_gate": "risk_review_after_pre_ai_eligibility_filter",
                }
            )
        if router_row["risk_flagged"] == "true":
            complex_rows.append(router_row)
        decision_requirement = router_row["decision_requirement"]
        recommended_next_gate = router_row["recommended_next_gate"]
        summary_rows.append(
            {
                "locus_id": locus_id,
                "seqid": router_row["seqid"],
                "start": router_row["start"],
                "end": router_row["end"],
                "candidate_model_count": router_row["candidate_model_count"],
                "source_count": router_row["source_count"],
                "complexity_class": router_row["complexity_class"],
                "risk_flagged": router_row["risk_flagged"],
                "risk_triggers": router_row["risk_triggers"],
                "decision_required": router_row["decision_required"],
                "decision_requirement": decision_requirement,
                "recommended_action_space": router_row["recommended_action_space"],
                "conflict_class": router_row["conflict_class"],
                "selectable_model_ids": ",".join(selectable_model_ids),
                "recommended_next_gate": recommended_next_gate,
                "risk_tags": ",".join(sorted(set(locus_risks))),
            }
        )
        cards.append(
            {
                "task_mode": args.task_mode,
                "locus_id": locus_id,
                "locus": {
                    "seqid": router_row["seqid"],
                    "start": router_row["start"],
                    "end": router_row["end"],
                },
                "risk_router": {
                    "risk_flagged": router_row["risk_flagged"] == "true",
                    "risk_triggers": [
                        item for item in str(router_row["risk_triggers"]).split(",") if item
                    ],
                    "decision_required": router_row["decision_required"] == "true",
                    "decision_requirement": decision_requirement,
                    "recommended_action_space": [
                        item for item in str(router_row["recommended_action_space"]).split(",") if item
                    ],
                    "conflict_class": router_row["conflict_class"],
                },
                "complexity": router_row,
                "arbitration_mode": "gene_set_arbitration",
                "unit_construction": {
                    "mode": "evidence_aware_arbitration_unit",
                    "edge_count": len(edges_by_locus.get(locus_id, [])),
                    "edge_types": sorted(set(str(edge.get("edge_type")) for edge in edges_by_locus.get(locus_id, []))),
                },
                "candidate_model_sets": candidate_model_sets,
                "candidate_models": selectable_models,
                "excluded_candidate_models": excluded_candidate_models,
                "allowed_decision_scope": [
                    item for item in str(router_row["recommended_action_space"]).split(",") if item
                ],
                "decision_requirement": decision_requirement,
                "scientific_boundaries": [
                    "AI must choose only listed set_id values for gene-set arbitration.",
                    "AI must not invent coordinates, exons, CDS, junctions, transcripts, gene models, or gene sets.",
                    "RNA none_observed or not_assessable is not negative evidence against current annotation.",
                    "In correction mode, current annotation is a candidate input, not truth.",
                    "In de_novo_annotation mode, no current annotation is defined.",
                    "risk_flagged loci are not ordinary AI-selection tasks and should be audited downstream.",
                ],
                "recommended_next_gate": recommended_next_gate,
            }
        )

    if args.output_units:
        write_tsv(args.output_units, UNIT_FIELDS, unit_rows_for_loci(loci, edges_by_pair))
    if args.output_unit_edges:
        edge_rows = sorted(
            [edge for edges in edges_by_pair.values() for edge in edges],
            key=lambda row: (str(row.get("locus_id", "")), str(row.get("model_id_a", "")), str(row.get("model_id_b", "")), str(row.get("edge_type", ""))),
        )
        write_tsv(args.output_unit_edges, EDGE_FIELDS, edge_rows)
    if args.output_model_sets:
        write_tsv(args.output_model_sets, SET_FIELDS, set_rows)
    write_tsv(args.output_models, MODEL_FIELDS, model_rows)
    write_tsv(args.output_complex_loci, COMPLEX_FIELDS, complex_rows)
    write_tsv(args.output_summary, SUMMARY_FIELDS, summary_rows)
    write_jsonl(args.output_jsonl, cards)


if __name__ == "__main__":
    main()
