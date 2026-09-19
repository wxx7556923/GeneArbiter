#!/usr/bin/env python3
# script_id_md5: 884cb623ed7a1b9069fbc137943943da
# created: 2026-06-24
# modified: 2026-07-10
# owner: project
# status: project_code
# purpose: 共享的 AI/LLM compact evidence card 构建核心函数；由 task-specific 入口调用。
# inputs: 命令行参数指定的 TSV/JSONL/GFF3/FASTA 输入。
# outputs: 命令行参数指定的 TSV/JSONL/GFF3/Markdown 输出。
# notes: 新增 run-level evidence normalization；缺失整类 evidence 时按运行级分母归一化，不按 locus/model 动态缩放。

"""共享的 compact evidence card 构建核心。

完整卡保留可追溯坐标，便于调试和后续导出。紧凑卡移除较长的
外显子/CDS/内含子坐标列表，只保留受约束 AI 模型选择需要的结构化证据。
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from gff_utils import write_jsonl, write_tsv


PROJECT_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_NIPPONBARE_TE3_TSV = PROJECT_ROOT / "reference" / "nipponbare" / "mh63_candidate_te3_support.tsv"


SUMMARY_FIELDS = [
    "card_id",
    "set_count",
    "component_count",
    "topology_class",
    "main_question",
    "possible_split_merge",
    "main_conflict",
    "discriminator_status",
    "recommended_default",
    "routing_decision",
    "auto_gate",
    "auto_selected_set_id",
    "auto_representative_real_set_id",
    "auto_risk_level",
    "review_router_class",
    "manual_review_novel_gene_candidate",
    "manual_review_deletion_candidate",
    "manual_review_unsupported_locus",
    "deletion_policy",
    "recommended_manual_action",
]

REVIEW_ROUTER_FIELDS = [
    "card_id",
    "router_class",
    "current_model_count",
    "tool_model_count",
    "current_selectable_model_count",
    "tool_selectable_model_count",
    "ai_set_count",
    "routing_decision",
    "auto_gate",
    "sent_to_ai",
    "manual_review_novel_gene_candidate",
    "manual_review_deletion_candidate",
    "manual_review_unsupported_locus",
    "deletion_allowed",
    "deletion_policy",
    "best_evidence_level",
    "unsupported_locus_reason_tags",
    "discovery_evidence_level",
    "review_priority",
    "recommended_manual_action",
    "review_reason_tags",
]

EVIDENCE_LEVEL_RANK = {
    "none": 0,
    "poor": 1,
    "weak": 2,
    "moderate": 3,
    "strong": 4,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-jsonl", required=True, help="Full model_arbitration_cards.jsonl")
    parser.add_argument("--output-jsonl", required=True, help="AI annotation card JSONL")
    parser.add_argument("--output-summary", required=True, help="AI annotation card summary TSV")
    parser.add_argument("--output-set-trace", required=True, help="Sidecar JSONL mapping AI set IDs to real member sets and deterministic representative real set IDs")
    parser.add_argument("--output-auto-decisions", required=True, help="JSONL decisions made deterministically before API calls by Gate A/B/C")
    parser.add_argument("--output-locus-review-router", required=True, help="TSV review router for deletion-disabled current-only loci and tool-only novel-gene candidates")
    parser.add_argument(
        "--include-risk-flagged",
        dest="include_risk_flagged",
        action="store_true",
        help="Include upstream risk_flagged / decision_required=false cards. Default skips them so ordinary AI selection only sees decisionable cards.",
    )
    parser.add_argument("--source-priority", action="append", default=[], help="Optional comma-separated source priority order used to choose representative sets after structure collapse.")
    parser.add_argument(
        "--scoring-normalization",
        choices=["legacy", "run_level"],
        default="run_level",
        help="Score normalization mode. run_level removes evidence classes not provided for this run from the score denominator.",
    )
    parser.add_argument(
        "--active-evidence-types",
        action="append",
        default=[],
        help="Comma-separated run-level evidence classes available to scoring: short_read,long_read,protein. Omitted means all are active for legacy direct use.",
    )
    return parser.parse_args()

def read_jsonl(path: str) -> List[Dict[str, object]]:
    rows = []
    with open(path, "r") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def parse_bool(value: object) -> bool | None:
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "supported", "present"}:
        return True
    if text in {"0", "false", "no", "n", "absent", "none", "none_observed", "no_support", "not_observed"}:
        return False
    return None


def numeric_positive(value: object) -> bool | None:
    try:
        return float(value or 0) > 0
    except (TypeError, ValueError):
        return None


def nipponbare_row_has_support(row: Dict[str, str]) -> bool | None:
    for field in ["has_rna_support", "rna_supported", "transcript_supported", "junction_supported"]:
        if field in row:
            parsed = parse_bool(row.get(field))
            if parsed is not None:
                return parsed
    for field in ["rna_support_level", "transcript_support_level", "junction_support_level"]:
        if field in row:
            text = str(row.get(field) or "").strip().lower()
            if text in {"strong", "moderate", "weak", "present", "supported"}:
                return True
            if text in {"absent", "none", "none_observed", "no_support", "not_observed"}:
                return False
    values = []
    for field in ["supported_junction_count", "junction_read_count", "expressed_sample_count", "nip_gene_count_sum"]:
        if field in row:
            value = numeric_positive(row.get(field))
            if value is not None:
                values.append(value)
    if values:
        return any(values)
    return None


def add_index_row(index: Dict[str, List[Dict[str, str]]], key: str, row: Dict[str, str]) -> None:
    if key and key != "NA":
        index.setdefault(key, []).append(row)


def nipponbare_support_paths(args: argparse.Namespace) -> List[Path]:
    paths: List[Path] = []
    if args.nipponbare_rna_support_tsv:
        paths.append(Path(args.nipponbare_rna_support_tsv))
    if not args.disable_default_nipponbare_te3_support and DEFAULT_NIPPONBARE_TE3_TSV.exists():
        paths.append(DEFAULT_NIPPONBARE_TE3_TSV)
    paths.extend(Path(item) for item in args.nipponbare_te3_support_tsv)
    seen = set()
    out: List[Path] = []
    for item in paths:
        key = str(item.resolve()) if item.exists() else str(item)
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return out


def load_nipponbare_support(paths: Sequence[Path]) -> Dict[str, object]:
    by_locus: Dict[str, List[Dict[str, str]]] = {}
    by_gene: Dict[str, List[Dict[str, str]]] = {}
    loaded_paths = []
    for path in paths:
        if not path.exists():
            continue
        loaded_paths.append(str(path))
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for raw_row in reader:
                row = dict(raw_row)
                row["_support_source_path"] = str(path)
                if "mh63_candidate_gene_id" in row or "nip_gene_count_sum" in row:
                    row["_support_table_type"] = "nipponbare_te3_projection"
                else:
                    row["_support_table_type"] = "nipponbare_rna_support"
                add_index_row(by_locus, row.get("locus_id") or "", row)
                for key in ["current_gene_id", "gene_id", "mh63_gene_id", "mh63_candidate_gene_id", "mh63_candidate_transcript_id"]:
                    add_index_row(by_gene, row.get(key) or "", row)
    return {"by_locus": by_locus, "by_gene": by_gene, "loaded_paths": loaded_paths}


def unique_nipponbare_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    seen = set()
    out = []
    for row in rows:
        key = (
            row.get("_support_source_path", ""),
            row.get("locus_id", ""),
            row.get("current_gene_id", ""),
            row.get("gene_id", ""),
            row.get("mh63_gene_id", ""),
            row.get("mh63_candidate_gene_id", ""),
            row.get("mh63_candidate_transcript_id", ""),
            row.get("irgs_gene_id", ""),
            row.get("irgs_transcript_id", ""),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def rows_have_support(rows: Sequence[Dict[str, str]]) -> bool | None:
    values = [nipponbare_row_has_support(row) for row in rows]
    known = [value for value in values if value is not None]
    if not known:
        return None
    return any(known)


def count_values(rows: Sequence[Dict[str, str]], field: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for row in rows:
        value = str(row.get(field) or "")
        if value:
            out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items()))


def numeric_max(rows: Sequence[Dict[str, str]], field: str) -> float | None:
    values = []
    for row in rows:
        try:
            values.append(float(row.get(field) or 0))
        except (TypeError, ValueError):
            pass
    return max(values) if values else None


def current_real_gene_ids(set_alias_map: Sequence[Dict[str, object]]) -> List[str]:
    out: List[str] = []
    for row in set_alias_map:
        if row.get("set_id") == "current_set" or row.get("real_source") == "current":
            out.extend(str(item) for item in row.get("real_gene_ids") or [])
    return sorted(set(out))


def summarize_nipponbare_support(locus_id: str, set_alias_map: Sequence[Dict[str, object]], support: Dict[str, object] | None, treat_missing_as_absent: bool) -> Tuple[bool | None, Dict[str, object]]:
    if not support or not support.get("loaded_paths"):
        return None, {"status": "not_configured"}
    by_locus = support.get("by_locus") or {}
    by_gene = support.get("by_gene") or {}
    rows: List[Dict[str, str]] = []
    rows.extend(by_locus.get(locus_id, []))
    matched_gene_ids = []
    for gene_id in current_real_gene_ids(set_alias_map):
        matched = by_gene.get(gene_id, [])
        if matched:
            matched_gene_ids.append(gene_id)
            rows.extend(matched)
    rows = unique_nipponbare_rows(rows)
    if not rows:
        return (True if treat_missing_as_absent else None), {
            "status": "missing_row",
            "matched_gene_ids": matched_gene_ids,
            "treated_missing_as_absent": bool(treat_missing_as_absent),
            "loaded_paths": support.get("loaded_paths") or [],
        }
    has_support = rows_have_support(rows)
    support_values = [nipponbare_row_has_support(row) for row in rows]
    counterpart_ids = sorted({str(row.get("irgs_gene_id") or "") for row in rows if row.get("irgs_gene_id")})
    return (None if has_support is None else not has_support), {
        "status": "available",
        "has_transcript_support": has_support,
        "matched_row_count": len(rows),
        "support_positive_row_count": sum(1 for value in support_values if value is True),
        "support_absent_row_count": sum(1 for value in support_values if value is False),
        "support_unknown_row_count": sum(1 for value in support_values if value is None),
        "support_table_types": count_values(rows, "_support_table_type"),
        "projection_relation_counts": count_values(rows, "mh63_vs_liftoff_relation"),
        "risk_tag_counts": count_values(rows, "risk_tags"),
        "max_nip_gene_count_sum": numeric_max(rows, "nip_gene_count_sum"),
        "counterpart_gene_ids_preview": counterpart_ids[:10],
        "matched_gene_ids": matched_gene_ids,
        "loaded_paths": support.get("loaded_paths") or [],
    }


def interval_len(items: Sequence[Sequence[object]]) -> int:
    total = 0
    for item in items:
        if len(item) >= 2:
            total += int(item[1]) - int(item[0]) + 1
    return total


def pair_intervals(items: Sequence[Sequence[object]]) -> List[Tuple[int, int]]:
    return sorted((int(item[0]), int(item[1])) for item in items if len(item) >= 2)


def current_models(models: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    return [model for model in models if model.get("source") == "current"]


def overlap_bp(a: Dict[str, object], b: Dict[str, object]) -> int:
    if a.get("seqid") != b.get("seqid"):
        return 0
    start = max(int(a.get("start") or 0), int(b.get("start") or 0))
    end = min(int(a.get("end") or 0), int(b.get("end") or 0))
    return max(0, end - start + 1)


def best_current_for_model(model: Dict[str, object], currents: Sequence[Dict[str, object]]) -> Dict[str, object]:
    if not currents:
        return {}
    return sorted(currents, key=lambda current: (-overlap_bp(model, current), str(current.get("model_id", ""))))[0]


def compare_to_current(model: Dict[str, object], currents: Sequence[Dict[str, object]]) -> Dict[str, object]:
    if model.get("source") == "current":
        return {
            "closest_current_model_id": model.get("model_id", ""),
            "relation_tags": ["current_model"],
            "span_delta_start": 0,
            "span_delta_end": 0,
            "exon_count_delta": 0,
            "intron_count_delta": 0,
            "cds_count_delta": 0,
            "cds_length_delta": 0,
        }
    current = best_current_for_model(model, currents)
    if not current:
        return {
            "closest_current_model_id": "",
            "relation_tags": ["no_current_model_in_locus"],
            "span_delta_start": "",
            "span_delta_end": "",
            "exon_count_delta": "",
            "intron_count_delta": "",
            "cds_count_delta": "",
            "cds_length_delta": "",
        }

    model_exons = pair_intervals(model.get("exons") or [])
    current_exons = pair_intervals(current.get("exons") or [])
    model_introns = pair_intervals(model.get("introns") or [])
    current_introns = pair_intervals(current.get("introns") or [])
    model_cds = pair_intervals(model.get("cds") or [])
    current_cds = pair_intervals(current.get("cds") or [])

    relation_tags = []
    if str(model.get("signature", "")) == str(current.get("signature", "")):
        relation_tags.append("same_structure_signature")
    if model_introns == current_introns:
        relation_tags.append("same_intron_chain")
    else:
        relation_tags.append("alternative_intron_chain")
    if model_exons == current_exons:
        relation_tags.append("same_exon_chain")
    if model_cds == current_cds:
        relation_tags.append("same_CDS_chain")
    elif len(model_cds) == len(current_cds):
        relation_tags.append("CDS_boundary_difference")
    elif len(model_cds) > len(current_cds):
        relation_tags.append("extra_CDS_blocks")
    elif len(model_cds) < len(current_cds):
        relation_tags.append("fewer_CDS_blocks")

    span_delta_start = int(model.get("start") or 0) - int(current.get("start") or 0)
    span_delta_end = int(model.get("end") or 0) - int(current.get("end") or 0)
    if model_introns == current_introns and (span_delta_start or span_delta_end):
        relation_tags.append("terminal_boundary_difference_only")
    if int(model.get("start") or 0) <= int(current.get("start") or 0) and int(model.get("end") or 0) >= int(current.get("end") or 0):
        relation_tags.append("candidate_contains_current_span")
    elif int(model.get("start") or 0) >= int(current.get("start") or 0) and int(model.get("end") or 0) <= int(current.get("end") or 0):
        relation_tags.append("candidate_within_current_span")

    return {
        "closest_current_model_id": current.get("model_id", ""),
        "relation_tags": sorted(set(relation_tags)),
        "span_delta_start": span_delta_start,
        "span_delta_end": span_delta_end,
        "exon_count_delta": len(model_exons) - len(current_exons),
        "intron_count_delta": len(model_introns) - len(current_introns),
        "cds_count_delta": len(model_cds) - len(current_cds),
        "cds_length_delta": interval_len(model_cds) - interval_len(current_cds),
    }


def _as_int(value: object) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def _block_pairs(model: Dict[str, object], field: str) -> List[Tuple[int, int]]:
    blocks: List[Tuple[int, int]] = []
    for item in model.get(field) or []:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            start = _as_int(item[0])
            end = _as_int(item[1])
            if start and end and start <= end:
                blocks.append((start, end))
    return sorted(set(blocks))


def _model_len(model: Dict[str, object]) -> int:
    start = _as_int(model.get("start"))
    end = _as_int(model.get("end"))
    return max(1, end - start + 1) if start and end else 1


def _span_overlap_fractions(left: Dict[str, object], right: Dict[str, object]) -> Tuple[float, float]:
    left_start = _as_int(left.get("start"))
    left_end = _as_int(left.get("end"))
    right_start = _as_int(right.get("start"))
    right_end = _as_int(right.get("end"))
    if not left_start or not left_end or not right_start or not right_end:
        return 0.0, 0.0
    overlap = min(left_end, right_end) - max(left_start, right_start) + 1
    if overlap <= 0:
        return 0.0, 0.0
    return overlap / float(_model_len(left)), overlap / float(_model_len(right))


def _jaccard(left: Sequence[Tuple[int, int]], right: Sequence[Tuple[int, int]]) -> float:
    left_set = set(left)
    right_set = set(right)
    if not left_set and not right_set:
        return 1.0
    if not left_set or not right_set:
        return 0.0
    return len(left_set & right_set) / float(len(left_set | right_set))


def _same_context(left: Dict[str, object], right: Dict[str, object]) -> bool:
    return (
        str(left.get("seqid", "")) == str(right.get("seqid", ""))
        and str(left.get("strand", "")) == str(right.get("strand", ""))
        and str(left.get("strand", "")) in {"+", "-"}
    )


def component_similarity_class(left: Dict[str, object], right: Dict[str, object]) -> str:
    """Classify whether two existing models describe the same biological component.

    This is compression only: it never creates coordinates or a rebuilt model.
    Split/merge alternatives should remain separate topology sets rather than be
    forced into one component.
    """
    if not _same_context(left, right):
        return "different_context"
    left_frac, right_frac = _span_overlap_fractions(left, right)
    if left_frac < 0.55 or right_frac < 0.55:
        return "separate_or_split_merge_component"

    if str(left.get("signature", "")) and str(left.get("signature", "")) == str(right.get("signature", "")):
        return "exact_same_structure"

    left_introns = _block_pairs(left, "introns")
    right_introns = _block_pairs(right, "introns")
    left_exons = _block_pairs(left, "exons")
    right_exons = _block_pairs(right, "exons")
    left_cds = _block_pairs(left, "cds")
    right_cds = _block_pairs(right, "cds")

    if left_introns and left_introns == right_introns:
        return "same_intron_chain_boundary_variant"
    if not left_introns and not right_introns and left_frac >= 0.70 and right_frac >= 0.70:
        return "single_exon_overlap_component"
    if _jaccard(left_introns, right_introns) >= 0.67 and left_frac >= 0.65 and right_frac >= 0.65:
        return "mostly_same_splice_component"
    if max(_jaccard(left_exons, right_exons), _jaccard(left_cds, right_cds)) >= 0.60 and left_frac >= 0.70 and right_frac >= 0.70:
        return "mostly_same_coding_component"
    return "distinct_component_or_topology_conflict"


def build_structure_groups(models: Sequence[Dict[str, object]]) -> Dict[str, str]:
    """Compress near-equivalent models into biological component groups.

    Exact identity is not required. Models with the same intron chain, strong
    reciprocal span overlap, or mostly matching exon/CDS structure are grouped so
    AI sees component/topology summaries instead of raw source-by-source models.
    """
    by_id = {str(model.get("model_id", "")): model for model in models if model.get("model_id")}
    graph: Dict[str, set[str]] = defaultdict(set)
    for model_id in by_id:
        graph[model_id].add(model_id)
    ids = sorted(by_id)
    for index, left_id in enumerate(ids):
        left = by_id[left_id]
        for right_id in ids[index + 1:]:
            relation = component_similarity_class(left, by_id[right_id])
            if relation in {
                "exact_same_structure",
                "same_intron_chain_boundary_variant",
                "single_exon_overlap_component",
                "mostly_same_splice_component",
                "mostly_same_coding_component",
            }:
                graph[left_id].add(right_id)
                graph[right_id].add(left_id)

    out: Dict[str, str] = {}
    seen: set[str] = set()
    component_index = 0
    for model_id in ids:
        if model_id in seen:
            continue
        component_index += 1
        stack = [model_id]
        members: List[str] = []
        seen.add(model_id)
        while stack:
            current = stack.pop()
            members.append(current)
            for neighbor in sorted(graph.get(current, [])):
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        component_id = "component_{0:03d}".format(component_index)
        for member in members:
            out[member] = component_id
    return out


def protein_support_status(model: Dict[str, object]) -> str:
    count = int((model.get("protein_summary") or {}).get("count") or 0)
    return "protein_overlap_present" if count else "protein_overlap_none_observed"


def compact_model(model: Dict[str, object], structure_groups: Dict[str, str], currents: Sequence[Dict[str, object]]) -> Dict[str, object]:
    validation = model.get("validation_summary") or {}
    rna = model.get("rna_junction_summary") or {}
    model_junction = model.get("model_level_junction_support_summary") or {}
    protein = model.get("protein_summary") or {}
    exons = model.get("exons") or []
    introns = model.get("introns") or []
    cds = model.get("cds") or []
    signature = str(model.get("signature", ""))
    return {
        "model_id": model.get("model_id", ""),
        "source": model.get("source", ""),
        "source_transcript_id": model.get("source_transcript_id", ""),
        "span": {
            "seqid": model.get("seqid", ""),
            "start": model.get("start", ""),
            "end": model.get("end", ""),
            "strand": model.get("strand", ""),
        },
        "structure_group_id": structure_groups.get(str(model.get("model_id", "")), "component_unknown"),
        "component_group_ids": [structure_groups.get(str(model.get("model_id", "")), "component_unknown")],
        "structure_summary": {
            "exon_count": len(exons),
            "intron_count": len(introns),
            "cds_count": len(cds),
            "cds_length": interval_len(cds),
        },
        "difference_from_current": compare_to_current(model, currents),
        "validation_status": validation.get("validation_status", "not_assessed"),
        "orf_status": validation.get("orf_status", ""),
        "phase_status": validation.get("phase_status", ""),
        "splice_motif_status": validation.get("splice_motif_status", ""),
        "protein_concordance_status": validation.get("protein_concordance_status", ""),
        "validation_fail_reasons": validation.get("validation_fail_reasons", ""),
        "rna_junction_status": rna.get("status", "not_assessable"),
        "junction_support_fraction": rna.get("fraction", "0.0000"),
        "conflicting_junction_count": rna.get("conflicting_count", 0),
        "model_level_junction_status": model_junction.get("status", "not_configured"),
        "model_level_junction_support_fraction": model_junction.get("junction_support_fraction", ""),
        "model_level_full_intron_chain_supported": model_junction.get("full_intron_chain_supported", ""),
        "model_level_unsupported_junction_count": model_junction.get("unsupported_junction_count", ""),
        "model_level_nearby_novel_junction_count": model_junction.get("nearby_novel_supported_junction_count", ""),
        "model_level_unsupported_junction_preview": model_junction.get("unsupported_junction_preview", []),
        "model_level_nearby_novel_junction_preview": model_junction.get("nearby_novel_supported_junction_preview", []),
        "protein_support_status": protein_support_status(model),
        "long_read_model_support": model.get("long_read_model_support") or {"status": "not_configured"},
        "protein_best_target": protein.get("best_target", ""),
        "protein_best_identity": protein.get("best_identity", ""),
        "evidence_levels": {
            "coding_validation": validation.get("validation_status", "not_assessed"),
            "short_read_junction": rna.get("status", "not_assessable"),
            "long_read_transcriptome": long_read_support_level(model.get("long_read_model_support") or {"status": "not_configured"}, 1),
            "homolog_protein": "present" if int(protein.get("count") or 0) else "none_observed",
        },
        "hard_flags": model.get("hard_flags") or [],
        "risk_tags": model.get("risk_tags") or [],
    }


def allowed_decisions(value: object) -> List[str]:
    if isinstance(value, str):
        decisions = [item.strip() for item in value.split(",") if item.strip()]
    elif isinstance(value, Sequence):
        decisions = [str(item).strip() for item in value if str(item).strip()]
    else:
        decisions = []
    return [decision for decision in decisions if decision in {"keep_current", "choose_candidate_model", "keep_current_set", "choose_candidate_set", "select_candidate_set", "propose_delete_current_set", "risk_flagged"}]


def ids_for(models: Sequence[Dict[str, object]], predicate) -> List[str]:
    return [str(model.get("model_id", "")) for model in models if predicate(model)]


def _float_or_none(value: object):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int_value(value: object) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def validation_support_level(validation: Dict[str, object], transcript_count: int) -> str:
    fail = _int_value(validation.get("fail"))
    warning = _int_value(validation.get("warning"))
    passed = _int_value(validation.get("pass"))
    total = max(1, transcript_count)
    if fail:
        return "fail"
    if warning:
        return "warning"
    if passed / total >= 0.95:
        return "strong"
    if passed / total >= 0.6:
        return "moderate"
    if passed:
        return "weak"
    return "not_assessed"


def short_read_support_level(summary: Dict[str, object], transcript_count: int) -> str:
    total = max(1, transcript_count)
    full = _int_value(summary.get("full_junction_support"))
    partial = _int_value(summary.get("partial_junction_support"))
    conflict = _int_value(summary.get("conflicting_junction"))
    none = _int_value(summary.get("none_observed"))
    na = _int_value(summary.get("not_assessable"))
    if conflict:
        return "conflict"
    if full / total >= 0.95:
        return "strong"
    if (full + partial) / total >= 0.8:
        return "moderate"
    if full + partial:
        return "weak"
    if na and not none:
        return "not_assessable"
    if none:
        return "none_observed"
    return "not_configured"


def long_read_support_level(summary: Dict[str, object], transcript_count: int) -> str:
    status = str(summary.get("status", "not_configured"))
    if status != "available":
        return status
    support_class = str(summary.get("support_class", ""))
    if support_class == "exact_chain_supported":
        return "strong_exact_chain"
    if support_class == "single_read_exact_chain":
        return "moderate_exact_chain_single_read"
    if support_class == "all_junctions_supported_no_exact_chain":
        return "moderate_junction"
    if support_class == "partial_junction_supported":
        return "weak_partial"
    if support_class == "mono_exon_not_assessable":
        return "not_assessable_mono_exon"
    if support_class == "unsupported_by_long_reads":
        return "unsupported"
    total = max(1, transcript_count)
    exact = _int_value(summary.get("exact_supported_transcripts"))
    partial = _int_value(summary.get("partial_supported_transcripts"))
    unsupported = _int_value(summary.get("unsupported_transcripts"))
    mono = _int_value(summary.get("mono_exon_not_assessable_count"))
    if exact / total >= 0.8:
        return "strong_exact_chain"
    if exact:
        return "moderate_exact_chain"
    if (exact + partial) / total >= 0.8:
        return "moderate_junction"
    if partial:
        return "weak_partial"
    if mono and not unsupported:
        return "not_assessable_mono_exon"
    if unsupported:
        return "unsupported"
    return "available_no_support"


def protein_support_level(supported_count: object, transcript_count: int) -> str:
    supported = _int_value(supported_count)
    total = max(1, transcript_count)
    if supported / total >= 0.8:
        return "strong"
    if supported:
        return "moderate"
    return "none_observed"


def evidence_conflict_summary(compact_models: Sequence[Dict[str, object]]) -> Dict[str, object]:
    current_ids = [
        model["model_id"]
        for model in compact_models
        if model.get("source") == "current" or model.get("model_role") == "current_model"
    ]
    validation_pass = ids_for(compact_models, lambda model: model.get("validation_status") == "pass")
    validation_warning = ids_for(compact_models, lambda model: model.get("validation_status") == "warning")
    validation_fail = ids_for(compact_models, lambda model: model.get("validation_status") == "fail")
    rna_supported = ids_for(
        compact_models,
        lambda model: model.get("rna_junction_status") in {"full_junction_support", "partial_junction_support"},
    )
    protein_supported = ids_for(compact_models, lambda model: model.get("protein_support_status") == "protein_overlap_present")
    conflicting_junction = ids_for(compact_models, lambda model: model.get("rna_junction_status") == "conflicting_junction")
    model_level_full_junction = ids_for(
        compact_models,
        lambda model: model.get("model_level_full_intron_chain_supported") is True,
    )

    rna_structures = {
        model["structure_group_id"]
        for model in compact_models
        if model["model_id"] in rna_supported
    }
    protein_structures = {
        model["structure_group_id"]
        for model in compact_models
        if model["model_id"] in protein_supported
    }
    supported_structures = rna_structures | protein_structures
    return {
        "current_model_ids": current_ids,
        "validation_pass_model_ids": validation_pass,
        "validation_warning_model_ids": validation_warning,
        "models_with_validation_fail": validation_fail,
        "rna_supported_model_ids": rna_supported,
        "protein_supported_model_ids": protein_supported,
        "models_with_conflicting_junction": conflicting_junction,
        "model_level_full_junction_supported_model_ids": model_level_full_junction,
        "RNA_protein_conflict": bool(rna_structures and protein_structures and rna_structures.isdisjoint(protein_structures)),
        "multiple_supported_distinct_structures": len(supported_structures) > 1,
        "current_supported_by_RNA": bool(set(current_ids) & set(rna_supported)),
        "current_supported_by_model_level_junction_summary": bool(set(current_ids) & set(model_level_full_junction)),
        "current_supported_by_protein": bool(set(current_ids) & set(protein_supported)),
    }


def compact_model_set(model_set: Dict[str, object]) -> Dict[str, object]:
    validation = model_set.get("validation_summary") or {}
    rna = model_set.get("rna_junction_summary") or {}
    transcript_count = _int_value(model_set.get("transcript_count"))
    validation_summary = {
        "pass": validation.get("pass", 0),
        "warning": validation.get("warning", 0),
        "fail": validation.get("fail", 0),
        "not_assessed": validation.get("not_assessed", 0),
    }
    short_read_summary = model_set.get("short_read_junction_summary") or {
        "full_junction_support": rna.get("full_junction_support", 0),
        "partial_junction_support": rna.get("partial_junction_support", 0),
        "conflicting_junction": rna.get("conflicting_junction", 0),
        "none_observed": rna.get("none_observed", 0),
        "not_assessable": rna.get("not_assessable", 0),
    }
    long_read_summary = model_set.get("long_read_set_support_summary") or {"status": "not_configured"}
    protein_count = model_set.get("protein_supported_model_count", 0)
    homolog_summary = model_set.get("homolog_protein_summary") or {
        "supported_model_count": protein_count,
        "status": "available" if int(protein_count or 0) else "none_observed",
    }
    record = {
        "set_id": model_set.get("set_id", ""),
        "source": model_set.get("source", ""),
        "set_role": model_set.get("set_role", ""),
        "gene_count": model_set.get("gene_count", 0),
        "transcript_count": model_set.get("transcript_count", 0),
        "model_ids": model_set.get("model_ids") or [],
        "gene_ids": model_set.get("gene_ids") or [],
        "span": {
            "seqid": model_set.get("seqid", ""),
            "start": model_set.get("start", ""),
            "end": model_set.get("end", ""),
            "strand": model_set.get("strand", ""),
        },
        "validation_summary": validation_summary,
        "rna_junction_summary": dict(short_read_summary),
        "protein_supported_model_count": protein_count,
        "homolog_protein_summary": homolog_summary,
        "short_read_junction_summary": dict(short_read_summary),
        "long_read_set_support_summary": long_read_summary,
        "optional_evidence_slots": model_set.get("optional_evidence_slots") or {
            "te_overlap_summary": {"status": "not_configured"},
            "synteny_summary": {"status": "not_configured"},
            "domain_summary": {"status": "not_configured"},
        },
        "evidence_levels": {
            "coding_validation": validation_support_level(validation_summary, transcript_count),
            "short_read_junction": short_read_support_level(short_read_summary, transcript_count),
            "long_read_transcriptome": long_read_support_level(long_read_summary, transcript_count),
            "homolog_protein": protein_support_level(protein_count, transcript_count),
        },
        "hard_flags": model_set.get("hard_flags") or [],
        "risk_tags": model_set.get("risk_tags") or [],
        "candidate_set_mode": model_set.get("candidate_set_mode", "source"),
        "supporting_sources": model_set.get("supporting_sources") or [],
        "representative_source": model_set.get("representative_source", model_set.get("source", "")),
        "equivalent_source_set_ids": model_set.get("equivalent_source_set_ids") or [],
        "structure_signature_key": model_set.get("structure_signature_key", ""),
    }
    return record



def _letters(index: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    out = ""
    value = index
    while value > 0:
        value -= 1
        out = letters[value % 26] + out
        value //= 26
    return out or "A"


def build_set_aliases(compact_sets: Sequence[Dict[str, object]], task_mode: str) -> Tuple[Dict[str, str], List[Dict[str, object]]]:
    aliases: Dict[str, str] = {}
    alias_rows: List[Dict[str, object]] = []
    candidate_index = 0
    for index, model_set in enumerate(compact_sets, start=1):
        real_set_id = str(model_set.get("set_id", ""))
        source = str(model_set.get("source", ""))
        if task_mode == "correction" and source == "current":
            alias = "current_set"
        elif task_mode == "correction":
            candidate_index += 1
            alias = "candidate_set_{0}".format(_letters(candidate_index))
        else:
            alias = "S{0:02d}".format(index)
        aliases[real_set_id] = alias
        alias_rows.append(
            {
                "set_id": alias,
                "real_set_id": real_set_id,
                "real_source": source,
                "real_model_ids": model_set.get("model_ids") or [],
                "real_gene_ids": model_set.get("gene_ids") or [],
            }
        )
    return aliases, alias_rows


def build_model_aliases(compact_models: Sequence[Dict[str, object]], task_mode: str) -> Tuple[Dict[str, str], List[Dict[str, object]]]:
    aliases: Dict[str, str] = {}
    alias_rows: List[Dict[str, object]] = []
    current_index = 0
    candidate_index = 0
    for index, model in enumerate(compact_models, start=1):
        real_model_id = str(model.get("model_id", ""))
        source = str(model.get("source", ""))
        if task_mode == "correction" and source == "current":
            current_index += 1
            alias = "current_model_{0:02d}".format(current_index)
        elif task_mode == "correction":
            candidate_index += 1
            alias = "candidate_model_{0:02d}".format(candidate_index)
        else:
            alias = "M{0:02d}".format(index)
        aliases[real_model_id] = alias
        alias_rows.append(
            {
                "model_id": alias,
                "real_model_id": real_model_id,
                "real_source": source,
                "real_source_transcript_id": model.get("source_transcript_id", ""),
            }
        )
    return aliases, alias_rows


def _alias_model_id(value: object, model_aliases: Dict[str, str]) -> object:
    if isinstance(value, str):
        return model_aliases.get(value, value)
    return value


def filter_ai_visible_risk_tags(tags: Sequence[object]) -> List[str]:
    return sorted({str(tag) for tag in tags or [] if str(tag)})


def sanitize_public_tags(tags: Sequence[object], task_mode: str) -> List[str]:
    out: List[str] = []
    for tag in tags or []:
        value = str(tag)
        if value == "not_available_for_de_novo_model":
            continue
        if task_mode == "de_novo_annotation" and "current" in value.lower():
            continue
        out.append(value)
    return sorted(set(out))




def public_long_read_set_summary(summary: Dict[str, object]) -> Dict[str, object]:
    if not isinstance(summary, dict):
        return {"status": "not_configured"}
    out = {
        "status": summary.get("status", "not_configured"),
        "exact_supported_transcripts": summary.get("exact_supported_transcripts", 0),
        "partial_supported_transcripts": summary.get("partial_supported_transcripts", 0),
        "unsupported_transcripts": summary.get("unsupported_transcripts", 0),
        "mono_exon_not_assessable_count": summary.get("mono_exon_not_assessable_count", 0),
        "conflicting_chain_count": summary.get("conflicting_chain_count", 0),
        "support_level_counts": summary.get("support_level_counts") or {},
        "model_support_class_counts": summary.get("model_support_class_counts") or {},
        "risk_tags": summary.get("risk_tags") or [],
    }
    return out


def public_model_for_task(model: Dict[str, object], task_mode: str, model_aliases: Dict[str, str]) -> Dict[str, object]:
    real_model_id = str(model.get("model_id", ""))
    source = str(model.get("source", ""))
    out = dict(model)
    out["model_id"] = model_aliases.get(real_model_id, real_model_id)
    out["model_role"] = "current_model" if task_mode == "correction" and source == "current" else "candidate_model"
    for key in ["source", "source_transcript_id", "span", "protein_best_target"]:
        out.pop(key, None)
    if task_mode == "de_novo_annotation":
        out.pop("difference_from_current", None)
    out["risk_tags"] = sanitize_public_tags(out.get("risk_tags") or [], task_mode)
    out["hard_flags"] = sanitize_public_tags(out.get("hard_flags") or [], task_mode)
    for key in ["model_level_unsupported_junction_preview", "model_level_nearby_novel_junction_preview"]:
        out.pop(key, None)
    lr_model = out.get("long_read_model_support")
    if isinstance(lr_model, dict):
        lr_model = dict(lr_model)
        lr_model.pop("best_chain_id", None)
        lr_model["risk_tags"] = sanitize_public_tags(lr_model.get("risk_tags") or [], task_mode)
        out["long_read_model_support"] = lr_model
    diff = out.get("difference_from_current")
    if isinstance(diff, dict):
        diff = dict(diff)
        diff.pop("span_delta_start", None)
        diff.pop("span_delta_end", None)
        diff["closest_current_model_id"] = _alias_model_id(diff.get("closest_current_model_id", ""), model_aliases)
        out["difference_from_current"] = diff
    return out


def public_set_for_task(model_set: Dict[str, object], task_mode: str, set_aliases: Dict[str, str], model_aliases: Dict[str, str], set_index: int) -> Dict[str, object]:
    real_set_id = str(model_set.get("set_id", ""))
    source = str(model_set.get("source", ""))
    out = dict(model_set)
    out["set_id"] = set_aliases.get(real_set_id, real_set_id)
    out["model_ids"] = [model_aliases.get(str(model_id), str(model_id)) for model_id in (model_set.get("model_ids") or [])]
    supporting_sources = [str(item) for item in (out.get("supporting_sources") or []) if str(item)]
    equivalent_source_sets = [str(item) for item in (out.get("equivalent_source_set_ids") or []) if str(item)]
    out["supporting_source_count"] = len(supporting_sources)
    out["supporting_sources"] = ["source_{0:02d}".format(i) for i in range(1, len(supporting_sources) + 1)]
    out["equivalent_source_set_count"] = len(equivalent_source_sets)
    out["equivalent_source_set_ids"] = ["equivalent_source_set_{0:02d}".format(i) for i in range(1, len(equivalent_source_sets) + 1)]
    out["representative_source"] = "source_01" if supporting_sources else ""
    out.pop("structure_signature_key", None)
    for key in ["source", "span", "gene_ids"]:
        out.pop(key, None)
    if task_mode == "de_novo_annotation":
        out["set_role"] = "candidate_gene_set"
        out["candidate_label"] = "candidate_{0:02d}".format(set_index)
    elif source != "current":
        out["set_role"] = "candidate_gene_set"
        out["candidate_label"] = out["set_id"].replace("candidate_set_", "candidate_")
    else:
        out["set_role"] = "current_gene_set"
        out["candidate_label"] = "current"
    out["risk_tags"] = sanitize_public_tags(out.get("risk_tags") or [], task_mode)
    out["hard_flags"] = sanitize_public_tags(out.get("hard_flags") or [], task_mode)
    out["long_read_set_support_summary"] = public_long_read_set_summary(out.get("long_read_set_support_summary") or {})
    if isinstance(out.get("evidence_summary"), dict):
        evidence = dict(out["evidence_summary"])
        evidence["long_read_transcriptome"] = public_long_read_set_summary(evidence.get("long_read_transcriptome") or {})
        out["evidence_summary"] = evidence
    return out


def task_policy_tags(task_mode: str) -> List[str]:
    common = [
        "choose_listed_set_ids_only",
        "choose_listed_ids_only",
        "no_new_coordinates",
        "risk_flagged_not_ordinary_ai_selection",
        "validation_fail_not_auto_apply",
        "model_only_arbitration",
        "rna_is_model_evidence_not_coordinate_source",
    ]
    if task_mode == "correction":
        return common + [
            "current_annotation_is_candidate_not_truth",
            "first_check_current_biological_plausibility",
            "replace_current_only_with_better_supported_candidate",
            "rna_none_observed_not_negative_for_current",
            "de_novo_alone_not_enough_to_override_current",
        ]
    return common + [
        "all_candidates_are_anonymous",
        "no_base_annotation_defined",
        "choose_best_supported_candidate_set",
    ]


def decision_values_for_task(task_mode: str) -> List[str]:
    if task_mode == "de_novo_annotation":
        return ["select_candidate_set"]
    return ["keep_current_set", "choose_candidate_set", "propose_delete_current_set"]


def candidate_set_summary(compact_sets: Sequence[Dict[str, object]]) -> Dict[str, object]:
    current_sets = [
        str(row.get("set_id", ""))
        for row in compact_sets
        if row.get("source") == "current" or row.get("set_role") == "current_gene_set"
    ]
    validation_clean = [
        str(row.get("set_id", ""))
        for row in compact_sets
        if int((row.get("validation_summary") or {}).get("fail") or 0) == 0
    ]
    rna_supported = [
        str(row.get("set_id", ""))
        for row in compact_sets
        if int((row.get("rna_junction_summary") or {}).get("full_junction_support") or 0)
        or int((row.get("rna_junction_summary") or {}).get("partial_junction_support") or 0)
    ]
    protein_supported = [
        str(row.get("set_id", ""))
        for row in compact_sets
        if int(row.get("protein_supported_model_count") or 0) > 0
    ]
    return {
        "current_set_ids": current_sets,
        "validation_clean_set_ids": validation_clean,
        "rna_supported_set_ids": rna_supported,
        "protein_supported_set_ids": protein_supported,
    }



def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _fraction(count: int, total: int) -> float:
    return count / max(1, total)


def _model_structure_lookup(compact_models: Sequence[Dict[str, object]]) -> Dict[str, Dict[str, object]]:
    return {str(model.get("model_id", "")): model for model in compact_models if model.get("model_id")}


def _structure_group_counts(compact_sets: Sequence[Dict[str, object]], compact_models: Sequence[Dict[str, object]]) -> Dict[str, int]:
    by_model = _model_structure_lookup(compact_models)
    counts: Dict[str, int] = {}
    for model_set in compact_sets:
        groups = {
            str((by_model.get(str(model_id), {}) or {}).get("structure_group_id", ""))
            for model_id in (model_set.get("model_ids") or [])
        }
        for group in groups:
            if group:
                counts[group] = counts.get(group, 0) + 1
    return counts


def candidate_score_for_set(
    model_set: Dict[str, object],
    compact_models: Sequence[Dict[str, object]],
    structure_counts: Dict[str, int],
) -> Dict[str, object]:
    """Build a deterministic 0-100 evidence-support score for prompt context.

    The score is not benchmark accuracy and is not a decision by itself. It
    summarizes listed evidence channels so the LLM can compare candidate gene
    sets without re-interpreting raw support-count fields.
    """
    transcript_count = max(1, _int_value(model_set.get("transcript_count")))
    validation = model_set.get("validation_summary") or {}
    rna = model_set.get("rna_junction_summary") or {}
    long_read = model_set.get("long_read_set_support_summary") or {}
    hard_flags = {str(flag) for flag in model_set.get("hard_flags") or []}
    raw_risk_tags = [str(tag) for tag in model_set.get("risk_tags") or []]

    supportive: List[str] = []
    not_assessable: List[str] = []
    risk_flags: List[str] = list(raw_risk_tags)

    pass_count = _int_value(validation.get("pass"))
    warning_count = _int_value(validation.get("warning"))
    fail_count = _int_value(validation.get("fail"))
    not_assessed_count = _int_value(validation.get("not_assessed"))
    pass_fraction = _fraction(pass_count, transcript_count)
    if fail_count:
        coding_score = max(0.0, 8.0 - 8.0 * _fraction(fail_count, transcript_count))
        risk_flags.append("set_contains_validation_fail")
    elif pass_fraction >= 0.95:
        coding_score = 25.0
        supportive.append("passes_basic_coding_qc")
    elif pass_fraction >= 0.6:
        coding_score = 17.0
        supportive.append("most_models_pass_basic_coding_qc")
    elif pass_count:
        coding_score = 8.0
        risk_flags.append("limited_coding_qc_support")
    else:
        coding_score = 0.0
        not_assessable.append("coding_qc_not_assessed")
    if warning_count:
        risk_flags.append("set_contains_validation_warning")
    if not_assessed_count:
        risk_flags.append("set_contains_not_assessed_models")

    full_count = _int_value(rna.get("full_junction_support"))
    partial_count = _int_value(rna.get("partial_junction_support"))
    conflict_count = _int_value(rna.get("conflicting_junction"))
    none_count = _int_value(rna.get("none_observed"))
    not_assessable_count = _int_value(rna.get("not_assessable"))
    full_fraction = _fraction(full_count, transcript_count)
    supported_fraction = _fraction(full_count + partial_count, transcript_count)
    if full_fraction >= 0.95:
        short_read_score = 25.0
        supportive.append("full_short_read_junction_support")
    elif supported_fraction >= 0.8:
        short_read_score = 18.0
        supportive.append("mostly_short_read_junction_supported")
        risk_flags.append("partial_short_read_junction_support")
    elif supported_fraction > 0:
        short_read_score = 9.0
        supportive.append("limited_short_read_junction_support")
        risk_flags.append("limited_short_read_junction_support")
    elif not_assessable_count and not none_count:
        short_read_score = 5.0
        not_assessable.append("short_read_junction_not_assessable")
    else:
        short_read_score = 0.0
        if none_count:
            risk_flags.append("short_read_junction_none_observed")
    if conflict_count:
        risk_flags.append("conflicting_short_read_junction")

    lr_status = str(long_read.get("status", "not_configured"))
    if lr_status == "available":
        lr_exact = _int_value(long_read.get("exact_supported_transcripts"))
        lr_partial = _int_value(long_read.get("partial_supported_transcripts"))
        lr_unsupported = _int_value(long_read.get("unsupported_transcripts"))
        lr_mono = _int_value(long_read.get("mono_exon_not_assessable_count"))
        lr_conflict = _int_value(long_read.get("conflicting_chain_count"))
        exact_fraction = _fraction(lr_exact, transcript_count)
        lr_supported_fraction = _fraction(lr_exact + lr_partial, transcript_count)
        if exact_fraction >= 0.8:
            long_read_score = 15.0
            supportive.append("long_read_exact_chain_support")
        elif lr_exact:
            long_read_score = 10.0
            supportive.append("long_read_exact_chain_support_partial_set")
        elif lr_supported_fraction >= 0.8:
            long_read_score = 10.0
            supportive.append("long_read_junction_chain_support")
        elif lr_partial:
            long_read_score = 5.0
            supportive.append("long_read_partial_junction_support")
        elif lr_mono and not lr_unsupported:
            long_read_score = 3.0
            not_assessable.append("long_read_mono_exon_not_assessable")
        else:
            long_read_score = 0.0
            risk_flags.append("long_read_support_absent")
        if lr_unsupported:
            risk_flags.append("long_read_unsupported_transcripts_present")
        if lr_conflict:
            risk_flags.append("long_read_conflicting_chains")
        if lr_mono:
            not_assessable.append("long_read_mono_exon_not_assessable")
    else:
        long_read_score = 0.0
        not_assessable.append("long_read_not_configured" if lr_status == "not_configured" else "long_read_" + lr_status)

    protein_supported = _int_value(model_set.get("protein_supported_model_count"))
    protein_fraction = _fraction(protein_supported, transcript_count)
    if not run_level_evidence_active("protein"):
        protein_score = 0.0
        not_assessable.append("homolog_protein_not_configured")
    elif protein_fraction >= 0.8:
        protein_score = 20.0
        supportive.append("homolog_protein_support")
    elif protein_supported:
        protein_score = 10.0
        supportive.append("partial_homolog_protein_support")
    else:
        protein_score = 0.0
        risk_flags.append("homolog_protein_not_observed")

    structure_score = 0.0
    if "CDS_present" in hard_flags:
        structure_score += 4.0
        supportive.append("CDS_present")
    if "all_introns_junction_supported" in hard_flags:
        structure_score += 4.0
        supportive.append("all_introns_junction_supported")
    if "model_level_full_intron_chain_supported" in hard_flags:
        structure_score += 2.0
        supportive.append("model_level_full_intron_chain_supported")
    structure_score = min(10.0, structure_score)

    by_model = _model_structure_lookup(compact_models)
    groups = {
        str((by_model.get(str(model_id), {}) or {}).get("structure_group_id", ""))
        for model_id in (model_set.get("model_ids") or [])
    }
    shared_groups = sorted(group for group in groups if group and structure_counts.get(group, 0) > 1)
    cross_source_score = 0.0

    risk_penalty = 0.0
    if fail_count:
        risk_penalty -= min(25.0, 12.0 + 5.0 * fail_count)
    if warning_count:
        risk_penalty -= min(8.0, 3.0 * warning_count)
    if conflict_count:
        risk_penalty -= min(15.0, 5.0 * conflict_count)
    lr_conflict_value = _int_value(long_read.get("conflicting_chain_count")) if isinstance(long_read, dict) else 0
    if lr_conflict_value:
        risk_penalty -= 5.0
    risk_penalty -= min(10.0, 2.0 * len(set(raw_risk_tags)))

    component_scores = {
        "coding_qc": coding_score,
        "short_read": short_read_score,
        "long_read": long_read_score,
        "protein": protein_score,
        "structure": structure_score,
    }
    raw_positive_score = sum(component_scores.values())
    raw_total = _clamp(raw_positive_score + risk_penalty, 0.0, 100.0)

    active_components = {"coding_qc", "structure"}
    active_components.update(kind for kind in RUN_LEVEL_EVIDENCE_TYPES if run_level_evidence_active(kind))
    active_weight = sum(SCORING_COMPONENT_WEIGHTS[kind] for kind in active_components)
    active_positive_score = sum(component_scores[kind] for kind in active_components)
    if SCORING_NORMALIZATION == "run_level" and active_weight > 0:
        positive_scale = SCORING_TOTAL_COMPONENT_WEIGHT / active_weight
        normalized_total = _clamp(active_positive_score * positive_scale + risk_penalty, 0.0, 100.0)
    else:
        positive_scale = 1.0
        normalized_total = raw_total
    missing_evidence = [kind for kind in RUN_LEVEL_EVIDENCE_TYPES if kind not in ACTIVE_RUN_LEVEL_EVIDENCE_TYPES]
    rank_score = normalized_total if SCORING_NORMALIZATION == "run_level" else raw_total

    score_breakdown = {
        "coding_qc": round(coding_score, 2),
        "short_read_junction_support": round(short_read_score, 2),
        "long_read_transcriptome_support": round(long_read_score, 2),
        "homolog_protein_support": round(protein_score, 2),
        "structure_plausibility": round(structure_score, 2),
        "candidate_generation_agreement_group_count": len(shared_groups),
        "cross_candidate_consistency": round(cross_source_score, 2),
        "risk_penalty": round(risk_penalty, 2),
        "raw_total": round(raw_total, 2),
        "normalized_total": round(normalized_total, 2),
        "total": round(rank_score, 2),
    }
    return {
        "candidate_absolute_score": round(raw_total, 2),
        "candidate_normalized_score": round(normalized_total, 2),
        "candidate_rank_score": round(rank_score, 2),
        "scoring_normalization": SCORING_NORMALIZATION,
        "active_run_level_evidence_types": sorted(ACTIVE_RUN_LEVEL_EVIDENCE_TYPES),
        "missing_run_level_evidence_types": missing_evidence,
        "active_scoring_weight": round(active_weight, 2),
        "total_scoring_weight": round(SCORING_TOTAL_COMPONENT_WEIGHT, 2),
        "positive_score_scale": round(positive_scale, 6),
        "score_breakdown": score_breakdown,
        "supportive_evidence": sorted(set(supportive)),
        "not_assessable_evidence": sorted(set(not_assessable)),
        "score_risk_flags": sorted(set(risk_flags)),
        "shared_structure_group_ids": shared_groups,
    }


def _support_tier(score: object) -> str:
    value = _float_or_none(score) or 0.0
    if value >= 75.0:
        return "high"
    if value >= 50.0:
        return "moderate"
    if value > 0.0:
        return "low"
    return "unsupported_or_not_assessed"


def _risk_level(risk_flags: Sequence[object]) -> str:
    risks = {str(flag) for flag in risk_flags or []}
    high_tokens = ("validation_fail", "conflicting_short_read_junction", "long_read_conflicting", "homolog_protein_not_observed")
    moderate_tokens = ("warning", "partial", "limited", "not_assessed", "not_configured")
    if any(any(token in risk for token in high_tokens) for risk in risks):
        return "high"
    if any(any(token in risk for token in moderate_tokens) for risk in risks):
        return "low_to_moderate"
    return "low"


def _support_profile_from_set(model_set: Dict[str, object], score: Dict[str, object]) -> Dict[str, object]:
    validation = model_set.get("validation_summary") or {}
    rna = model_set.get("rna_junction_summary") or {}
    long_read = model_set.get("long_read_set_support_summary") or {}
    transcript_count = max(1, _int_value(model_set.get("transcript_count")))

    fail_count = _int_value(validation.get("fail"))
    warning_count = _int_value(validation.get("warning"))
    pass_count = _int_value(validation.get("pass"))
    if fail_count:
        coding_qc = "fail"
    elif warning_count:
        coding_qc = "warning"
    elif pass_count / transcript_count >= 0.95:
        coding_qc = "pass"
    elif pass_count:
        coding_qc = "partial_pass"
    else:
        coding_qc = "not_assessed"

    if not run_level_evidence_active("short_read"):
        short_read = "not_configured"
    else:
        full_count = _int_value(rna.get("full_junction_support"))
        partial_count = _int_value(rna.get("partial_junction_support"))
        conflict_count = _int_value(rna.get("conflicting_junction"))
        none_count = _int_value(rna.get("none_observed"))
        not_assessable_count = _int_value(rna.get("not_assessable"))
        if conflict_count:
            short_read = "conflicting_support"
        elif full_count / transcript_count >= 0.95:
            short_read = "full_support"
        elif (full_count + partial_count) / transcript_count >= 0.8:
            short_read = "mostly_supported"
        elif full_count + partial_count:
            short_read = "partial_support"
        elif not_assessable_count and not none_count:
            short_read = "not_assessable"
        elif none_count:
            short_read = "none_observed"
        else:
            short_read = "not_configured"

    lr_status = str(long_read.get("status", "not_configured"))
    if not run_level_evidence_active("long_read"):
        long_read_profile = "not_configured"
    elif lr_status == "available":
        lr_exact = _int_value(long_read.get("exact_supported_transcripts"))
        lr_partial = _int_value(long_read.get("partial_supported_transcripts"))
        lr_unsupported = _int_value(long_read.get("unsupported_transcripts"))
        lr_mono = _int_value(long_read.get("mono_exon_not_assessable_count"))
        lr_conflict = _int_value(long_read.get("conflicting_chain_count"))
        if lr_conflict:
            long_read_profile = "conflicting_chains_present"
        elif lr_exact / transcript_count >= 0.8:
            long_read_profile = "exact_chain_support"
        elif (lr_exact + lr_partial) / transcript_count >= 0.8:
            long_read_profile = "mostly_supported"
        elif lr_exact + lr_partial:
            long_read_profile = "partial_support"
        elif lr_mono and not lr_unsupported:
            long_read_profile = "not_assessable_mono_exon"
        elif lr_unsupported:
            long_read_profile = "unsupported"
        else:
            long_read_profile = "available_no_support"
    else:
        long_read_profile = lr_status

    protein_supported = _int_value(model_set.get("protein_supported_model_count"))
    if not run_level_evidence_active("protein"):
        protein = "not_configured"
    elif protein_supported / transcript_count >= 0.8:
        protein = "supported"
    elif protein_supported:
        protein = "partially_supported"
    else:
        protein = "not_observed"

    shared_groups = score.get("shared_structure_group_ids") or []
    return {
        "coding_qc": coding_qc,
        "short_read_junction_support": short_read,
        "long_read_transcriptome_support": long_read_profile,
        "protein_support": protein,
        "candidate_generation_agreement": "present" if shared_groups else "not_detected",
        "cross_candidate_consistency": "candidate_generation_agreement_only" if shared_groups else "not_detected",
        "risk_level": _risk_level(score.get("score_risk_flags") or []),
    }


def add_interpreted_scores(
    compact_sets: Sequence[Dict[str, object]],
    compact_models: Sequence[Dict[str, object]],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    structure_counts = _structure_group_counts(compact_sets, compact_models)
    public_sets: List[Dict[str, object]] = []
    score_trace: List[Dict[str, object]] = []
    for model_set in compact_sets:
        score = candidate_score_for_set(model_set, compact_models, structure_counts)
        public = dict(model_set)
        public["support_tier"] = _support_tier(score["candidate_rank_score"])
        public["support_profile"] = _support_profile_from_set(model_set, score)
        public["supportive_evidence"] = score["supportive_evidence"]
        public["not_assessable_evidence"] = score["not_assessable_evidence"]
        public["interpreted_risk_flags"] = filter_ai_visible_risk_tags(score["score_risk_flags"])
        public["shared_structure_group_ids"] = score["shared_structure_group_ids"]
        public_sets.append(public)
        score_trace.append(
            {
                "set_id": public.get("set_id", ""),
                "set_role": public.get("set_role", ""),
                "candidate_absolute_score": score["candidate_absolute_score"],
                "candidate_normalized_score": score["candidate_normalized_score"],
                "candidate_rank_score": score["candidate_rank_score"],
                "scoring_normalization": score["scoring_normalization"],
                "active_run_level_evidence_types": score["active_run_level_evidence_types"],
                "missing_run_level_evidence_types": score["missing_run_level_evidence_types"],
                "active_scoring_weight": score["active_scoring_weight"],
                "total_scoring_weight": score["total_scoring_weight"],
                "positive_score_scale": score["positive_score_scale"],
                "score_breakdown": score["score_breakdown"],
                "support_tier": public["support_tier"],
                "support_profile": public["support_profile"],
                "supportive_evidence": score["supportive_evidence"],
                "not_assessable_evidence": score["not_assessable_evidence"],
                "risk_flags": filter_ai_visible_risk_tags(score["score_risk_flags"]),
                "shared_structure_group_ids": score["shared_structure_group_ids"],
            }
        )
    ranked = sorted(score_trace, key=lambda row: (-float(row.get("candidate_rank_score") or 0), str(row.get("set_id", ""))))
    previous_score = None
    previous_rank = 0
    for index, row in enumerate(ranked, start=1):
        score = float(row.get("candidate_rank_score") or 0)
        if previous_score is None or score < previous_score:
            previous_rank = index
            previous_score = score
        row["rank"] = previous_rank
    return public_sets, ranked


def _set_structure_keys(model_set: Dict[str, object], compact_models: Sequence[Dict[str, object]]) -> List[str]:
    by_model = _model_structure_lookup(compact_models)
    keys = sorted(
        {
            str((by_model.get(str(model_id), {}) or {}).get("structure_group_id", ""))
            for model_id in (model_set.get("model_ids") or [])
            if str((by_model.get(str(model_id), {}) or {}).get("structure_group_id", ""))
        }
    )
    return keys or [str(model_set.get("set_id", ""))]


def _support_tier_rank(value: object) -> int:
    return {"high": 3, "moderate": 2, "low": 1, "unsupported_or_not_assessed": 0}.get(str(value), 0)



def _model_component_support_summary(model: Dict[str, object]) -> List[str]:
    summary: List[str] = []
    if str(model.get("orf_status", "")) == "complete_ORF":
        summary.append("complete_ORF")
    if str(model.get("phase_status", "")) == "phase_consistent":
        summary.append("phase_consistent")
    if str(model.get("splice_motif_status", "")) == "canonical_splice_site":
        summary.append("canonical_splice_site")
    rna_status = str(model.get("rna_junction_status", ""))
    if rna_status == "full_junction_support":
        summary.append("full_short_read_junction_support")
    elif rna_status == "partial_junction_support":
        summary.append("partial_short_read_junction_support")
    elif rna_status == "conflicting_junction":
        summary.append("conflicting_short_read_junction")
    elif rna_status in {"not_assessable", "none_observed"}:
        summary.append("short_read_junction_" + rna_status)
    if str(model.get("protein_support_status", "")) == "protein_overlap_present":
        summary.append("homolog_protein_support")
    long_read = model.get("long_read_model_support") or {}
    if isinstance(long_read, dict):
        support_class = str(long_read.get("support_class") or long_read.get("status") or "")
        if support_class in {"exact_chain_supported", "single_read_exact_chain"}:
            summary.append("long_read_exact_chain_support")
        elif support_class in {"partial_junction_supported", "all_junctions_supported_no_exact_chain"}:
            summary.append("long_read_junction_support")
        elif support_class == "mono_exon_not_assessable":
            summary.append("long_read_mono_exon_not_assessable")
    return sorted(set(summary))


def _component_type_for_models(models: Sequence[Dict[str, object]], member_set_count: int, total_set_count: int) -> str:
    intron_counts = [_int_value((model.get("structure_summary") or {}).get("intron_count")) for model in models]
    if intron_counts and all(count == 0 for count in intron_counts):
        if member_set_count < total_set_count:
            return "possible_extra_or_mono_exon_component"
        return "mono_exon_protein_coding_component"
    return "protein_coding_transcript_component"


def _component_claim(component_type: str, member_set_count: int, total_set_count: int) -> str:
    if component_type == "possible_extra_or_mono_exon_component":
        return "possible extra mono-exon gene or candidate-specific coding unit"
    if member_set_count == total_set_count:
        return "shared supported protein-coding transcript component"
    if member_set_count > 1:
        return "supported coding component shared by a subset of candidate sets"
    return "candidate-specific protein-coding component"


def _component_role(component_type: str, support_tier: str, member_set_count: int, total_set_count: int) -> str:
    if component_type == "possible_extra_or_mono_exon_component":
        return "uncertain_extra_component"
    if member_set_count == total_set_count:
        return "shared_core_component"
    if support_tier == "high":
        return "strong_component"
    return "candidate_specific_component"


def _component_uncertainty(component_type: str, member_set_count: int, total_set_count: int) -> str:
    if component_type == "possible_extra_or_mono_exon_component":
        return "mono_exon_or_extra_gene_status"
    if member_set_count == total_set_count:
        return "terminal_boundary_or_transcript_span"
    return "split_merge_or_component_presence"



def derived_component_fit_summary(
    component_mapping: Dict[str, str],
    component_by_id: Dict[str, Dict[str, object]],
) -> Dict[str, List[str]]:
    strong_components = [
        component_id
        for component_id, component in component_by_id.items()
        if component.get("component_role_in_decision") in {"shared_core_component", "strong_component"}
    ]
    uncertain_components = [
        component_id
        for component_id, component in component_by_id.items()
        if component.get("component_role_in_decision") == "uncertain_extra_component"
    ]
    unsupported_components = [
        component_id
        for component_id, component in component_by_id.items()
        if component.get("support_tier") in {"low", "unsupported_or_not_assessed"}
        and component.get("component_role_in_decision") not in {"shared_core_component", "strong_component", "uncertain_extra_component"}
    ]
    covered = {component_id for component_id, state in component_mapping.items() if state != "not_included"}
    return {
        "strong_components_covered": sorted(component_id for component_id in strong_components if component_id in covered),
        "strong_components_missing": sorted(component_id for component_id in strong_components if component_id not in covered),
        "uncertain_components_included": sorted(component_id for component_id in uncertain_components if component_id in covered),
        "unsupported_components_included": sorted(component_id for component_id in unsupported_components if component_id in covered),
    }


def derived_topology_risk_summary(
    topology_claim: str,
    component_mapping: Dict[str, str],
    component_fit_summary: Dict[str, List[str]],
    set_level_risks: Sequence[object],
) -> Dict[str, str]:
    risks = {str(item) for item in set_level_risks or []}
    mapping_values = {str(value) for value in component_mapping.values()}
    fusion = "possible" if "possible_fusion_if_components_are_independent" in risks or "covered_within_same_gene" in mapping_values else "none"
    split = "possible" if "possible_over_split_if_transcriptional_bridge_supports_fusion" in risks or "independent_genes" in topology_claim else "none"
    extra = "possible" if "may_accept_uncertain_extra_gene" in risks or component_fit_summary.get("uncertain_components_included") else "none"
    missing = "possible" if component_fit_summary.get("strong_components_missing") else "none"
    return {
        "fusion_risk": fusion,
        "split_risk": split,
        "extra_gene_risk": extra,
        "missing_component_risk": missing,
    }


def validate_candidate_set_view_derivations(view: Dict[str, object], component_by_id: Dict[str, Dict[str, object]]) -> None:
    component_mapping = view.get("component_mapping") or {}
    if not isinstance(component_mapping, dict):
        raise ValueError("candidate_set_view component_mapping must be a dict for {0}".format(view.get("set_id", "")))
    set_level_risks = view.get("set_level_risks") or []
    expected_fit = derived_component_fit_summary(component_mapping, component_by_id)
    if view.get("component_fit_summary") != expected_fit:
        raise ValueError(
            "component_fit_summary is not derived from component_mapping/component_role for {0}".format(view.get("set_id", ""))
        )
    expected_risk = derived_topology_risk_summary(
        str(view.get("topology_claim", "")),
        component_mapping,
        expected_fit,
        set_level_risks,
    )
    if view.get("topology_risk_summary") != expected_risk:
        raise ValueError(
            "topology_risk_summary is not derived from topology_claim/component_mapping/risk flags for {0}".format(view.get("set_id", ""))
        )


def _public_candidate_set_views(
    locus_topology_summary: Dict[str, object],
    internal_views: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    if locus_topology_summary.get("topology_class") == "single_shared_component":
        return []
    rows: List[Dict[str, object]] = []
    for view in internal_views:
        fit = view.get("component_fit_summary") or {}
        public_fit = {
            "strong_component_coverage": "all_covered" if not fit.get("strong_components_missing") else "missing_strong_components",
            "uncertain_extra_component": "included" if fit.get("uncertain_components_included") else "not_included",
            "unsupported_component": "included" if fit.get("unsupported_components_included") else "not_included",
        }
        row: Dict[str, object] = {
            "set_id": view.get("set_id", ""),
            "gene_count": view.get("gene_count", 0),
            "transcript_count": view.get("transcript_count", 0),
            "component_fit_summary": public_fit,
        }
        set_level_risks = view.get("set_level_risks") or []
        if set_level_risks:
            row["set_level_risks"] = set_level_risks
        rows.append(row)
    return rows


def build_arbitration_component_layers(
    compact_sets: Sequence[Dict[str, object]],
    compact_models: Sequence[Dict[str, object]],
) -> Tuple[Dict[str, object], List[Dict[str, object]], List[Dict[str, object]]]:
    by_model = _model_structure_lookup(compact_models)
    group_to_model_ids: Dict[str, List[str]] = {}
    model_to_set_ids: Dict[str, List[str]] = {}
    for model_set in compact_sets:
        set_id = str(model_set.get("set_id", ""))
        for model_id in model_set.get("model_ids") or []:
            model_id = str(model_id)
            model = by_model.get(model_id, {})
            groups = [str(item) for item in (model.get("component_group_ids") or []) if str(item)] or [str(model.get("structure_group_id") or model_id)]
            for group in groups:
                group_to_model_ids.setdefault(group, []).append(model_id)
            model_to_set_ids.setdefault(model_id, []).append(set_id)

    component_id_by_group: Dict[str, str] = {}
    components: List[Dict[str, object]] = []
    total_set_count = len([row for row in compact_sets if row.get("set_id")])
    ordered_groups = sorted(group_to_model_ids)
    for index, group in enumerate(ordered_groups, start=1):
        component_id = "C{0:02d}".format(index)
        component_id_by_group[group] = component_id
        member_model_ids = sorted(set(group_to_model_ids[group]))
        member_models = [by_model[model_id] for model_id in member_model_ids if model_id in by_model]
        member_set_ids = sorted({set_id for model_id in member_model_ids for set_id in model_to_set_ids.get(model_id, [])})
        member_sets = [row for row in compact_sets if str(row.get("set_id", "")) in set(member_set_ids)]
        tiers = [str(row.get("support_tier", "unsupported_or_not_assessed")) for row in member_sets]
        support_tier = sorted(tiers, key=lambda tier: -_support_tier_rank(tier))[0] if tiers else "unsupported_or_not_assessed"
        component_type = _component_type_for_models(member_models, len(member_set_ids), total_set_count)
        support_summary = sorted({item for model in member_models for item in _model_component_support_summary(model)})
        risk_flags = sorted({item for row in member_sets for item in (row.get("interpreted_risk_flags") or [])})
        components.append(
            {
                "component_id": component_id,
                "component_type": component_type,
                "biological_claim": _component_claim(component_type, len(member_set_ids), total_set_count),
                "support_tier": "uncertain" if component_type == "possible_extra_or_mono_exon_component" and support_tier == "low" else support_tier,
                "support_summary": support_summary,
                "main_uncertainty": _component_uncertainty(component_type, len(member_set_ids), total_set_count),
                "component_role_in_decision": _component_role(component_type, support_tier, len(member_set_ids), total_set_count),
                "component_risk_flags": risk_flags,
            }
        )

    all_component_ids = [row["component_id"] for row in components]
    component_by_id = {str(row["component_id"]): row for row in components}
    views: List[Dict[str, object]] = []
    topology_claims: List[str] = []
    for model_set in compact_sets:
        set_id = str(model_set.get("set_id", ""))
        set_groups = _set_structure_keys(model_set, compact_models)
        included = [component_id_by_group.get(group, group) for group in set_groups]
        included = [component_id for component_id in included if component_id in component_by_id]
        gene_count = _int_value(model_set.get("gene_count"))
        mapping: Dict[str, str] = {}
        includes_uncertain_extra = False
        for component_id in all_component_ids:
            component = component_by_id[component_id]
            if component_id not in included:
                mapping[component_id] = "not_included"
                continue
            if component.get("component_role_in_decision") == "uncertain_extra_component":
                mapping[component_id] = "covered_as_extra_mono_exon_gene"
                includes_uncertain_extra = True
            elif gene_count <= 1 and len(included) > 1:
                mapping[component_id] = "covered_within_same_gene"
            elif len(included) == 1:
                mapping[component_id] = "covered_as_single_gene"
            else:
                mapping[component_id] = "covered_as_independent_gene"

        if not included:
            topology_claim = "no_component_mapping"
        elif gene_count <= 1 and len(included) == 1:
            topology_claim = "single_gene_covering_{0}".format(included[0])
        elif gene_count <= 1 and len(included) > 1:
            topology_claim = "single_gene_covering_{0}".format("_and_".join(included))
        elif gene_count >= len(included) and len(included) > 1:
            topology_claim = "split_into_{0}_independent_genes".format(len(included))
        else:
            topology_claim = "multi_gene_component_arrangement"
        if includes_uncertain_extra:
            topology_claim += "_plus_uncertain_extra_component"
        topology_claims.append(topology_claim)

        risks: List[str] = []
        if gene_count <= 1 and len(included) > 1:
            risks.append("possible_fusion_if_components_are_independent")
        if gene_count >= len(included) and len(included) > 1:
            risks.append("possible_over_split_if_transcriptional_bridge_supports_fusion")
        if includes_uncertain_extra:
            risks.append("may_accept_uncertain_extra_gene")
            risks.append("mono_exon_component_not_junction_assessable")
        risks.extend(str(item) for item in model_set.get("interpreted_risk_flags") or [])
        set_level_risks = sorted(set(risks))
        component_fit_summary = derived_component_fit_summary(mapping, component_by_id)
        topology_risk_summary = derived_topology_risk_summary(topology_claim, mapping, component_fit_summary, set_level_risks)
        internal_view = {
            "set_id": set_id,
            "gene_count": model_set.get("gene_count", 0),
            "transcript_count": model_set.get("transcript_count", 0),
            "topology_claim": topology_claim,
            "component_mapping": mapping,
            "component_fit_summary": component_fit_summary,
            "topology_risk_summary": topology_risk_summary,
            "set_level_risks": set_level_risks,
        }
        validate_candidate_set_view_derivations(internal_view, component_by_id)
        views.append(internal_view)

    component_count = len(components)
    uncertain_count = sum(1 for row in components if row.get("component_role_in_decision") == "uncertain_extra_component")
    has_fusion_risk = any("possible_fusion_if_components_are_independent" in (row.get("set_level_risks") or []) for row in views)
    has_split_risk = any("possible_over_split_if_transcriptional_bridge_supports_fusion" in (row.get("set_level_risks") or []) for row in views)
    if component_count <= 1:
        topology_class = "single_shared_component"
        main_uncertainty = "terminal_boundary_or_transcript_span"
    elif has_fusion_risk and has_split_risk:
        topology_class = "split_merge_conflict"
        main_uncertainty = "whether_components_are_independent_genes_or_one_merged_gene"
    elif uncertain_count:
        topology_class = "extra_component_conflict"
        main_uncertainty = "whether_uncertain_extra_components_should_be_accepted"
    else:
        topology_class = "multi_component_candidate_locus"
        main_uncertainty = "component_topology"
    summary = {
        "topology_class": topology_class,
        "component_count": component_count,
        "candidate_set_count": total_set_count,
        "candidate_set_topology_pattern_count": len(set(topology_claims)),
        "has_possible_split_merge": bool(has_fusion_risk or has_split_risk),
        "has_uncertain_extra_component": bool(uncertain_count),
        "main_uncertainty": main_uncertainty,
    }
    return summary, components, views


def _candidate_equivalence_key(model_set: Dict[str, object], compact_models: Sequence[Dict[str, object]], collapse_terminal_boundary: bool) -> str:
    if not collapse_terminal_boundary:
        return "+".join(_set_structure_keys(model_set, compact_models))
    profile = model_set.get("support_profile") or {}
    support_key = (
        str(model_set.get("support_tier", "")),
        str(profile.get("coding_qc", "")),
        str(profile.get("short_read_junction_support", "")),
        str(profile.get("long_read_transcriptome_support", "")),
        str(profile.get("protein_support", "")),
        str(profile.get("risk_level", "")),
        ",".join(sorted(str(item) for item in (model_set.get("interpreted_risk_flags") or []))),
    )
    return "biological_equivalence:" + "|".join(support_key)


def candidate_structure_groups(
    compact_sets: Sequence[Dict[str, object]],
    compact_models: Sequence[Dict[str, object]],
    locus_topology_summary: Dict[str, object] | None = None,
) -> List[Dict[str, object]]:
    collapse_terminal_boundary = bool(
        locus_topology_summary
        and locus_topology_summary.get("topology_class") == "single_shared_component"
        and not locus_topology_summary.get("has_possible_split_merge")
        and not locus_topology_summary.get("has_uncertain_extra_component")
    )
    grouped: Dict[str, List[Dict[str, object]]] = {}
    for model_set in compact_sets:
        key = _candidate_equivalence_key(model_set, compact_models, collapse_terminal_boundary)
        grouped.setdefault(key, []).append(model_set)
    rows: List[Dict[str, object]] = []
    ordered_groups = sorted(
        grouped.items(),
        key=lambda item: (
            -max(_support_tier_rank(row.get("support_tier")) for row in item[1]),
            str(item[0]),
        ),
    )
    top_key = ordered_groups[0][0] if ordered_groups else ""
    for index, (key, members) in enumerate(ordered_groups, start=1):
        member_ids = sorted(str(row.get("set_id", "")) for row in members if row.get("set_id"))
        canonical = member_ids[0] if member_ids else ""
        tiers = [str(row.get("support_tier", "")) for row in members]
        group_tier = sorted(tiers, key=lambda tier: -_support_tier_rank(tier))[0] if tiers else "unsupported_or_not_assessed"
        supportive = sorted({item for row in members for item in (row.get("supportive_evidence") or [])})
        risks = sorted({item for row in members for item in (row.get("interpreted_risk_flags") or [])})
        equivalent_members = len(member_ids) > 1
        confidence_interpretation = {
            "group_support_confidence": group_tier,
            "specific_set_selection_confidence": "low" if equivalent_members else group_tier,
            "reason": "{0} is selected as canonical representative; available evidence does not show it is biologically superior to other member sets".format(canonical)
            if equivalent_members
            else "single-member group; selected set confidence follows group support confidence",
        }
        row: Dict[str, object] = {
            "group_id": "G{0:02d}".format(index),
            "member_set_ids": member_ids,
            "canonical_representative_set_id": canonical,
            "group_claim": "biologically_equivalent_supported_protein_coding_candidate_sets"
            if collapse_terminal_boundary and equivalent_members and key == top_key
            else ("top_supported_coding_intron_structure_group" if key == top_key else "similar_supported_structure_group"),
            "support_tier": group_tier,
            "supportive_evidence": supportive,
            "risk_flags": risks,
            "biological_equivalence": "member_sets_are_equivalent_under_available_evidence"
            if equivalent_members
            else "single_member_group",
            "shared_biological_component": "same_supported_protein_coding_component" if collapse_terminal_boundary else "same_candidate_component",
            "within_group_discriminator": "none_detected" if equivalent_members else "single_member_group",
            "evidence_not_discriminating_members": ["terminal_boundary_or_transcript_span"] if collapse_terminal_boundary and equivalent_members else [],
            "specific_set_selection_basis": "canonical_representative_only" if equivalent_members else "single_member_group",
            "confidence_interpretation": confidence_interpretation,
        }
        if key != top_key:
            row["difference_from_top_group"] = "alternative_structure_or_terminal_span_difference_without_clear_boundary_specific_support"
        rows.append(row)
    return rows


def between_candidate_evidence_interpretation(compact_sets: Sequence[Dict[str, object]]) -> Dict[str, object]:
    profiles = [row.get("support_profile") or {} for row in compact_sets]
    def varies(key: str) -> bool:
        return len({str(profile.get(key, "")) for profile in profiles}) > 1
    return {
        "coding_evidence_discriminates_between_candidates": varies("coding_qc"),
        "junction_evidence_discriminates_between_candidates": varies("short_read_junction_support") or varies("long_read_transcriptome_support"),
        "protein_evidence_discriminates_between_candidates": varies("protein_support"),
        "terminal_boundary_evidence_discriminates_between_candidates": False,
        "main_uncertainty": "terminal_boundary_or_transcript_span",
        "decision_meaning": "select a representative of the top supported equivalence group when listed evidence does not distinguish candidates",
    }


def _support_level_is_reliable(value: object) -> bool:
    text = str(value or "")
    return text.startswith("strong") or text.startswith("moderate") or text in {"present", "pass"}


def _current_model_has_intrinsic_risk(model: Dict[str, object], min_cds_length: int = 300) -> bool:
    validation = str(model.get("validation_status", ""))
    if validation in {"fail", "warning"}:
        return True
    if str(model.get("splice_motif_status", "")) == "noncanonical_splice_site":
        return True
    risk_tags = {str(tag) for tag in model.get("risk_tags") or []}
    risk_tokens = (
        "validation_warning",
        "validation_fail",
        "noncanonical_splice",
        "internal_stop",
        "TE_overlap_high",
        "high_TE_overlap",
        "weak_fragmented_structure",
        "very_short_ORF",
    )
    if any(any(token in tag for token in risk_tokens) for tag in risk_tags):
        return True
    cds_length = _int_value((model.get("structure_summary") or {}).get("cds_length"))
    return 0 < cds_length < min_cds_length


def _candidate_set_is_reliable_replacement(model_set: Dict[str, object]) -> bool:
    levels = model_set.get("evidence_levels") or {}
    validation = str(levels.get("coding_validation", ""))
    if validation in {"fail", "warning"}:
        return False
    return any(
        _support_level_is_reliable(levels.get(key))
        for key in ["short_read_junction", "long_read_transcriptome", "homolog_protein"]
    )


def deletion_evidence_summary(compact_sets: Sequence[Dict[str, object]], compact_models: Sequence[Dict[str, object]], task_mode: str, locus_id: str = "", set_alias_map: Sequence[Dict[str, object]] = (), nipponbare_support: Dict[str, object] | None = None, treat_missing_nipponbare_as_absent: bool = False) -> Dict[str, object]:
    """Summarize evidence for a correction-mode delete proposal.

    This is an LLM-visible summary and a deterministic audit input. It does not
    delete anything by itself; deletion acceptance is handled by the audit script.
    """
    if task_mode != "correction":
        return {"status": "not_applicable_for_de_novo"}
    current_sets = [row for row in compact_sets if row.get("set_role") == "current_gene_set"]
    if not current_sets:
        return {"status": "no_current_set"}
    current_set = current_sets[0]
    current_model_ids = {str(model_id) for model_id in current_set.get("model_ids") or []}
    current_models = [model for model in compact_models if str(model.get("model_id", "")) in current_model_ids]
    candidate_sets = [row for row in compact_sets if row.get("set_role") == "candidate_gene_set"]
    reliable_candidates = [row for row in candidate_sets if _candidate_set_is_reliable_replacement(row)]

    sr = current_set.get("short_read_junction_summary") or {}
    lr = current_set.get("long_read_set_support_summary") or {}
    short_read_supported = _int_value(sr.get("full_junction_support")) + _int_value(sr.get("partial_junction_support"))
    long_read_supported = _int_value(lr.get("exact_supported_transcripts")) + _int_value(lr.get("partial_supported_transcripts"))
    protein_supported = _int_value(current_set.get("protein_supported_model_count"))
    current_intrinsic_risk = any(_current_model_has_intrinsic_risk(model) for model in current_models)
    validation = current_set.get("validation_summary") or {}
    if _int_value(validation.get("fail")) or _int_value(validation.get("warning")):
        current_intrinsic_risk = True

    nipponbare_absent, nipponbare_status = summarize_nipponbare_support(
        locus_id,
        set_alias_map,
        nipponbare_support,
        treat_missing_nipponbare_as_absent,
    )
    criteria = {
        "no_reliable_replacement_candidate": not reliable_candidates,
        "mh63_rna_support_absent": short_read_supported == 0 and long_read_supported == 0,
        "mh63_homolog_protein_absent": protein_supported == 0,
        "nipponbare_transcript_support_absent": nipponbare_absent,
        "current_intrinsic_risk": current_intrinsic_risk,
    }
    return {
        "status": "available",
        "criteria": criteria,
        "candidate_replacement_count": len(candidate_sets),
        "reliable_replacement_candidate_count": len(reliable_candidates),
        "mh63_short_read_supported_transcripts": short_read_supported,
        "mh63_long_read_supported_transcripts": long_read_supported,
        "mh63_homolog_protein_supported_models": protein_supported,
        "nipponbare_transcript_support": nipponbare_status,
        "current_intrinsic_risk_tags": sorted({tag for model in current_models for tag in (model.get("risk_tags") or [])}),
        "interpretation": "Deletion is only for weak current sets with no reliable replacement and absent MH63/Nipponbare protection evidence; deterministic deletion audit decides whether deletion is accepted. Uncertain cases should be kept with flags rather than broadly routed to manual review.",
    }


def _support_word(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return "not_assessed"
    return text


def _validation_state(model_set: Dict[str, object]) -> str:
    validation = model_set.get("validation_summary") or {}
    fail = _int_value(validation.get("fail"))
    warning = _int_value(validation.get("warning"))
    passed = _int_value(validation.get("pass"))
    not_assessed = _int_value(validation.get("not_assessed"))
    if fail:
        return "fail"
    if warning:
        return "warning"
    if passed and not not_assessed:
        return "pass"
    if passed:
        return "partial_pass"
    return "not_assessed"


def _short_read_state(model_set: Dict[str, object]) -> str:
    if not run_level_evidence_active("short_read"):
        return "not_configured"
    rna = model_set.get("rna_junction_summary") or {}
    if _int_value(rna.get("full_junction_support")):
        return "full"
    if _int_value(rna.get("partial_junction_support")):
        return "partial"
    if _int_value(rna.get("conflicting_junction")):
        return "conflict"
    if _int_value(rna.get("not_assessable")) and not _int_value(rna.get("none_observed")):
        return "not_assessable"
    if _int_value(rna.get("none_observed")):
        return "not_observed"
    return "not_configured"


def _long_read_state(model_set: Dict[str, object]) -> str:
    if not run_level_evidence_active("long_read"):
        return "not_configured"
    lr = model_set.get("long_read_set_support_summary") or {}
    status = str(lr.get("status", "not_configured"))
    if status != "available":
        return "not_available"
    if _int_value(lr.get("exact_supported_transcripts")):
        return "exact"
    if _int_value(lr.get("partial_supported_transcripts")):
        return "partial"
    if _int_value(lr.get("conflicting_chain_count")):
        return "conflict"
    if _int_value(lr.get("mono_exon_not_assessable_count")) and not _int_value(lr.get("unsupported_transcripts")):
        return "not_assessable"
    return "not_observed"


def _protein_state(model_set: Dict[str, object]) -> str:
    if not run_level_evidence_active("protein"):
        return "not_configured"
    supported = _int_value(model_set.get("protein_supported_model_count"))
    transcripts = max(1, _int_value(model_set.get("transcript_count")))
    if supported >= transcripts:
        return "supported"
    if supported:
        return "partial"
    return "not_observed"


def _component_fit_label(view: Dict[str, object]) -> str:
    fit = view.get("component_fit_summary") or {}
    if fit.get("strong_components_missing"):
        return "misses_core_component"
    if fit.get("uncertain_components_included"):
        return "includes_uncertain_extra_component"
    if fit.get("unsupported_components_included"):
        return "includes_unsupported_component"
    mapping_values = {str(value) for value in (view.get("component_mapping") or {}).values()}
    if "covered_within_same_gene" in mapping_values:
        return "fusion_of_components"
    return "covers_core_components"


def _set_interpretation(model_set: Dict[str, object], view_by_set: Dict[str, Dict[str, object]]) -> str:
    set_id = str(model_set.get("set_id", ""))
    fit = _component_fit_label(view_by_set.get(set_id, {}))
    validation = _validation_state(model_set)
    short_read = _short_read_state(model_set)
    long_read = _long_read_state(model_set)
    protein = _protein_state(model_set)
    if validation == "fail":
        return "biologically_problematic"
    if fit in {"misses_core_component", "includes_unsupported_component"}:
        return "topology_risk"
    if short_read == "conflict" or long_read == "conflict":
        return "evidence_conflict"
    if short_read in {"full", "partial"} or long_read in {"exact", "partial"} or protein in {"supported", "partial"}:
        return "biologically_supported"
    return "plausible_low_evidence"


def _clean_tags(values: object, limit: int = 8) -> List[str]:
    if not isinstance(values, list):
        return []
    blocked_tokens = (
        "not_observed",
        "none_observed",
        "not_assessable",
        "not_configured",
        "not_computed",
        "support_absent",
        "no_exact_chain",
        "missing_junction_support",
        "unsupported_by_long_reads",
        "unsupported_transcripts_present",
        "no_splice_junction_to_assess",
        "mono_exon_no_chain",
        "partial_junction_supported",
        "partial_junction_support",
        "partial_short_read_junction_support",
        "limited_short_read_junction_support",
        "model_level_no_junction_support",
        "missing_for_current",
        "long_read_none",
        "long_read_not_available",
        "single_read_exact",
    )
    out = []
    for value in values:
        text = str(value)
        lower = text.lower()
        if not text or text in out:
            continue
        if any(token in lower for token in blocked_tokens):
            continue
        out.append(text)
        if len(out) >= limit:
            break
    return out


def _set_quality_rank(row: Dict[str, object]) -> int:
    interpretation = str(row.get("interpretation", ""))
    if interpretation == "biologically_supported":
        rank = 4
    elif interpretation == "evidence_conflict":
        rank = 2
    elif interpretation == "plausible_low_evidence":
        rank = 1
    else:
        rank = 0
    priority = str(row.get("candidate_priority", "normal"))
    if priority.startswith("low_priority"):
        return min(rank, 1)
    return rank


DEFAULT_REPRESENTATIVE_SOURCE_ORDER = ["current"]
REPRESENTATIVE_SOURCE_ORDER = list(DEFAULT_REPRESENTATIVE_SOURCE_ORDER)

RUN_LEVEL_EVIDENCE_TYPES = ("short_read", "long_read", "protein")
SCORING_COMPONENT_WEIGHTS = {
    "coding_qc": 25.0,
    "short_read": 25.0,
    "long_read": 15.0,
    "protein": 20.0,
    "structure": 10.0,
}
SCORING_TOTAL_COMPONENT_WEIGHT = sum(SCORING_COMPONENT_WEIGHTS.values())
SCORING_NORMALIZATION = "run_level"
ACTIVE_RUN_LEVEL_EVIDENCE_TYPES = set(RUN_LEVEL_EVIDENCE_TYPES)


def parse_source_priority(values: Sequence[str]) -> List[str]:
    order: List[str] = []
    for value in values or []:
        for item in str(value).split(","):
            source = item.strip()
            if source and source not in order:
                order.append(source)
    if "current" not in order:
        order.insert(0, "current")
    return order or list(DEFAULT_REPRESENTATIVE_SOURCE_ORDER)


def parse_active_evidence_types(values: Sequence[str]) -> set[str]:
    if not values:
        return set(RUN_LEVEL_EVIDENCE_TYPES)
    aliases = {
        "rna": "short_read",
        "short_reads": "short_read",
        "short-read": "short_read",
        "short_read_junction": "short_read",
        "short_read_junctions": "short_read",
        "long_reads": "long_read",
        "long-read": "long_read",
        "long_read_transcriptome": "long_read",
        "homolog": "protein",
        "homolog_protein": "protein",
        "protein_gff": "protein",
    }
    active: set[str] = set()
    for value in values:
        for item in str(value).split(","):
            key = item.strip().lower().replace("-", "_")
            if not key:
                continue
            key = aliases.get(key, key)
            if key not in RUN_LEVEL_EVIDENCE_TYPES:
                raise SystemExit(f"Unknown --active-evidence-types value: {item}")
            active.add(key)
    return active


def run_level_evidence_active(kind: str) -> bool:
    if SCORING_NORMALIZATION != "run_level":
        return True
    return kind in ACTIVE_RUN_LEVEL_EVIDENCE_TYPES


def _model_feature_key(model: Dict[str, object], fields: Sequence[str]) -> Tuple[object, ...]:
    parts: List[object] = [str(model.get("seqid", "")), str(model.get("strand", ""))]
    for field in fields:
        blocks = []
        for item in model.get(field) or []:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                block = [int(item[0]), int(item[1])]
                if field == "cds" and len(item) >= 3:
                    block.append(str(item[2]))
                blocks.append(tuple(block))
        parts.append(tuple(sorted(blocks)))
    return tuple(parts)


def _set_structure_key(model_set: Dict[str, object], raw_model_by_id: Dict[str, Dict[str, object]], fields: Sequence[str]) -> Tuple[object, ...]:
    model_keys = []
    for model_id in model_set.get("model_ids") or []:
        model = raw_model_by_id.get(str(model_id), {})
        if model:
            model_keys.append(_model_feature_key(model, fields))
    if not model_keys:
        return ("no_models", str(model_set.get("set_id", "")))
    return tuple(sorted(model_keys))


def _source_priority(source: object) -> int:
    try:
        return REPRESENTATIVE_SOURCE_ORDER.index(str(source))
    except ValueError:
        return len(REPRESENTATIVE_SOURCE_ORDER)


def _representative_member(members: Sequence[Dict[str, object]]) -> Dict[str, object]:
    return sorted(
        members,
        key=lambda row: (_source_priority(row.get("source", "")), str(row.get("set_id", ""))),
    )[0]


def build_set_collapse_trace(
    locus_id: str,
    real_compact_sets: Sequence[Dict[str, object]],
    raw_models: Sequence[Dict[str, object]],
    set_aliases: Dict[str, str],
) -> List[Dict[str, object]]:
    raw_model_by_id = {str(model.get("model_id", "")): model for model in raw_models if model.get("model_id")}
    groups: Dict[Tuple[object, ...], List[Dict[str, object]]] = {}
    for model_set in real_compact_sets:
        key = _set_structure_key(model_set, raw_model_by_id, ["cds", "introns"])
        groups.setdefault(key, []).append(model_set)

    trace_rows: List[Dict[str, object]] = []
    for members in groups.values():
        representative = _representative_member(members)
        representative_real_set_id = str(representative.get("set_id", ""))
        representative_ai_set_id = set_aliases.get(representative_real_set_id, representative_real_set_id)
        member_real_set_ids = [str(row.get("set_id", "")) for row in members]
        member_ai_set_ids = [set_aliases.get(real_id, real_id) for real_id in member_real_set_ids]
        exon_span_keys = {
            _set_structure_key(row, raw_model_by_id, ["exons"])
            + tuple(str(raw_model_by_id.get(str(model_id), {}).get(key, "")) for model_id in row.get("model_ids") or [] for key in ["start", "end"])
            for row in members
        }
        if len(members) == 1:
            equivalence_scope = "single_set"
            non_decision_difference = "none"
            representative_basis = "single_member"
        elif len(exon_span_keys) == 1:
            equivalence_scope = "CDS_intron_exon_identical"
            non_decision_difference = "none"
            representative_basis = "current_if_CDS_intron_equivalent_else_configured_source_order"
        else:
            equivalence_scope = "CDS_intron_identical_UTR_or_span_differs"
            non_decision_difference = "UTR_or_transcript_span"
            representative_basis = "current_if_CDS_intron_equivalent_else_configured_source_order"
        trace_rows.append(
            {
                "card_id": locus_id,
                "ai_set_id": representative_ai_set_id,
                "representative_real_set_id": representative_real_set_id,
                "representative_ai_set_id": representative_ai_set_id,
                "representative_basis": representative_basis,
                "equivalence_scope": equivalence_scope,
                "non_decision_difference": non_decision_difference,
                "member_real_set_ids": sorted(member_real_set_ids),
                "member_ai_set_ids": sorted(member_ai_set_ids),
                "member_real_sources": sorted({str(row.get("source", "")) for row in members}),
                "representative_real_model_ids": representative.get("model_ids") or [],
                "member_real_model_ids_by_set": {
                    str(row.get("set_id", "")): [str(model_id) for model_id in (row.get("model_ids") or [])]
                    for row in members
                },
            }
        )
    return sorted(trace_rows, key=lambda row: str(row.get("ai_set_id", "")))


def _best_evidence_state(values: Sequence[object], order: Sequence[str]) -> str:
    ranks = {value: index for index, value in enumerate(order)}
    known = [str(value) for value in values]
    return sorted(known, key=lambda value: ranks.get(value, len(order)))[0] if known else "not_configured"


def _interpretation_for_ai_row(row: Dict[str, object]) -> str:
    evidence = row.get("evidence") or {}
    fit = str(row.get("component_fit", ""))
    if str(evidence.get("coding_qc", "")) == "fail":
        return "biologically_problematic"
    if fit in {"misses_core_component", "includes_unsupported_component"}:
        return "topology_risk"
    if evidence.get("short_read_junction") in {"full", "partial"} or evidence.get("long_read") in {"exact", "partial"} or evidence.get("protein") in {"supported", "partial"}:
        return "biologically_supported"
    if evidence.get("short_read_junction") == "conflict" or evidence.get("long_read") == "conflict":
        return "evidence_conflict"
    return "plausible_low_evidence"


def _merge_equivalent_ai_rows(members: Sequence[Dict[str, object]], trace: Dict[str, object]) -> Dict[str, object]:
    by_set = {str(row.get("set_id", "")): row for row in members}
    representative = dict(by_set.get(str(trace.get("ai_set_id", "")), sorted(members, key=lambda row: str(row.get("set_id", "")))[0]))
    evidence_rows = [row.get("evidence") or {} for row in members]
    representative["evidence"] = {
        "coding_qc": _best_evidence_state([row.get("coding_qc") for row in evidence_rows], ["pass", "partial_pass", "warning", "not_assessed", "fail"]),
        "short_read_junction": _best_evidence_state([row.get("short_read_junction") for row in evidence_rows], ["full", "partial", "conflict", "not_assessable", "not_observed", "not_configured"]),
        "long_read": _best_evidence_state([row.get("long_read") for row in evidence_rows], ["exact", "partial", "conflict", "not_assessable", "not_observed", "not_available", "not_configured"]),
        "protein": _best_evidence_state([row.get("protein") for row in evidence_rows], ["supported", "partial", "not_observed", "not_configured"]),
    }
    representative["risks"] = _clean_tags([risk for row in members for risk in (row.get("risks") or [])])
    representative["interpretation"] = _interpretation_for_ai_row(representative)
    return representative


def _collapse_ai_sets_by_structure(rows: Sequence[Dict[str, object]], set_collapse_trace: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    by_set = {str(row.get("set_id", "")): row for row in rows}
    out: List[Dict[str, object]] = []
    for trace in set_collapse_trace:
        member_ids = [str(item) for item in (trace.get("member_ai_set_ids") or [])]
        members = [by_set[set_id] for set_id in member_ids if set_id in by_set]
        if not members:
            continue
        out.append(_merge_equivalent_ai_rows(members, trace))
    return sorted(out, key=lambda row: str(row.get("set_id", "")))


def _trace_by_ai_set(set_trace: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    return {str(row.get("ai_set_id", "")): row for row in (set_trace.get("sets") or []) if row.get("ai_set_id")}


def _strong_transcript_support_present(evidence: Dict[str, object]) -> bool:
    return evidence.get("short_read_junction") == "full" or evidence.get("long_read") == "exact"


def _has_strong_risk(row: Dict[str, object]) -> bool:
    risks = {str(item) for item in (row.get("risks") or [])}
    strong_tokens = (
        "conflicting_RNA_junction",
        "conflicting_short_read_junction",
        "possible_fusion_if_components_are_independent",
        "misses_core_component",
        "includes_unsupported_component",
    )
    return any(token in risks for token in strong_tokens)


def _has_strong_competing_evidence(row: Dict[str, object]) -> bool:
    evidence = row.get("evidence") or {}
    return (
        row.get("component_fit") == "covers_core_components"
        and row.get("interpretation") == "biologically_supported"
        and evidence.get("coding_qc") == "pass"
        and _strong_transcript_support_present(evidence)
    )


def _higher_auto_risk(left: str, right: str) -> str:
    ranks = {"low": 0, "low_to_moderate": 1, "moderate": 2, "high": 3}
    return left if ranks.get(left, 0) >= ranks.get(right, 0) else right


def _auto_risk_annotation(selected: Dict[str, object], gate: str) -> Tuple[str, List[str]]:
    interpretation = str(selected.get("interpretation", ""))
    risks = {str(item) for item in (selected.get("risks") or [])}
    reasons: List[str] = []

    if interpretation == "biologically_supported":
        level = "low"
        reasons.append("selected_set_has_direct_biological_support")
    elif interpretation == "plausible_low_evidence":
        level = "moderate"
        reasons.append("single_selected_set_has_limited_direct_evidence")
    elif interpretation == "evidence_conflict":
        level = "moderate"
        reasons.append("selected_set_has_conflicting_evidence")
    else:
        level = "high"
        reasons.append("selected_set_has_topology_or_coding_risk")

    if _has_strong_risk(selected):
        level = "high"
        reasons.append("strong_risk_tag_present")
    elif risks:
        level = _higher_auto_risk(level, "low_to_moderate")
        reasons.append("nonblocking_risk_tags_present")

    if gate == "Gate_A_single_set_after_CDS_intron_collapse":
        reasons.append("no_alternative_CDS_intron_set_after_collapse")
    elif gate == "Gate_B_dominant_supported_set":
        reasons.append("dominant_supported_set_passed_strict_gate")
    elif gate == "Gate_C_terminal_boundary_equivalent":
        reasons.append("terminal_boundary_only_without_visible_evidence_discriminator")

    return level, sorted(set(reasons))


def _terminal_gate_visible_signature(row: Dict[str, object]) -> Tuple[object, ...]:
    evidence = row.get("evidence") or {}
    return (
        str(row.get("component_fit", "")),
        str(row.get("topology_claim", "")),
        str(row.get("interpretation", "")),
        tuple((field, str(evidence.get(field, ""))) for field in EVIDENCE_FIELD_ORDER),
        tuple(sorted(str(item) for item in (row.get("risks") or []))),
    )


def _terminal_boundary_sets_are_ai_equivalent(card: Dict[str, object]) -> bool:
    sets = card.get("sets") or []
    if len(sets) <= 1:
        return False
    differences = card.get("set_differences") or {}
    structure = differences.get("structure_differences") or {}
    matrix = card.get("component_matrix") or {}
    if structure.get("overall") != "terminal_boundary_only":
        return False
    if structure.get("intron_chain") != "identical":
        return False
    if structure.get("gene_count") != "same" or structure.get("transcript_count") != "same":
        return False
    if matrix.get("variable_components"):
        return False
    if differences.get("discriminating_evidence"):
        return False
    return len({_terminal_gate_visible_signature(row) for row in sets}) == 1


def _terminal_gate_rank(set_id: str, trace_map: Dict[str, Dict[str, object]]) -> Tuple[int, str]:
    trace = trace_map.get(set_id, {})
    source_ranks = [_source_priority(src) for src in (trace.get("member_real_sources") or [])]
    return (min(source_ranks) if source_ranks else len(REPRESENTATIVE_SOURCE_ORDER), set_id)


def _ai_component_uncertainty(value: object) -> str:
    text = str(value or "")
    if text == "terminal_boundary_or_transcript_span":
        return "shared_component_context"
    return text


def _ai_topology_main_question(set_differences: Dict[str, object], locus_topology_summary: Dict[str, object]) -> str:
    structure = set_differences.get("structure_differences") or {}
    overall = str(structure.get("overall", ""))
    if overall == "split_merge_or_gene_count_difference":
        return "component_topology_or_gene_count"
    if overall == "intron_chain_difference":
        return "which_intron_chain_is_best_supported"
    if overall == "internal_exon_or_CDS_difference":
        return "which_internal_CDS_or_exon_structure_is_best_supported"
    if overall == "terminal_boundary_only":
        return "terminal_boundary_tie_break"
    return str(locus_topology_summary.get("main_uncertainty", "which listed annotation set is best supported"))


def _decision_focus_from_differences(set_differences: Dict[str, object]) -> str:
    structure = set_differences.get("structure_differences") or {}
    overall = str(structure.get("overall", ""))
    evidence_fields = [
        str(item.get("field", ""))
        for item in (set_differences.get("discriminating_evidence") or [])
        if isinstance(item, dict) and item.get("type") == "evidence_difference" and item.get("field")
    ]
    evidence_suffix = ""
    if evidence_fields:
        evidence_suffix = " Evidence differs in: {0}.".format(", ".join(sorted(set(evidence_fields))))
    if overall == "split_merge_or_gene_count_difference":
        return "Compare component topology: independent genes versus fused, extra, or missing components." + evidence_suffix
    if overall == "intron_chain_difference":
        return "Compare alternative intron chains using transcript junction and long-read chain evidence." + evidence_suffix
    if overall == "internal_exon_or_CDS_difference":
        return "Compare internal CDS/exon structure using coding QC, transcript support, and protein support." + evidence_suffix
    if overall == "terminal_boundary_only":
        return "Candidate sets mainly differ at terminal CDS/exon boundaries; listed evidence may not distinguish them." + evidence_suffix
    if overall == "single_set":
        return "Only one annotation set remains after deterministic collapse."
    return "Select the listed annotation set best supported by coding QC, transcript evidence, protein evidence, and component topology."


def auto_decision_for_card(card: Dict[str, object], set_trace: Dict[str, object]) -> Dict[str, object] | None:
    sets = card.get("sets") or []
    trace_map = _trace_by_ai_set(set_trace)
    comparison = card.get("set_comparison") or {}
    topology = card.get("topology") or {}
    if not sets:
        return None

    if len(sets) == 1:
        selected = sets[0]
        selected_id = str(selected.get("set_id", ""))
        trace = trace_map.get(selected_id, {})
        gate = "Gate_A_single_set_after_CDS_intron_collapse"
        risk_level, risk_reasons = _auto_risk_annotation(selected, gate)
        return {
            "card_id": card.get("card_id", ""),
            "decision": "select_annotation_set",
            "selected_set_id": selected_id,
            "representative_real_set_id": trace.get("representative_real_set_id", ""),
            "auto_gate": gate,
            "confidence": "high" if risk_level == "low" else "medium" if risk_level in {"low_to_moderate", "moderate"} else "low",
            "auto_risk_level": risk_level,
            "auto_risk_reason_tags": risk_reasons,
            "reason_tags": ["single_candidate_after_CDS_intron_collapse"],
            "risk_tags": selected.get("risks") or [],
            "evidence_summary": "Only one CDS/intron annotation set remains after deterministic structural collapse; no API arbitration is needed.",
            "requires_validation": True,
        }

    if _terminal_boundary_sets_are_ai_equivalent(card):
        selected_id = sorted([str(row.get("set_id", "")) for row in sets], key=lambda set_id: _terminal_gate_rank(set_id, trace_map))[0]
        selected = next(row for row in sets if str(row.get("set_id", "")) == selected_id)
        trace = trace_map.get(selected_id, {})
        gate = "Gate_C_terminal_boundary_equivalent"
        risk_level, risk_reasons = _auto_risk_annotation(selected, gate)
        return {
            "card_id": card.get("card_id", ""),
            "decision": "select_annotation_set",
            "selected_set_id": selected_id,
            "representative_real_set_id": trace.get("representative_real_set_id", ""),
            "auto_gate": gate,
            "confidence": "medium",
            "auto_risk_level": risk_level,
            "auto_risk_reason_tags": risk_reasons,
            "reason_tags": [
                "terminal_boundary_only",
                "same_ai_visible_evidence",
                "same_component_topology",
                "deterministic_representative_selected",
            ],
            "risk_tags": selected.get("risks") or [],
            "evidence_summary": "Multiple annotation sets differ only at terminal CDS/exon boundaries and have identical AI-visible evidence, topology, component fit, and risk tags; API arbitration is skipped and the deterministic representative is selected.",
            "requires_validation": True,
        }

    if comparison.get("main_conflict") != "none" or topology.get("possible_split_merge") is True:
        return None

    all_sources = {str(src) for row in (set_trace.get("sets") or []) for src in (row.get("member_real_sources") or []) if str(src)}
    if len(all_sources) < 3:
        return None

    candidates = []
    for row in sets:
        set_id = str(row.get("set_id", ""))
        evidence = row.get("evidence") or {}
        trace = trace_map.get(set_id, {})
        source_count = len({str(src) for src in (trace.get("member_real_sources") or []) if str(src)})
        if (
            source_count > len(all_sources) / 2.0
            and row.get("component_fit") == "covers_core_components"
            and row.get("interpretation") == "biologically_supported"
            and evidence.get("coding_qc") == "pass"
            and _strong_transcript_support_present(evidence)
            and not _has_strong_risk(row)
        ):
            candidates.append((source_count, set_id, row, trace))
    if len(candidates) != 1:
        return None
    source_count, selected_id, selected, trace = sorted(candidates, key=lambda item: (-item[0], item[1]))[0]
    strong_competitors = [
        str(row.get("set_id", ""))
        for row in sets
        if str(row.get("set_id", "")) != selected_id and _has_strong_competing_evidence(row)
    ]
    if strong_competitors:
        return None
    gate = "Gate_B_dominant_supported_set"
    risk_level, risk_reasons = _auto_risk_annotation(selected, gate)
    return {
        "card_id": card.get("card_id", ""),
        "decision": "select_annotation_set",
        "selected_set_id": selected_id,
        "representative_real_set_id": trace.get("representative_real_set_id", ""),
        "auto_gate": gate,
        "confidence": "high",
        "auto_risk_level": risk_level,
        "auto_risk_reason_tags": risk_reasons,
        "reason_tags": ["dominant_source_support", "coding_qc_pass", "strong_transcript_support_present", "no_strong_competing_set"],
        "risk_tags": selected.get("risks") or [],
        "evidence_summary": "One CDS/intron set has strict majority source support, clean coding QC, full short-read junction support or exact long-read support, no strong topology conflict, and no competing set with comparable strong evidence; API arbitration is skipped.",
        "requires_validation": True,
        "dominant_source_count": source_count,
        "total_source_count": len(all_sources),
    }


def _pairs_from_blocks(values: Sequence[object]) -> Tuple[Tuple[int, int], ...]:
    blocks = []
    for item in values or []:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            try:
                blocks.append((int(item[0]), int(item[1])))
            except (TypeError, ValueError):
                continue
    return tuple(sorted(set(blocks)))


def _cds_blocks_from_model(model: Dict[str, object]) -> Tuple[Tuple[int, int, str], ...]:
    blocks = []
    for item in model.get("cds") or []:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            try:
                blocks.append((int(item[0]), int(item[1]), str(item[2]) if len(item) >= 3 else "."))
            except (TypeError, ValueError):
                continue
    return tuple(sorted(set(blocks)))


def _set_model_ids_for_ai_set(ai_set_id: str, set_collapse_trace: Sequence[Dict[str, object]]) -> List[str]:
    for row in set_collapse_trace:
        if str(row.get("ai_set_id", "")) == ai_set_id:
            return [str(item) for item in (row.get("representative_real_model_ids") or []) if str(item)]
    return []


def _raw_models_for_ai_set(ai_set_id: str, raw_model_by_id: Dict[str, Dict[str, object]], set_collapse_trace: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    return [raw_model_by_id[model_id] for model_id in _set_model_ids_for_ai_set(ai_set_id, set_collapse_trace) if model_id in raw_model_by_id]


def _feature_signature_for_models(models: Sequence[Dict[str, object]], feature: str) -> Tuple[object, ...]:
    signatures = []
    for model in models:
        if feature == "cds":
            value = _cds_blocks_from_model(model)
        elif feature == "introns":
            value = _pairs_from_blocks(model.get("introns") or [])
        elif feature == "exons":
            value = _pairs_from_blocks(model.get("exons") or [])
        else:
            value = ()
        signatures.append((str(model.get("seqid", "")), str(model.get("strand", "")), value))
    return tuple(sorted(signatures))


def _terminal_only_difference(signatures: Sequence[Tuple[object, ...]]) -> bool:
    if len(signatures) <= 1:
        return False
    normalized = []
    for signature in signatures:
        feature_sets = []
        for _seqid, _strand, blocks in signature:
            block_pairs = [(int(item[0]), int(item[1])) for item in blocks if len(item) >= 2]
            if len(block_pairs) <= 2:
                feature_sets.append(tuple())
            else:
                feature_sets.append(tuple(block_pairs[1:-1]))
        normalized.append(tuple(feature_sets))
    return len(set(normalized)) == 1


def _feature_difference_state(signatures: Sequence[Tuple[object, ...]], feature: str) -> str:
    unique = {signature for signature in signatures}
    if len(unique) <= 1:
        return "identical"
    if feature in {"cds", "exons"} and _terminal_only_difference(signatures):
        return "terminal_only_difference"
    return "different"


def _structure_difference_summary(
    sets: Sequence[Dict[str, object]],
    raw_models: Sequence[Dict[str, object]],
    set_collapse_trace: Sequence[Dict[str, object]],
    component_matrix: Dict[str, object],
) -> Dict[str, object]:
    raw_model_by_id = {str(model.get("model_id", "")): model for model in raw_models if model.get("model_id")}
    ai_sets = [str(row.get("set_id", "")) for row in sets if row.get("set_id")]
    models_by_set = {set_id: _raw_models_for_ai_set(set_id, raw_model_by_id, set_collapse_trace) for set_id in ai_sets}
    cds_sigs = [_feature_signature_for_models(models_by_set.get(set_id, []), "cds") for set_id in ai_sets]
    intron_sigs = [_feature_signature_for_models(models_by_set.get(set_id, []), "introns") for set_id in ai_sets]
    exon_sigs = [_feature_signature_for_models(models_by_set.get(set_id, []), "exons") for set_id in ai_sets]
    cds_state = _feature_difference_state(cds_sigs, "cds")
    intron_state = _feature_difference_state(intron_sigs, "introns")
    exon_state = _feature_difference_state(exon_sigs, "exons")
    gene_counts = {str(row.get("set_id", "")): str(row.get("gene_count", "")) for row in sets}
    transcript_counts = {str(row.get("set_id", "")): str(row.get("transcript_count", "")) for row in sets}
    gene_count_state = "same" if len({value for value in gene_counts.values() if value}) <= 1 else "different"
    transcript_count_state = "same" if len({value for value in transcript_counts.values() if value}) <= 1 else "different"
    if len(ai_sets) <= 1:
        overall = "single_set"
        relevance = "low"
    elif component_matrix.get("variable_components") or gene_count_state == "different":
        overall = "split_merge_or_gene_count_difference"
        relevance = "high"
    elif intron_state == "different":
        overall = "intron_chain_difference"
        relevance = "high"
    elif cds_state == "different" or exon_state == "different":
        overall = "internal_exon_or_CDS_difference"
        relevance = "medium"
    elif cds_state == "terminal_only_difference" or exon_state == "terminal_only_difference":
        overall = "terminal_boundary_only"
        relevance = "low"
    else:
        overall = "no_visible_structure_difference_after_collapse"
        relevance = "low"
    return {
        "overall": overall,
        "cds": cds_state,
        "intron_chain": intron_state,
        "exon_structure": exon_state,
        "gene_count": gene_count_state,
        "transcript_count": transcript_count_state,
        "decision_relevance": relevance,
    }


EVIDENCE_FIELD_ORDER = ["coding_qc", "short_read_junction", "long_read", "protein"]


def _collapsed_member_set_ids(set_id: str, set_collapse_trace: Sequence[Dict[str, object]]) -> List[str]:
    for row in set_collapse_trace:
        if str(row.get("ai_set_id", "")) == set_id:
            members = [str(item) for item in (row.get("member_ai_set_ids") or []) if str(item)]
            return members or [set_id]
    return [set_id]


def _component_state_for_values(values: Sequence[object]) -> str:
    states = {str(value or "not_included") for value in values}
    present = {state for state in states if state != "not_included"}
    if not present:
        return "absent"
    if len(present) > 1:
        return "mixed"
    state = next(iter(present))
    if state == "covered_within_same_gene":
        return "present_as_fused_component"
    if state == "covered_as_independent_gene":
        return "present_as_independent_gene"
    if state == "covered_as_extra_mono_exon_gene":
        return "present_as_uncertain_extra"
    return "present"


def _component_matrix(
    sets: Sequence[Dict[str, object]],
    internal_views: Sequence[Dict[str, object]],
    set_collapse_trace: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    view_by_set = {str(view.get("set_id", "")): view for view in internal_views}
    component_ids = sorted({
        str(component_id)
        for view in internal_views
        for component_id in (view.get("component_mapping") or {})
        if str(component_id)
    })
    rows: List[Dict[str, object]] = []
    variable_components = []
    for row in sets:
        set_id = str(row.get("set_id", ""))
        member_ids = _collapsed_member_set_ids(set_id, set_collapse_trace)
        component_states: Dict[str, str] = {}
        for component_id in component_ids:
            values = []
            for member_id in member_ids:
                mapping = (view_by_set.get(member_id, {}) or {}).get("component_mapping") or {}
                values.append(mapping.get(component_id, "not_included"))
            component_states[component_id] = _component_state_for_values(values)
        rows.append({"set_id": set_id, "components": component_states})
    for component_id in component_ids:
        states = {str(row["components"].get(component_id, "absent")) for row in rows}
        if len(states) > 1:
            variable_components.append(component_id)
    return {
        "component_ids": component_ids,
        "rows": rows,
        "variable_components": variable_components,
    }


def _evidence_difference_summary(sets: Sequence[Dict[str, object]]) -> Tuple[List[str], List[Dict[str, object]]]:
    shared: List[str] = []
    discriminating: List[Dict[str, object]] = []
    if not sets:
        return shared, discriminating
    for field in EVIDENCE_FIELD_ORDER:
        values_by_set = {str(row.get("set_id", "")): str((row.get("evidence") or {}).get(field, "")) for row in sets}
        values = {value for value in values_by_set.values() if value}
        if len(values) == 1:
            shared.append("{0}={1}".format(field, next(iter(values))))
        elif len(values) > 1:
            discriminating.append({
                "type": "evidence_difference",
                "field": field,
                "values_by_set": values_by_set,
            })
    return shared, discriminating


def _component_difference_rows(component_matrix: Dict[str, object]) -> List[Dict[str, object]]:
    rows = component_matrix.get("rows") or []
    out: List[Dict[str, object]] = []
    for component_id in component_matrix.get("variable_components") or []:
        out.append({
            "type": "component_presence_difference",
            "component_id": component_id,
            "values_by_set": {str(row.get("set_id", "")): str((row.get("components") or {}).get(component_id, "absent")) for row in rows},
        })
    return out


def _direct_support_label_from_evidence(evidence: Dict[str, object]) -> str:
    supports: List[str] = []
    if str(evidence.get("short_read_junction", "")) == "full":
        supports.append("splice")
    if str(evidence.get("long_read", "")) == "exact":
        supports.append("long_read")
    if str(evidence.get("protein", "")) == "supported":
        supports.append("protein")
    if "splice" in supports and "long_read" in supports:
        return "splice_and_long_read"
    if "splice" in supports:
        return "splice"
    if "long_read" in supports:
        return "long_read"
    if "protein" in supports:
        return "protein_only"
    return "none"


def _direct_support_label_from_tags(tags: Sequence[object]) -> str:
    values = {str(tag) for tag in tags}
    splice = "full_short_read_junction_support" in values
    long_read = "long_read_exact_chain_support" in values
    protein = "homolog_protein_support" in values
    if splice and long_read:
        return "splice_and_long_read"
    if splice:
        return "splice"
    if long_read:
        return "long_read"
    if protein:
        return "protein_only"
    return "none"


def _max_direct_support_label(labels: Sequence[str]) -> str:
    rank = {"none": 0, "protein_only": 1, "splice": 2, "long_read": 2, "splice_and_long_read": 3}
    values = [label for label in labels if label]
    if not values:
        return "none"
    return sorted(values, key=lambda label: (-rank.get(label, 0), label))[0]


def _member_views_for_ai_set(
    set_id: str,
    view_by_set: Dict[str, Dict[str, object]],
    set_collapse_trace: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    member_ids = _collapsed_member_set_ids(set_id, set_collapse_trace)
    return [view_by_set[member_id] for member_id in member_ids if member_id in view_by_set]


def _component_ids_matching_roles(component_ids: Sequence[str], component_by_id: Dict[str, Dict[str, object]], roles: set[str]) -> List[str]:
    return sorted(
        component_id
        for component_id in component_ids
        if str((component_by_id.get(component_id) or {}).get("component_role_in_decision", "")) in roles
    )


def _component_decision_summary(
    row: Dict[str, object],
    view_by_set: Dict[str, Dict[str, object]],
    component_by_id: Dict[str, Dict[str, object]],
    set_collapse_trace: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    set_id = str(row.get("set_id", ""))
    member_views = _member_views_for_ai_set(set_id, view_by_set, set_collapse_trace)
    included_components: set[str] = set()
    missing_strong_components: set[str] = set()
    extra_or_unsupported_components: set[str] = set()
    split_or_fusion_possible = False
    for view in member_views:
        mapping = view.get("component_mapping") or {}
        for component_id, state in mapping.items():
            component_id = str(component_id)
            if str(state) != "not_included":
                included_components.add(component_id)
        fit = view.get("component_fit_summary") or {}
        missing_strong_components.update(str(item) for item in (fit.get("strong_components_missing") or []) if str(item))
        extra_or_unsupported_components.update(str(item) for item in (fit.get("uncertain_components_included") or []) if str(item))
        extra_or_unsupported_components.update(str(item) for item in (fit.get("unsupported_components_included") or []) if str(item))
        risk = view.get("topology_risk_summary") or {}
        if risk.get("fusion_risk") == "possible" or risk.get("split_risk") == "possible":
            split_or_fusion_possible = True

    supported_components = _component_ids_matching_roles(
        sorted(included_components),
        component_by_id,
        {"shared_core_component", "strong_component", "candidate_specific_component", "uncertain_extra_component"},
    )
    extra_support_labels = [
        _direct_support_label_from_tags((component_by_id.get(component_id) or {}).get("support_summary") or [])
        for component_id in sorted(extra_or_unsupported_components)
    ]
    evidence = row.get("evidence") or {}
    split_support = _direct_support_label_from_evidence(evidence) if split_or_fusion_possible else "not_applicable"
    extra_support = _max_direct_support_label(extra_support_labels) if extra_or_unsupported_components else "not_applicable"

    reasons: List[str] = []
    if missing_strong_components:
        reasons.append("missing_strong_component")
    if extra_or_unsupported_components and extra_support in {"none", "protein_only"}:
        reasons.append("unsupported_extra_component")
    if split_or_fusion_possible and split_support in {"none", "protein_only"}:
        reasons.append("unsupported_split_or_fusion_topology")
    if str(evidence.get("coding_qc", "")) == "fail":
        reasons.append("coding_qc_fail")

    if "coding_qc_fail" in reasons:
        priority = "low_priority_coding_problem"
    elif "missing_strong_component" in reasons:
        priority = "low_priority_missing_strong_component"
    elif "unsupported_extra_component" in reasons:
        priority = "low_priority_extra_component"
    elif "unsupported_split_or_fusion_topology" in reasons:
        priority = "low_priority_unsupported_split_or_fusion"
    else:
        priority = "normal"

    return {
        "supported_component_count": len(supported_components),
        "unsupported_extra_component_count": len(extra_or_unsupported_components),
        "missing_strong_component_count": len(missing_strong_components),
        "split_or_fusion_direct_support": split_support,
        "extra_component_direct_support": extra_support,
        "priority_hint": priority,
        "priority_reasons": sorted(set(reasons)),
    }


def _annotate_component_decision_summaries(
    sets: Sequence[Dict[str, object]],
    internal_views: Sequence[Dict[str, object]],
    arbitration_components: Sequence[Dict[str, object]],
    set_collapse_trace: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    view_by_set = {str(view.get("set_id", "")): view for view in internal_views}
    component_by_id = {str(component.get("component_id", "")): component for component in arbitration_components if component.get("component_id")}
    out: List[Dict[str, object]] = []
    for row in sets:
        annotated = dict(row)
        summary = _component_decision_summary(annotated, view_by_set, component_by_id, set_collapse_trace)
        annotated["component_decision_summary"] = summary
        annotated["candidate_priority"] = summary["priority_hint"]
        annotated["candidate_priority_reasons"] = summary["priority_reasons"]
        out.append(annotated)
    return out


def _difference_class(
    sets: Sequence[Dict[str, object]],
    topology: Dict[str, object],
    component_matrix: Dict[str, object],
    discriminating: Sequence[Dict[str, object]],
    structure_differences: Dict[str, object],
) -> str:
    if len(sets) <= 1:
        return "single_set"
    if topology.get("has_possible_split_merge") or structure_differences.get("overall") == "split_merge_or_gene_count_difference":
        return "split_merge_or_component_topology"
    if component_matrix.get("variable_components"):
        return "component_presence_difference"
    if structure_differences.get("overall") in {"intron_chain_difference", "internal_exon_or_CDS_difference", "terminal_boundary_only"}:
        return str(structure_differences.get("overall"))
    if discriminating:
        return "evidence_difference_between_structures"
    return "CDS_or_intron_structure_difference_without_visible_evidence_discriminator"


def _set_difference_summary(
    sets: Sequence[Dict[str, object]],
    locus_topology_summary: Dict[str, object],
    component_matrix: Dict[str, object],
    between_evidence: Dict[str, object],
    structure_differences: Dict[str, object],
) -> Dict[str, object]:
    shared_evidence, evidence_differences = _evidence_difference_summary(sets)
    component_differences = _component_difference_rows(component_matrix)
    discriminating = component_differences + evidence_differences
    unresolved: List[str] = []
    if len(sets) > 1 and not discriminating:
        unresolved.append("listed_sets_differ_in_CDS_or_intron_structure_but_visible_evidence_does_not_discriminate")
    note = str(between_evidence.get("main_uncertainty", ""))
    if note == "terminal_boundary_or_transcript_span":
        note = ""
    if note and note not in unresolved:
        unresolved.append(note)
    return {
        "difference_class": _difference_class(sets, locus_topology_summary, component_matrix, discriminating, structure_differences),
        "structure_differences": structure_differences,
        "shared_evidence": shared_evidence,
        "discriminating_evidence": discriminating,
        "unresolved_differences": unresolved,
    }


def _compare_ai_sets(rows: Sequence[Dict[str, object]], locus_topology_summary: Dict[str, object]) -> Dict[str, object]:
    if not rows:
        return {
            "discriminator_status": "none",
            "main_conflict": "none",
            "recommended_default": "risk_if_unresolved_conflict",
            "recommended_set_id": "",
        }
    ranked = sorted(rows, key=lambda row: (-_set_quality_rank(row), str(row.get("set_id", ""))))
    best_rank = _set_quality_rank(ranked[0])
    best_rows = [row for row in ranked if _set_quality_rank(row) == best_rank]
    main_conflict = "split_merge" if locus_topology_summary.get("has_possible_split_merge") else "none"
    if best_rank <= 1 and main_conflict != "none" and len(rows) > 1:
        return {
            "discriminator_status": "conflicting",
            "main_conflict": main_conflict,
            "recommended_default": "risk_if_unresolved_conflict",
            "recommended_set_id": str(ranked[0].get("set_id", "")),
        }
    if len(rows) == 1:
        status = "none"
    elif len(best_rows) == 1 and best_rank >= 3:
        status = "clear"
    else:
        status = "weak"
    return {
        "discriminator_status": status,
        "main_conflict": main_conflict,
        "recommended_default": "select_best_listed_set",
        "recommended_set_id": str(ranked[0].get("set_id", "")),
    }


def build_ai_annotation_card(
    locus_id: str,
    compact_sets: Sequence[Dict[str, object]],
    locus_topology_summary: Dict[str, object],
    arbitration_components: Sequence[Dict[str, object]],
    internal_candidate_set_views: Sequence[Dict[str, object]],
    candidate_groups: Sequence[Dict[str, object]],
    between_evidence: Dict[str, object],
    set_collapse_trace: Sequence[Dict[str, object]],
    raw_models: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    view_by_set = {str(view.get("set_id", "")): view for view in internal_candidate_set_views}
    raw_sets: List[Dict[str, object]] = []
    for model_set in compact_sets:
        set_id = str(model_set.get("set_id", ""))
        view = view_by_set.get(set_id, {})
        raw_sets.append(
            {
                "set_id": set_id,
                "component_fit": _component_fit_label(view),
                "topology_claim": view.get("topology_claim", "single_or_unmapped_annotation_hypothesis"),
                "evidence": {
                    "coding_qc": _validation_state(model_set),
                    "short_read_junction": _short_read_state(model_set),
                    "long_read": _long_read_state(model_set),
                    "protein": _protein_state(model_set),
                },
                "interpretation": _set_interpretation(model_set, view_by_set),
                "risks": _clean_tags((view.get("set_level_risks") or []) + (model_set.get("interpreted_risk_flags") or [])),
            }
        )
    sets = _collapse_ai_sets_by_structure(raw_sets, set_collapse_trace)
    sets = _annotate_component_decision_summaries(sets, internal_candidate_set_views, arbitration_components, set_collapse_trace)
    set_comparison = _compare_ai_sets(sets, locus_topology_summary)
    component_matrix = _component_matrix(sets, internal_candidate_set_views, set_collapse_trace)
    structure_differences = _structure_difference_summary(sets, raw_models, set_collapse_trace, component_matrix)
    set_differences = _set_difference_summary(sets, locus_topology_summary, component_matrix, between_evidence, structure_differences)

    components = []
    for component in arbitration_components:
        components.append(
            {
                "id": component.get("component_id", ""),
                "claim": component.get("biological_claim", ""),
                "main_uncertainty": _ai_component_uncertainty(component.get("main_uncertainty", "")),
                "evidence": _clean_tags(component.get("support_summary") or [], limit=6),
                "risks": _clean_tags(component.get("component_risk_flags") or [], limit=6),
            }
        )

    return {
        "card_id": locus_id,
        "decision_rule": {
            "allowed_decisions": ["select_annotation_set", "risk_flagged"],
            "default_action": "select_annotation_set",
            "risk_flagged_only_for": "mutually exclusive candidate sets with no card evidence discriminator",
            "must_select_listed_set_unless_unresolved_conflict": True,
            "do_not_mix_sets": True,
            "do_not_create_coordinates": True,
        },
        "topology": {
            "class": locus_topology_summary.get("topology_class", "single_or_unclassified"),
            "main_question": _ai_topology_main_question(set_differences, locus_topology_summary),
            "possible_split_merge": bool(locus_topology_summary.get("has_possible_split_merge")),
        },
        "components": components,
        "sets": sets,
        "component_matrix": component_matrix,
        "set_differences": set_differences,
        "decision_focus": _decision_focus_from_differences(set_differences),
        "set_comparison": {
            "discriminator_status": set_comparison.get("discriminator_status", ""),
            "main_conflict": set_comparison.get("main_conflict", ""),
            "recommended_default": set_comparison.get("recommended_default", ""),
        },
    }


def compact_card_with_trace(card: Dict[str, object]) -> Tuple[Dict[str, object], Dict[str, object]]:
    """Build the AI-facing annotation card.

    The output is intentionally a single-layer biological summary for model
    selection. It hides coordinates, source/tool identity, raw per-model tables,
    numeric score traces, and correction/deletion task labels. Full-card trace
    fields remain only in the upstream JSONL.
    """
    models = card.get("candidate_models") or []
    currents = current_models(models)
    structure_groups = build_structure_groups(models)

    real_compact_models = [compact_model(model, structure_groups, currents) for model in models]
    real_compact_sets = [compact_model_set(model_set) for model_set in (card.get("candidate_model_sets") or [])]
    alias_mode = "de_novo_annotation"
    set_aliases, _set_alias_map = build_set_aliases(real_compact_sets, alias_mode)
    set_collapse_rows = build_set_collapse_trace(str(card.get("locus_id") or ""), real_compact_sets, models, set_aliases)
    model_aliases, _model_alias_map = build_model_aliases(real_compact_models, alias_mode)

    compact_models = [public_model_for_task(model, alias_mode, model_aliases) for model in real_compact_models]
    compact_sets = [
        public_set_for_task(model_set, alias_mode, set_aliases, model_aliases, index)
        for index, model_set in enumerate(real_compact_sets, start=1)
    ]
    compact_sets, score_trace = add_interpreted_scores(compact_sets, compact_models)

    locus_topology_summary, arbitration_components, candidate_set_views = build_arbitration_component_layers(
        compact_sets,
        compact_models,
    )
    candidate_groups = candidate_structure_groups(compact_sets, compact_models, locus_topology_summary)
    between_evidence = between_candidate_evidence_interpretation(compact_sets)
    ai_card = build_ai_annotation_card(
        str(card.get("locus_id") or ""),
        compact_sets,
        locus_topology_summary,
        arbitration_components,
        candidate_set_views,
        candidate_groups,
        between_evidence,
        set_collapse_rows,
        models,
    )
    return ai_card, {"card_id": str(card.get("locus_id") or ""), "sets": set_collapse_rows, "score_trace": score_trace}


def compact_card(card: Dict[str, object]) -> Dict[str, object]:
    return compact_card_with_trace(card)[0]


def _model_source(model: Dict[str, object]) -> str:
    return str(model.get("source", ""))


def _all_review_models(full_card: Dict[str, object]) -> List[Dict[str, object]]:
    models: List[Dict[str, object]] = []
    for key in ["candidate_models", "excluded_candidate_models"]:
        for model in full_card.get(key) or []:
            if isinstance(model, dict):
                models.append(model)
    return models


def _count_by_current(models: Sequence[Dict[str, object]]) -> Tuple[int, int]:
    current_count = 0
    tool_count = 0
    for model in models:
        if _model_source(model) == "current":
            current_count += 1
        else:
            tool_count += 1
    return current_count, tool_count


def _evidence_level_for_set(row: Dict[str, object]) -> str:
    evidence = row.get("evidence") or {}
    coding = str(evidence.get("coding_qc", ""))
    short_read = str(evidence.get("short_read_junction", ""))
    long_read = str(evidence.get("long_read", ""))
    protein = str(evidence.get("protein", ""))
    interpretation = str(row.get("interpretation", ""))
    if coding == "fail" or interpretation == "biologically_problematic":
        return "poor"
    coding_ok = coding in {"pass", "partial_pass", "warning", "not_assessed"}
    if coding_ok and (short_read == "full" or long_read == "exact"):
        return "strong"
    if coding_ok and (short_read == "partial" or long_read == "partial" or protein in {"supported", "partial"}):
        return "moderate"
    if coding_ok:
        return "weak"
    return "none"


def _best_discovery_evidence_level(ai_card: Dict[str, object]) -> str:
    levels = [_evidence_level_for_set(row) for row in (ai_card.get("sets") or []) if isinstance(row, dict)]
    if not levels:
        return "none"
    return sorted(levels, key=lambda level: -EVIDENCE_LEVEL_RANK.get(level, 0))[0]


def _unsupported_locus_reason_tags(ai_card: Dict[str, object], evidence_level: str) -> List[str]:
    if EVIDENCE_LEVEL_RANK.get(evidence_level, 0) >= EVIDENCE_LEVEL_RANK["moderate"]:
        return []
    sets = [row for row in (ai_card.get("sets") or []) if isinstance(row, dict)]
    evidence_rows = [row.get("evidence") or {} for row in sets]
    tags = ["no_listed_set_has_moderate_or_strong_biological_support"]
    if not any(row.get("short_read_junction") == "full" for row in evidence_rows):
        tags.append("no_full_short_read_junction_support")
    if not any(row.get("long_read") == "exact" for row in evidence_rows):
        tags.append("no_exact_long_read_chain_support")
    if not any(row.get("short_read_junction") in {"full", "partial"} or row.get("long_read") in {"exact", "partial"} for row in evidence_rows):
        tags.append("no_transcript_structure_support")
    if not any(row.get("protein") in {"supported", "partial"} for row in evidence_rows):
        tags.append("no_homolog_protein_support")
    if any(row.get("coding_qc") == "fail" for row in evidence_rows):
        tags.append("one_or_more_sets_fail_coding_qc")
    if not sets:
        tags.append("no_ai_visible_annotation_set")
    return sorted(set(tags))


def _review_priority(router_class: str, evidence_level: str, ai_set_count: int) -> str:
    if router_class != "tool_only_novel_gene_candidate":
        return "not_applicable"
    if evidence_level == "strong" and ai_set_count <= 1:
        return "high"
    if evidence_level == "strong":
        return "high_after_arbitration"
    if evidence_level == "moderate":
        return "medium"
    if evidence_level == "weak":
        return "low"
    return "qc_review"


def build_locus_review_router_row(full_card: Dict[str, object], ai_card: Dict[str, object]) -> Dict[str, object]:
    all_models = _all_review_models(full_card)
    selectable_models = [model for model in (full_card.get("candidate_models") or []) if isinstance(model, dict)]
    current_count, tool_count = _count_by_current(all_models)
    current_selectable_count, tool_selectable_count = _count_by_current(selectable_models)
    ai_set_count = len(ai_card.get("sets") or [])
    evidence_level = _best_discovery_evidence_level(ai_card)
    unsupported_reason_tags = _unsupported_locus_reason_tags(ai_card, evidence_level)
    manual_review_unsupported = "true" if unsupported_reason_tags else "false"

    deletion_allowed = "false"
    deletion_policy = "no_deletion_default_retain_current"
    manual_review_deletion = "false"
    manual_review_novel = "false"
    recommended_action = "annotation_set_arbitration"
    reason_tags: List[str] = []

    if current_count == 0 and tool_count > 0:
        router_class = "tool_only_novel_gene_candidate"
        manual_review_novel = "true" if tool_selectable_count > 0 else "false"
        recommended_action = "manual_review_novel_gene_candidate_after_ai_arbitration" if ai_set_count > 1 else "manual_review_novel_gene_candidate"
        reason_tags.extend(["no_current_annotation_overlap", "tool_prediction_present", "novel_gene_discovery_only"])
        if tool_selectable_count == 0:
            reason_tags.append("no_ai_selectable_tool_model")
        if manual_review_unsupported == "true":
            reason_tags.append("unsupported_novel_candidate_requires_manual_review")
    elif current_count > 0 and tool_count == 0:
        router_class = "current_only_retain_no_deletion"
        recommended_action = "retain_current_no_deletion"
        reason_tags.extend(["current_only_locus", "deletion_disabled_by_policy", "insufficient_evidence_is_not_absence"])
        if manual_review_unsupported == "true":
            reason_tags.append("unsupported_current_only_locus_retained_no_deletion")
    elif current_count > 0 and tool_count > 0:
        router_class = "mixed_current_tool_arbitration"
        recommended_action = "annotation_set_arbitration"
        reason_tags.extend(["current_and_tool_models_present", "send_to_annotation_set_arbitration"])
        if manual_review_unsupported == "true":
            reason_tags.append("unsupported_mixed_locus_after_arbitration_review")
    else:
        router_class = "no_selectable_annotation_model"
        recommended_action = "manual_review_no_model"
        reason_tags.append("no_current_or_tool_model_in_card")

    return {
        "card_id": ai_card.get("card_id", full_card.get("locus_id", "")),
        "router_class": router_class,
        "current_model_count": current_count,
        "tool_model_count": tool_count,
        "current_selectable_model_count": current_selectable_count,
        "tool_selectable_model_count": tool_selectable_count,
        "ai_set_count": ai_set_count,
        "routing_decision": "",
        "auto_gate": "",
        "sent_to_ai": "",
        "manual_review_novel_gene_candidate": manual_review_novel,
        "manual_review_deletion_candidate": manual_review_deletion,
        "manual_review_unsupported_locus": manual_review_unsupported,
        "deletion_allowed": deletion_allowed,
        "deletion_policy": deletion_policy,
        "best_evidence_level": evidence_level,
        "unsupported_locus_reason_tags": ",".join(unsupported_reason_tags),
        "discovery_evidence_level": evidence_level,
        "review_priority": _review_priority(router_class, evidence_level, ai_set_count),
        "recommended_manual_action": recommended_action,
        "review_reason_tags": ",".join(sorted(set(reason_tags))),
    }


def finalize_review_router_row(row: Dict[str, object], auto_decision: Dict[str, object] | None) -> Dict[str, object]:
    out = dict(row)
    sent_to_ai = auto_decision is None
    out["routing_decision"] = "send_to_ai" if sent_to_ai else "auto_select"
    out["auto_gate"] = (auto_decision or {}).get("auto_gate", "")
    out["sent_to_ai"] = str(sent_to_ai).lower()
    return out


def annotate_auto_decision_with_review(auto_decision: Dict[str, object] | None, review_row: Dict[str, object]) -> Dict[str, object] | None:
    if not auto_decision:
        return None
    out = dict(auto_decision)
    out["review_router_class"] = review_row.get("router_class", "")
    out["manual_review_novel_gene_candidate"] = review_row.get("manual_review_novel_gene_candidate", "false")
    out["manual_review_deletion_candidate"] = review_row.get("manual_review_deletion_candidate", "false")
    out["manual_review_unsupported_locus"] = review_row.get("manual_review_unsupported_locus", "false")
    out["deletion_policy"] = review_row.get("deletion_policy", "")
    out["recommended_manual_action"] = review_row.get("recommended_manual_action", "")
    tags = [str(tag) for tag in (out.get("reason_tags") or []) if str(tag)]
    if review_row.get("manual_review_novel_gene_candidate") == "true":
        tags.append("novel_gene_candidate_locus")
    if review_row.get("manual_review_unsupported_locus") == "true":
        tags.append("unsupported_locus_manual_review")
    if tags:
        out["reason_tags"] = sorted(set(tags))
    return out


def evidence_tags_for_model(model: Dict[str, object]) -> List[str]:
    tags: List[str] = []
    validation = str(model.get("validation_status", "not_assessed"))
    rna_status = str(model.get("rna_junction_status", "not_assessable"))
    protein_status = str(model.get("protein_support_status", "protein_overlap_none_observed"))
    splice_status = str(model.get("splice_motif_status", ""))
    phase_status = str(model.get("phase_status", ""))
    if validation:
        tags.append("validation_" + validation)
    if rna_status:
        tags.append("rna_" + rna_status)
    if protein_status == "protein_overlap_present":
        tags.append("protein_supported")
    else:
        tags.append("protein_not_observed")
    if splice_status:
        tags.append("splice_" + splice_status)
    if phase_status:
        tags.append("phase_" + phase_status)
    return sorted(set(tags))


def blind_risk_tags(model: Dict[str, object]) -> List[str]:
    """移除会泄漏来源、当前注释身份或仅当前注释证据通道的标签。"""
    blocked_tokens = (
        "current",
        "de_novo",
        "denovo",
        "source",
        "protein_chain",
        "model_level",
        "same_structure_as_current",
    )
    out = []
    for tag in model.get("risk_tags") or []:
        value = str(tag)
        lower = value.lower()
        if any(token in lower for token in blocked_tokens):
            continue
        out.append(value)
    return sorted(set(out))


def blind_model(alias: str, model: Dict[str, object]) -> Dict[str, object]:
    structure = model.get("structure_summary") or {}
    return {
        "alias": alias,
        "model_type": "full_model",
        "structure": {
            "exon_count": structure.get("exon_count", 0),
            "intron_count": structure.get("intron_count", 0),
            "cds_count": structure.get("cds_count", 0),
            "cds_length": structure.get("cds_length", 0),
        },
        "coding": {
            "status": model.get("validation_status", "not_assessed"),
            "orf": model.get("orf_status", ""),
            "phase": model.get("phase_status", ""),
            "splice_motif": model.get("splice_motif_status", ""),
            "protein_concordance": model.get("protein_concordance_status", ""),
            "validation_fail_reasons": model.get("validation_fail_reasons", ""),
        },
        "rna": {
            "junction_status": model.get("rna_junction_status", "not_assessable"),
            "support_fraction": model.get("junction_support_fraction", "0.0000"),
            "conflicting_junctions": model.get("conflicting_junction_count", 0),
        },
        "protein": {
            "supported": model.get("protein_support_status") == "protein_overlap_present",
            "best_identity": model.get("protein_best_identity", ""),
        },
        "evidence_tags": evidence_tags_for_model(model),
        "risk_tags": blind_risk_tags(model),
        "_structure_group_id": model.get("structure_group_id", ""),
    }


def blind_structure_relations(alias_models: Sequence[Dict[str, object]]) -> Dict[str, Dict[str, List[str]]]:
    by_alias: Dict[str, Dict[str, List[str]]] = {
        str(model["alias"]): {
            "same_structure_as": [],
            "same_intron_count_as": [],
            "alternative_structure_vs": [],
            "similar_cds_length_as": [],
        }
        for model in alias_models
    }
    for i, left in enumerate(alias_models):
        left_alias = str(left["alias"])
        left_structure = left.get("_structure_group_id", "")
        left_struct = left.get("structure") or {}
        left_intron_count = int(left_struct.get("intron_count") or 0)
        left_cds_length = int(left_struct.get("cds_length") or 0)
        for right in alias_models[i + 1 :]:
            right_alias = str(right["alias"])
            right_structure = right.get("_structure_group_id", "")
            right_struct = right.get("structure") or {}
            right_intron_count = int(right_struct.get("intron_count") or 0)
            right_cds_length = int(right_struct.get("cds_length") or 0)
            if left_structure and left_structure == right_structure:
                by_alias[left_alias]["same_structure_as"].append(right_alias)
                by_alias[right_alias]["same_structure_as"].append(left_alias)
            elif left_intron_count and right_intron_count:
                by_alias[left_alias]["alternative_structure_vs"].append(right_alias)
                by_alias[right_alias]["alternative_structure_vs"].append(left_alias)
            if left_intron_count == right_intron_count:
                by_alias[left_alias]["same_intron_count_as"].append(right_alias)
                by_alias[right_alias]["same_intron_count_as"].append(left_alias)
            if left_cds_length and right_cds_length and abs(left_cds_length - right_cds_length) <= 30:
                by_alias[left_alias]["similar_cds_length_as"].append(right_alias)
                by_alias[right_alias]["similar_cds_length_as"].append(left_alias)
    for relation in by_alias.values():
        for key, aliases in list(relation.items()):
            relation[key] = sorted(set(aliases))
    return by_alias


def blind_rank_card(card: Dict[str, object], task_mode: str = "correction") -> Dict[str, object]:
    standard = compact_card(card, task_mode)
    models = standard.get("candidate_model_table") or []
    alias_models: List[Dict[str, object]] = []
    alias_map: List[Dict[str, object]] = []
    for index, model in enumerate(models, start=1):
        alias = "M{0:02d}".format(index)
        alias_model = blind_model(alias, model)
        alias_models.append(alias_model)
        alias_map.append(
            {
                "alias": alias,
                "object_type": "full_model",
                "model_id": model.get("model_id", ""),
                "source": model.get("source", ""),
            }
        )
    relations = blind_structure_relations(alias_models)
    for model in alias_models:
        model["relations"] = relations.get(str(model["alias"]), {})
        model.pop("_structure_group_id", None)
    policy = standard.get("_decision_policy") or standard.get("decision_policy") or {}
    return {
        "card_schema": "blind_evidence_ranking",
        "locus_id": standard.get("locus_id", ""),
        "locus_summary": {
            "candidate_model_count": len(alias_models),
        },
        "decision_policy": policy,
        "anonymous_candidate_table": alias_models,
        "ranking_task": {
            "objective": "rank anonymous full gene-model aliases by evidence support",
            "tie_rule_for_ai": "if evidence cannot distinguish multiple top aliases, return all tied aliases as best_aliases",
        },
        "required_output_schema": {
            "must_decide_required_fields": [
                "locus_id",
                "best_aliases",
                "model_scores",
                "tie_break_needed",
                "reason_tags",
                "risk_tags",
                "evidence_summary",
            ],
            "risk_pool_extra_fields": [
                "automatic_decision_ready",
                "risk_reasons",
            ],
            "score_rule": "model_scores must include every listed alias with a 0-100 within-card evidence support score; scores are not benchmark accuracy",
            "id_rules": "aliases must be listed in anonymous_candidate_table",
        },
        "_alias_map": alias_map,
    }


def summary_row(
    card: Dict[str, object],
    auto_decision: Dict[str, object] | None = None,
    review_row: Dict[str, object] | None = None,
) -> Dict[str, object]:
    topology = card.get("topology") or {}
    comparison = card.get("set_comparison") or {}
    review = review_row or {}
    return {
        "card_id": card.get("card_id", ""),
        "set_count": len(card.get("sets") or []),
        "component_count": len(card.get("components") or []),
        "topology_class": topology.get("class", ""),
        "main_question": topology.get("main_question", ""),
        "possible_split_merge": str(bool(topology.get("possible_split_merge"))).lower(),
        "main_conflict": comparison.get("main_conflict", ""),
        "discriminator_status": comparison.get("discriminator_status", ""),
        "recommended_default": comparison.get("recommended_default", ""),
        "routing_decision": "auto_select" if auto_decision else "send_to_ai",
        "auto_gate": (auto_decision or {}).get("auto_gate", ""),
        "auto_selected_set_id": (auto_decision or {}).get("selected_set_id", ""),
        "auto_representative_real_set_id": (auto_decision or {}).get("representative_real_set_id", ""),
        "auto_risk_level": (auto_decision or {}).get("auto_risk_level", ""),
        "review_router_class": review.get("router_class", ""),
        "manual_review_novel_gene_candidate": review.get("manual_review_novel_gene_candidate", ""),
        "manual_review_deletion_candidate": review.get("manual_review_deletion_candidate", ""),
        "manual_review_unsupported_locus": review.get("manual_review_unsupported_locus", ""),
        "deletion_policy": review.get("deletion_policy", ""),
        "recommended_manual_action": review.get("recommended_manual_action", ""),
    }


def is_risk_flagged_card(card: Dict[str, object]) -> bool:
    router = card.get("risk_router") or card.get("manual_review_router") or {}
    decision_required = router.get("decision_required")
    decision_requirement = str(router.get("decision_requirement") or card.get("decision_requirement") or "")
    action_space = {str(item) for item in (router.get("recommended_action_space") or card.get("allowed_decision_scope") or [])}
    return (
        decision_required is False
        or decision_requirement in {"risk_flagged_no_ai_decision", "risk_flagged_no_ai_eligible_candidate_models", "human_review_only"}
        or action_space in [{"risk_flagged"}, {"human_review"}]
    )


def main() -> None:
    args = parse_args()
    global REPRESENTATIVE_SOURCE_ORDER, SCORING_NORMALIZATION, ACTIVE_RUN_LEVEL_EVIDENCE_TYPES
    REPRESENTATIVE_SOURCE_ORDER = parse_source_priority(args.source_priority)
    SCORING_NORMALIZATION = args.scoring_normalization
    ACTIVE_RUN_LEVEL_EVIDENCE_TYPES = parse_active_evidence_types(args.active_evidence_types)
    full_cards = read_jsonl(args.input_jsonl)
    if not args.include_risk_flagged:
        full_cards = [card for card in full_cards if not is_risk_flagged_card(card)]
    records = [compact_card_with_trace(card) for card in full_cards]
    compact_cards = [card for card, _trace in records]
    set_traces = [trace for _card, trace in records]
    review_rows_raw = [build_locus_review_router_row(full_card, compact_card) for full_card, compact_card in zip(full_cards, compact_cards)]
    auto_decisions_raw = [auto_decision_for_card(card, trace) for card, trace in records]
    review_rows = [finalize_review_router_row(row, decision) for row, decision in zip(review_rows_raw, auto_decisions_raw)]
    auto_decisions = [annotate_auto_decision_with_review(decision, row) for decision, row in zip(auto_decisions_raw, review_rows)]
    auto_records = [decision for decision in auto_decisions if decision]
    ai_cards = [card for card, decision in zip(compact_cards, auto_decisions_raw) if not decision]
    write_jsonl(args.output_jsonl, ai_cards)
    write_jsonl(args.output_set_trace, set_traces)
    write_jsonl(args.output_auto_decisions, auto_records)
    write_tsv(args.output_locus_review_router, REVIEW_ROUTER_FIELDS, review_rows)
    write_tsv(args.output_summary, SUMMARY_FIELDS, [summary_row(card, decision, review) for card, decision, review in zip(compact_cards, auto_decisions, review_rows)])


if __name__ == "__main__":
    main()
