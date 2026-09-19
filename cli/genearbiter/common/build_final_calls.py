#!/usr/bin/env python3
# script_id_md5: d15f90b6604ea0c47bab183776fb45b6
# created: 2026-07-04
# modified: 2026-07-05
# owner: project
# status: project_code
# purpose: 合并自动和模型 annotation-set decisions，生成冻结 final calls 表。
# inputs: set trace JSONL、auto decisions JSONL、model decisions JSONL、locus review router TSV。
# outputs: final_annotation_calls.tsv、final_annotation_call_summary.tsv、final_annotation_call_report.md。
# notes: 不读取 manual truth，不生成 GFF，不允许删除；novel/unsupported 只作为审核标签，不阻断 annotation-set proposal。

"""Build final annotation calls from routed annotation-set decisions."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence


DETAIL_FIELDS = [
    "card_id",
    "locus_id",
    "decision_source",
    "profile",
    "model",
    "api_status",
    "decision_valid",
    "validation_errors",
    "decision",
    "selected_ai_set_id",
    "selected_id",
    "selected_real_set_id",
    "selected_source",
    "selected_real_sources",
    "selected_model_id",
    "selected_real_model_ids",
    "gene_set_selected_set_id",
    "member_real_set_ids",
    "member_real_sources",
    "representative_basis",
    "equivalence_scope",
    "non_decision_difference",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "confidence",
    "requires_validation",
    "review_router_class",
    "manual_review_novel_gene_candidate",
    "manual_review_deletion_candidate",
    "manual_review_unsupported_locus",
    "deletion_allowed",
    "deletion_policy",
    "best_evidence_level",
    "discovery_evidence_level",
    "review_priority",
    "recommended_manual_action",
    "routing_decision",
    "auto_gate",
    "auto_risk_level",
    "auto_risk_reason_tags",
    "reason_tags",
    "risk_tags",
    "evidence_summary",
    "blocking_reasons",
]

SUMMARY_FIELDS = ["metric", "value"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set-trace", required=True, help="GeneArbiter set-trace JSONL")
    parser.add_argument("--auto-decisions", required=True, help="GeneArbiter automatic-decision JSONL")
    parser.add_argument(
        "--api-decisions",
        action="append",
        default=[],
        help="AI decisions.jsonl, or a profile/out directory containing decisions.jsonl. Repeatable.",
    )
    parser.add_argument("--locus-review-router", required=True, help="GeneArbiter locus review-routing TSV")
    parser.add_argument("--out-dir", required=True, help="Output directory for final annotation calls")
    parser.add_argument(
        "--allow-duplicate-api-decisions",
        action="store_true",
        help="Keep the first API decision for a card when repeated decision files contain duplicates.",
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if not path.exists():
        raise SystemExit(f"Missing JSONL input: {path}")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"Missing TSV input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Dict[str, Any]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_report(path: Path, summary_rows: Sequence[Dict[str, Any]]) -> None:
    lines = [
        "# GeneArbiter Final Annotation Calls",
        "",
        f"Generated: {utc_now()}",
        "",
        "Scope: frozen GeneArbiter annotation-set call table for export and audit.",
        "This table does not generate coordinates or GFF records and does not permit deletion.",
        "",
        "## Summary",
        "",
    ]
    for row in summary_rows:
        lines.append(f"- {row['metric']}: {row['value']}")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `keep_current` / `keep_current_flagged` retain the current representative set under the no-deletion policy.",
            "- `annotation_set_selected_pending_validation` selects a non-current listed set but still requires downstream validation/export handling.",
            "- `novel_gene_candidate_manual_review` marks tool-only discovery candidates; it remains a selected annotation-set proposal and should be counted in proposal-level benchmark.",
            "- `annotation_set_selected_needs_manual_review` marks selected non-current sets with unsupported-locus review flags; it remains a selected annotation-set proposal.",
            "- `risk_flagged_or_blocked`, `invalid_decision_blocked`, and `missing_decision` are not export-ready calls; risk decisions with a selected set remain provisional proposals with review/export-risk labels.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def as_list_text(values: Any) -> str:
    if values is None:
        return ""
    if isinstance(values, list):
        return ";".join(str(item) for item in values if str(item))
    return str(values)


def is_true(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def card_id(row: Dict[str, Any]) -> str:
    return str(row.get("card_id") or row.get("locus_id") or "")


def trace_by_card(path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for row in read_jsonl(path):
        cid = card_id(row)
        if cid:
            out[cid] = row
    return out


def router_by_card(path: Path) -> Dict[str, Dict[str, str]]:
    return {row.get("card_id", ""): row for row in read_tsv(path) if row.get("card_id")}


def auto_by_card(path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for row in read_jsonl(path):
        cid = card_id(row)
        if cid:
            out[cid] = row
    return out


def decision_files(paths: Sequence[str]) -> List[Path]:
    files: List[Path] = []
    for value in paths:
        path = Path(value)
        if path.is_file():
            files.append(path)
        elif path.is_dir() and (path / "decisions.jsonl").exists():
            files.append(path / "decisions.jsonl")
        elif path.is_dir():
            files.extend(sorted(child / "decisions.jsonl" for child in path.iterdir() if (child / "decisions.jsonl").exists()))
        else:
            raise SystemExit(f"Missing API decision input: {path}")
    return files


def api_by_card(paths: Sequence[str], allow_duplicate: bool) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    duplicates: List[str] = []
    for path in decision_files(paths):
        for row in read_jsonl(path):
            cid = card_id(row)
            if not cid:
                continue
            if cid in out:
                duplicates.append(cid)
                if allow_duplicate:
                    continue
            out[cid] = row
    if duplicates and not allow_duplicate:
        preview = ",".join(sorted(set(duplicates))[:10])
        raise SystemExit(f"Duplicate API decisions for card_id(s): {preview}. Use --allow-duplicate-api-decisions to keep the first occurrence.")
    return out


def trace_set_map(trace: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {str(row.get("ai_set_id", "")): row for row in (trace.get("sets") or []) if row.get("ai_set_id")}


def selected_trace(decision: Dict[str, Any], trace: Dict[str, Any]) -> Dict[str, Any]:
    selected_id = str(decision.get("selected_set_id") or "")
    if not selected_id:
        return {}
    return trace_set_map(trace).get(selected_id, {})


def representative_source(real_set_id: str, member_sources: Sequence[Any]) -> str:
    if real_set_id == "current_set":
        return "current"
    if real_set_id.endswith("_set"):
        return real_set_id[: -len("_set")]
    values = [str(item) for item in member_sources if str(item)]
    return values[0] if values else ""


def decision_is_valid(decision: Dict[str, Any]) -> bool:
    if not decision:
        return False
    if "valid" in decision:
        return is_true(decision.get("valid"))
    status = str(decision.get("api_status", ""))
    return status in {"ok", "local_ok"} or not status


def classify_call(decision: Dict[str, Any], trace_row: Dict[str, Any], router: Dict[str, str]) -> Dict[str, str]:
    blocking: List[str] = []
    decision_value = str(decision.get("decision", ""))
    selected_ai = str(decision.get("selected_set_id") or "")
    selected_real_set = str(trace_row.get("representative_real_set_id", ""))
    selected_sources = [str(item) for item in (trace_row.get("member_real_sources") or []) if str(item)]
    selected_is_current = selected_real_set == "current_set" or "current" in selected_sources
    novel = is_true(router.get("manual_review_novel_gene_candidate"))
    unsupported = is_true(router.get("manual_review_unsupported_locus"))

    if not decision:
        blocking.append("missing_decision")
        return {
            "final_call": "missing_decision",
            "call_group": "blocked",
            "relation_to_current": "unknown",
            "export_readiness": "blocked_or_review",
            "blocking_reasons": ";".join(blocking),
        }
    if not decision_is_valid(decision):
        blocking.extend(str(item) for item in (decision.get("validation_errors") or []) if str(item))
        if not blocking:
            blocking.append("invalid_decision")
        return {
            "final_call": "invalid_decision_blocked",
            "call_group": "blocked",
            "relation_to_current": "unknown",
            "export_readiness": "blocked_or_review",
            "blocking_reasons": ";".join(blocking),
        }
    if decision_value not in {"select_annotation_set", "risk_flagged"}:
        blocking.append("unsupported_decision_value:{0}".format(decision_value or "missing"))
        return {
            "final_call": "invalid_decision_blocked",
            "call_group": "blocked",
            "relation_to_current": "unknown",
            "export_readiness": "blocked_or_review",
            "blocking_reasons": ";".join(blocking),
        }
    if not selected_ai:
        blocking.append("missing_selected_set_id")
        if decision_value == "risk_flagged":
            blocking.append("legacy_risk_flagged_without_provisional_selection")
        return {
            "final_call": "risk_flagged_or_blocked" if decision_value == "risk_flagged" else "invalid_selected_set_blocked",
            "call_group": "human_review_or_blocked" if decision_value == "risk_flagged" else "blocked",
            "relation_to_current": "unresolved_annotation_set_conflict" if decision_value == "risk_flagged" else "unknown",
            "export_readiness": "blocked_or_review",
            "blocking_reasons": ";".join(blocking),
        }
    if selected_ai and not trace_row:
        blocking.append("selected_ai_set_id_not_found_in_trace")
        return {
            "final_call": "invalid_selected_set_blocked",
            "call_group": "blocked",
            "relation_to_current": "unknown",
            "export_readiness": "blocked_or_review",
            "blocking_reasons": ";".join(blocking),
        }
    risk_review = decision_value == "risk_flagged"
    if risk_review:
        blocking.append("risk_flagged_manual_review")
    if novel:
        blocking.append("manual_review_novel_gene_candidate")
        return {
            "final_call": "novel_gene_candidate_manual_review",
            "call_group": "annotation_set_change_pending_validation",
            "relation_to_current": "novel_gene_candidate_no_current_overlap",
            "export_readiness": "manual_review_required_before_export",
            "blocking_reasons": ";".join(blocking),
        }
    if selected_is_current:
        if unsupported:
            blocking.append("manual_review_unsupported_locus")
        return {
            "final_call": "keep_current_flagged" if (unsupported or risk_review) else "keep_current",
            "call_group": "keep_current",
            "relation_to_current": "current_set",
            "export_readiness": "manual_review_required_before_export" if risk_review else "no_change",
            "blocking_reasons": ";".join(blocking),
        }
    if unsupported:
        blocking.append("manual_review_unsupported_locus")
        return {
            "final_call": "annotation_set_selected_needs_manual_review",
            "call_group": "annotation_set_change_pending_validation",
            "relation_to_current": "candidate_set_selected_but_unsupported_locus",
            "export_readiness": "manual_review_required_before_export",
            "blocking_reasons": ";".join(blocking),
        }
    return {
        "final_call": "annotation_set_selected_needs_manual_review" if risk_review else "annotation_set_selected_pending_validation",
        "call_group": "annotation_set_change_pending_validation",
        "relation_to_current": "candidate_set_selected_with_risk_flag" if risk_review else "candidate_set_selected",
        "export_readiness": "manual_review_required_before_export" if risk_review else "pending_validation_before_export",
        "blocking_reasons": ";".join(blocking),
    }


def build_detail_row(cid: str, decision: Dict[str, Any], decision_source: str, trace: Dict[str, Any], router: Dict[str, str]) -> Dict[str, Any]:
    trace_row = selected_trace(decision, trace) if decision else {}
    call = classify_call(decision, trace_row, router)
    selected_real_set = str(trace_row.get("representative_real_set_id", ""))
    selected_models = trace_row.get("representative_real_model_ids") or []
    selected_sources = trace_row.get("member_real_sources") or []
    selected_source = representative_source(selected_real_set, selected_sources)
    return {
        "card_id": cid,
        "locus_id": cid,
        "decision_source": decision_source,
        "profile": decision.get("profile", "") if decision else "",
        "model": decision.get("model", "") if decision else "",
        "api_status": decision.get("api_status", "") if decision else "",
        "decision_valid": str(decision_is_valid(decision)).lower() if decision else "false",
        "validation_errors": as_list_text(decision.get("validation_errors", "") if decision else "missing_decision"),
        "decision": decision.get("decision", "") if decision else "",
        "selected_ai_set_id": decision.get("selected_set_id", "") if decision else "",
        "selected_id": selected_real_set,
        "selected_real_set_id": selected_real_set,
        "selected_source": selected_source,
        "selected_real_sources": ";".join(str(item) for item in selected_sources if str(item)),
        "selected_model_id": as_list_text(selected_models),
        "selected_real_model_ids": as_list_text(selected_models),
        "gene_set_selected_set_id": selected_real_set,
        "member_real_set_ids": as_list_text(trace_row.get("member_real_set_ids") or []),
        "member_real_sources": as_list_text(trace_row.get("member_real_sources") or []),
        "representative_basis": trace_row.get("representative_basis", ""),
        "equivalence_scope": trace_row.get("equivalence_scope", ""),
        "non_decision_difference": trace_row.get("non_decision_difference", ""),
        "final_call": call["final_call"],
        "call_group": call["call_group"],
        "relation_to_current": call["relation_to_current"],
        "export_readiness": call["export_readiness"],
        "confidence": decision.get("confidence", "") if decision else "",
        "requires_validation": str(is_true(decision.get("requires_validation"))).lower() if decision else "false",
        "review_router_class": router.get("router_class", ""),
        "manual_review_novel_gene_candidate": router.get("manual_review_novel_gene_candidate", ""),
        "manual_review_deletion_candidate": router.get("manual_review_deletion_candidate", ""),
        "manual_review_unsupported_locus": router.get("manual_review_unsupported_locus", ""),
        "deletion_allowed": router.get("deletion_allowed", ""),
        "deletion_policy": router.get("deletion_policy", ""),
        "best_evidence_level": router.get("best_evidence_level", ""),
        "discovery_evidence_level": router.get("discovery_evidence_level", ""),
        "review_priority": router.get("review_priority", ""),
        "recommended_manual_action": router.get("recommended_manual_action", ""),
        "routing_decision": router.get("routing_decision", ""),
        "auto_gate": decision.get("auto_gate", router.get("auto_gate", "")) if decision else router.get("auto_gate", ""),
        "auto_risk_level": decision.get("auto_risk_level", "") if decision else "",
        "auto_risk_reason_tags": as_list_text(decision.get("auto_risk_reason_tags", []) if decision else []),
        "reason_tags": as_list_text(decision.get("reason_tags", []) if decision else []),
        "risk_tags": as_list_text(decision.get("risk_tags", []) if decision else []),
        "evidence_summary": decision.get("evidence_summary", "") if decision else "",
        "blocking_reasons": call["blocking_reasons"],
    }


def append_counter(rows: List[Dict[str, Any]], prefix: str, counter: Counter) -> None:
    for key, value in sorted(counter.items()):
        rows.append({"metric": f"{prefix}:{key}", "value": value})


def main() -> None:
    args = parse_args()
    traces = trace_by_card(Path(args.set_trace))
    routers = router_by_card(Path(args.locus_review_router))
    autos = auto_by_card(Path(args.auto_decisions))
    apis = api_by_card(args.api_decisions, args.allow_duplicate_api_decisions)

    all_card_ids = sorted(set(traces) | set(routers) | set(autos) | set(apis))
    detail_rows: List[Dict[str, Any]] = []
    for cid in all_card_ids:
        if cid in autos:
            decision = autos[cid]
            source = "auto"
        else:
            decision = apis.get(cid, {})
            source = "api" if decision else "missing"
        detail_rows.append(build_detail_row(cid, decision, source, traces.get(cid, {}), routers.get(cid, {})))

    summary_rows: List[Dict[str, Any]] = [
        {"metric": "generated_at_utc", "value": utc_now()},
        {"metric": "final_call_rows", "value": len(detail_rows)},
        {"metric": "set_trace_rows", "value": len(traces)},
        {"metric": "auto_decision_rows", "value": len(autos)},
        {"metric": "api_decision_rows", "value": len(apis)},
        {"metric": "router_rows", "value": len(routers)},
    ]
    append_counter(summary_rows, "decision_source", Counter(row["decision_source"] for row in detail_rows))
    append_counter(summary_rows, "decision", Counter(row["decision"] or "missing" for row in detail_rows))
    append_counter(summary_rows, "final_call", Counter(row["final_call"] for row in detail_rows))
    append_counter(summary_rows, "call_group", Counter(row["call_group"] for row in detail_rows))
    append_counter(summary_rows, "review_router_class", Counter(row["review_router_class"] or "missing" for row in detail_rows))
    append_counter(summary_rows, "manual_review_novel_gene_candidate", Counter(row["manual_review_novel_gene_candidate"] or "missing" for row in detail_rows))
    append_counter(summary_rows, "manual_review_unsupported_locus", Counter(row["manual_review_unsupported_locus"] or "missing" for row in detail_rows))
    append_counter(summary_rows, "deletion_allowed", Counter(row["deletion_allowed"] or "missing" for row in detail_rows))

    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)
    detail_path = out_dir / "final_annotation_calls.tsv"
    summary_path = out_dir / "final_annotation_call_summary.tsv"
    report_path = out_dir / "final_annotation_call_report.md"
    write_tsv(detail_path, DETAIL_FIELDS, detail_rows)
    write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    write_report(report_path, summary_rows)
    print(f"final_call_rows\t{len(detail_rows)}")
    print(f"detail\t{detail_path}")
    print(f"summary\t{summary_path}")


if __name__ == "__main__":
    main()
