#!/usr/bin/env python3
# script_id_md5: add87160922af22f9f1ac8193820774f
# created: 2026-07-10
# modified: 2026-07-10
# owner: project
# status: project_code
# purpose: 从 GeneArbiter final_calls/full_cards 生成供人工审核的轻量 review bundle。
# inputs: final_annotation_calls.tsv; model_arbitration_full_cards.jsonl。
# outputs: review_queue.tsv; manual_review_decisions.template.tsv; per-locus Markdown review cards。
# notes: 第一版只生成审核材料，不把人工审核结果接回 final GFF。

"""Build a lightweight manual-review bundle from GeneArbiter final calls."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


QUEUE_FIELDS = [
    "review_bucket",
    "review_rank",
    "locus_id",
    "seqid",
    "start",
    "end",
    "igv_region",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "selected_set_id",
    "selected_source",
    "selected_real_sources",
    "selected_model_id",
    "selected_real_model_ids",
    "confidence",
    "requires_validation",
    "review_router_class",
    "best_evidence_level",
    "discovery_evidence_level",
    "risk_tags",
    "reason_tags",
    "blocking_reasons",
    "why_review",
    "recommended_manual_action",
    "review_card",
]

TEMPLATE_FIELDS = ["locus_id", "decision", "selected_set_id", "comment", "reviewer"]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_jsonl(path: Path) -> Iterable[Dict[str, object]]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def load_card_index(path: Path) -> Dict[str, Dict[str, object]]:
    if not path.exists():
        return {}
    cards: Dict[str, Dict[str, object]] = {}
    for card in read_jsonl(path):
        locus_id = str(card.get("locus_id") or "")
        if locus_id:
            cards[locus_id] = card
    return cards


def split_values(text: object) -> List[str]:
    values: List[str] = []
    for item in str(text or "").replace(";", ",").split(","):
        item = item.strip()
        if item:
            values.append(item)
    return values


def truthy(text: object) -> bool:
    return str(text or "").strip().lower() in {"1", "true", "yes", "y"}


def card_locus_fields(card: Mapping[str, object]) -> Dict[str, object]:
    locus = card.get("locus") if isinstance(card.get("locus"), Mapping) else {}
    complexity = card.get("complexity") if isinstance(card.get("complexity"), Mapping) else {}
    seqid = str(locus.get("seqid") or complexity.get("seqid") or "")
    start = locus.get("start") or complexity.get("start") or ""
    end = locus.get("end") or complexity.get("end") or ""
    try:
        start_i = int(start)
        end_i = int(end)
        igv = f"{seqid}:{max(1, start_i - 500)}-{end_i + 500}" if seqid else ""
    except (TypeError, ValueError):
        igv = ""
    return {"seqid": seqid, "start": start, "end": end, "igv_region": igv}


def classify_review(row: Mapping[str, str]) -> tuple[str, int, str]:
    reasons: List[str] = []
    manual_flags = [
        "manual_review_novel_gene_candidate",
        "manual_review_deletion_candidate",
        "manual_review_unsupported_locus",
    ]
    if any(truthy(row.get(field)) for field in manual_flags):
        reasons.append("manual_review_router_flag")
    if split_values(row.get("blocking_reasons")):
        reasons.append("blocking_reasons_present")
    export_readiness = str(row.get("export_readiness") or "").lower()
    if export_readiness in {"blocked", "manual_review", "manual_review_required", "not_ready", "not_exportable"}:
        reasons.append(f"export_readiness={export_readiness}")
    final_call = str(row.get("final_call") or "").lower()
    call_group = str(row.get("call_group") or "").lower()
    if "blocked" in final_call or "blocked" in call_group:
        reasons.append("blocked_final_call")
    if reasons:
        return "manual_required", 0, ";".join(reasons)

    risk_tags = split_values(row.get("risk_tags"))
    reason_tags = split_values(row.get("reason_tags"))
    selected_source = str(row.get("selected_source") or "")
    relation = str(row.get("relation_to_current") or "").lower()
    if risk_tags:
        reasons.append("risk_tags_present")
    if reason_tags:
        reasons.append("reason_tags_present")
    if selected_source and selected_source != "current":
        reasons.append("selected_non_current_source")
    if relation and relation not in {"current_set", "no_change", "keep_current"}:
        reasons.append(f"relation_to_current={relation}")
    if str(row.get("confidence") or "").lower() == "low":
        reasons.append("low_confidence")
    if reasons:
        return "review_recommended", 1, ";".join(reasons)

    if truthy(row.get("requires_validation")) or str(row.get("confidence") or "").lower() == "medium":
        return "optional_review", 2, "requires_validation_or_medium_confidence"
    return "auto_accept", 3, "low_review_priority"


def make_markdown(row: Mapping[str, object]) -> str:
    lines = [
        f"# {row.get('locus_id', '')}",
        "",
        "## Review Summary",
        f"- Review bucket: {row.get('review_bucket', '')}",
        f"- Why review: {row.get('why_review', '')}",
        f"- Recommended manual action: {row.get('recommended_manual_action', '')}",
        "",
        "## Locus",
        f"- Region: {row.get('igv_region', '')}",
        f"- Coordinates: {row.get('seqid', '')}:{row.get('start', '')}-{row.get('end', '')}",
        "",
        "## AI Call",
        f"- Final call: {row.get('final_call', '')}",
        f"- Call group: {row.get('call_group', '')}",
        f"- Relation to current: {row.get('relation_to_current', '')}",
        f"- Export readiness: {row.get('export_readiness', '')}",
        f"- Confidence: {row.get('confidence', '')}",
        "",
        "## Selected Set",
        f"- Selected set ID: {row.get('selected_set_id', '')}",
        f"- Selected source: {row.get('selected_source', '')}",
        f"- Real sources: {row.get('selected_real_sources', '')}",
        f"- Selected model ID: {row.get('selected_model_id', '')}",
        f"- Real model IDs: {row.get('selected_real_model_ids', '')}",
        "",
        "## Evidence And Risk",
        f"- Best evidence level: {row.get('best_evidence_level', '')}",
        f"- Discovery evidence level: {row.get('discovery_evidence_level', '')}",
        f"- Risk tags: {row.get('risk_tags', '')}",
        f"- Reason tags: {row.get('reason_tags', '')}",
        f"- Blocking reasons: {row.get('blocking_reasons', '')}",
        "",
        "## Manual Decision Template",
        "",
        "```text",
        "decision: accept_ai | keep_current | choose_alternative | manual_block",
        f"locus_id: {row.get('locus_id', '')}",
        f"selected_set_id: {row.get('selected_set_id', '')}",
        "comment:",
        "reviewer:",
        "```",
    ]
    return "\n".join(lines) + "\n"


def build_review_bundle(args: argparse.Namespace) -> int:
    final_calls = Path(args.final_calls).expanduser().resolve()
    full_cards = Path(args.full_cards).expanduser().resolve() if args.full_cards else Path("")
    out_dir = Path(args.out_dir).expanduser().resolve()
    cards_dir = out_dir / "review_cards"
    rows = read_tsv(final_calls)
    card_index = load_card_index(full_cards) if full_cards else {}

    queue_rows: List[Dict[str, object]] = []
    for row in rows:
        locus_id = row.get("locus_id") or row.get("card_id") or ""
        card = card_index.get(locus_id, {})
        locus_fields = card_locus_fields(card) if card else {"seqid": "", "start": "", "end": "", "igv_region": ""}
        bucket, rank, why = classify_review(row)
        card_name = f"{locus_id}.md" if locus_id else "unknown_locus.md"
        selected_set = row.get("selected_real_set_id") or row.get("selected_set_id") or row.get("gene_set_selected_set_id") or ""
        queue_row: Dict[str, object] = {
            **{field: row.get(field, "") for field in QUEUE_FIELDS},
            **locus_fields,
            "review_bucket": bucket,
            "review_rank": rank,
            "locus_id": locus_id,
            "selected_set_id": selected_set,
            "why_review": why,
            "review_card": str(Path("review_cards") / card_name),
        }
        queue_rows.append(queue_row)
        cards_dir.mkdir(parents=True, exist_ok=True)
        (cards_dir / card_name).write_text(make_markdown(queue_row), encoding="utf-8")

    queue_rows.sort(key=lambda item: (int(item.get("review_rank", 9)), str(item.get("seqid", "")), int(item.get("start") or 0), str(item.get("locus_id", ""))))
    write_tsv(out_dir / "review_queue.tsv", queue_rows, QUEUE_FIELDS)
    template_rows = [
        {
            "locus_id": row.get("locus_id", ""),
            "decision": "accept_ai",
            "selected_set_id": row.get("selected_set_id", ""),
            "comment": "",
            "reviewer": "",
        }
        for row in queue_rows
        if row.get("review_bucket") in {"manual_required", "review_recommended"}
    ]
    write_tsv(out_dir / "manual_review_decisions.template.tsv", template_rows, TEMPLATE_FIELDS)
    summary = {
        "final_calls": str(final_calls),
        "full_cards": str(full_cards) if full_cards else "",
        "review_queue": str(out_dir / "review_queue.tsv"),
        "manual_review_decisions_template": str(out_dir / "manual_review_decisions.template.tsv"),
        "review_cards_dir": str(cards_dir),
        "locus_count": len(queue_rows),
        "bucket_counts": {},
        "note": "This bundle is for human review only. apply-review is not implemented in this first version.",
    }
    for row in queue_rows:
        bucket = str(row.get("review_bucket", ""))
        summary["bucket_counts"][bucket] = int(summary["bucket_counts"].get(bucket, 0)) + 1
    (out_dir / "review_bundle_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"review_queue\t{out_dir / 'review_queue.tsv'}")
    print(f"manual_review_decisions_template\t{out_dir / 'manual_review_decisions.template.tsv'}")
    print(f"review_cards_dir\t{cards_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-calls", required=True)
    parser.add_argument("--full-cards", default="")
    parser.add_argument("--out-dir", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    return build_review_bundle(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
