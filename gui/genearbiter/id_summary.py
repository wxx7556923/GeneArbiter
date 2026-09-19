#!/usr/bin/env python3
# script_id_md5: 4783ad45747c4553cd9fd87ce125b0bb
# created: 2026-07-08
# modified: 2026-07-08
# owner: project
# status: project_code
# purpose: 汇总 GeneArbiter correction/arbitration 输出中的 current replacement 和 source provenance ID 映射。
# inputs: full cards JSONL; final_annotation_calls.tsv; optional proposal_catalog_records.tsv; optional gene/transcript mapping TSV。
# outputs: genearbiter_id_mapping_summary.tsv。
# notes: 不读取 truth annotation；correction mode 记录 current replacement，arbitration mode 记录 selected source provenance。

"""Build a human-readable ID replacement/provenance summary table."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence


FIELDS = [
    "workflow_mode",
    "locus_id",
    "change_class",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "current_gene_ids",
    "retired_current_gene_ids",
    "reuse_current_gene_ids",
    "final_gene_ids",
    "final_transcript_ids",
    "proposal_gene_ids",
    "proposal_transcript_ids",
    "selected_set_id",
    "selected_source",
    "selected_model_ids",
    "source_sources",
    "source_gene_ids",
    "source_transcript_ids",
    "source_model_ids",
    "replacement_interpretation",
    "review_flags",
    "risk_tags",
    "blocking_reasons",
    "notes",
]


def split_list(text: Any) -> List[str]:
    return [item.strip() for item in str(text or "").replace(";", ",").split(",") if item.strip()]


def public_proposal_id(value: Any) -> str:
    text = str(value or "")
    if text.startswith("V2PROP_"):
        return "GeneArbiterProposal_" + text[len("V2PROP_") :]
    return text


def join_items(items: Iterable[Any]) -> str:
    seen = set()
    out = []
    for item in items:
        text = str(item or "").strip()
        if text and text not in seen:
            seen.add(text)
            out.append(text)
    return ";".join(out)


def read_tsv(path: str) -> List[Dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_jsonl(path: str) -> List[Dict[str, Any]]:
    rows = []
    if not path or not Path(path).exists():
        return rows
    with open(path, "r") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def write_tsv(path: str, rows: Sequence[Mapping[str, Any]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def model_source_info(model: Mapping[str, Any]) -> Dict[str, str]:
    model_id = str(model.get("model_id", ""))
    parts = model_id.split(":", 2)
    source = str(model.get("source") or (parts[0] if parts else ""))
    gene_id = str(model.get("source_gene_id") or model.get("gene_id") or (parts[1] if len(parts) > 1 else ""))
    tx_id = str(model.get("source_transcript_id") or model.get("transcript_id") or (parts[2] if len(parts) > 2 else ""))
    return {"source": source, "source_gene_id": gene_id, "source_transcript_id": tx_id, "source_model_id": model_id}


def card_indexes(full_cards: Sequence[Mapping[str, Any]]) -> tuple[Dict[str, Mapping[str, Any]], Dict[str, Dict[str, Mapping[str, Any]]]]:
    by_locus = {}
    models_by_locus = {}
    for card in full_cards:
        locus_id = str(card.get("locus_id") or card.get("card_id") or "")
        if not locus_id:
            continue
        by_locus[locus_id] = card
        models_by_locus[locus_id] = {str(model.get("model_id", "")): model for model in card.get("candidate_models", []) if model.get("model_id")}
    return by_locus, models_by_locus


def rows_by_key(rows: Sequence[Mapping[str, str]], key: str) -> Dict[str, List[Mapping[str, str]]]:
    out: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        out[str(row.get(key, ""))].append(row)
    return out


def selected_model_ids(call: Mapping[str, str]) -> List[str]:
    return split_list(call.get("selected_real_model_ids") or call.get("selected_model_id"))


def source_info_for_call(call: Mapping[str, str], model_index: Mapping[str, Mapping[str, Any]]) -> List[Dict[str, str]]:
    infos = []
    for model_id in selected_model_ids(call):
        model = model_index.get(model_id)
        if model:
            infos.append(model_source_info(model))
        else:
            parts = model_id.split(":", 2)
            infos.append(
                {
                    "source": parts[0] if parts else "",
                    "source_gene_id": parts[1] if len(parts) > 1 else "",
                    "source_transcript_id": parts[2] if len(parts) > 2 else "",
                    "source_model_id": model_id,
                }
            )
    return infos


def infer_change_class(mode: str, call: Mapping[str, str], prop_rows: Sequence[Mapping[str, str]], gene_rows: Sequence[Mapping[str, str]]) -> str:
    if mode == "arbitration":
        return "no_current_arbitration_selection"
    relations = [row.get("final_relation", "") for row in gene_rows if row.get("final_relation")]
    if relations:
        return join_items(relations)
    actions = [row.get("proposal_action", "") for row in prop_rows if row.get("proposal_action")]
    if actions:
        action = join_items(actions)
        if action == "keep_current":
            return "keep_current"
        if action == "add_novel_selected_set":
            return "novel_gene"
        if action == "replace_current_with_selected_set":
            return "replacement_without_final_gene_mapping"
        return action
    final_call = str(call.get("final_call", ""))
    if final_call.startswith("keep_current"):
        return "keep_current"
    return str(call.get("relation_to_current") or final_call or "unknown")


def interpretation(mode: str, change_class: str) -> str:
    if mode == "arbitration":
        return "No base/current GFF was used; final IDs trace selected source models but do not replace current IDs."
    if change_class == "keep_current":
        return "Current gene model is retained."
    if change_class == "one_to_one_replacement":
        return "One current gene is replaced by one proposal gene; current gene ID may be reused according to policy."
    if change_class == "split_replacement":
        return "One current gene is split into multiple final proposal genes."
    if change_class == "merge_replacement":
        return "Multiple current genes are merged into one final proposal gene."
    if change_class == "complex_replacement":
        return "Many-to-many or complex current/proposal replacement."
    if change_class == "novel_gene":
        return "Selected proposal is a novel/tool-only gene candidate with no current gene replacement."
    return "See mapping/proposal records for detailed provenance."


def build_summary(
    mode: str,
    full_cards: Sequence[Mapping[str, Any]],
    final_calls: Sequence[Mapping[str, str]],
    proposal_records: Sequence[Mapping[str, str]],
    gene_mapping: Sequence[Mapping[str, str]],
    tx_mapping: Sequence[Mapping[str, str]],
) -> List[Dict[str, str]]:
    _cards_by_locus, models_by_locus = card_indexes(full_cards)
    prop_by_locus = rows_by_key(proposal_records, "locus_id")
    gene_by_locus = rows_by_key(gene_mapping, "locus_id")
    tx_by_locus = rows_by_key(tx_mapping, "locus_id")
    rows = []
    for call in final_calls:
        locus_id = str(call.get("locus_id") or call.get("card_id") or "")
        prop_rows = prop_by_locus.get(locus_id, [])
        gene_rows = gene_by_locus.get(locus_id, [])
        tx_rows = tx_by_locus.get(locus_id, [])
        infos = source_info_for_call(call, models_by_locus.get(locus_id, {}))
        change = infer_change_class(mode, call, prop_rows, gene_rows)
        current_gene_ids = join_items(
            item
            for row in list(prop_rows) + list(gene_rows)
            for item in split_list(row.get("current_gene_ids") or row.get("removed_current_gene_ids"))
        )
        retired_current_gene_ids = join_items(item for row in gene_rows for item in split_list(row.get("retired_current_gene_ids")))
        final_gene_ids = join_items(row.get("final_gene_id", "") for row in gene_rows)
        if not final_gene_ids and change == "keep_current":
            final_gene_ids = current_gene_ids
        rows.append(
            {
                "workflow_mode": mode,
                "locus_id": locus_id,
                "change_class": change,
                "final_call": call.get("final_call", ""),
                "call_group": call.get("call_group", ""),
                "relation_to_current": call.get("relation_to_current", ""),
                "export_readiness": call.get("export_readiness", ""),
                "current_gene_ids": current_gene_ids,
                "retired_current_gene_ids": retired_current_gene_ids,
                "reuse_current_gene_ids": join_items(row.get("final_gene_id", "") for row in gene_rows if row.get("reuse_current_gene_id") == "true"),
                "final_gene_ids": final_gene_ids,
                "final_transcript_ids": join_items(row.get("final_transcript_id", "") for row in tx_rows),
                "proposal_gene_ids": join_items(public_proposal_id(row.get("proposal_gene_id", "") or row.get("v2_proposal_gene_id", "") or row.get("export_gene_id", "")) for row in list(gene_rows) + list(prop_rows)),
                "proposal_transcript_ids": join_items(public_proposal_id(row.get("proposal_transcript_id", "") or row.get("v2_proposal_transcript_id", "") or row.get("export_transcript_id", "")) for row in list(tx_rows) + list(prop_rows)),
                "selected_set_id": call.get("selected_real_set_id") or call.get("gene_set_selected_set_id") or call.get("selected_set_id", ""),
                "selected_source": call.get("selected_source", ""),
                "selected_model_ids": join_items(selected_model_ids(call)),
                "source_sources": join_items(info["source"] for info in infos) or join_items(row.get("source", "") for row in list(prop_rows) + list(gene_rows)),
                "source_gene_ids": join_items(info["source_gene_id"] for info in infos) or join_items(row.get("source_gene_id", "") for row in list(prop_rows) + list(gene_rows)),
                "source_transcript_ids": join_items(info["source_transcript_id"] for info in infos) or join_items(row.get("source_transcript_id", "") for row in list(prop_rows) + list(tx_rows)),
                "source_model_ids": join_items(info["source_model_id"] for info in infos) or join_items(row.get("source_model_id", "") for row in list(prop_rows) + list(gene_rows)),
                "replacement_interpretation": interpretation(mode, change),
                "review_flags": join_items(row.get("review_flags", "") for row in prop_rows),
                "risk_tags": join_items([call.get("risk_tags", "")] + [row.get("risk_tags", "") for row in prop_rows]),
                "blocking_reasons": join_items([call.get("blocking_reasons", "")] + [row.get("blocking_reasons", "") for row in prop_rows]),
                "notes": "one row per locus; semicolon-separated fields may contain split/merge/provenance multiplicity",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-mode", choices=["correction", "arbitration"], required=True)
    parser.add_argument("--full-cards", required=True)
    parser.add_argument("--final-calls", required=True)
    parser.add_argument("--proposal-records", default="")
    parser.add_argument("--gene-id-mapping", default="")
    parser.add_argument("--transcript-id-mapping", default="")
    parser.add_argument("--out-tsv", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_summary(
        args.workflow_mode,
        read_jsonl(args.full_cards),
        read_tsv(args.final_calls),
        read_tsv(args.proposal_records),
        read_tsv(args.gene_id_mapping),
        read_tsv(args.transcript_id_mapping),
    )
    write_tsv(args.out_tsv, rows)
    print("id_mapping_summary\t{0}".format(args.out_tsv))
    print("rows\t{0}".format(len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
