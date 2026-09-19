#!/usr/bin/env python3
# script_id_md5: 3dcebab33be256a5406c06f1fcc8a1e0
# created: 2026-07-05
# modified: 2026-07-05
# owner: project
# status: project_code
# purpose: 导出 GeneArbiter annotation-set proposal catalog 的完整 GFF3 草案。
# inputs: current/base GFF3；full cards JSONL；final_annotation_calls.tsv。
# outputs: proposal_catalog.gff3；proposal_catalog_records.tsv；proposal_catalog_review_flags.tsv；proposal_catalog_summary.tsv；proposal_catalog_report.md。
# notes: 不读取 truth annotation；不生成新坐标；输出是 proposal catalog，不是 final annotation。

"""Export a complete GeneArbiter annotation-set proposal catalog GFF3.

The exporter keeps the current annotation as the whole-genome backbone. For each
locus with a selected non-current annotation set, it removes the current gene(s)
represented in that card and writes the selected candidate model(s) as draft
proposal gene records. Novel/tool-only candidates are appended without deleting
current genes. Deletion is never applied here; deletion-like or unsupported loci
are carried as review flags.

This script does not read final-evaluation truth and does not invent coordinates.
All proposal coordinates are copied from selected models in full cards; CDS
phase is recomputed only for GFF3 consistency.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple
from urllib.parse import quote

from gff_utils import ensure_dir, iter_gff_rows, open_text, parse_attrs, write_tsv


RECORD_FIELDS = [
    "locus_id",
    "proposal_action",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "selected_set_id",
    "selected_source",
    "selected_model_id",
    "current_gene_ids",
    "removed_current_gene_ids",
    "export_gene_id",
    "export_transcript_id",
    "source_model_id",
    "source",
    "source_gene_id",
    "source_transcript_id",
    "seqid",
    "start",
    "end",
    "strand",
    "exon_count",
    "cds_count",
    "cds_length",
    "review_flags",
    "risk_tags",
    "blocking_reasons",
    "notes",
]

REVIEW_FIELDS = [
    "locus_id",
    "proposal_action",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "selected_set_id",
    "selected_source",
    "current_gene_ids",
    "review_flags",
    "risk_tags",
    "blocking_reasons",
    "recommended_manual_action",
]

SUMMARY_FIELDS = ["metric", "value"]

FEATURE_RANK = {
    "gene": 0,
    "mrna": 1,
    "transcript": 1,
    "rna": 1,
    "exon": 2,
    "cds": 3,
}

TRUTH_LIKE_RE = re.compile(r"truth|manual", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-gff", required=True, help="Current/base annotation GFF3 or GFF3.gz used as whole-genome backbone.")
    parser.add_argument("--full-cards", required=True, help="GeneArbiter full decision cards JSONL.")
    parser.add_argument("--final-calls", required=True, help="GeneArbiter final annotation calls TSV.")
    parser.add_argument("--out-dir", required=True, help="Output directory for proposal_catalog.* files.")
    parser.add_argument("--source-label", default="GeneArbiter", help="GFF source column for newly exported proposal records.")
    parser.add_argument("--sort-records", action="store_true", help="Sort output records by seqid/start instead of preserving current GFF order and appending proposals.")
    parser.add_argument("--fail-on-truth-path", action="store_true", help="Fail if an input path looks like a truth/manual evaluation file.")
    return parser.parse_args()


def read_jsonl(path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(path, "r") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError("{0}:{1}: invalid JSON: {2}".format(path, line_no, exc)) from exc
    return rows


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_list(text: Any) -> List[str]:
    return [item.strip() for item in str(text or "").replace(";", ",").split(",") if item.strip()]


def safe_token(text: Any) -> str:
    token = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(text or ""))
    return token or "unknown"


def attr_value(text: Any) -> str:
    return quote(str(text or ""), safe="._:-|,")


def attrs(items: Iterable[Tuple[str, Any]]) -> str:
    return ";".join("{0}={1}".format(key, attr_value(value)) for key, value in items if value not in {"", None})


def as_int(value: Any) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def truth_path_guard(paths: Sequence[str]) -> None:
    for path in paths:
        if TRUTH_LIKE_RE.search(path):
            raise SystemExit("Refusing truth/manual-like input path for proposal export: {0}".format(path))


def parse_model_id(model_id: Any) -> Tuple[str, str, str]:
    parts = str(model_id or "").split(":", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    if len(parts) == 2:
        return parts[0], parts[1], parts[1]
    return "", str(model_id or ""), str(model_id or "")


def model_source_info(model: Dict[str, Any]) -> Tuple[str, str, str]:
    source, gene_id, tx_id = parse_model_id(model.get("model_id", ""))
    return (
        str(model.get("source") or source),
        str(model.get("source_gene_id") or model.get("gene_id") or gene_id),
        str(model.get("source_transcript_id") or model.get("transcript_id") or tx_id),
    )


def model_cds_length(model: Dict[str, Any]) -> int:
    total = 0
    for block in model.get("cds", []):
        if len(block) >= 2:
            total += as_int(block[1]) - as_int(block[0]) + 1
    return total


def recalc_cds_phase(cds: Sequence[Sequence[Any]], strand: str) -> List[Tuple[int, int, str]]:
    blocks = [(as_int(block[0]), as_int(block[1])) for block in cds if len(block) >= 2]
    ordered = sorted(blocks, key=lambda item: (item[0], item[1]), reverse=(strand == "-"))
    phased: List[Tuple[int, int, str]] = []
    cumulative = 0
    for index, (start, end) in enumerate(ordered):
        phase = 0 if index == 0 else (3 - (cumulative % 3)) % 3
        phased.append((start, end, str(phase)))
        cumulative += end - start + 1
    return sorted(phased, key=lambda item: (item[0], item[1], item[2]))


def gff_line(seqid: Any, source: str, feature: str, start: Any, end: Any, strand: Any, phase: Any, attr_text: str) -> str:
    return "\t".join([str(seqid), source, feature, str(start), str(end), ".", str(strand or "."), str(phase or "."), attr_text])


def collect_current_gene_maps(current_gff: str) -> Tuple[Dict[str, Set[str]], Dict[str, str], Set[str]]:
    gene_to_tx: Dict[str, Set[str]] = defaultdict(set)
    tx_to_gene: Dict[str, str] = {}
    gene_ids: Set[str] = set()
    for _seqid, _source, feature, _start, _end, _score, _strand, _phase, row_attrs in iter_gff_rows(current_gff):
        feature_l = feature.lower()
        if feature_l == "gene":
            gene_id = row_attrs.get("ID") or row_attrs.get("Name")
            if gene_id:
                gene_ids.add(gene_id)
        elif feature_l in {"mrna", "transcript", "rna"}:
            tx_id = row_attrs.get("ID") or row_attrs.get("Name")
            parents = split_list(row_attrs.get("Parent", ""))
            if tx_id and parents:
                gene_id = parents[0]
                tx_to_gene[tx_id] = gene_id
                gene_to_tx[gene_id].add(tx_id)
                gene_ids.add(gene_id)
    return gene_to_tx, tx_to_gene, gene_ids


def current_gene_ids_from_card(card: Dict[str, Any]) -> List[str]:
    genes: Set[str] = set()
    for model in card.get("candidate_models", []):
        source, gene_id, _tx_id = model_source_info(model)
        if source == "current" and gene_id:
            genes.add(gene_id)
    return sorted(genes)


def model_ids_from_call(call: Dict[str, str]) -> List[str]:
    return split_list(call.get("selected_model_id") or call.get("selected_real_model_ids") or "")


def selected_set_id_from_call(call: Dict[str, str]) -> str:
    return call.get("selected_real_set_id") or call.get("gene_set_selected_set_id") or call.get("selected_id") or ""


def review_flags_from_call(call: Dict[str, str]) -> List[str]:
    flags: List[str] = []
    final_call = str(call.get("final_call", ""))
    call_group = str(call.get("call_group", ""))
    decision = str(call.get("decision", ""))
    export_readiness = str(call.get("export_readiness", ""))
    for field in [
        "manual_review_novel_gene_candidate",
        "manual_review_deletion_candidate",
        "manual_review_unsupported_locus",
    ]:
        if str(call.get(field, "")).lower() == "true":
            flags.append(field)
    if "manual_review" in final_call or "manual_review" in call_group:
        flags.append("manual_review_required")
    if final_call == "keep_current_flagged" or decision == "risk_flagged":
        flags.append("risk_flagged")
    if "blocked" in call_group or call.get("blocking_reasons"):
        flags.append("blocking_reason_present")
    if export_readiness == "manual_review_required_before_export":
        flags.append(export_readiness)
    if final_call.startswith("delete") or str(call.get("manual_review_deletion_candidate", "")).lower() == "true":
        flags.append("deletion_not_applied_by_proposal_exporter")
    return sorted(set(flag for flag in flags if flag))


def group_models_by_source_gene(models: Sequence[Dict[str, Any]]) -> List[Tuple[str, List[Dict[str, Any]]]]:
    grouped: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = defaultdict(list)
    for model in models:
        _source, gene_id, _tx_id = model_source_info(model)
        key = (gene_id or str(model.get("model_id", "")), str(model.get("seqid", "")), str(model.get("strand", ".")))
        grouped[key].append(model)
    return [(key[0], rows) for key, rows in sorted(grouped.items(), key=lambda item: (min(as_int(m.get("start")) for m in item[1]), item[0]))]


def build_indexes(cards: Sequence[Dict[str, Any]], calls: Sequence[Dict[str, str]]) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, str]]]:
    cards_by_locus = {str(card.get("locus_id", "")): card for card in cards if card.get("locus_id")}
    calls_by_locus = {str(row.get("locus_id", "")): row for row in calls if row.get("locus_id")}
    return cards_by_locus, calls_by_locus


def determine_action(call: Dict[str, str], selected_models: Sequence[Dict[str, Any]], current_gene_ids: Sequence[str]) -> str:
    final_call = call.get("final_call", "")
    selected_sources = {model_source_info(model)[0] for model in selected_models}
    if final_call in {"keep_current", "keep_current_flagged"}:
        return "keep_current"
    if str(final_call).startswith("delete"):
        return "keep_current_deletion_not_applied"
    if selected_models and selected_sources and selected_sources <= {"current"}:
        return "keep_current"
    if not selected_models:
        return "keep_current_no_selected_model"
    if current_gene_ids:
        return "replace_current_with_selected_set"
    return "add_novel_selected_set"


def proposal_lines_for_call(
    card: Dict[str, Any],
    call: Dict[str, str],
    selected_models: Sequence[Dict[str, Any]],
    action: str,
    current_gene_ids: Sequence[str],
    source_label: str,
) -> Tuple[List[str], List[Dict[str, Any]]]:
    locus_id = str(card.get("locus_id", ""))
    selected_set_id = selected_set_id_from_call(call)
    selected_source = call.get("selected_source", "")
    relation = call.get("relation_to_current", "")
    review_flags = ";".join(review_flags_from_call(call))
    risk_tags = call.get("risk_tags", "")
    blocking_reasons = call.get("blocking_reasons", "")
    removed = ";".join(current_gene_ids) if action == "replace_current_with_selected_set" else ""
    lines: List[str] = []
    records: List[Dict[str, Any]] = []

    for gene_index, (_group_gene_id, gene_models) in enumerate(group_models_by_source_gene(selected_models), start=1):
        seqids = sorted({str(model.get("seqid", "")) for model in gene_models if model.get("seqid", "")})
        strands = sorted({str(model.get("strand", ".")) for model in gene_models})
        seqid = seqids[0] if seqids else ""
        strand = strands[0] if strands else "."
        starts = [as_int(model.get("start")) for model in gene_models if as_int(model.get("start"))]
        ends = [as_int(model.get("end")) for model in gene_models if as_int(model.get("end"))]
        for model in gene_models:
            starts.extend(as_int(block[0]) for block in model.get("exons", []) if len(block) >= 2)
            ends.extend(as_int(block[1]) for block in model.get("exons", []) if len(block) >= 2)
            starts.extend(as_int(block[0]) for block in model.get("cds", []) if len(block) >= 2)
            ends.extend(as_int(block[1]) for block in model.get("cds", []) if len(block) >= 2)
        start = min(starts) if starts else 0
        end = max(ends) if ends else 0
        source_gene_ids = sorted({model_source_info(model)[1] for model in gene_models if model_source_info(model)[1]})
        source_model_ids = [str(model.get("model_id", "")) for model in gene_models]
        export_gene_id = "GeneArbiterProposal_{0}_{1}_g{2:02d}".format(safe_token(locus_id), safe_token(selected_set_id), gene_index)
        common_gene_attrs = [
            ("ID", export_gene_id),
            ("Name", export_gene_id),
            ("proposal_catalog", "true"),
            ("proposal_status", "not_final_annotation"),
            ("locus_id", locus_id),
            ("proposal_action", action),
            ("selected_set_id", selected_set_id),
            ("selected_source", selected_source),
            ("source_gene_id", ",".join(source_gene_ids)),
            ("source_model_id", ",".join(source_model_ids)),
            ("replaces_current_genes", removed),
            ("final_call", call.get("final_call", "")),
            ("call_group", call.get("call_group", "")),
            ("relation_to_current", relation),
            ("export_readiness", call.get("export_readiness", "")),
            ("review_flags", review_flags),
            ("risk_tags", risk_tags),
            ("blocking_reasons", blocking_reasons),
        ]
        lines.append(gff_line(seqid, source_label, "gene", start, end, strand, ".", attrs(common_gene_attrs)))

        for tx_index, model in enumerate(sorted(gene_models, key=lambda item: (as_int(item.get("start")), str(item.get("model_id", "")))), start=1):
            source, source_gene_id, source_tx_id = model_source_info(model)
            model_id = str(model.get("model_id", ""))
            exons = [(as_int(block[0]), as_int(block[1])) for block in model.get("exons", []) if len(block) >= 2]
            cds = recalc_cds_phase(model.get("cds", []), str(model.get("strand", strand)))
            if not exons and cds:
                exons = [(start_i, end_i) for start_i, end_i, _phase in cds]
            tx_starts = [as_int(model.get("start"))] + [start_i for start_i, _end_i in exons] + [start_i for start_i, _end_i, _phase in cds]
            tx_ends = [as_int(model.get("end"))] + [end_i for _start_i, end_i in exons] + [end_i for _start_i, end_i, _phase in cds]
            tx_starts = [item for item in tx_starts if item]
            tx_ends = [item for item in tx_ends if item]
            tx_start = min(tx_starts) if tx_starts else start
            tx_end = max(tx_ends) if tx_ends else end
            export_tx_id = "{0}.t{1:02d}".format(export_gene_id, tx_index)
            tx_attrs = [
                ("ID", export_tx_id),
                ("Parent", export_gene_id),
                ("Name", export_tx_id),
                ("proposal_catalog", "true"),
                ("proposal_status", "not_final_annotation"),
                ("locus_id", locus_id),
                ("selected_set_id", selected_set_id),
                ("source_model_id", model_id),
                ("source", source),
                ("source_gene_id", source_gene_id),
                ("source_transcript_id", source_tx_id),
                ("review_flags", review_flags),
            ]
            lines.append(gff_line(model.get("seqid", seqid), source_label, "mRNA", tx_start, tx_end, model.get("strand", strand), ".", attrs(tx_attrs)))
            for exon_index, (exon_start, exon_end) in enumerate(sorted(set(exons)), start=1):
                lines.append(
                    gff_line(
                        model.get("seqid", seqid),
                        source_label,
                        "exon",
                        exon_start,
                        exon_end,
                        model.get("strand", strand),
                        ".",
                        attrs([("ID", "{0}.exon{1}".format(export_tx_id, exon_index)), ("Parent", export_tx_id)]),
                    )
                )
            for cds_index, (cds_start, cds_end, phase) in enumerate(cds, start=1):
                lines.append(
                    gff_line(
                        model.get("seqid", seqid),
                        source_label,
                        "CDS",
                        cds_start,
                        cds_end,
                        model.get("strand", strand),
                        phase,
                        attrs([("ID", "{0}.cds{1}".format(export_tx_id, cds_index)), ("Parent", export_tx_id)]),
                    )
                )
            records.append(
                {
                    "locus_id": locus_id,
                    "proposal_action": action,
                    "final_call": call.get("final_call", ""),
                    "call_group": call.get("call_group", ""),
                    "relation_to_current": relation,
                    "export_readiness": call.get("export_readiness", ""),
                    "selected_set_id": selected_set_id,
                    "selected_source": selected_source,
                    "selected_model_id": call.get("selected_model_id", ""),
                    "current_gene_ids": ";".join(current_gene_ids),
                    "removed_current_gene_ids": removed,
                    "export_gene_id": export_gene_id,
                    "export_transcript_id": export_tx_id,
                    "source_model_id": model_id,
                    "source": source,
                    "source_gene_id": source_gene_id,
                    "source_transcript_id": source_tx_id,
                    "seqid": model.get("seqid", seqid),
                    "start": tx_start,
                    "end": tx_end,
                    "strand": model.get("strand", strand),
                    "exon_count": len(exons),
                    "cds_count": len(cds),
                    "cds_length": model_cds_length(model),
                    "review_flags": review_flags,
                    "risk_tags": risk_tags,
                    "blocking_reasons": blocking_reasons,
                    "notes": "candidate coordinates copied from full card; CDS phase recomputed",
                }
            )
    return lines, records


def current_lines_excluding(current_gff: str, remove_gene_ids: Set[str], gene_to_tx: Dict[str, Set[str]]) -> Tuple[List[str], List[str], Counter]:
    remove_tx_ids = {tx_id for gene_id in remove_gene_ids for tx_id in gene_to_tx.get(gene_id, set())}
    header_lines: List[str] = []
    kept_lines: List[str] = []
    counts: Counter = Counter()
    with open_text(current_gff) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if not line.startswith("##gff-version"):
                    header_lines.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                kept_lines.append(line)
                counts["kept_malformed_lines"] += 1
                continue
            feature_l = fields[2].lower()
            row_attrs = parse_attrs(fields[8])
            skip = False
            if feature_l == "gene":
                gene_id = row_attrs.get("ID") or row_attrs.get("Name")
                skip = bool(gene_id and gene_id in remove_gene_ids)
            elif feature_l in {"mrna", "transcript", "rna"}:
                tx_id = row_attrs.get("ID") or row_attrs.get("Name")
                parents = split_list(row_attrs.get("Parent", ""))
                skip = bool((tx_id and tx_id in remove_tx_ids) or any(parent in remove_gene_ids for parent in parents))
            else:
                parents = split_list(row_attrs.get("Parent", ""))
                skip = any(parent in remove_tx_ids for parent in parents)
            if skip:
                counts["removed_current_gff_lines"] += 1
                counts["removed_current_feature:{0}".format(feature_l)] += 1
            else:
                kept_lines.append(line)
                counts["kept_current_gff_lines"] += 1
    return header_lines, kept_lines, counts


def sort_gff_records(lines: Sequence[str]) -> List[str]:
    def key(line: str) -> Tuple[str, int, int, int, str]:
        fields = line.split("\t")
        if len(fields) < 9:
            return ("", 0, 0, 99, line)
        return (fields[0], as_int(fields[3]), as_int(fields[4]), FEATURE_RANK.get(fields[2].lower(), 50), fields[8])

    return sorted(lines, key=key)


def write_report(path: str, args: argparse.Namespace, summary_rows: Sequence[Dict[str, Any]]) -> None:
    with open(path, "w") as handle:
        handle.write("# GeneArbiter Proposal Catalog GFF Export Report\n\n")
        handle.write("Generated: {0}\n\n".format(datetime.now(timezone.utc).isoformat()))
        handle.write("Scope: complete proposal catalog GFF, not final annotation. Truth/manual annotation is not used.\n\n")
        handle.write("Inputs:\n\n")
        handle.write("- current_gff: `{0}`\n".format(args.current_gff))
        handle.write("- full_cards: `{0}`\n".format(args.full_cards))
        handle.write("- final_calls: `{0}`\n\n".format(args.final_calls))
        handle.write("Outputs:\n\n")
        handle.write("- proposal_catalog.gff3\n")
        handle.write("- proposal_catalog_records.tsv\n")
        handle.write("- proposal_catalog_review_flags.tsv\n")
        handle.write("- proposal_catalog_summary.tsv\n\n")
        handle.write("## Summary\n\n")
        for row in summary_rows:
            handle.write("- {0}: {1}\n".format(row["metric"], row["value"]))
        handle.write("\n## Boundary\n\n")
        handle.write("- Current GFF records outside replaced current genes are copied from the input current GFF.\n")
        handle.write("- Candidate proposal coordinates are copied from selected models in full cards.\n")
        handle.write("- Novel selected sets are appended; deletion is not applied by this exporter.\n")
        handle.write("- Review/risk labels are written as attributes and TSV flags, but they do not block proposal output.\n")
        handle.write("- Draft IDs are proposal IDs, not final MH manual-style gene names.\n")


def main() -> None:
    args = parse_args()
    for path in [args.current_gff, args.full_cards, args.final_calls]:
        if not os.path.exists(path):
            raise SystemExit("Missing required input: {0}".format(path))
    if args.fail_on_truth_path:
        truth_path_guard([args.current_gff, args.full_cards, args.final_calls])

    cards = read_jsonl(args.full_cards)
    calls = read_tsv(args.final_calls)
    cards_by_locus, calls_by_locus = build_indexes(cards, calls)
    gene_to_tx, _tx_to_gene, current_gff_gene_ids = collect_current_gene_maps(args.current_gff)

    remove_gene_ids: Set[str] = set()
    proposal_lines: List[str] = []
    records: List[Dict[str, Any]] = []
    review_rows: List[Dict[str, Any]] = []
    summary = Counter()
    missing_call_loci = 0
    missing_selected_models = 0
    missing_current_genes_in_gff: Set[str] = set()

    for locus_id in sorted(cards_by_locus):
        card = cards_by_locus[locus_id]
        call = calls_by_locus.get(locus_id, {})
        if not call:
            missing_call_loci += 1
        model_by_id = {str(model.get("model_id", "")): model for model in card.get("candidate_models", []) if model.get("model_id")}
        selected_ids = model_ids_from_call(call)
        if call.get("final_call") in {"keep_current", "keep_current_flagged"} and not selected_ids:
            selected_ids = [model_id for model_id, model in model_by_id.items() if model_source_info(model)[0] == "current"]
        selected_models = [model_by_id[model_id] for model_id in selected_ids if model_id in model_by_id]
        missing_selected_models += len([model_id for model_id in selected_ids if model_id not in model_by_id])
        current_gene_ids = current_gene_ids_from_card(card)
        action = determine_action(call, selected_models, current_gene_ids)
        summary["action:{0}".format(action)] += 1
        summary["final_call:{0}".format(call.get("final_call", "missing_final_call"))] += 1

        flags = review_flags_from_call(call)
        if flags:
            review_rows.append(
                {
                    "locus_id": locus_id,
                    "proposal_action": action,
                    "final_call": call.get("final_call", ""),
                    "call_group": call.get("call_group", ""),
                    "relation_to_current": call.get("relation_to_current", ""),
                    "export_readiness": call.get("export_readiness", ""),
                    "selected_set_id": selected_set_id_from_call(call),
                    "selected_source": call.get("selected_source", ""),
                    "current_gene_ids": ";".join(current_gene_ids),
                    "review_flags": ";".join(flags),
                    "risk_tags": call.get("risk_tags", ""),
                    "blocking_reasons": call.get("blocking_reasons", ""),
                    "recommended_manual_action": call.get("recommended_manual_action", ""),
                }
            )

        if action == "replace_current_with_selected_set":
            remove_gene_ids.update(current_gene_ids)
            missing_current_genes_in_gff.update(gene_id for gene_id in current_gene_ids if gene_id not in current_gff_gene_ids)
        if action in {"replace_current_with_selected_set", "add_novel_selected_set"}:
            lines, proposal_records = proposal_lines_for_call(card, call, selected_models, action, current_gene_ids, args.source_label)
            proposal_lines.extend(lines)
            records.extend(proposal_records)
        else:
            records.append(
                {
                    "locus_id": locus_id,
                    "proposal_action": action,
                    "final_call": call.get("final_call", ""),
                    "call_group": call.get("call_group", ""),
                    "relation_to_current": call.get("relation_to_current", ""),
                    "export_readiness": call.get("export_readiness", ""),
                    "selected_set_id": selected_set_id_from_call(call),
                    "selected_source": call.get("selected_source", ""),
                    "selected_model_id": call.get("selected_model_id", ""),
                    "current_gene_ids": ";".join(current_gene_ids),
                    "removed_current_gene_ids": "",
                    "review_flags": ";".join(flags),
                    "risk_tags": call.get("risk_tags", ""),
                    "blocking_reasons": call.get("blocking_reasons", ""),
                    "notes": "current backbone retained; no candidate GFF record exported for this locus",
                }
            )

    header_lines, kept_current_lines, current_counts = current_lines_excluding(args.current_gff, remove_gene_ids, gene_to_tx)
    summary_rows: List[Dict[str, Any]] = [
        {"metric": "generated_at_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"metric": "card_loci", "value": len(cards_by_locus)},
        {"metric": "final_call_rows", "value": len(calls_by_locus)},
        {"metric": "missing_call_loci", "value": missing_call_loci},
        {"metric": "proposal_gff_new_records", "value": len(proposal_lines)},
        {"metric": "proposal_transcript_records", "value": sum(1 for row in records if row.get("export_transcript_id"))},
        {"metric": "removed_current_genes", "value": len(remove_gene_ids)},
        {"metric": "missing_selected_models", "value": missing_selected_models},
        {"metric": "missing_current_genes_in_gff", "value": len(missing_current_genes_in_gff)},
        {"metric": "review_flag_loci", "value": len(review_rows)},
    ]
    for key, value in sorted(summary.items()):
        summary_rows.append({"metric": key, "value": value})
    for key, value in sorted(current_counts.items()):
        summary_rows.append({"metric": key, "value": value})

    ensure_dir(args.out_dir)
    gff_path = os.path.join(args.out_dir, "proposal_catalog.gff3")
    records_path = os.path.join(args.out_dir, "proposal_catalog_records.tsv")
    review_path = os.path.join(args.out_dir, "proposal_catalog_review_flags.tsv")
    summary_path = os.path.join(args.out_dir, "proposal_catalog_summary.tsv")
    report_path = os.path.join(args.out_dir, "proposal_catalog_report.md")

    metadata = [
        "##gff-version 3",
        "# proposal_catalog=true",
        "# proposal_status=not_final_annotation",
        "# generated_at_utc={0}".format(datetime.now(timezone.utc).isoformat()),
        "# current_gff={0}".format(args.current_gff),
        "# full_cards={0}".format(args.full_cards),
        "# final_calls={0}".format(args.final_calls),
    ]
    if header_lines:
        metadata.extend("# source_current_header: {0}".format(line.lstrip("# ")) for line in header_lines[:100])
    data_lines = kept_current_lines + proposal_lines
    if args.sort_records:
        data_lines = sort_gff_records(data_lines)
    with open(gff_path, "w") as handle:
        handle.write("\n".join(metadata) + "\n")
        handle.write("\n".join(data_lines) + "\n")

    write_tsv(records_path, RECORD_FIELDS, records)
    write_tsv(review_path, REVIEW_FIELDS, review_rows)
    write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    write_report(report_path, args, summary_rows)

    print("gff\t{0}".format(gff_path))
    print("records\t{0}".format(records_path))
    print("review_flags\t{0}".format(review_path))
    print("summary\t{0}".format(summary_path))
    print("removed_current_genes\t{0}".format(len(remove_gene_ids)))
    print("proposal_transcripts\t{0}".format(sum(1 for row in records if row.get("export_transcript_id"))))


if __name__ == "__main__":
    main()
