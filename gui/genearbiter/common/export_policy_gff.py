#!/usr/bin/env python3
# script_id_md5: 7c1425f17c1467da296301a9f97bfaa8
# created: 2026-07-12
# modified: 2026-07-13
# owner: project
# status: project_code
# purpose: 从 GeneArbiter proposal catalog 导出保守 policy GFF，并生成匹配的 proposal records 供 ID reconciliation 使用。
# inputs: proposal_catalog.gff3；proposal_catalog_records.tsv；final_annotation_calls.tsv；model_arbitration_models.tsv；model_arbitration_full_cards.jsonl；current/base GFF。
# outputs: policy_catalog.gff3；policy_catalog_records.tsv；policy_catalog_summary.tsv；policy_catalog_report.md。
# notes: 不读取 truth annotation；不重跑 AI；只过滤或恢复已有 current/proposal 记录，不生成新坐标。
#        strict_* policies apply novel-specific filters to all no-current/add_novel proposals before export.

"""Export a conservative GeneArbiter policy GFF from a proposal catalog.

The raw proposal catalog intentionally keeps a broad review-prioritized union of
selected replacement and novel candidates. This post-export policy layer keeps
the current backbone, restores current genes when excluded replacement proposals
are removed, and can apply conservative rules to current-uncovered novel genes.

The default policy, supported_novel_singletool_junction_supported, keeps
GeneArbiter-supported novel candidates only when the final-call support layer marks
them strong/moderate and not unsupported; single-tool novel candidates additionally
need full or partial short-read junction support. Coordinates are copied from
current/proposal GFF records. Truth/manual evaluation files are refused.

Strict multi-tool support is evaluated against the exact CDS structure selected
for each exported gene. A supporting source must reproduce the same chromosome,
strand, ordered CDS coordinates, and CDS phases. Multiple transcripts from one
source count as one prediction tool.

The strict_* policies differ from legacy policy names by evaluating all
current-uncovered/add_novel proposals with the novel filter before considering
pending-validation export readiness. This prevents single-tool no-current novel
candidates from being exported merely because AI/auto selected them.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple
from urllib.parse import unquote

from gff_utils import ensure_dir, write_tsv


POLICY_CHOICES = {
    "supported_novel_only",
    "supported_novel_plus_risk_replacement",
    "supported_novel_multitool_only",
    "supported_novel_singletool_junction_supported",
    "strict_supported_novel_multitool_only",
    "strict_supported_novel_singletool_junction_supported",
}

KEEP_JUNCTION_STATUS = {"full_junction_support", "partial_junction_support"}
SUPPORTED_EVIDENCE_LEVELS = {"strong", "moderate"}
TRUTH_LIKE_RE = re.compile(r"truth|manual", re.IGNORECASE)

FEATURE_RANK = {
    "gene": 0,
    "mRNA": 1,
    "transcript": 1,
    "rna": 1,
    "exon": 2,
    "CDS": 2,
    "five_prime_UTR": 2,
    "three_prime_UTR": 2,
    "five_prime_utr": 2,
    "three_prime_utr": 2,
}

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
    "selected_cds_support_sources",
    "selected_cds_support_count",
    "policy_decision_reason",
    "notes",
]

SUMMARY_FIELDS = ["metric", "value"]


@dataclass
class Record:
    raw: str
    parts: List[str]
    attrs: Dict[str, str]
    idx: int


@dataclass
class GeneBlock:
    gene_id: str
    gene_record: Record
    records: List[Record] = field(default_factory=list)

    @property
    def sort_key(self) -> Tuple[object, int, int, str]:
        fields = self.gene_record.parts
        return (natural_seq_key(fields[0]), int_or_big(fields[3]), int_or_big(fields[4]), self.gene_id)

    def ordered_records(self) -> List[Record]:
        return sorted(self.records, key=lambda rec: (FEATURE_RANK.get(rec.parts[2], 5), rec.idx))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposal-gff", required=True, help="Raw GeneArbiter proposal_catalog.gff3.")
    parser.add_argument("--proposal-records", required=True, help="Raw proposal_catalog_records.tsv.")
    parser.add_argument("--final-calls", required=True, help="final_annotation_calls.tsv.")
    parser.add_argument("--models", required=True, help="model_arbitration_models.tsv with RNA junction support status.")
    parser.add_argument("--full-cards", default="", help="model_arbitration_full_cards.jsonl with exact candidate CDS structures.")
    parser.add_argument("--exact-support-cache", default="", help="Reusable TSV cache of exact CDS support; regenerated when --full-cards is supplied.")
    parser.add_argument("--out-dir", required=True, help="Output directory for policy_catalog.* files.")
    parser.add_argument("--policy", default="supported_novel_singletool_junction_supported", choices=sorted(POLICY_CHOICES))
    parser.add_argument("--current-gff", default="", help="Current/base GFF. Defaults to # current_gff in proposal GFF header.")
    parser.add_argument("--fail-on-truth-path", action="store_true", help="Fail if input/output paths look like truth/manual files.")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open(encoding="utf-8")


def natural_seq_key(seqid: str) -> Tuple[object, ...]:
    return tuple(int(token) if token.isdigit() else token for token in re.split(r"(\d+)", seqid))


def int_or_big(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 10**18


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def split_list(value: object) -> List[str]:
    decoded = unquote(str(value or ""))
    return [item.strip() for item in decoded.replace(";", ",").split(",") if item.strip()]


def attr_values(attrs: Dict[str, str], key: str) -> List[str]:
    value = attrs.get(key, "")
    return [item for item in re.split(r"[, ]+", value) if item] if value else []


def fail_if_truth_like(paths: Sequence[str]) -> None:
    for path in paths:
        if TRUTH_LIKE_RE.search(str(path)):
            raise SystemExit("Refusing truth/manual-like path for policy export: {0}".format(path))


def current_gff_from_header(proposal_gff: Path) -> Path:
    with open_text(proposal_gff) as handle:
        for line in handle:
            if line.startswith("# current_gff="):
                return Path(line.strip().split("=", 1)[1])
            if not line.startswith("#"):
                break
    raise SystemExit("Missing # current_gff header in proposal GFF: {0}".format(proposal_gff))


def parse_gene_blocks(path: Path) -> Tuple[List[str], Dict[str, GeneBlock], Counter]:
    headers: List[str] = []
    records: List[Record] = []
    counts: Counter = Counter()
    with open_text(path) as handle:
        for idx, line in enumerate(handle):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                headers.append(line)
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                counts["malformed_rows"] += 1
                continue
            records.append(Record(line, parts, parse_attrs(parts[8]), idx))

    gene_records: List[Tuple[str, str, Record]] = []
    gene_keys_by_id: Dict[str, List[str]] = {}
    gene_record_key_by_idx: Dict[int, str] = {}
    tx_to_gene_id: Dict[str, str] = {}
    gene_seen: Counter = Counter()
    for rec in records:
        row_id = rec.attrs.get("ID", "")
        feature = rec.parts[2].lower()
        if feature == "gene" and row_id:
            gene_seen[row_id] += 1
            if gene_seen[row_id] > 1:
                counts["duplicate_gene_ids"] += 1
            block_key = row_id if gene_seen[row_id] == 1 else "{0}__duplicate_{1}".format(row_id, gene_seen[row_id])
            gene_records.append((block_key, row_id, rec))
            gene_keys_by_id.setdefault(row_id, []).append(block_key)
            gene_record_key_by_idx[rec.idx] = block_key
        elif feature in {"mrna", "transcript", "rna"} and row_id:
            parent = ""
            parents = attr_values(rec.attrs, "Parent")
            if parents:
                parent = parents[0]
            parent = rec.attrs.get("geneID") or rec.attrs.get("gene_id") or parent
            if parent:
                tx_to_gene_id[row_id] = parent

    assigned: Dict[str, List[Record]] = {block_key: [] for block_key, _gene_id, _record in gene_records}
    for rec in records:
        row_id = rec.attrs.get("ID", "")
        target_gene_id = ""
        if rec.parts[2].lower() == "gene" and rec.idx in gene_record_key_by_idx:
            assigned[gene_record_key_by_idx[rec.idx]].append(rec)
            continue
        for parent in attr_values(rec.attrs, "Parent"):
            if parent in gene_keys_by_id:
                target_gene_id = parent
                break
            if parent in tx_to_gene_id:
                target_gene_id = tx_to_gene_id[parent]
                break
        if not target_gene_id:
            for key in ("geneID", "gene_id"):
                value = rec.attrs.get(key, "")
                if value in gene_keys_by_id:
                    target_gene_id = value
                    break
        if target_gene_id:
            keys = gene_keys_by_id.get(target_gene_id, [])
            if len(keys) == 1:
                assigned[keys[0]].append(rec)
            elif keys:
                assigned[keys[0]].append(rec)
                counts["ambiguous_duplicate_gene_child_rows"] += 1

    blocks: Dict[str, GeneBlock] = {}
    for block_key, gene_id, gene_record in gene_records:
        rows = assigned.get(block_key) or [gene_record]
        if gene_record not in rows:
            rows.insert(0, gene_record)
        blocks[block_key] = GeneBlock(gene_id, gene_record, rows)
    counts["gene_blocks"] = len(blocks)
    return headers, blocks, counts


def read_tsv_by_key(path: Path, key: str) -> Dict[str, Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row.get(key, ""): row for row in csv.DictReader(handle, delimiter="\t") if row.get(key, "")}


def read_records_by_gene(path: Path) -> Dict[str, List[Dict[str, str]]]:
    by_gene: Dict[str, List[Dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row.get("export_gene_id", "")
            if gene_id:
                by_gene.setdefault(gene_id, []).append(row)
    return by_gene


def cds_signature(seqid: object, strand: object, cds_rows: Iterable[Sequence[object]]) -> str:
    """Return a stable exact CDS signature including phase."""
    normalized: List[Tuple[int, int, str]] = []
    for item in cds_rows:
        if len(item) < 2:
            continue
        try:
            start = int(item[0])
            end = int(item[1])
        except (TypeError, ValueError):
            continue
        phase = str(item[2]) if len(item) >= 3 else "."
        normalized.append((min(start, end), max(start, end), phase))
    if not normalized:
        return ""
    payload = [str(seqid or ""), str(strand or "."), sorted(set(normalized))]
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def block_cds_signatures(block: GeneBlock) -> List[str]:
    """Collect exact CDS signatures for every coding transcript in a proposal gene."""
    transcript_ids: Set[str] = set()
    for rec in block.records:
        if rec.parts[2].lower() in {"mrna", "transcript", "rna"} and rec.attrs.get("ID"):
            transcript_ids.add(rec.attrs["ID"])

    cds_by_transcript: Dict[str, List[Tuple[int, int, str]]] = {}
    for rec in block.records:
        if rec.parts[2].lower() != "cds":
            continue
        parents = attr_values(rec.attrs, "Parent")
        targets = [parent for parent in parents if parent in transcript_ids]
        if not targets and (not transcript_ids or block.gene_id in parents):
            targets = ["__gene_cds__"]
        for transcript_id in targets:
            cds_by_transcript.setdefault(transcript_id, []).append(
                (int(rec.parts[3]), int(rec.parts[4]), rec.parts[7])
            )

    seqid = block.gene_record.parts[0]
    strand = block.gene_record.parts[6]
    return sorted(
        {
            signature
            for signature in (
                cds_signature(seqid, strand, rows) for rows in cds_by_transcript.values()
            )
            if signature
        }
    )


def load_exact_cds_support_from_full_cards(
    path: Path,
    target_loci: Set[str],
) -> Dict[str, Dict[str, Set[str]]]:
    """Stream full cards and index exact CDS signatures only for relevant loci."""
    support: Dict[str, Dict[str, Set[str]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                card = json.loads(line)
            except json.JSONDecodeError as exc:
                raise SystemExit("Invalid full-card JSON at {0}:{1}: {2}".format(path, line_number, exc))
            locus_id = str(card.get("locus_id", ""))
            if locus_id not in target_loci:
                continue
            locus_support = support.setdefault(locus_id, {})
            for model in card.get("candidate_models") or []:
                source = str(model.get("source", ""))
                if not source or source.lower() == "current":
                    continue
                signature = cds_signature(model.get("seqid", ""), model.get("strand", "."), model.get("cds") or [])
                if signature:
                    locus_support.setdefault(signature, set()).add(source)
    return support


def write_exact_support_cache(path: Path, support: Dict[str, Dict[str, Set[str]]]) -> None:
    rows: List[Dict[str, object]] = []
    for locus_id in sorted(support):
        for signature, sources in sorted(support[locus_id].items()):
            rows.append(
                {
                    "locus_id": locus_id,
                    "cds_signature": signature,
                    "supporting_sources": ",".join(sorted(sources)),
                    "support_count": len(sources),
                }
            )
    ensure_dir(str(path.parent))
    write_tsv(str(path), ["locus_id", "cds_signature", "supporting_sources", "support_count"], rows)


def read_exact_support_cache(path: Path) -> Dict[str, Dict[str, Set[str]]]:
    support: Dict[str, Dict[str, Set[str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus_id = row.get("locus_id", "")
            signature = row.get("cds_signature", "")
            if locus_id and signature:
                support.setdefault(locus_id, {})[signature] = {
                    source for source in split_list(row.get("supporting_sources", ""))
                    if source.lower() != "current"
                }
    return support


def exact_selected_cds_sources(
    block: GeneBlock,
    locus_support: Dict[str, Set[str]],
) -> Set[str]:
    """Require the same independent sources to support every selected coding transcript."""
    signatures = block_cds_signatures(block)
    if not signatures:
        return set()
    source_sets = [set(locus_support.get(signature, set())) for signature in signatures]
    return set.intersection(*source_sets) if source_sets else set()


def call_supported_novel(call: Dict[str, str]) -> bool:
    return (
        call.get("final_call") == "novel_gene_candidate_manual_review"
        and call.get("best_evidence_level") in SUPPORTED_EVIDENCE_LEVELS
        and call.get("manual_review_unsupported_locus") == "false"
    )


def call_has_supported_novel_evidence(call: Dict[str, str]) -> bool:
    return (
        call.get("best_evidence_level") in SUPPORTED_EVIDENCE_LEVELS
        and call.get("manual_review_unsupported_locus") == "false"
    )


def member_sources(call: Dict[str, str]) -> Set[str]:
    return {item for item in split_list(call.get("member_real_sources", "")) if item and item != "current"}


def is_current_uncovered_novel(attrs: Dict[str, str]) -> bool:
    return (
        attrs.get("proposal_action") == "add_novel_selected_set"
        or attrs.get("relation_to_current") == "novel_gene_candidate_no_current_overlap"
    )


def include_strict_novel_block(
    block: GeneBlock,
    call: Dict[str, str],
    records_by_gene: Dict[str, List[Dict[str, str]]],
    models: Dict[str, Dict[str, str]],
    exact_support: Dict[str, Dict[str, Set[str]]],
    policy: str,
) -> Tuple[bool, str]:
    if not call_has_supported_novel_evidence(call):
        return False, "strict_drop_novel_not_strong_moderate_or_unsupported"
    locus_id = block.gene_record.attrs.get("locus_id", "")
    if len(exact_selected_cds_sources(block, exact_support.get(locus_id, {}))) >= 2:
        return True, "strict_keep_supported_exact_cds_multitool_novel"
    if policy == "strict_supported_novel_singletool_junction_supported":
        if has_full_or_partial_junction(block, records_by_gene, models):
            return True, "strict_keep_supported_singletool_junction_novel"
        return False, "strict_drop_supported_singletool_no_full_or_partial_junction"
    return False, "strict_drop_supported_singletool_or_unknown_novel"


def source_model_ids(block: GeneBlock, records_by_gene: Dict[str, List[Dict[str, str]]]) -> List[str]:
    ids: List[str] = []
    ids.extend(split_list(block.gene_record.attrs.get("source_model_id", "")))
    for row in records_by_gene.get(block.gene_id, []):
        ids.extend(split_list(row.get("source_model_id", "")))
    seen: Set[str] = set()
    out: List[str] = []
    for model_id in ids:
        if model_id not in seen:
            seen.add(model_id)
            out.append(model_id)
    return out


def has_full_or_partial_junction(block: GeneBlock, records_by_gene: Dict[str, List[Dict[str, str]]], models: Dict[str, Dict[str, str]]) -> bool:
    for model_id in source_model_ids(block, records_by_gene):
        if models.get(model_id, {}).get("rna_junction_status", "") in KEEP_JUNCTION_STATUS:
            return True
    return False


def include_policy_block(
    block: GeneBlock,
    call: Dict[str, str],
    records_by_gene: Dict[str, List[Dict[str, str]]],
    models: Dict[str, Dict[str, str]],
    exact_support: Dict[str, Dict[str, Set[str]]],
    policy: str,
) -> Tuple[bool, str]:
    attrs = block.gene_record.attrs
    export_readiness = attrs.get("export_readiness", "")
    if policy.startswith("strict_") and is_current_uncovered_novel(attrs):
        return include_strict_novel_block(block, call, records_by_gene, models, exact_support, policy)

    if not export_readiness:
        return True, "keep_current_or_base_backbone"
    if export_readiness == "pending_validation_before_export":
        return True, "keep_pending_validation_replacement"

    final_call = attrs.get("final_call", "")
    relation = attrs.get("relation_to_current", "")
    if final_call == "novel_gene_candidate_manual_review":
        if not call_supported_novel(call):
            return False, "drop_novel_not_strong_moderate_or_unsupported"
        if policy in {"supported_novel_only", "supported_novel_plus_risk_replacement"}:
            return True, "keep_supported_novel"
        if policy == "supported_novel_multitool_only":
            if len(member_sources(call)) >= 2:
                return True, "keep_supported_multitool_novel"
            return False, "drop_supported_singletool_or_unknown_novel"
        if policy == "supported_novel_singletool_junction_supported":
            if len(member_sources(call)) >= 2:
                return True, "keep_supported_multitool_novel"
            if has_full_or_partial_junction(block, records_by_gene, models):
                return True, "keep_supported_singletool_junction_novel"
            return False, "drop_supported_singletool_no_full_or_partial_junction"

    if policy == "supported_novel_plus_risk_replacement" and relation == "candidate_set_selected_with_risk_flag":
        return True, "keep_risk_flagged_replacement"
    return False, "drop_manual_review_or_blocked_proposal"


def count_gene_tx(blocks: Iterable[GeneBlock]) -> Counter:
    counts: Counter = Counter()
    seen_gene_ids: Set[str] = set()
    for block in blocks:
        counts["genes"] += 1
        if block.gene_id in seen_gene_ids:
            counts["duplicate_gene_ids"] += 1
        seen_gene_ids.add(block.gene_id)
        for rec in block.records:
            feature = rec.parts[2].lower()
            if feature in {"mrna", "transcript", "rna"}:
                counts["transcripts"] += 1
    return counts


def write_report(path: Path, args: argparse.Namespace, current_gff: Path, summary_rows: Sequence[Dict[str, object]]) -> None:
    lines = [
        "# GeneArbiter Policy GFF Export Report",
        "",
        "Generated: {0}".format(datetime.now(timezone.utc).isoformat()),
        "",
        "- policy: `{0}`".format(args.policy),
        "- proposal_gff: `{0}`".format(args.proposal_gff),
        "- proposal_records: `{0}`".format(args.proposal_records),
        "- final_calls: `{0}`".format(args.final_calls),
        "- models: `{0}`".format(args.models),
        "- full_cards: `{0}`".format(args.full_cards or "not_supplied"),
        "- exact_support_cache: `{0}`".format(args.exact_support_cache or "not_supplied"),
        "- current_gff: `{0}`".format(current_gff),
        "",
        "## Summary",
        "",
    ]
    for row in summary_rows:
        lines.append("- {0}: {1}".format(row["metric"], row["value"]))
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- This is a post-export policy filter over frozen GeneArbiter proposal/final-call artifacts.",
            "- Truth/manual evaluation annotation is not used.",
            "- Coordinates are copied from current/proposal GFF records; excluded replacement proposals restore current gene blocks.",
            "- Strict multi-tool novel support means at least two independent non-current sources match every selected transcript's exact CDS coordinates, strand, and phase.",
            "- ID reconciliation is performed by the GeneArbiter reconcile_ids step when policy_reconcile_ids is run.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    proposal_gff = Path(args.proposal_gff)
    proposal_records = Path(args.proposal_records)
    final_calls = Path(args.final_calls)
    models_path = Path(args.models)
    out_dir = Path(args.out_dir)
    current_gff = Path(args.current_gff) if args.current_gff else current_gff_from_header(proposal_gff)
    full_cards = Path(args.full_cards) if args.full_cards else None
    exact_support_cache = Path(args.exact_support_cache) if args.exact_support_cache else None

    inputs = [proposal_gff, proposal_records, final_calls, models_path, current_gff]
    if full_cards:
        inputs.append(full_cards)
    elif exact_support_cache:
        inputs.append(exact_support_cache)
    for path in inputs:
        if not path.exists():
            raise SystemExit("Missing required input: {0}".format(path))
    if args.fail_on_truth_path:
        fail_if_truth_like([str(path) for path in inputs] + [str(out_dir)])

    _proposal_headers, proposal_blocks, proposal_parse_counts = parse_gene_blocks(proposal_gff)
    _current_headers, current_blocks, current_parse_counts = parse_gene_blocks(current_gff)
    calls = read_tsv_by_key(final_calls, "locus_id")
    models = read_tsv_by_key(models_path, "model_id")
    records_by_gene = read_records_by_gene(proposal_records)
    novel_loci = {
        block.gene_record.attrs.get("locus_id", "")
        for block in proposal_blocks.values()
        if is_current_uncovered_novel(block.gene_record.attrs)
    }
    novel_loci.discard("")
    if full_cards:
        exact_support = load_exact_cds_support_from_full_cards(full_cards, novel_loci)
        if exact_support_cache:
            write_exact_support_cache(exact_support_cache, exact_support)
        exact_support_source = "full_cards"
    elif exact_support_cache:
        exact_support = read_exact_support_cache(exact_support_cache)
        exact_support_source = "cache"
    elif args.policy.startswith("strict_"):
        raise SystemExit(
            "Strict policy requires --full-cards or --exact-support-cache so multi-tool support can be checked against exact CDS structures"
        )
    else:
        exact_support = {}
        exact_support_source = "not_supplied"

    selected: Dict[str, GeneBlock] = {}
    restored_current_gene_ids: Set[str] = set()
    included_proposal_gene_ids: Set[str] = set()
    decision_counts: Counter = Counter()
    final_call_counts: Counter = Counter()
    selected_policy_counts: Counter = Counter()
    policy_decision_by_gene: Dict[str, str] = {}
    exact_sources_by_gene: Dict[str, Set[str]] = {}

    for gene_id, block in proposal_blocks.items():
        attrs = block.gene_record.attrs
        locus_id = attrs.get("locus_id", "")
        call = calls.get(locus_id, {})
        exact_sources_by_gene[gene_id] = exact_selected_cds_sources(block, exact_support.get(locus_id, {}))
        include, reason = include_policy_block(block, call, records_by_gene, models, exact_support, args.policy)
        policy_decision_by_gene[gene_id] = reason
        decision_counts[reason] += 1
        final_call_counts[attrs.get("final_call", "current_or_base")] += 1
        if include:
            selected["proposal::{0}".format(gene_id)] = block
            if attrs.get("proposal_catalog", "").lower() == "true":
                included_proposal_gene_ids.add(gene_id)
                selected_policy_counts[reason] += 1
            continue

        for current_gene_id in split_list(attrs.get("replaces_current_genes", "")):
            if current_gene_id:
                restored_current_gene_ids.add(current_gene_id)

    current_blocks_by_gene_id: Dict[str, List[GeneBlock]] = {}
    for block in current_blocks.values():
        current_blocks_by_gene_id.setdefault(block.gene_id, []).append(block)
    selected_current_raw = {
        block.gene_record.raw
        for block in selected.values()
        if block.gene_record.attrs.get("proposal_catalog", "").lower() != "true"
    }
    restored_count = 0
    for current_gene_id in sorted(restored_current_gene_ids):
        for block in current_blocks_by_gene_id.get(current_gene_id, []):
            if block.gene_record.raw in selected_current_raw:
                continue
            selected["restore::{0}::{1}".format(current_gene_id, block.gene_record.idx)] = block
            selected_current_raw.add(block.gene_record.raw)
            restored_count += 1

    selected_blocks = sorted(selected.values(), key=lambda block: block.sort_key)
    policy_records_rows: List[Dict[str, str]] = []
    for gene_id in sorted(included_proposal_gene_ids):
        sources = sorted(exact_sources_by_gene.get(gene_id, set()))
        for source_row in records_by_gene.get(gene_id, []):
            row = dict(source_row)
            row["selected_cds_support_sources"] = ",".join(sources)
            row["selected_cds_support_count"] = str(len(sources))
            row["policy_decision_reason"] = policy_decision_by_gene.get(gene_id, "")
            policy_records_rows.append(row)

    ensure_dir(str(out_dir))
    gff_path = out_dir / "policy_catalog.gff3"
    records_path = out_dir / "policy_catalog_records.tsv"
    summary_path = out_dir / "policy_catalog_summary.tsv"
    report_path = out_dir / "policy_catalog_report.md"

    with gff_path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write("# generated_by=GeneArbiter export_policy_gff\n")
        handle.write("# generated_at_utc={0}\n".format(datetime.now(timezone.utc).isoformat()))
        handle.write("# policy={0}\n".format(args.policy))
        handle.write("# policy_status=candidate_correction_not_final_annotation\n")
        handle.write("# proposal_gff={0}\n".format(proposal_gff))
        handle.write("# final_calls={0}\n".format(final_calls))
        handle.write("# current_gff={0}\n".format(current_gff))
        handle.write("# exact_cds_support_source={0}\n".format(exact_support_source))
        handle.write("# truth_usage=none\n")
        for block in selected_blocks:
            for rec in block.ordered_records():
                handle.write(rec.raw + "\n")

    output_counts = count_gene_tx(selected_blocks)
    summary_rows: List[Dict[str, object]] = [
        {"metric": "generated_at_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"metric": "policy", "value": args.policy},
        {"metric": "exact_cds_support_source", "value": exact_support_source},
        {"metric": "exact_cds_support_loci", "value": len(exact_support)},
        {"metric": "proposal_gene_blocks", "value": proposal_parse_counts.get("gene_blocks", 0)},
        {"metric": "current_gene_blocks", "value": current_parse_counts.get("gene_blocks", 0)},
        {"metric": "output_gene_blocks", "value": output_counts.get("genes", 0)},
        {"metric": "output_transcripts", "value": output_counts.get("transcripts", 0)},
        {"metric": "included_proposal_gene_blocks", "value": len(included_proposal_gene_ids)},
        {"metric": "policy_records", "value": len(policy_records_rows)},
        {"metric": "restored_current_gene_blocks", "value": restored_count},
        {"metric": "duplicate_gene_ids_in_output_blocks", "value": output_counts.get("duplicate_gene_ids", 0)},
    ]
    for key, value in sorted(decision_counts.items()):
        summary_rows.append({"metric": "decision:" + key, "value": value})
    for key, value in sorted(selected_policy_counts.items()):
        summary_rows.append({"metric": "included:" + key, "value": value})
    for key, value in sorted(final_call_counts.items()):
        summary_rows.append({"metric": "final_call:" + (key or "missing"), "value": value})

    write_tsv(str(records_path), RECORD_FIELDS, policy_records_rows)
    write_tsv(str(summary_path), SUMMARY_FIELDS, summary_rows)
    write_report(report_path, args, current_gff, summary_rows)

    print("policy_gff\t{0}".format(gff_path))
    print("policy_records\t{0}".format(records_path))
    print("summary\t{0}".format(summary_path))
    print("policy\t{0}".format(args.policy))
    print("output_gene_blocks\t{0}".format(output_counts.get("genes", 0)))
    print("included_proposal_gene_blocks\t{0}".format(len(included_proposal_gene_ids)))
    print("restored_current_gene_blocks\t{0}".format(restored_count))


if __name__ == "__main__":
    main()
