#!/usr/bin/env python3
# script_id_md5: f9606eede378895d8a589985025a58d5
# created: 2026-07-07
# modified: 2026-07-12
# owner: project
# status: project_code
# purpose: 为 GeneArbiter proposal catalog 生成 final naming GFF 和 gene/transcript ID 映射表。
# inputs: proposal_catalog.gff3；proposal_catalog_records.tsv。
# outputs: final_named_annotation.gff3；final_named_annotation.auto_pass.clean.gff3；final_named_annotation.auto_pass.hisat.gtf；gene_id_mapping.tsv；transcript_id_mapping.tsv；final_naming_summary.tsv；final_naming_report.md。
# notes: 不读取 truth annotation；不改变坐标；proposal catalog 仍需人工/规则验证后才能作为 final annotation release。

"""Reconcile temporary proposal IDs into a release-like naming layer.

The upstream proposal catalog deliberately uses temporary IDs so that
source/tool IDs are not mistaken for final annotation IDs. This script creates
a second naming layer:

* kept current genes remain unchanged;
* one-to-one replacement can reuse the retired current gene ID by default;
* split/merge proposals preserve one structurally closest current ID when possible;
* additional/novel proposals can receive position-based IDs in the current namespace;
* transcript/exon/CDS IDs are regenerated under the final gene ID;
* explicit mapping tables preserve current, proposal and source provenance.

This is a naming and mapping step only. It does not decide whether a proposal is
biologically accepted, and it does not read final-evaluation truth.
"""

from __future__ import annotations

import argparse
import csv
import re
from bisect import bisect_right
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple
from urllib.parse import quote, unquote

from gff_utils import ensure_dir, iter_gff_rows, open_text, parse_attrs, write_tsv


GENE_MAP_FIELDS = [
    "locus_id",
    "proposal_action",
    "final_relation",
    "reuse_current_gene_id",
    "final_gene_id",
    "proposal_gene_id",
    "current_gene_ids",
    "retired_current_gene_ids",
    "source",
    "source_gene_id",
    "source_model_id",
    "selected_set_id",
    "selected_source",
    "seqid",
    "start",
    "end",
    "strand",
    "id_assignment_basis",
    "id_assignment_status",
    "anchor_current_gene_id",
    "notes",
]

TX_MAP_FIELDS = [
    "locus_id",
    "final_relation",
    "final_gene_id",
    "final_transcript_id",
    "proposal_gene_id",
    "proposal_transcript_id",
    "current_gene_ids",
    "source",
    "source_gene_id",
    "source_transcript_id",
    "source_model_id",
    "selected_set_id",
    "selected_source",
    "seqid",
    "start",
    "end",
    "strand",
]

CURRENT_DISAMBIG_FIELDS = [
    "original_id",
    "final_id",
    "feature",
    "seqid",
    "start",
    "end",
    "strand",
    "duplicate_index",
    "original_parent",
    "final_parent",
    "reason",
]

SUMMARY_FIELDS = ["metric", "value"]

LOCUS_MAP_FIELDS = [
    "locus_id",
    "proposal_action",
    "final_relation",
    "seqid",
    "locus_start",
    "locus_end",
    "current_gene_ids",
    "corrected_gene_ids",
    "retained_current_gene_ids",
    "retired_current_gene_ids",
    "added_gene_ids",
    "mapping_status",
    "mapping_basis",
]

TRUTH_LIKE_RE = re.compile(r"truth|manual", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposal-gff", required=True, help="GeneArbiter proposal_catalog.gff3 from export_proposal_gff.py.")
    parser.add_argument("--proposal-records", required=True, help="proposal_catalog_records.tsv with current/source/proposal ID mapping.")
    parser.add_argument("--out-dir", required=True, help="Output directory for final naming artifacts.")
    parser.add_argument("--output-gff-name", default="final_named_annotation.gff3", help="Output GFF3 filename.")
    parser.add_argument("--clean-gff-name", default="final_named_annotation.auto_pass.clean.gff3", help="Clean auto-pass GFF3 filename; empty disables clean GFF output.")
    parser.add_argument("--clean-source-label", default="GeneArbiter", help="Source label used for proposal records in the clean auto-pass GFF.")
    parser.add_argument("--proposal-source-label", default="", help="Optional source label applied to proposal records in the named GFF (for example GeneArbiter).")
    parser.add_argument("--gtf-name", default="final_named_annotation.auto_pass.hisat.gtf", help="RNA-seq/HISAT-friendly GTF filename; empty disables GTF output.")
    parser.add_argument("--gtf-exclude-organelle", action="store_true", default=True, help="Exclude seqids whose region has genome=chloroplast/mitochondrion from the GTF.")
    parser.add_argument("--new-gene-prefix", default="auto", help="Fallback prefix for IDs whose source style cannot be inferred; 'auto' uses GeneArbiterG.")
    parser.add_argument("--new-gene-start", type=int, default=1, help="Starting integer for newly allocated gene IDs.")
    parser.add_argument("--new-gene-width", type=int, default=6, help="Zero-padding width for newly allocated gene IDs.")
    parser.add_argument(
        "--id-strategy",
        choices=["sequential", "current_style_positional", "current_style_or_sequential"],
        default="sequential",
        help="Allocate new IDs sequentially, positionally, or by source-style auto inference with a safe sequential fallback.",
    )
    parser.add_argument(
        "--current-gff",
        default="",
        help="Original current GFF used only for ID namespace/structural lineage. Required by current_style_positional.",
    )
    parser.add_argument(
        "--transcript-template",
        default="{gene_id}.t{index:02d}",
        help="Python format template for final transcript IDs; available keys: gene_id, index.",
    )
    parser.add_argument(
        "--transcript-id-strategy",
        choices=["template", "current_style"],
        default="template",
        help="Use --transcript-template or convert the final gene ID's terminal G to T and append .<isoform>.",
    )
    parser.add_argument(
        "--one-to-one-policy",
        choices=["reuse_current_gene_id", "new_id"],
        default="reuse_current_gene_id",
        help="How to name a single proposal gene replacing exactly one current gene.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite output files if they already exist.")
    parser.add_argument("--fail-on-truth-path", action="store_true", help="Fail if input paths look like truth/manual evaluation files.")
    return parser.parse_args()


def split_list(text: object) -> List[str]:
    return [item.strip() for item in str(text or "").replace(";", ",").split(",") if item.strip()]


def public_proposal_id(value: object) -> str:
    """Return a stable public proposal ID without an internal workflow-version label."""
    text = str(value or "")
    if text.startswith("V2PROP_"):
        return "GeneArbiterProposal_" + text[len("V2PROP_") :]
    return text


def attr_value(text: object) -> str:
    return quote(str(text or ""), safe="._:-|,")


def attrs(items: Iterable[Tuple[str, object]]) -> str:
    return ";".join("{0}={1}".format(key, attr_value(value)) for key, value in items if value not in {"", None})


def gff_line(seqid: object, source: object, feature: object, start: object, end: object, score: object, strand: object, phase: object, attr_text: str) -> str:
    return "\t".join([str(seqid), str(source), str(feature), str(start), str(end), str(score or "."), str(strand or "."), str(phase or "."), attr_text])


def fail_if_truth_like(paths: Sequence[str]) -> None:
    for path in paths:
        if TRUTH_LIKE_RE.search(path):
            raise SystemExit("Refusing truth/manual-like input path for ID reconciliation: {0}".format(path))


def read_records(path: str) -> List[Dict[str, str]]:
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_new_gene_id(prefix: str, index: int, width: int) -> str:
    return "{0}{1:0{2}d}".format(prefix, index, width)


def proposal_gene_order(records: Sequence[Dict[str, str]]) -> List[str]:
    seen: Set[str] = set()
    ordered = sorted(
        records,
        key=lambda row: (
            str(row.get("seqid", "")),
            int(row.get("start") or 0),
            int(row.get("end") or 0),
            str(row.get("export_gene_id", "")),
        ),
    )
    result: List[str] = []
    for row in ordered:
        gene_id = row.get("export_gene_id", "")
        if gene_id and gene_id not in seen:
            seen.add(gene_id)
            result.append(gene_id)
    return result


def relation_for_locus(action: str, current_count: int, proposal_gene_count: int) -> str:
    if action == "add_novel_selected_set":
        return "novel_gene"
    if action != "replace_current_with_selected_set":
        return "other_proposal"
    if current_count == 1 and proposal_gene_count == 1:
        return "one_to_one_replacement"
    if current_count == 1 and proposal_gene_count > 1:
        return "split_replacement"
    if current_count > 1 and proposal_gene_count == 1:
        return "merge_replacement"
    return "complex_replacement"


def collect_existing_ids(proposal_gff: str) -> Tuple[Set[str], Set[str]]:
    gene_ids: Set[str] = set()
    tx_ids: Set[str] = set()
    for _seqid, _source, feature, _start, _end, _score, _strand, _phase, row_attrs in iter_gff_rows(proposal_gff):
        feature_l = feature.lower()
        row_id = row_attrs.get("ID")
        is_proposal = row_attrs.get("proposal_catalog", "").lower() == "true"
        if not row_id or is_proposal:
            continue
        if feature_l == "gene":
            gene_ids.add(row_id)
        elif feature_l in {"mrna", "transcript", "rna"}:
            tx_ids.add(row_id)
    return gene_ids, tx_ids


NUMERIC_TOKEN_RE = re.compile(r"\d+")


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    merged: List[List[int]] = []
    for start, end in sorted(set(intervals)):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def collect_current_genes(path: str) -> Dict[str, Dict[str, object]]:
    """Read current gene/CDS models for namespace reservation and lineage matching.

    Some published annotations (for example the tested cotton annotation) use
    parentless mRNA rows as top-level gene models. Those rows are accepted as
    gene anchors instead of inventing a different biological model here.
    """
    genes: Dict[str, Dict[str, object]] = {}
    transcripts: Dict[str, Dict[str, object]] = {}
    cds_by_parent: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for seqid, _source, feature, start, end, _score, strand, _phase, row_attrs in iter_gff_rows(path):
        feature_l = feature.lower()
        row_id = unquote(row_attrs.get("ID", ""))
        if feature_l == "gene" and row_id:
            genes[row_id] = {
                "gene_id": row_id,
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "cds": (),
            }
        elif feature_l in {"mrna", "transcript", "rna"} and row_id:
            parents = [unquote(parent) for parent in split_list(row_attrs.get("Parent", ""))]
            transcripts[row_id] = {
                "parents": parents,
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            }
        elif feature_l == "cds":
            for parent in split_list(row_attrs.get("Parent", "")):
                cds_by_parent[unquote(parent)].append((int(start), int(end)))

    gene_cds: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for tx_id, tx in transcripts.items():
        parents = list(tx["parents"])
        if not parents:
            # A parentless transcript is the top-level stable model in several
            # source annotations. Preserve its ID as the gene-lineage ID.
            genes.setdefault(
                tx_id,
                {
                    "gene_id": tx_id,
                    "seqid": tx["seqid"],
                    "start": tx["start"],
                    "end": tx["end"],
                    "strand": tx["strand"],
                    "cds": (),
                },
            )
            parents = [tx_id]
        for gene_id in parents:
            gene_cds[gene_id].extend(cds_by_parent.get(tx_id, []))
    for gene_id in genes:
        # Also support CDS rows whose Parent directly names a gene.
        gene_cds[gene_id].extend(cds_by_parent.get(gene_id, []))
        genes[gene_id]["cds"] = merge_intervals(gene_cds[gene_id])
    return genes


def collect_proposal_cds(
    proposal_gff: str, records: Sequence[Dict[str, str]]
) -> Dict[str, Tuple[Tuple[int, int], ...]]:
    """Collect unioned CDS bases for each proposal gene from the policy GFF."""
    tx_to_gene = {
        row.get("export_transcript_id", ""): row.get("export_gene_id", "")
        for row in records
        if row.get("export_transcript_id", "") and row.get("export_gene_id", "")
    }
    cds_by_gene: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for _seqid, _source, feature, start, end, _score, _strand, _phase, row_attrs in iter_gff_rows(proposal_gff):
        if feature.lower() != "cds":
            continue
        for parent in split_list(row_attrs.get("Parent", "")):
            gene_id = tx_to_gene.get(unquote(parent), "")
            if gene_id:
                cds_by_gene[gene_id].append((int(start), int(end)))
    return {gene_id: merge_intervals(intervals) for gene_id, intervals in cds_by_gene.items()}


def proposal_spans(records: Sequence[Dict[str, str]]) -> Dict[str, Dict[str, object]]:
    spans: Dict[str, Dict[str, object]] = {}
    for row in records:
        gene_id = row.get("export_gene_id", "")
        if not gene_id:
            continue
        start = int(row.get("start") or 0)
        end = int(row.get("end") or 0)
        if gene_id not in spans:
            spans[gene_id] = {
                "gene_id": gene_id,
                "seqid": row.get("seqid", ""),
                "start": start,
                "end": end,
                "strand": row.get("strand", ""),
            }
        else:
            spans[gene_id]["start"] = min(int(spans[gene_id]["start"]), start)
            spans[gene_id]["end"] = max(int(spans[gene_id]["end"]), end)
    return spans


def interval_set_overlap(left: Sequence[Tuple[int, int]], right: Sequence[Tuple[int, int]]) -> int:
    total = 0
    left_index = right_index = 0
    while left_index < len(left) and right_index < len(right):
        lstart, lend = left[left_index]
        rstart, rend = right[right_index]
        total += interval_overlap(lstart, lend, rstart, rend)
        if lend < rend:
            left_index += 1
        else:
            right_index += 1
    return total


def structural_score(
    proposal: Mapping[str, object], current: Mapping[str, object]
) -> Tuple[int, int, int, int, int]:
    """Rank lineage using CDS contribution, then coordinate-only fallbacks."""
    pstart = int(proposal.get("start", 0))
    pend = int(proposal.get("end", 0))
    cstart = int(current.get("start", 0))
    cend = int(current.get("end", 0))
    same_context = bool(current) and proposal.get("seqid") == current.get("seqid") and (
        proposal.get("strand") in {"", ".", current.get("strand")} or current.get("strand") in {"", "."}
    )
    if not same_context:
        return (0, 0, 0, 0, -10**18)
    proposal_cds = tuple(proposal.get("cds", ()))
    current_cds = tuple(current.get("cds", ()))
    cds_overlap = interval_set_overlap(proposal_cds, current_cds) if proposal_cds and current_cds else 0
    exact_cds = len(set(proposal_cds) & set(current_cds)) if proposal_cds and current_cds else 0
    cds_union = sum(end - start + 1 for start, end in proposal_cds) + sum(
        end - start + 1 for start, end in current_cds
    ) - cds_overlap
    cds_jaccard_ppm = int(cds_overlap * 1_000_000 / cds_union) if cds_union else 0
    span_overlap = interval_overlap(pstart, pend, cstart, cend)
    distance = abs(((pstart + pend) // 2) - ((cstart + cend) // 2))
    return (cds_overlap, exact_cds, cds_jaccard_ppm, span_overlap, -distance)


def add_scores(left: Tuple[int, ...], right: Tuple[int, ...]) -> Tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right))


def optimal_structural_pairs(
    proposal_ids: Sequence[str],
    current_ids: Sequence[str],
    proposals: Mapping[str, Mapping[str, object]],
    current_genes: Mapping[str, Mapping[str, object]],
) -> List[Tuple[str, str, Tuple[int, int, int, int, int]]]:
    """Return an exact maximum-weight one-to-one assignment for one locus.

    All members of the smaller side are assigned. Zero-CDS cases therefore
    receive a deterministic span/position fallback instead of a review state.
    """
    proposal_ids = sorted(proposal_ids, key=lambda item: (int(proposals[item].get("start", 0)), item))
    current_ids = sorted(
        current_ids,
        key=lambda item: (int(current_genes.get(item, {}).get("start", 0)), item),
    )
    proposal_is_small = len(proposal_ids) <= len(current_ids)
    small = proposal_ids if proposal_is_small else current_ids
    large = current_ids if proposal_is_small else proposal_ids
    if not small or not large:
        return []
    if len(large) > 20:
        # This is not expected for accepted locus windows; keep deterministic
        # behavior without exponential memory if an external dataset violates it.
        ranked = []
        for proposal_id in proposal_ids:
            for current_id in current_ids:
                ranked.append((structural_score(proposals[proposal_id], current_genes.get(current_id, {})), proposal_id, current_id))
        used_p: Set[str] = set()
        used_c: Set[str] = set()
        result = []
        for score, proposal_id, current_id in sorted(ranked, key=lambda item: (item[0], item[1], item[2]), reverse=True):
            if proposal_id not in used_p and current_id not in used_c:
                used_p.add(proposal_id)
                used_c.add(current_id)
                result.append((proposal_id, current_id, score))
                if len(result) == min(len(proposal_ids), len(current_ids)):
                    break
        return result

    score_matrix: List[List[Tuple[int, int, int, int, int]]] = []
    for small_id in small:
        row = []
        for large_id in large:
            proposal_id, current_id = (small_id, large_id) if proposal_is_small else (large_id, small_id)
            row.append(structural_score(proposals[proposal_id], current_genes.get(current_id, {})))
        score_matrix.append(row)

    zero = (0, 0, 0, 0, 0)

    @lru_cache(maxsize=None)
    def solve(index: int, used_mask: int) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
        if index == len(small):
            return zero, ()
        best_score: Optional[Tuple[int, ...]] = None
        best_choice: Optional[Tuple[int, ...]] = None
        for large_index in range(len(large)):
            if used_mask & (1 << large_index):
                continue
            suffix_score, suffix_choice = solve(index + 1, used_mask | (1 << large_index))
            total_score = add_scores(score_matrix[index][large_index], suffix_score)
            choice = (large_index,) + suffix_choice
            if best_score is None or total_score > best_score or (total_score == best_score and choice < best_choice):
                best_score, best_choice = total_score, choice
        assert best_score is not None and best_choice is not None
        return best_score, best_choice

    _score, choices = solve(0, 0)
    result = []
    for index, large_index in enumerate(choices):
        proposal_id, current_id = (small[index], large[large_index]) if proposal_is_small else (large[large_index], small[index])
        result.append((proposal_id, current_id, score_matrix[index][large_index]))
    return result


def assign_current_anchors(
    proposal_ids: Sequence[str],
    current_ids: Sequence[str],
    proposals: Mapping[str, Mapping[str, object]],
    current_genes: Mapping[str, Mapping[str, object]],
    relation: str,
    reuse_one_to_one: bool,
) -> Dict[str, Dict[str, str]]:
    anchors: Dict[str, Dict[str, str]] = {}
    if relation == "one_to_one_replacement" and reuse_one_to_one and current_ids:
        anchors[proposal_ids[0]] = {
            "current_id": current_ids[0],
            "basis": "one_to_one_lineage",
            "status": "deterministic_one_to_one_inheritance",
        }
        return anchors
    if relation not in {"split_replacement", "merge_replacement", "complex_replacement"}:
        return anchors

    for proposal_id, current_id, score in optimal_structural_pairs(proposal_ids, current_ids, proposals, current_genes):
        if score[0] > 0:
            basis = relation.replace("_replacement", "") + "_maximum_cds_contribution"
            status = "deterministic_cds_anchor"
        elif score[3] > 0:
            basis = relation.replace("_replacement", "") + "_maximum_span_fallback"
            status = "deterministic_span_fallback"
        else:
            basis = relation.replace("_replacement", "") + "_genomic_order_fallback"
            status = "deterministic_order_fallback"
        anchors[proposal_id] = {
            "current_id": current_id,
            "basis": basis,
            "status": status,
        }
    return anchors


def parse_current_style_gene_id(gene_id: str) -> Tuple[str, int, int, str]:
    """Split an ID around its gene-number token.

    A numeric token immediately following ``G/g`` is preferred because scaffold
    numbers can be wider than gene counters (for example
    ``GH_scaffold49195_objG0001``). Otherwise the longest numeric token captures
    static forms such as CSS and maize IDs while retaining version suffixes.
    """
    matches = list(NUMERIC_TOKEN_RE.finditer(gene_id))
    if not matches:
        raise ReleaseExportStyleError("Current gene ID has no numeric token: {0}".format(gene_id))
    gene_number_matches = [match for match in matches if match.start() > 0 and gene_id[match.start() - 1] in {"G", "g"}]
    match = gene_number_matches[-1] if gene_number_matches else min(matches, key=lambda item: (-len(item.group(0)), item.start()))
    number_text = match.group(0)
    return gene_id[: match.start()], int(number_text), len(number_text), gene_id[match.end() :]


class ReleaseExportStyleError(ValueError):
    pass


def source_namespace_is_reserved(prefix: str) -> bool:
    """Avoid inventing identifiers in externally assigned registry namespaces."""
    return prefix.upper().endswith("LOC") or "GENE-LOC" in prefix.upper()


def numeric_step(numbers: Sequence[int]) -> int:
    unique = sorted(set(numbers))
    differences = [right - left for left, right in zip(unique, unique[1:]) if right > left]
    if not differences:
        return 1
    step = differences[0]
    for difference in differences[1:]:
        step = gcd(step, difference)
        if step == 1:
            break
    return max(1, step)


def allocate_source_style_ids(
    unassigned: Sequence[str],
    proposals: Mapping[str, Mapping[str, object]],
    current_genes: Mapping[str, Mapping[str, object]],
    reserved_gene_ids: Set[str],
) -> Dict[str, Dict[str, str]]:
    """Infer source gene-ID templates without species-specific configuration.

    Local seqid styles are preferred. If a source embeds the exact seqid in its
    gene ID (for example ``Bna{seqid}G#######ZS``), that relationship is learned
    across the annotation and can be instantiated for a previously unannotated
    seqid. A single dominant static style (for example ``CSS#######``) is also
    reusable. NCBI LOC identifiers are externally assigned and are never minted.
    """
    parsed: List[Dict[str, object]] = []
    local_counts: Dict[str, Counter] = defaultdict(Counter)
    global_counts: Counter = Counter()
    dynamic_counts: Counter = Counter()
    dynamic_numbers: Dict[Tuple[str, int, str], List[int]] = defaultdict(list)
    dynamic_starts: Dict[Tuple[str, int, str], List[int]] = defaultdict(list)
    for gene_id, model in current_genes.items():
        try:
            prefix, number, width, suffix = parse_current_style_gene_id(gene_id)
        except ReleaseExportStyleError:
            continue
        seqid = str(model.get("seqid", ""))
        schema = (prefix, width, suffix)
        row = {"gene_id": gene_id, "seqid": seqid, "prefix": prefix, "number": number, "width": width, "suffix": suffix}
        parsed.append(row)
        local_counts[seqid][schema] += 1
        global_counts[schema] += 1
        if seqid and seqid in prefix:
            dynamic = (prefix.replace(seqid, "{seqid}", 1), width, suffix)
            dynamic_counts[dynamic] += 1
            dynamic_numbers[dynamic].append(number)
            dynamic_starts[dynamic].append(number)
        elif seqid and seqid in suffix:
            dynamic = (prefix, width, suffix.replace(seqid, "{seqid}", 1))
            dynamic_counts[dynamic] += 1
            dynamic_numbers[dynamic].append(number)
            dynamic_starts[dynamic].append(number)

    dominant_static: Optional[Tuple[str, int, str]] = None
    if global_counts:
        candidate, count = global_counts.most_common(1)[0]
        if count / sum(global_counts.values()) >= 0.80:
            dominant_static = candidate
    dominant_dynamic: Optional[Tuple[str, int, str]] = None
    if dynamic_counts:
        candidate, count = dynamic_counts.most_common(1)[0]
        if count >= 10:
            dominant_dynamic = candidate

    parsed_reserved: Dict[Tuple[str, int, str], Set[int]] = defaultdict(set)
    for gene_id in reserved_gene_ids:
        try:
            prefix, number, width, suffix = parse_current_style_gene_id(gene_id)
        except ReleaseExportStyleError:
            continue
        parsed_reserved[(prefix, width, suffix)].add(number)

    grouped: Dict[Tuple[str, int, str, str], List[str]] = defaultdict(list)
    for proposal_id in unassigned:
        seqid = str(proposals[proposal_id].get("seqid", ""))
        schema: Optional[Tuple[str, int, str]] = None
        basis = ""
        if local_counts.get(seqid):
            schema = local_counts[seqid].most_common(1)[0][0]
            basis = "inferred_local_source_id_style"
        elif dominant_dynamic is not None:
            pattern_prefix, width, pattern_suffix = dominant_dynamic
            schema = (pattern_prefix.format(seqid=seqid), width, pattern_suffix.format(seqid=seqid))
            basis = "inferred_seqid_source_id_style"
        elif dominant_static is not None:
            schema = dominant_static
            basis = "inferred_dominant_source_id_style"
        if schema is None or source_namespace_is_reserved(schema[0]):
            continue
        grouped[(schema[0], schema[1], schema[2], basis)].append(proposal_id)

    result: Dict[str, Dict[str, str]] = {}
    for group_key, proposal_ids in sorted(grouped.items()):
        prefix, width, suffix, basis = group_key
        schema = (prefix, width, suffix)
        used = parsed_reserved[schema]
        current_numbers = [int(row["number"]) for row in parsed if (row["prefix"], row["width"], row["suffix"]) == schema]
        if current_numbers:
            step = numeric_step(current_numbers)
            candidate = max(current_numbers) + step
        else:
            if basis == "inferred_seqid_source_id_style" and dominant_dynamic is not None:
                step = numeric_step(dynamic_numbers[dominant_dynamic])
                candidate = min(dynamic_starts[dominant_dynamic])
            else:
                step = 1
                candidate = 1
        limit = 10**width - 1
        for proposal_id in sorted(proposal_ids, key=lambda item: (str(proposals[item].get("seqid", "")), int(proposals[item].get("start", 0)), item)):
            while candidate in used and candidate <= limit:
                candidate += step
            if candidate > limit:
                break
            final_id = "{0}{1:0{2}d}{3}".format(prefix, candidate, width, suffix)
            used.add(candidate)
            result[proposal_id] = {"final_gene_id": final_id, "basis": basis}
            candidate += step
    return result


def nearest_free_number(target: float, available: Set[int]) -> int:
    if not available:
        raise ReleaseExportStyleError("No free numeric ID remains between neighboring current genes")
    for modulus in (10, 5, 1):
        candidates = [number for number in available if number % modulus == 0]
        if candidates:
            return min(candidates, key=lambda number: (abs(number - target), number))
    raise AssertionError("unreachable")


def allocate_positional_ids(
    unassigned: Sequence[str],
    proposals: Mapping[str, Mapping[str, object]],
    current_genes: Mapping[str, Mapping[str, object]],
    reserved_gene_ids: Set[str],
    strict: bool = True,
) -> Dict[str, str]:
    """Allocate current-style IDs in genomic gaps without recycling old IDs."""
    current_by_seqid: Dict[str, List[Tuple[int, int, str, str, int, str]]] = defaultdict(list)
    for current_id, model in current_genes.items():
        try:
            prefix, number, width, suffix = parse_current_style_gene_id(current_id)
        except ReleaseExportStyleError:
            continue
        position = int(model["start"])
        current_by_seqid[str(model["seqid"])].append((position, number, prefix, current_id, width, suffix))
    for seqid in current_by_seqid:
        current_by_seqid[seqid].sort()
    positions_by_seqid = {seqid: [anchor[0] for anchor in anchors] for seqid, anchors in current_by_seqid.items()}

    namespace_quality: Dict[Tuple[str, str, str], Tuple[int, int]] = {}
    for seqid, anchors in current_by_seqid.items():
        by_namespace: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
        for position, number, prefix, _current_id, _width, suffix in anchors:
            by_namespace[(prefix, suffix)].append((position, number))
        for namespace, values in by_namespace.items():
            values.sort()
            adjacent = list(zip(values, values[1:]))
            increasing = sum(1 for (left_pos, left_num), (right_pos, right_num) in adjacent if right_pos >= left_pos and right_num > left_num)
            namespace_quality[(seqid, namespace[0], namespace[1])] = (increasing, len(adjacent))

    parsed_reserved: Dict[Tuple[str, str], Set[int]] = defaultdict(set)
    for gene_id in reserved_gene_ids:
        try:
            prefix, number, _width, suffix = parse_current_style_gene_id(gene_id)
        except ReleaseExportStyleError:
            continue
        parsed_reserved[(prefix, suffix)].add(number)

    gap_groups: Dict[Tuple[str, int, int], List[str]] = defaultdict(list)
    gap_meta: Dict[Tuple[str, int, int], Tuple[str, int, str, int, int]] = {}
    for proposal_id in unassigned:
        proposal = proposals[proposal_id]
        seqid = str(proposal["seqid"])
        anchors = current_by_seqid.get(seqid, [])
        if not anchors:
            if strict:
                raise ReleaseExportStyleError("No parseable current-style gene IDs on seqid {0}".format(seqid))
            continue
        proposal_position = int(proposal["start"])
        insertion = bisect_right(positions_by_seqid[seqid], proposal_position)
        left_index = insertion - 1
        right_index = insertion
        key = (seqid, left_index, right_index)
        gap_groups[key].append(proposal_id)
        reference = anchors[left_index] if left_index >= 0 else anchors[right_index]
        gap_meta[key] = (reference[2], reference[4], reference[5], left_index, right_index)

    result: Dict[str, str] = {}
    for key, proposal_ids in sorted(gap_groups.items()):
        try:
            seqid, _left_key, _right_key = key
            anchors = current_by_seqid[seqid]
            prefix, width, suffix, left_index, right_index = gap_meta[key]
            proposal_ids = sorted(proposal_ids, key=lambda item: (int(proposals[item]["start"]), int(proposals[item]["end"]), item))
            namespace = (prefix, suffix)
            increasing, adjacent = namespace_quality.get((seqid, prefix, suffix), (0, 0))
            if not strict and (adjacent < 2 or increasing / adjacent < 0.98):
                raise ReleaseExportStyleError("Current numeric namespace is not reliably positional on {0}".format(seqid))
            if left_index >= 0 and (anchors[left_index][2], anchors[left_index][5]) != namespace:
                raise ReleaseExportStyleError("Mixed current ID prefixes around {0}".format(seqid))
            if right_index < len(anchors) and (anchors[right_index][2], anchors[right_index][5]) != namespace:
                raise ReleaseExportStyleError("Mixed current ID prefixes around {0}".format(seqid))
            used_numbers = parsed_reserved[namespace]
            if left_index >= 0 and right_index < len(anchors):
                lower = anchors[left_index][1]
                upper = anchors[right_index][1]
                if upper <= lower:
                    raise ReleaseExportStyleError("Non-increasing current numeric IDs on {0}: {1}, {2}".format(seqid, lower, upper))
                if upper - lower > 1_000_000:
                    raise ReleaseExportStyleError("Implausibly large numeric namespace gap on {0}: {1}, {2}".format(seqid, lower, upper))
                available = set(range(lower + 1, upper)) - used_numbers
                if len(available) < len(proposal_ids):
                    raise ReleaseExportStyleError("Insufficient unused IDs between {0} and {1}".format(anchors[left_index][3], anchors[right_index][3]))
                candidates = []
                for index, _proposal_id in enumerate(proposal_ids, start=1):
                    target = lower + (upper - lower) * index / (len(proposal_ids) + 1)
                    number = nearest_free_number(target, available)
                    available.remove(number)
                    candidates.append(number)
            else:
                numeric_ids = sorted(
                    number
                    for _mid, number, item_prefix, _gid, _width, item_suffix in anchors
                    if (item_prefix, item_suffix) == namespace
                )
                differences = [right - left for left, right in zip(numeric_ids, numeric_ids[1:]) if right > left]
                step = sorted(differences)[len(differences) // 2] if differences else 100
                if right_index < len(anchors):
                    base = anchors[right_index][1]
                    available = set(range(1, base)) - used_numbers
                    if len(available) < len(proposal_ids):
                        raise ReleaseExportStyleError("Insufficient positive IDs before the first current gene on {0}".format(seqid))
                    candidates = []
                    for index in range(1, len(proposal_ids) + 1):
                        target = base * index / (len(proposal_ids) + 1)
                        number = nearest_free_number(target, available)
                        available.remove(number)
                        candidates.append(number)
                else:
                    base = anchors[left_index][1]
                    candidates = []
                    candidate = base
                    limit = 10**width - 1
                    for _proposal_id in proposal_ids:
                        candidate += step
                        while candidate in used_numbers and candidate <= limit:
                            candidate += 10
                        if candidate > limit:
                            raise ReleaseExportStyleError("Numeric ID width exhausted after the last current gene on {0}".format(seqid))
                        candidates.append(candidate)
            for proposal_id, number in zip(proposal_ids, candidates):
                used_numbers.add(number)
                result[proposal_id] = "{0}{1:0{2}d}{3}".format(prefix, number, width, suffix)
        except ReleaseExportStyleError:
            if strict:
                raise
    return result


def build_gene_mapping(
    records: Sequence[Dict[str, str]],
    existing_gene_ids: Set[str],
    one_to_one_policy: str,
    new_gene_prefix: str,
    new_gene_start: int,
    new_gene_width: int,
    id_strategy: str = "sequential",
    current_genes: Optional[Mapping[str, Mapping[str, object]]] = None,
    proposal_cds: Optional[Mapping[str, Sequence[Tuple[int, int]]]] = None,
) -> Tuple[Dict[str, Dict[str, str]], List[Dict[str, str]], int]:
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_locus: Dict[str, List[str]] = defaultdict(list)
    for row in records:
        export_gene_id = row.get("export_gene_id", "")
        if not export_gene_id:
            continue
        by_gene[export_gene_id].append(row)
        locus_id = row.get("locus_id", "")
        if export_gene_id not in by_locus[locus_id]:
            by_locus[locus_id].append(export_gene_id)

    current_by_locus: Dict[str, List[str]] = {}
    action_by_locus: Dict[str, str] = {}
    for locus_id, gene_ids in by_locus.items():
        current_ids: Set[str] = set()
        actions: Set[str] = set()
        for gene_id in gene_ids:
            for row in by_gene[gene_id]:
                current_ids.update(split_list(row.get("removed_current_gene_ids") or row.get("current_gene_ids")))
                actions.add(row.get("proposal_action", ""))
        current_by_locus[locus_id] = sorted(current_ids)
        if "replace_current_with_selected_set" in actions:
            action_by_locus[locus_id] = "replace_current_with_selected_set"
        elif "add_novel_selected_set" in actions:
            action_by_locus[locus_id] = "add_novel_selected_set"
        else:
            action_by_locus[locus_id] = sorted(actions)[0] if actions else ""

    current_genes = dict(current_genes or {})
    proposals = proposal_spans(records)
    for proposal_id, intervals in (proposal_cds or {}).items():
        if proposal_id in proposals:
            proposals[proposal_id]["cds"] = tuple(intervals)
    ordered_gene_ids = proposal_gene_order(records)
    anchor_by_gene: Dict[str, Dict[str, str]] = {}
    for locus_id, locus_gene_ids in by_locus.items():
        current_ids = current_by_locus.get(locus_id, [])
        relation = relation_for_locus(action_by_locus.get(locus_id, ""), len(current_ids), len(locus_gene_ids))
        anchor_by_gene.update(
            assign_current_anchors(
                locus_gene_ids,
                current_ids,
                proposals,
                current_genes,
                relation,
                one_to_one_policy == "reuse_current_gene_id",
            )
        )

    reserved_gene_ids = set(existing_gene_ids) | set(current_genes)
    unassigned = [gene_id for gene_id in ordered_gene_ids if gene_id not in anchor_by_gene]
    positional_ids: Dict[str, str] = {}
    source_style_ids: Dict[str, Dict[str, str]] = {}
    if id_strategy in {"current_style_positional", "current_style_or_sequential"} and unassigned:
        if not current_genes:
            raise SystemExit("--id-strategy current_style_positional requires --current-gff with gene records")
        try:
            positional_ids = allocate_positional_ids(
                unassigned,
                proposals,
                current_genes,
                reserved_gene_ids,
                strict=id_strategy == "current_style_positional",
            )
        except ReleaseExportStyleError as exc:
            raise SystemExit("Current-style positional ID allocation failed: {0}".format(exc))
    if id_strategy == "current_style_or_sequential":
        source_style_ids = allocate_source_style_ids(
            [gene_id for gene_id in unassigned if gene_id not in positional_ids],
            proposals,
            current_genes,
            reserved_gene_ids | set(positional_ids.values()),
        )

    mapping: Dict[str, Dict[str, str]] = {}
    rows: List[Dict[str, str]] = []
    assigned_final_ids: Set[str] = set()
    next_index = new_gene_start

    for export_gene_id in ordered_gene_ids:
        gene_rows = by_gene[export_gene_id]
        first = gene_rows[0]
        locus_id = first.get("locus_id", "")
        current_ids = current_by_locus.get(locus_id, [])
        proposal_count = len(by_locus.get(locus_id, []))
        action = action_by_locus.get(locus_id, first.get("proposal_action", ""))
        relation = relation_for_locus(action, len(current_ids), proposal_count)

        anchor = anchor_by_gene.get(export_gene_id, {})
        reuse_current = bool(anchor)
        if reuse_current:
            final_gene_id = anchor["current_id"]
            assignment_basis = anchor["basis"]
            assignment_status = anchor["status"]
        elif export_gene_id in positional_ids:
            final_gene_id = positional_ids[export_gene_id]
            assignment_basis = "current_namespace_genomic_gap"
            assignment_status = "deterministic_new_position_id"
        elif export_gene_id in source_style_ids:
            final_gene_id = source_style_ids[export_gene_id]["final_gene_id"]
            assignment_basis = source_style_ids[export_gene_id]["basis"]
            assignment_status = "deterministic_new_source_style_id"
        else:
            fallback_prefix = "GeneArbiterG" if new_gene_prefix.strip().lower() == "auto" else new_gene_prefix
            while True:
                final_gene_id = make_new_gene_id(fallback_prefix, next_index, new_gene_width)
                next_index += 1
                if final_gene_id not in reserved_gene_ids and final_gene_id not in assigned_final_ids:
                    break
            assignment_basis = "sequential_new_namespace"
            assignment_status = "deterministic_new_namespace_id"

        if final_gene_id in assigned_final_ids or (not reuse_current and final_gene_id in reserved_gene_ids):
            raise SystemExit(
                "Final gene ID collision for {0}: {1}. "
                "Review lineage anchors or ID allocation settings.".format(export_gene_id, final_gene_id)
            )
        assigned_final_ids.add(final_gene_id)

        retained_current = {item.get("current_id", "") for item in anchor_by_gene.values() if item.get("current_id", "") in current_ids}
        retired_current = [current_id for current_id in current_ids if current_id not in retained_current]
        source_gene_ids = sorted({item for row in gene_rows for item in split_list(row.get("source_gene_id"))})
        source_model_ids = sorted({item for row in gene_rows for item in split_list(row.get("source_model_id"))})
        sources = sorted({row.get("source", "") for row in gene_rows if row.get("source", "")})
        seqid = first.get("seqid", "")
        start = min(int(row.get("start") or 0) for row in gene_rows)
        end = max(int(row.get("end") or 0) for row in gene_rows)

        mapping[export_gene_id] = {
            "final_gene_id": final_gene_id,
            "final_relation": relation,
            "reuse_current_gene_id": str(reuse_current).lower(),
            "current_gene_ids": ";".join(current_ids),
            "retired_current_gene_ids": ";".join(retired_current),
            "id_assignment_basis": assignment_basis,
            "id_assignment_status": assignment_status,
            "anchor_current_gene_id": anchor.get("current_id", ""),
        }
        rows.append(
            {
                "locus_id": locus_id,
                "proposal_action": action,
                "final_relation": relation,
                "reuse_current_gene_id": str(reuse_current).lower(),
                "final_gene_id": final_gene_id,
                "proposal_gene_id": public_proposal_id(export_gene_id),
                "current_gene_ids": ";".join(current_ids),
                "retired_current_gene_ids": ";".join(retired_current),
                "source": ",".join(sources),
                "source_gene_id": ",".join(source_gene_ids),
                "source_model_id": ",".join(source_model_ids),
                "selected_set_id": first.get("selected_set_id", ""),
                "selected_source": first.get("selected_source", ""),
                "seqid": seqid,
                "start": str(start),
                "end": str(end),
                "strand": first.get("strand", ""),
                "id_assignment_basis": assignment_basis,
                "id_assignment_status": assignment_status,
                "anchor_current_gene_id": anchor.get("current_id", ""),
                "notes": "current_gene_id_reused" if reuse_current else "new_final_gene_id_allocated",
            }
        )
    return mapping, rows, next_index


def build_tx_mapping(
    records: Sequence[Dict[str, str]],
    gene_mapping: Dict[str, Dict[str, str]],
    existing_tx_ids: Set[str],
    transcript_template: str,
    transcript_id_strategy: str = "template",
) -> Tuple[Dict[str, Dict[str, str]], List[Dict[str, str]]]:
    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in records:
        export_gene_id = row.get("export_gene_id", "")
        export_tx_id = row.get("export_transcript_id", "")
        if export_gene_id and export_tx_id:
            by_gene[export_gene_id].append(row)

    tx_mapping: Dict[str, Dict[str, str]] = {}
    rows: List[Dict[str, str]] = []
    used_tx_ids = set(existing_tx_ids)
    for export_gene_id, gene_rows in sorted(by_gene.items(), key=lambda item: gene_mapping[item[0]]["final_gene_id"]):
        final_gene_id = gene_mapping[export_gene_id]["final_gene_id"]
        seen_tx: Set[str] = set()
        tx_rows = []
        for row in sorted(gene_rows, key=lambda item: (int(item.get("start") or 0), str(item.get("export_transcript_id", "")))):
            export_tx_id = row.get("export_transcript_id", "")
            if export_tx_id and export_tx_id not in seen_tx:
                seen_tx.add(export_tx_id)
                tx_rows.append(row)
        for tx_index, row in enumerate(tx_rows, start=1):
            export_tx_id = row.get("export_transcript_id", "")
            if transcript_id_strategy == "current_style":
                if not re.search(r"G\d+$", final_gene_id):
                    raise SystemExit("Cannot derive current-style transcript ID from gene ID: {0}".format(final_gene_id))
                transcript_stem = re.sub(r"G(?=\d+$)", "T", final_gene_id)
                final_tx_id = "{0}.{1}".format(transcript_stem, tx_index)
            else:
                final_tx_id = transcript_template.format(gene_id=final_gene_id, index=tx_index)
            if final_tx_id in used_tx_ids:
                raise SystemExit(
                    "Final transcript ID collision for {0}: {1}. Adjust --transcript-template.".format(export_tx_id, final_tx_id)
                )
            used_tx_ids.add(final_tx_id)
            tx_mapping[export_tx_id] = {
                "final_transcript_id": final_tx_id,
                "final_gene_id": final_gene_id,
                "final_relation": gene_mapping[export_gene_id]["final_relation"],
                "proposal_gene_id": public_proposal_id(export_gene_id),
            }
            rows.append(
                {
                    "locus_id": row.get("locus_id", ""),
                    "final_relation": gene_mapping[export_gene_id]["final_relation"],
                    "final_gene_id": final_gene_id,
                    "final_transcript_id": final_tx_id,
                    "proposal_gene_id": public_proposal_id(export_gene_id),
                    "proposal_transcript_id": public_proposal_id(export_tx_id),
                    "current_gene_ids": gene_mapping[export_gene_id]["current_gene_ids"],
                    "source": row.get("source", ""),
                    "source_gene_id": row.get("source_gene_id", ""),
                    "source_transcript_id": row.get("source_transcript_id", ""),
                    "source_model_id": row.get("source_model_id", ""),
                    "selected_set_id": row.get("selected_set_id", ""),
                    "selected_source": row.get("selected_source", ""),
                    "seqid": row.get("seqid", ""),
                    "start": row.get("start", ""),
                    "end": row.get("end", ""),
                    "strand": row.get("strand", ""),
                }
            )
    return tx_mapping, rows


def build_locus_mapping(gene_rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    by_locus: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in gene_rows:
        by_locus[row.get("locus_id", "")].append(row)
    result: List[Dict[str, str]] = []
    for locus_id, rows in sorted(by_locus.items(), key=lambda item: (item[1][0].get("seqid", ""), min(int(row.get("start") or 0) for row in item[1]), item[0])):
        current_ids = sorted({item for row in rows for item in split_list(row.get("current_gene_ids"))})
        final_ids = [row.get("final_gene_id", "") for row in sorted(rows, key=lambda row: (int(row.get("start") or 0), row.get("final_gene_id", "")))]
        retained = sorted({row.get("anchor_current_gene_id", "") for row in rows if row.get("anchor_current_gene_id", "")})
        retired = sorted(set(current_ids) - set(retained))
        added = [gene_id for gene_id in final_ids if gene_id not in current_ids]
        statuses = sorted({row.get("id_assignment_status", "") for row in rows if row.get("id_assignment_status", "")})
        bases = sorted({row.get("id_assignment_basis", "") for row in rows if row.get("id_assignment_basis", "")})
        result.append(
            {
                "locus_id": locus_id,
                "proposal_action": rows[0].get("proposal_action", ""),
                "final_relation": rows[0].get("final_relation", ""),
                "seqid": rows[0].get("seqid", ""),
                "locus_start": str(min(int(row.get("start") or 0) for row in rows)),
                "locus_end": str(max(int(row.get("end") or 0) for row in rows)),
                "current_gene_ids": ";".join(current_ids),
                "corrected_gene_ids": ";".join(final_ids),
                "retained_current_gene_ids": ";".join(retained),
                "retired_current_gene_ids": ";".join(retired),
                "added_gene_ids": ";".join(added),
                "mapping_status": ";".join(statuses),
                "mapping_basis": ";".join(bases),
            }
        )
    return result


def update_attr_id(row_attrs: MutableMapping[str, str], new_id: str, old_key: str) -> None:
    old_id = row_attrs.get("ID", "")
    if old_id:
        row_attrs[old_key] = public_proposal_id(old_id)
    row_attrs["ID"] = new_id
    if row_attrs.get("Name", "") == old_id:
        row_attrs["Name"] = new_id


CLEAN_REMOVE_ATTRS = {
    "proposal_catalog",
    "proposal_status",
    "review_flags",
    "proposal_action",
    "final_call",
    "call_group",
    "relation_to_current",
    "export_readiness",
    "proposal_gene_id",
    "proposal_transcript_id",
    "proposal_feature_id",
    "risk_tags",
    "blocking_reasons",
}


def write_clean_auto_pass_gff(source_gff: Path, clean_gff: Path, clean_source_label: str) -> Counter:
    counts: Counter = Counter()
    with source_gff.open(encoding="utf-8") as in_handle, clean_gff.open("w", encoding="utf-8") as out_handle:
        out_handle.write("##gff-version 3\n")
        for line in in_handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                counts["malformed_gff_rows"] += 1
                out_handle.write(line)
                continue
            if parts[1] in {"AI_annotation_proposal", "GeneArbiter"}:
                parts[1] = clean_source_label
                counts["proposal_source_relabelled"] += 1
            row_attrs = parse_attrs(parts[8])
            for key in list(row_attrs):
                if key in CLEAN_REMOVE_ATTRS:
                    del row_attrs[key]
                    counts[f"removed_attr:{key}"] += 1
            parts[8] = attrs(row_attrs.items())
            out_handle.write("\t".join(parts) + "\n")
            counts["written_gff_rows"] += 1
    return counts


GTF_TRANSCRIPT_FEATURES = {"mrna", "transcript", "lnc_rna"}
GTF_CHILD_FEATURES = {"exon", "cds"}
GTF_ORGANELLE_GENOMES = {"chloroplast", "mitochondrion"}


def gtf_escape(value: object) -> str:
    return str(value or "").replace("\\", "\\\\").replace('"', '\"')


def gtf_attr_text(items: Iterable[Tuple[str, object]]) -> str:
    return " ".join('{0} "{1}";'.format(key, gtf_escape(value)) for key, value in items if value not in {"", None})


def collect_gtf_context(source_gff: Path, exclude_organelle: bool) -> Tuple[Set[str], Dict[str, List[str]], Set[str]]:
    excluded_seqids: Set[str] = set()
    transcript_to_gene: Dict[str, List[str]] = {}
    allowed_genes: Set[str] = set()
    with source_gff.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature, _start, _end, _score, _strand, _phase, attr_text = parts
            row_attrs = parse_attrs(attr_text)
            if exclude_organelle and feature.lower() == "region" and row_attrs.get("genome", "").lower() in GTF_ORGANELLE_GENOMES:
                excluded_seqids.add(seqid)
                continue
            if seqid in excluded_seqids:
                continue
            if feature.lower() in GTF_TRANSCRIPT_FEATURES:
                tx_id = row_attrs.get("ID", "")
                parents = split_list(row_attrs.get("Parent", ""))
                gene_ids = [parent for parent in parents if parent]
                if tx_id and gene_ids:
                    transcript_to_gene[tx_id] = gene_ids
                    allowed_genes.update(gene_ids)
    return excluded_seqids, transcript_to_gene, allowed_genes


def write_hisat_gtf(source_gff: Path, gtf_path: Path, exclude_organelle: bool) -> Counter:
    excluded_seqids, transcript_to_gene, allowed_genes = collect_gtf_context(source_gff, exclude_organelle)
    counts: Counter = Counter()
    with source_gff.open(encoding="utf-8") as in_handle, gtf_path.open("w", encoding="utf-8") as out_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                counts["malformed_gff_rows"] += 1
                continue
            seqid, source, feature, start, end, score, strand, phase, attr_text = parts
            if seqid in excluded_seqids:
                counts["excluded_organelle_rows"] += 1
                continue
            feature_l = feature.lower()
            row_attrs = parse_attrs(attr_text)
            row_id = row_attrs.get("ID", "")
            out_feature = feature
            attr_items: List[Tuple[str, object]] = []
            if feature_l == "gene":
                if row_id not in allowed_genes:
                    counts["skipped_gene_without_selected_transcript"] += 1
                    continue
                attr_items = [("gene_id", row_id), ("gene_name", row_attrs.get("Name") or row_attrs.get("gene") or row_id)]
            elif feature_l in GTF_TRANSCRIPT_FEATURES:
                gene_ids = transcript_to_gene.get(row_id, [])
                if not gene_ids:
                    counts["skipped_transcript_without_gene"] += 1
                    continue
                out_feature = "transcript"
                attr_items = [
                    ("gene_id", gene_ids[0]),
                    ("transcript_id", row_id),
                    ("gene_name", row_attrs.get("gene") or gene_ids[0]),
                ]
            elif feature_l in GTF_CHILD_FEATURES:
                parents = [parent for parent in split_list(row_attrs.get("Parent", "")) if parent in transcript_to_gene]
                if not parents:
                    counts["skipped_child_without_selected_transcript"] += 1
                    continue
                for parent in parents:
                    gene_ids = transcript_to_gene[parent]
                    child_attrs: List[Tuple[str, object]] = [
                        ("gene_id", gene_ids[0]),
                        ("transcript_id", parent),
                    ]
                    exon_number = row_attrs.get("exon_number") or row_attrs.get("number")
                    if exon_number:
                        child_attrs.append(("exon_number", exon_number))
                    out_handle.write(
                        "\t".join([seqid, source, feature, start, end, score or ".", strand or ".", phase or ".", gtf_attr_text(child_attrs)])
                        + "\n"
                    )
                    counts[f"written_{feature_l}_rows"] += 1
                continue
            else:
                counts[f"skipped_feature:{feature_l}"] += 1
                continue
            out_handle.write("\t".join([seqid, source, out_feature, start, end, score or ".", strand or ".", phase or ".", gtf_attr_text(attr_items)]) + "\n")
            counts[f"written_{out_feature.lower()}_rows"] += 1
    counts["excluded_organelle_seqids"] = len(excluded_seqids)
    counts["selected_transcripts"] = len(transcript_to_gene)
    counts["selected_genes"] = len(allowed_genes)
    return counts


def rewrite_current_backbone_attrs(
    row_attrs: MutableMapping[str, str],
    feature: str,
    seqid: str,
    start: str,
    end: str,
    strand: str,
    current_id_counts: Counter,
    current_latest_id: Dict[str, str],
    current_used_ids: Set[str],
    current_disambig_rows: List[Dict[str, str]],
) -> Counter:
    counts: Counter = Counter()
    feature_l = feature.lower()
    original_parents = split_list(row_attrs.get("Parent", ""))
    final_parents = [current_latest_id.get(parent, parent) for parent in original_parents]
    if final_parents != original_parents:
        row_attrs["Parent"] = ",".join(final_parents)
        counts["renamed_current_parent_refs"] += 1

    original_id = row_attrs.get("ID", "")
    if not original_id:
        return counts

    if feature_l not in {"gene", "mrna", "transcript", "rna"}:
        return counts

    current_id_counts[original_id] += 1
    duplicate_index = current_id_counts[original_id]
    final_id = original_id
    reason = ""
    if duplicate_index > 1:
        final_id = "{0}.dup{1:02d}".format(original_id, duplicate_index)
        while final_id in current_used_ids:
            duplicate_index += 1
            final_id = "{0}.dup{1:02d}".format(original_id, duplicate_index)
        row_attrs["original_current_feature_id"] = original_id
        row_attrs["current_duplicate_index"] = str(duplicate_index)
        row_attrs["ID"] = final_id
        if row_attrs.get("Name", "") == original_id:
            row_attrs["Name"] = final_id
        counts["renamed_current_duplicate_ids"] += 1
        if feature_l == "gene":
            counts["renamed_current_duplicate_gene_ids"] += 1
        reason = "duplicate_current_gene_or_transcript_id"

    current_latest_id[original_id] = final_id
    current_used_ids.add(final_id)
    if reason:
        current_disambig_rows.append(
            {
                "original_id": original_id,
                "final_id": final_id,
                "feature": feature,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "duplicate_index": str(duplicate_index),
                "original_parent": ";".join(original_parents),
                "final_parent": ";".join(final_parents),
                "reason": reason,
            }
        )
    return counts


def rewrite_gff(
    proposal_gff: str,
    out_gff: str,
    gene_mapping: Dict[str, Dict[str, str]],
    tx_mapping: Dict[str, Dict[str, str]],
    proposal_source_label: str = "",
) -> Tuple[Counter, List[Dict[str, str]]]:
    counts: Counter = Counter()
    child_index: Dict[Tuple[str, str], int] = defaultdict(int)
    current_id_counts: Counter = Counter()
    current_latest_id: Dict[str, str] = {}
    current_used_ids: Set[str] = set()
    current_disambig_rows: List[Dict[str, str]] = []
    with open_text(proposal_gff) as in_handle, open(out_gff, "w") as out_handle:
        out_handle.write("##gff-version 3\n")
        out_handle.write("# generated_by=GeneArbiter reconcile_ids\n")
        out_handle.write("# generated_at_utc={0}\n".format(datetime.now(timezone.utc).isoformat()))
        for line in in_handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    continue
                out_handle.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                counts["malformed_gff_rows"] += 1
                out_handle.write(line)
                continue
            seqid, source, feature, start, end, score, strand, phase, attr_text = parts
            row_attrs = parse_attrs(attr_text)
            feature_l = feature.lower()
            parents = split_list(row_attrs.get("Parent", ""))
            is_proposal = row_attrs.get("proposal_catalog", "").lower() == "true"
            if feature_l in {"exon", "cds"} and any(parent in tx_mapping for parent in parents):
                is_proposal = True
            if not is_proposal:
                counts.update(
                    rewrite_current_backbone_attrs(
                        row_attrs,
                        feature,
                        seqid,
                        start,
                        end,
                        strand,
                        current_id_counts,
                        current_latest_id,
                        current_used_ids,
                        current_disambig_rows,
                    )
                )
                counts["kept_current_or_nonproposal_records"] += 1
                out_handle.write(gff_line(seqid, source, feature, start, end, score, strand, phase, attrs(row_attrs.items())) + "\n")
                continue

            if proposal_source_label:
                source = proposal_source_label
                counts["relabeled_proposal_source"] += 1

            if feature_l == "gene":
                export_gene_id = row_attrs.get("ID", "")
                mapping = gene_mapping.get(export_gene_id)
                if not mapping:
                    counts["proposal_gene_without_mapping"] += 1
                else:
                    update_attr_id(row_attrs, mapping["final_gene_id"], "proposal_gene_id")
                    row_attrs["final_id_relation"] = mapping["final_relation"]
                    row_attrs["reuse_current_gene_id"] = mapping["reuse_current_gene_id"]
                    row_attrs["current_gene_ids"] = mapping["current_gene_ids"]
                    row_attrs["retired_current_gene_ids"] = mapping["retired_current_gene_ids"]
                    counts["renamed_proposal_gene_records"] += 1
            elif feature_l in {"mrna", "transcript", "rna"}:
                export_tx_id = row_attrs.get("ID", "")
                mapping = tx_mapping.get(export_tx_id)
                if not mapping:
                    counts["proposal_transcript_without_mapping"] += 1
                else:
                    update_attr_id(row_attrs, mapping["final_transcript_id"], "proposal_transcript_id")
                    row_attrs["Parent"] = mapping["final_gene_id"]
                    row_attrs["final_id_relation"] = mapping["final_relation"]
                    row_attrs["proposal_gene_id"] = mapping["proposal_gene_id"]
                    counts["renamed_proposal_transcript_records"] += 1
            elif feature_l in {"exon", "cds"}:
                if len(parents) != 1 or parents[0] not in tx_mapping:
                    counts["proposal_child_without_mapping"] += 1
                else:
                    old_id = row_attrs.get("ID", "")
                    final_tx_id = tx_mapping[parents[0]]["final_transcript_id"]
                    child_index[(final_tx_id, feature_l)] += 1
                    row_attrs["Parent"] = final_tx_id
                    if old_id:
                        row_attrs["proposal_feature_id"] = public_proposal_id(old_id)
                    row_attrs["ID"] = "{0}.{1}{2}".format(final_tx_id, feature_l, child_index[(final_tx_id, feature_l)])
                    counts["renamed_proposal_child_records"] += 1
            else:
                counts["other_proposal_records"] += 1

            out_handle.write(gff_line(seqid, source, feature, start, end, score, strand, phase, attrs(row_attrs.items())) + "\n")
    return counts, current_disambig_rows


def guard_outputs(paths: Sequence[Path], force: bool) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing and not force:
        raise SystemExit("Output exists; use --force to overwrite: {0}".format(", ".join(existing)))


def write_report(path: Path, summary_rows: Sequence[Dict[str, str]], args: argparse.Namespace) -> None:
    metrics = {row["metric"]: row["value"] for row in summary_rows}
    lines = [
        "# GeneArbiter Final Naming Report",
        "",
        "- proposal_gff: `{0}`".format(args.proposal_gff),
        "- proposal_records: `{0}`".format(args.proposal_records),
        "- output_gff: `{0}`".format(path.parent / args.output_gff_name),
        "- clean_auto_pass_gff: `{0}`".format(path.parent / args.clean_gff_name) if args.clean_gff_name else "- clean_auto_pass_gff: disabled",
        "- hisat_gtf: `{0}`".format(path.parent / args.gtf_name) if args.gtf_name else "- hisat_gtf: disabled",
        "- one_to_one_policy: `{0}`".format(args.one_to_one_policy),
        "- id_strategy: `{0}`".format(args.id_strategy),
        "- transcript_id_strategy: `{0}`".format(args.transcript_id_strategy),
        "- new_gene_prefix: `{0}`".format(args.new_gene_prefix),
        "",
        "## Summary",
        "",
    ]
    for key in [
        "proposal_genes",
        "one_to_one_replacement",
        "split_replacement",
        "merge_replacement",
        "complex_replacement",
        "novel_gene",
        "reused_current_gene_ids",
        "new_final_gene_ids",
        "new_position_ids",
        "new_source_style_ids",
        "new_fallback_namespace_ids",
        "proposal_transcripts",
        "renamed_proposal_child_records",
    ]:
        if key in metrics:
            lines.append("- {0}: {1}".format(key, metrics[key]))
    lines.extend(
        [
            "",
            "Notes:",
            "- This is an ID reconciliation layer, not a biological acceptance step.",
            "- Kept current genes retain their existing IDs unless the source current GFF reuses the same ID for multiple blocks; those duplicate current IDs are disambiguated with `.dupNN` and recorded in `current_id_disambiguation.tsv`.",
            "- One-to-one replacements retain the current gene lineage ID. Split, merge and complex loci use deterministic maximum-weight one-to-one matching: CDS overlap is primary, exact CDS blocks and reciprocal CDS similarity are secondary, and gene span/genomic order are fallbacks.",
            "- Under current_style_positional, extra/novel genes are inserted between neighboring current IDs only when that source namespace is demonstrably positional.",
            "- Under current_style_or_sequential, source ID prefix/width/suffix, embedded seqid templates and numeric step are inferred from the current GFF. No species name or species prefix is hard-coded.",
            "- Externally assigned NCBI LOC namespaces are never minted; only those unsafe/uninferable cases use the configured fallback namespace.",
            "- `locus_id_mapping.tsv` records current, corrected, retained, retired and added IDs per changed locus. Accepted structures always receive a deterministic naming result; mapping statuses describe the evidence used and are not review states.",
            "- Retired current IDs remain reserved and are never recycled for unrelated genes.",
            "- The clean auto-pass GFF removes proposal/review/risk display attributes from the release-like copy; provenance remains in mapping tables and the unclean GFF.",
            "- The HISAT GTF is a derived RNA-seq helper file: it excludes organelle seqids by default and emits gene/transcript/exon/CDS rows with gene_id/transcript_id attributes.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.fail_on_truth_path:
        fail_if_truth_like([args.proposal_gff, args.proposal_records, args.current_gff, args.out_dir])

    out_dir = Path(args.out_dir)
    ensure_dir(str(out_dir))
    out_gff = out_dir / args.output_gff_name
    clean_gff = out_dir / args.clean_gff_name if args.clean_gff_name else None
    hisat_gtf = out_dir / args.gtf_name if args.gtf_name else None
    gene_map_path = out_dir / "gene_id_mapping.tsv"
    tx_map_path = out_dir / "transcript_id_mapping.tsv"
    current_disambig_path = out_dir / "current_id_disambiguation.tsv"
    locus_map_path = out_dir / "locus_id_mapping.tsv"
    summary_path = out_dir / "final_naming_summary.tsv"
    report_path = out_dir / "final_naming_report.md"
    guarded_outputs = [out_gff, gene_map_path, tx_map_path, current_disambig_path, locus_map_path, summary_path, report_path]
    if clean_gff is not None:
        guarded_outputs.append(clean_gff)
    if hisat_gtf is not None:
        guarded_outputs.append(hisat_gtf)
    guard_outputs(guarded_outputs, args.force)

    records = read_records(args.proposal_records)
    existing_gene_ids, existing_tx_ids = collect_existing_ids(args.proposal_gff)
    current_genes = collect_current_genes(args.current_gff) if args.current_gff else {}
    proposal_cds = collect_proposal_cds(args.proposal_gff, records)
    gene_mapping, gene_rows, next_index = build_gene_mapping(
        records,
        existing_gene_ids,
        args.one_to_one_policy,
        args.new_gene_prefix,
        args.new_gene_start,
        args.new_gene_width,
        args.id_strategy,
        current_genes,
        proposal_cds,
    )
    tx_mapping, tx_rows = build_tx_mapping(
        records,
        gene_mapping,
        existing_tx_ids,
        args.transcript_template,
        args.transcript_id_strategy,
    )
    locus_rows = build_locus_mapping(gene_rows)
    rewrite_counts, current_disambig_rows = rewrite_gff(
        args.proposal_gff,
        str(out_gff),
        gene_mapping,
        tx_mapping,
        args.proposal_source_label,
    )
    clean_counts: Counter = Counter()
    if clean_gff is not None:
        clean_counts = write_clean_auto_pass_gff(out_gff, clean_gff, args.clean_source_label)
    gtf_counts: Counter = Counter()
    if hisat_gtf is not None:
        gtf_counts = write_hisat_gtf(clean_gff or out_gff, hisat_gtf, args.gtf_exclude_organelle)

    relation_counts = Counter(row["final_relation"] for row in gene_rows)
    summary_rows = [
        {"metric": "generated_at_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"metric": "proposal_genes", "value": str(len(gene_rows))},
        {"metric": "proposal_transcripts", "value": str(len(tx_rows))},
        {"metric": "changed_loci", "value": str(len(locus_rows))},
        {"metric": "id_strategy", "value": args.id_strategy},
        {"metric": "transcript_id_strategy", "value": args.transcript_id_strategy},
        {"metric": "deterministic_cds_anchors", "value": str(sum(1 for row in gene_rows if row["id_assignment_status"] == "deterministic_cds_anchor"))},
        {"metric": "deterministic_fallback_anchors", "value": str(sum(1 for row in gene_rows if row["id_assignment_status"] in {"deterministic_span_fallback", "deterministic_order_fallback"}))},
        {"metric": "kept_current_gene_ids_in_backbone", "value": str(len(existing_gene_ids))},
        {"metric": "one_to_one_replacement", "value": str(relation_counts.get("one_to_one_replacement", 0))},
        {"metric": "split_replacement", "value": str(relation_counts.get("split_replacement", 0))},
        {"metric": "merge_replacement", "value": str(relation_counts.get("merge_replacement", 0))},
        {"metric": "complex_replacement", "value": str(relation_counts.get("complex_replacement", 0))},
        {"metric": "novel_gene", "value": str(relation_counts.get("novel_gene", 0))},
        {"metric": "reused_current_gene_ids", "value": str(sum(1 for row in gene_rows if row["reuse_current_gene_id"] == "true"))},
        {"metric": "new_final_gene_ids", "value": str(sum(1 for row in gene_rows if row["reuse_current_gene_id"] != "true"))},
        {"metric": "new_position_ids", "value": str(sum(1 for row in gene_rows if row["id_assignment_status"] == "deterministic_new_position_id"))},
        {"metric": "new_source_style_ids", "value": str(sum(1 for row in gene_rows if row["id_assignment_status"] == "deterministic_new_source_style_id"))},
        {"metric": "new_fallback_namespace_ids", "value": str(sum(1 for row in gene_rows if row["id_assignment_status"] == "deterministic_new_namespace_id"))},
        {"metric": "next_new_gene_index", "value": str(next_index)},
    ]
    for key, value in sorted(rewrite_counts.items()):
        summary_rows.append({"metric": key, "value": str(value)})
    if clean_gff is not None:
        summary_rows.append({"metric": "clean_gff_path", "value": str(clean_gff)})
        for key, value in sorted(clean_counts.items()):
            summary_rows.append({"metric": "clean_" + key, "value": str(value)})
    if hisat_gtf is not None:
        summary_rows.append({"metric": "hisat_gtf_path", "value": str(hisat_gtf)})
        for key, value in sorted(gtf_counts.items()):
            summary_rows.append({"metric": "hisat_gtf_" + key, "value": str(value)})

    summary_rows.append({"metric": "current_disambiguated_ids", "value": str(len(current_disambig_rows))})
    write_tsv(str(gene_map_path), GENE_MAP_FIELDS, gene_rows)
    write_tsv(str(tx_map_path), TX_MAP_FIELDS, tx_rows)
    write_tsv(str(locus_map_path), LOCUS_MAP_FIELDS, locus_rows)
    write_tsv(str(current_disambig_path), CURRENT_DISAMBIG_FIELDS, current_disambig_rows)
    write_tsv(str(summary_path), SUMMARY_FIELDS, summary_rows)
    write_report(report_path, summary_rows, args)

    print("wrote\t{0}".format(out_gff))
    if clean_gff is not None:
        print("clean_gff\t{0}".format(clean_gff))
    if hisat_gtf is not None:
        print("hisat_gtf\t{0}".format(hisat_gtf))
    print("gene_mappings\t{0}".format(len(gene_rows)))
    print("transcript_mappings\t{0}".format(len(tx_rows)))


if __name__ == "__main__":
    main()
