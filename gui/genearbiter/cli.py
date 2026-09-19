#!/usr/bin/env python3
# script_id_md5: c94f70a33071ffe38d3049a55c046fe8
# created: 2026-07-07
# modified: 2026-07-10
# owner: project
# status: project_code
# purpose: 提供 GeneArbiter 配置驱动的科研命令行入口。
# inputs: JSON/YAML workflow config; ordered GFF list with current/base first; optional GFF normalization; optional RNA/protein evidence; model API key env。
# outputs: full cards, AI cards, AI decisions, final calls, optional review bundle, proposal/final-naming GFF, ID replacement summary。
# notes: v0.4 支持有序 GFF 配置和 run-files 简化入口；scoring.normalization=run_level 按实际 evidence 类型归一化。

"""Run the GeneArbiter annotation-set workflow from a config file.

The CLI keeps one implementation of the biological logic while making the
pipeline straightforward to install, run and audit.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import runpy
import shutil
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence

try:
    from genearbiter import __version__
    from genearbiter.short_read import add_extract_short_read_parser
    from genearbiter.evidence_validation import add_validate_evidence_parser
except ModuleNotFoundError:  # direct `python genearbiter/cli.py` execution from the source tree
    from __init__ import __version__
    from short_read import add_extract_short_read_parser
    from evidence_validation import add_validate_evidence_parser


STEP_ORDER = [
    "normalize_gffs",
    "validate_candidates",
    "audit_eligibility",
    "cards",
    "ai_cards",
    "api",
    "final_calls",
    "review_bundle",
    "export_gff",
    "reconcile_ids",
    "id_summary",
    "policy_export",
    "policy_reconcile_ids",
    "release_gff_qc",
]

SUPPORTED_EVIDENCE_TYPES = {
    "splice_junctions",
    "model_junction_support",
    "long_read_support_dir",
    "protein_gff",
    "validation_summary",
    "eligibility_summary",
    "gene_list",
    "region_table",
}
CORE_SCRIPT_NAMES = {
    "build_cards": "build_arbitration_cards.py",
    "build_ai_cards": "build_decision_cards.py",
    "run_ai": "run_decisions.py",
    "final_calls": "build_final_calls.py",
    "export_gff": "export_proposal_gff.py",
    "reconcile_ids": "reconcile_ids.py",
    "policy_export": "export_policy_gff.py",
    "release_gff_qc": "release_gff.py",
}
CANDIDATE_SCRIPT_NAMES = {
    "validate_candidates": "validate_models.py",
    "audit_eligibility": "audit_eligibility.py",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def load_config(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        try:
            import yaml  # type: ignore
        except ImportError as exc:
            raise SystemExit(
                "YAML config requires PyYAML. Install with `pip install -e .[yaml]` "
                "or use the JSON example config."
            ) from exc
        data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise SystemExit("Config root must be a mapping/object.")
    return data


def is_empty_value(value: Any) -> bool:
    return value is None or value == ""


def as_list(value: Any) -> List[Any]:
    if is_empty_value(value):
        return []
    if isinstance(value, list):
        return value
    return [value]


def config_root(config: Mapping[str, Any], config_path: Path) -> Path:
    project = config.get("project") or {}
    root = project.get("root") if isinstance(project, Mapping) else ""
    if root:
        path = Path(str(root)).expanduser()
        if path.is_absolute():
            return path.resolve()
        return (config_path.parent / path).resolve()
    return config_path.parent.resolve()


def resolve_path(value: Any, root: Path) -> str:
    if is_empty_value(value):
        return ""
    path = Path(str(value)).expanduser()
    if path.is_absolute():
        return str(path)
    return str((root / path).resolve())


def resolve_optional_dir(value: Any, root: Path) -> str:
    return resolve_path(value, root) if value else ""


def script_root() -> Path:
    """Locate the workflow modules bundled with GeneArbiter."""
    override = os.environ.get("GENEARBITER_SCRIPT_ROOT", "").strip()
    if override:
        root = Path(override).expanduser().resolve()
        if not (root / "common").exists():
            raise SystemExit(f"GENEARBITER_SCRIPT_ROOT lacks common/: {root}")
        return root

    package_root = Path(__file__).resolve().parent
    if (package_root / "common").exists():
        return package_root
    raise SystemExit(
        "Cannot locate GeneArbiter script root. Set GENEARBITER_SCRIPT_ROOT to a directory "
        "containing common/."
    )


def common_script(name: str) -> Path:
    path = script_root() / "common" / CORE_SCRIPT_NAMES[name]
    if not path.exists():
        raise SystemExit(f"Missing GeneArbiter workflow module: {path}")
    return path


def candidate_script(name: str) -> Path:
    path = script_root() / "candidates" / CANDIDATE_SCRIPT_NAMES[name]
    if not path.exists():
        raise SystemExit(f"Missing GeneArbiter candidate-QC module: {path}")
    return path


def run_dir(config: Mapping[str, Any], root: Path) -> Path:
    outputs = config.get("outputs") or {}
    if not isinstance(outputs, Mapping) or not outputs.get("run_dir"):
        raise SystemExit("Config must set outputs.run_dir.")
    return Path(resolve_path(outputs["run_dir"], root))


def strategy_name(config: Mapping[str, Any]) -> str:
    outputs = config.get("outputs") or {}
    return str(outputs.get("strategy") or "genearbiter_default")


def workflow_mode(config: Mapping[str, Any]) -> str:
    mode = str(config.get("mode") or "").strip()
    params = config.get("params") or {}
    if not mode and isinstance(params, Mapping):
        task_mode = str(params.get("task_mode") or "correction")
        mode = "arbitration" if task_mode == "de_novo_annotation" else "correction"
    if not mode:
        mode = "correction"
    aliases = {"de_novo_annotation": "arbitration", "multi_gff_arbitration": "arbitration", "refinement": "correction"}
    mode = aliases.get(mode, mode)
    if mode not in {"correction", "arbitration"}:
        raise SystemExit("Config mode must be correction or arbitration.")
    return mode


def card_subdir(config: Mapping[str, Any]) -> str:
    return "correction" if workflow_mode(config) == "correction" else "arbitration"


def input_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    inputs = config.get("inputs") or {}
    if not isinstance(inputs, Mapping):
        raise SystemExit("Config inputs must be a mapping/object.")
    return inputs


def params_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    params = config.get("params") or {}
    if not isinstance(params, Mapping):
        raise SystemExit("Config params must be a mapping/object.")
    return params


def api_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    api = config.get("api") or {}
    if not isinstance(api, Mapping):
        raise SystemExit("Config api must be a mapping/object.")
    return api


def export_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    export = config.get("export") or {}
    if not isinstance(export, Mapping):
        raise SystemExit("Config export must be a mapping/object.")
    return export


def policy_export_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    section = config.get("policy_export")
    if section is None:
        export = export_section(config)
        section = export.get("policy_export") or {}
    if not isinstance(section, Mapping):
        raise SystemExit("Config policy_export must be a mapping/object.")
    return section


def policy_export_enabled(config: Mapping[str, Any]) -> bool:
    section = policy_export_section(config)
    return bool(section.get("enabled", False))


def policy_name(config: Mapping[str, Any]) -> str:
    section = policy_export_section(config)
    return str(section.get("policy") or "supported_novel_singletool_junction_supported")


def policy_reconcile_enabled(config: Mapping[str, Any]) -> bool:
    section = policy_export_section(config)
    return bool(section.get("reconcile_ids", True))


def release_gff_qc_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    section = config.get("release_gff_qc")
    if section is None:
        policy = policy_export_section(config)
        section = policy.get("release_gff_qc") or {}
    if not isinstance(section, Mapping):
        raise SystemExit("Config release_gff_qc must be a mapping/object.")
    return section


def release_gff_qc_enabled(config: Mapping[str, Any]) -> bool:
    section = release_gff_qc_section(config)
    return bool(section.get("enabled", True))


def review_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    review = config.get("review") or {}
    if not isinstance(review, Mapping):
        raise SystemExit("Config review must be a mapping/object.")
    return review


def preprocess_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    preprocess = config.get("preprocess") or {}
    if not isinstance(preprocess, Mapping):
        raise SystemExit("Config preprocess must be a mapping/object.")
    return preprocess


def scoring_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    scoring = config.get("scoring") or {}
    if not isinstance(scoring, Mapping):
        raise SystemExit("Config scoring must be a mapping/object.")
    return scoring


def scoring_normalization(config: Mapping[str, Any]) -> str:
    scoring = scoring_section(config)
    mode = str(scoring.get("normalization") or "run_level").strip().lower().replace("-", "_")
    aliases = {"none": "legacy", "raw": "legacy", "runlevel": "run_level", "run_level_available_evidence": "run_level"}
    mode = aliases.get(mode, mode)
    if mode not in {"legacy", "run_level"}:
        raise SystemExit("Config scoring.normalization must be legacy or run_level.")
    return mode


def candidate_qc_section(config: Mapping[str, Any]) -> Mapping[str, Any]:
    candidate_qc = config.get("candidate_qc") or {}
    if not isinstance(candidate_qc, Mapping):
        raise SystemExit("Config candidate_qc must be a mapping/object.")
    return candidate_qc


def genome_fasta_path(config: Mapping[str, Any], root: Path) -> str:
    value = config.get("genome_fasta")
    if not value:
        inputs = input_section(config)
        value = inputs.get("genome_fasta") or inputs.get("genome")
    return resolve_path(value, root) if value else ""


def candidate_qc_enabled(config: Mapping[str, Any], root: Path) -> bool:
    section = candidate_qc_section(config)
    if "enabled" in section:
        return bool(section.get("enabled"))
    return bool(genome_fasta_path(config, root))


def normalize_gffs_enabled(config: Mapping[str, Any]) -> bool:
    preprocess = preprocess_section(config)
    return bool(preprocess.get("normalize_gffs", True))


def safe_label(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in text) or "source"


def review_mode(config: Mapping[str, Any], override: str = "") -> str:
    value = str(override or review_section(config).get("mode") or "auto").strip().lower()
    aliases = {"none": "auto", "off": "auto", "false": "auto", "manual": "review", "human": "review"}
    value = aliases.get(value, value)
    if value not in {"auto", "review"}:
        raise SystemExit("review.mode must be auto or review.")
    return value


def gff_entries(config: Mapping[str, Any], root: Path) -> List[Dict[str, str]]:
    """Return normalized GFF entries. Top-level gffs are ordered: first=current, rest=candidate."""
    raw = config.get("gffs")
    entries: List[Dict[str, Any]] = []
    ordered_top_level = raw is not None
    if raw is None:
        inputs = input_section(config)
        if inputs.get("current_gff"):
            entries.append({"name": "current", "role": "current", "path": inputs.get("current_gff")})
        candidates = inputs.get("candidate_gffs") or []
        if isinstance(candidates, Mapping):
            candidates = [{"name": name, "path": path} for name, path in candidates.items()]
        for row in as_list(candidates):
            if not isinstance(row, Mapping):
                raise SystemExit("Each inputs.candidate_gffs entry must be a mapping with name/path.")
            item = dict(row)
            item.setdefault("role", "candidate")
            entries.append(item)
    else:
        if isinstance(raw, Mapping):
            raw = [{"name": name, "path": value} for name, value in raw.items()]
        for row in as_list(raw):
            if not isinstance(row, Mapping):
                raise SystemExit("Each gffs entry must be a mapping with name/path.")
            entries.append(dict(row))

    normalized: List[Dict[str, str]] = []
    for idx, row in enumerate(entries, start=1):
        name = str(row.get("name") or row.get("label") or row.get("source") or f"GFF{idx}").strip()
        path_value = row.get("path") or row.get("gff") or row.get("file")
        if not name or not path_value:
            raise SystemExit("Each GFF entry requires at least name and path.")
        if ordered_top_level and workflow_mode(config) == "correction":
            role = "current" if idx == 1 else "candidate"
            if idx == 1:
                name = "current"
        else:
            role = str(row.get("role") or row.get("kind") or "candidate").strip().lower()
            if role in {"base", "backbone", "reference", "current_gff"}:
                role = "current"
            elif role in {"tool", "source", "prediction", "candidate_gff"}:
                role = "candidate"
            if role not in {"current", "candidate"}:
                raise SystemExit(f"Unsupported GFF role for {name}: {role}")
        normalized.append({"name": name, "role": role, "path": resolve_path(path_value, root)})
    return normalized


def current_gff_path(config: Mapping[str, Any], root: Path) -> str:
    currents = [entry for entry in gff_entries(config, root) if entry["role"] == "current"]
    if not currents:
        return ""
    if len(currents) > 1:
        names = ",".join(entry["name"] for entry in currents)
        raise SystemExit(f"Correction mode accepts exactly one current/base GFF; found: {names}")
    return currents[0]["path"]


def normalized_candidate_path(paths: Mapping[str, Path], index: int, name: str) -> str:
    return str(paths["normalized_gff_dir"] / f"{index:02d}_{safe_label(name)}.normalized.gff3")


def candidate_gff_args(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path] | None = None) -> List[str]:
    args: List[str] = []
    candidate_index = 0
    for entry in gff_entries(config, root):
        if entry["role"] == "candidate":
            candidate_index += 1
            path = entry["path"]
            if paths is not None and normalize_gffs_enabled(config):
                path = normalized_candidate_path(paths, candidate_index, entry["name"])
            args.extend(["--tool-gff", f"{entry['name']}={path}"])
    return args


def source_priority(config: Mapping[str, Any], root: Path) -> List[str]:
    priority: List[str] = []
    for entry in gff_entries(config, root):
        name = "current" if entry["role"] == "current" else entry["name"]
        if name not in priority:
            priority.append(name)
    if "current" not in priority:
        priority.insert(0, "current")
    return priority


def infer_evidence_type(path: str) -> str:
    p = Path(path)
    name = p.name.lower()
    suffixes = [item.lower() for item in p.suffixes]
    if p.is_dir():
        long_read_markers = {
            "long_read_model_support_by_model.tsv",
            "long_read_set_support_by_locus.tsv",
            "long_read_merged_junctions.tsv",
        }
        child_names = {child.name for child in p.glob("**/*.tsv")} if p.exists() else set()
        if long_read_markers & child_names or "long_read" in name:
            return "long_read_support_dir"
        return "unknown"
    if suffixes and any(suffix in {".gff", ".gff3", ".gtf"} for suffix in suffixes):
        return "protein_gff"
    if not p.exists():
        return "unknown"
    try:
        opener = open
        if suffixes and suffixes[-1] == ".gz":
            import gzip

            opener = gzip.open  # type: ignore[assignment]
        with opener(p, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
            header = ""
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    header = line
                    break
    except OSError:
        return "unknown"
    columns = {col.strip() for col in header.replace(",", "\t").split("\t") if col.strip()}
    if {"chrom", "intron_start_1based", "intron_end_1based"} <= columns:
        return "splice_junctions"
    if {"seqid", "intron_start", "intron_end"} <= columns:
        return "splice_junctions"
    if {"transcript_id", "junction_support_fraction"} <= columns:
        return "model_junction_support"
    if {"model_set", "model_id", "support_class"} <= columns:
        return "long_read_model_support"
    if {"model_set", "locus_id", "support_level"} <= columns:
        return "long_read_locus_support"
    return "unknown"


def evidence_entries(config: Mapping[str, Any], root: Path) -> List[Dict[str, str]]:
    """Return normalized evidence records from top-level evidence or legacy inputs."""
    raw = config.get("evidence")
    entries: List[Dict[str, Any]] = []
    if raw is None:
        inputs = input_section(config)
        for item in as_list(inputs.get("junctions")):
            entries.append({"name": "splice_junctions", "type": "splice_junctions", "path": item})
        legacy_single = [
            ("junction_support_summary", "model_junction_support"),
            ("long_read_support_dir", "long_read_support_dir"),
            ("validation_summary", "validation_summary"),
            ("eligibility_summary", "eligibility_summary"),
            ("gene_list", "gene_list"),
            ("region_table", "region_table"),
        ]
        for key, kind in legacy_single:
            if inputs.get(key):
                entries.append({"name": key, "type": kind, "path": inputs.get(key)})
        for item in as_list(inputs.get("protein_gffs")):
            if isinstance(item, Mapping):
                entry = dict(item)
                entry.setdefault("type", "protein_gff")
                entries.append(entry)
            else:
                entries.append({"name": "protein_gff", "type": "protein_gff", "path": item})
    else:
        if isinstance(raw, Mapping):
            raw = [{"name": name, **(value if isinstance(value, Mapping) else {"path": value})} for name, value in raw.items()]
        for row in as_list(raw):
            if isinstance(row, Mapping):
                entries.append(dict(row))
            else:
                entries.append({"path": row})

    normalized: List[Dict[str, str]] = []
    for idx, row in enumerate(entries, start=1):
        path_value = row.get("path") or row.get("file") or row.get("dir")
        if not path_value:
            raise SystemExit("Each evidence entry requires path.")
        path = resolve_path(path_value, root)
        kind = str(row.get("type") or row.get("kind") or "auto").strip().lower()
        if kind in {"", "auto"}:
            kind = infer_evidence_type(path)
        aliases = {
            "junctions": "splice_junctions",
            "splice_junction": "splice_junctions",
            "rna_junctions": "splice_junctions",
            "junction_support": "model_junction_support",
            "junction_support_summary": "model_junction_support",
            "long_read": "long_read_support_dir",
            "long_read_dir": "long_read_support_dir",
            "protein": "protein_gff",
            "protein_alignment_gff": "protein_gff",
            "miniprot_gff": "protein_gff",
        }
        kind = aliases.get(kind, kind)
        if kind == "unknown":
            raise SystemExit(
                f"Cannot infer evidence type from {path}. Set an explicit supported type; "
                "unrecognized evidence is not silently ignored."
            )
        if kind in {"long_read_model_support", "long_read_locus_support"}:
            raise SystemExit(
                f"Standalone {kind} is not a complete long-read input. "
                "Provide the containing evidence directory as type: long_read_support_dir."
            )
        if kind not in SUPPORTED_EVIDENCE_TYPES:
            allowed = ",".join(sorted(SUPPORTED_EVIDENCE_TYPES))
            raise SystemExit(f"Unsupported evidence type {kind!r} for {path}. Supported types: {allowed}")
        name = str(row.get("name") or row.get("label") or Path(path).stem or f"evidence{idx}").strip()
        normalized.append({"name": name, "type": kind, "path": path})
    return normalized


def evidence_args(config: Mapping[str, Any], root: Path) -> List[str]:
    args: List[str] = []
    singles: Dict[str, tuple[str, str]] = {
        "model_junction_support": ("--junction-support-summary", ""),
        "long_read_support_dir": ("--long-read-support-dir", ""),
        "validation_summary": ("--validation-summary", ""),
        "eligibility_summary": ("--eligibility-summary", ""),
        "gene_list": ("--gene-list", ""),
        "region_table": ("--region-table", ""),
    }
    seen_single: set[str] = set()
    for entry in evidence_entries(config, root):
        kind = entry["type"]
        path = entry["path"]
        if kind == "splice_junctions":
            args.extend(["--junctions", path])
        elif kind == "protein_gff":
            args.extend(["--protein-gff", f"{entry['name']}={path}"])
        elif kind in singles and kind not in seen_single:
            flag, _ = singles[kind]
            args.extend([flag, path])
            seen_single.add(kind)
        else:
            raise SystemExit(f"Evidence type {kind!r} is validated but not connected to the workflow.")
    return args


def active_run_level_evidence_types(config: Mapping[str, Any], root: Path) -> List[str]:
    kinds = {entry["type"] for entry in evidence_entries(config, root)}
    active: List[str] = []
    if kinds & {"splice_junctions", "model_junction_support"}:
        active.append("short_read")
    if kinds & {"long_read_support_dir", "long_read_model_support", "long_read_locus_support"}:
        active.append("long_read")
    if "protein_gff" in kinds:
        active.append("protein")
    return active


def output_paths(config: Mapping[str, Any], root: Path) -> Dict[str, Path]:
    run = run_dir(config, root)
    strategy = strategy_name(config)
    profile = str(api_section(config).get("profile") or "deepseek_flash")
    api_label = str(api_section(config).get("run_label") or profile)
    mode_dir = card_subdir(config)
    policy = policy_name(config)
    policy_dir = run / "05_proposal_catalog" / strategy / "policy_exports" / policy
    policy_final_naming_dir = policy_dir / "final_naming_v1"
    release_source_dir = (
        policy_final_naming_dir
        if policy_export_enabled(config) and policy_reconcile_enabled(config)
        else run / "05_proposal_catalog" / strategy / "final_naming_v1"
    )
    release_dir = release_source_dir / "release_gff_qc"
    release_section = release_gff_qc_section(config)
    release_gff_name = str(release_section.get("output_gff_name") or "annotation.trace.gff3")
    release_clean_name = str(release_section.get("clean_gff_name") or "annotation.clean.gff3")
    release_summary_name = str(release_section.get("summary_name") or "release_gff_qc_summary.tsv")
    release_mapping_name = str(release_section.get("mapping_name") or "release_gff_id_mapping.tsv")
    release_report_name = str(release_section.get("report_name") or "release_gff_report.md")
    return {
        "run_dir": run,
        "logs": run / "logs",
        "normalized_gff_dir": run / "00_normalized_gffs",
        "normalized_gff_manifest": run / "00_normalized_gffs" / "normalized_gff_manifest.tsv",
        "normalized_gff_report": run / "00_normalized_gffs" / "normalized_gff_report.tsv",
        "candidate_qc_dir": run / "01_candidate_qc" / mode_dir,
        "validation_summary": run / "01_candidate_qc" / mode_dir / "candidate_model_validation.tsv",
        "eligibility_summary": run / "01_candidate_qc" / mode_dir / "candidate_model_eligibility.tsv",
        "eligibility_reason_summary": run / "01_candidate_qc" / mode_dir / "candidate_model_eligibility_reason_summary.tsv",
        "models": run / "01_windows" / mode_dir / "model_arbitration_models.tsv",
        "complex_loci": run / "01_windows" / mode_dir / "model_arbitration_complex_loci.tsv",
        "summary": run / "01_windows" / mode_dir / "model_arbitration_summary.tsv",
        "units": run / "01_windows" / mode_dir / "model_arbitration_units.tsv",
        "unit_edges": run / "01_windows" / mode_dir / "model_arbitration_unit_edges.tsv",
        "model_sets": run / "01_windows" / mode_dir / "model_arbitration_sets.tsv",
        "full_cards": run / "02_cards" / mode_dir / "model_arbitration_full_cards.jsonl",
        "ai_cards": run / "02_cards" / strategy / "genearbiter_decision_cards.jsonl",
        "ai_summary": run / "02_cards" / strategy / "genearbiter_decision_summary.tsv",
        "set_trace": run / "02_cards" / strategy / "genearbiter_set_trace.jsonl",
        "auto_decisions": run / "02_cards" / strategy / "genearbiter_auto_decisions.jsonl",
        "review_router": run / "02_cards" / strategy / "locus_review_router.tsv",
        "api_run_dir": run / "03_decisions" / strategy / api_label,
        "api_decisions": run / "03_decisions" / strategy / api_label / profile / "decisions.jsonl",
        "final_calls_dir": run / "04_final_calls" / strategy,
        "final_calls": run / "04_final_calls" / strategy / "final_annotation_calls.tsv",
        "review_bundle_dir": run / "06_review_bundle" / strategy,
        "review_queue": run / "06_review_bundle" / strategy / "review_queue.tsv",
        "manual_review_template": run / "06_review_bundle" / strategy / "manual_review_decisions.template.tsv",
        "review_cards_dir": run / "06_review_bundle" / strategy / "review_cards",
        "proposal_dir": run / "05_proposal_catalog" / strategy,
        "proposal_gff": run / "05_proposal_catalog" / strategy / "proposal_catalog.gff3",
        "proposal_records": run / "05_proposal_catalog" / strategy / "proposal_catalog_records.tsv",
        "final_naming_dir": run / "05_proposal_catalog" / strategy / "final_naming_v1",
        "gene_id_mapping": run / "05_proposal_catalog" / strategy / "final_naming_v1" / "gene_id_mapping.tsv",
        "transcript_id_mapping": run / "05_proposal_catalog" / strategy / "final_naming_v1" / "transcript_id_mapping.tsv",
        "clean_final_named_gff": run / "05_proposal_catalog" / strategy / "final_naming_v1" / "final_named_annotation.auto_pass.clean.gff3",
        "hisat_gtf": run / "05_proposal_catalog" / strategy / "final_naming_v1" / "final_named_annotation.auto_pass.hisat.gtf",
        "policy_dir": policy_dir,
        "policy_gff": policy_dir / "policy_catalog.gff3",
        "policy_records": policy_dir / "policy_catalog_records.tsv",
        "policy_summary": policy_dir / "policy_catalog_summary.tsv",
        "policy_final_naming_dir": policy_final_naming_dir,
        "policy_gene_id_mapping": policy_final_naming_dir / "gene_id_mapping.tsv",
        "policy_transcript_id_mapping": policy_final_naming_dir / "transcript_id_mapping.tsv",
        "policy_named_gff": policy_final_naming_dir / "policy_named_annotation.gff3",
        "policy_clean_final_named_gff": policy_final_naming_dir / "policy_named_annotation.auto_pass.clean.gff3",
        "policy_hisat_gtf": policy_final_naming_dir / "policy_named_annotation.auto_pass.hisat.gtf",
        "release_gff_qc_dir": release_dir,
        "release_source_gff": (
            policy_final_naming_dir / "policy_named_annotation.gff3"
            if policy_export_enabled(config) and policy_reconcile_enabled(config)
            else run / "05_proposal_catalog" / strategy / "final_naming_v1" / "final_named_annotation.gff3"
        ),
        "release_gff": release_dir / release_gff_name,
        "release_clean_gff": release_dir / release_clean_name,
        "release_gff_summary": release_dir / release_summary_name,
        "release_gff_mapping": release_dir / release_mapping_name,
        "release_gff_report": release_dir / release_report_name,
        "id_summary": (
            run / "05_proposal_catalog" / strategy / "final_naming_v1" / "genearbiter_id_mapping_summary.tsv"
            if workflow_mode(config) == "correction"
            else run / "04_final_calls" / strategy / "genearbiter_id_mapping_summary.tsv"
        ),
    }


def normalize_gffs_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    cmd = [
        sys.executable,
        str(Path(__file__).with_name("gff_normalize.py")),
        "--out-dir",
        str(paths["normalized_gff_dir"]),
        "--manifest",
        str(paths["normalized_gff_manifest"]),
        "--report",
        str(paths["normalized_gff_report"]),
    ]
    for entry in gff_entries(config, root):
        if entry["role"] == "candidate":
            cmd.extend(["--gff", f"{entry['name']}={entry['path']}"])
    return cmd


def qc_common_args(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    genome_fasta = genome_fasta_path(config, root)
    if not genome_fasta:
        raise SystemExit("candidate_qc requires genome_fasta at config top level or inputs.genome_fasta.")
    args = [
        "--genome-fasta",
        genome_fasta,
        "--current-gff",
        current_gff_path(config, root),
    ]
    args.extend(candidate_gff_args(config, root, paths))
    params = params_section(config)
    for cfg_key, flag in [("gene_list", "--gene-list"), ("region_table", "--region-table"), ("region_flank", "--region-flank")]:
        if cfg_key in params:
            args.extend([flag, str(params[cfg_key])])
    return args


def validate_candidates_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    return [
        sys.executable,
        str(candidate_script("validate_candidates")),
        *qc_common_args(config, root, paths),
        "--output-validation",
        str(paths["validation_summary"]),
    ]


def audit_eligibility_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    params = params_section(config)
    cmd = [
        sys.executable,
        str(candidate_script("audit_eligibility")),
        *qc_common_args(config, root, paths),
    ]
    for entry in evidence_entries(config, root):
        if entry["type"] == "splice_junctions":
            cmd.extend(["--junctions", entry["path"]])
    for cfg_key, flag in [
        ("min_junction_support_count", "--min-junction-support-count"),
        ("min_junction_sample_count", "--min-junction-sample-count"),
        ("same_source_conflict_fraction", "--same-source-conflict-fraction"),
    ]:
        if cfg_key in params:
            cmd.extend([flag, str(params[cfg_key])])
    cmd.extend([
        "--output-eligibility",
        str(paths["eligibility_summary"]),
        "--output-summary",
        str(paths["eligibility_reason_summary"]),
    ])
    return cmd


def build_cards_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    params = params_section(config)
    mode = workflow_mode(config)
    task_mode = "correction" if mode == "correction" else "de_novo_annotation"
    cmd = [
        sys.executable,
        str(common_script("build_cards")),
        "--task-mode",
        task_mode,
        "--candidate-set-mode",
        str(params.get("candidate_set_mode") or "structure_consensus"),
        "--source-priority",
        ",".join(source_priority(config, root)),
    ]
    if mode == "correction":
        current_gff = current_gff_path(config, root)
        if not current_gff:
            raise SystemExit("Correction mode requires one GFF with role: current.")
        cmd.extend(["--current-gff", current_gff])
    cmd.extend(candidate_gff_args(config, root, paths))
    cmd.extend(evidence_args(config, root))
    if candidate_qc_enabled(config, root):
        cmd.extend(["--validation-summary", str(paths["validation_summary"])])
        cmd.extend(["--eligibility-summary", str(paths["eligibility_summary"])])
    for cfg_key, flag in [
        ("signature_mode", "--signature-mode"),
        ("min_overlap_fraction", "--min-overlap-fraction"),
        ("min_junction_support_count", "--min-junction-support-count"),
        ("min_junction_sample_count", "--min-junction-sample-count"),
        ("boundary_window", "--boundary-window"),
        ("nearby_window", "--nearby-window"),
        ("protein_flank", "--protein-flank"),
        ("max_candidate_models", "--max-candidate-models"),
        ("max_unique_structures", "--max-unique-structures"),
        ("max_transcripts_per_source", "--max-transcripts-per-source"),
        ("max_locus_bp", "--max-locus-bp"),
        ("ai_soft_locus_bp", "--ai-soft-locus-bp"),
        ("ai_hard_locus_bp", "--ai-hard-locus-bp"),
        ("max_ai_current_genes_in_long_window", "--max-ai-current-genes-in-long-window"),
        ("max_ai_models_per_source_in_long_window", "--max-ai-models-per-source-in-long-window"),
        ("region_flank", "--region-flank"),
    ]:
        if cfg_key in params:
            cmd.extend([flag, str(params[cfg_key])])
    cmd.extend(
        [
            "--output-units",
            str(paths["units"]),
            "--output-unit-edges",
            str(paths["unit_edges"]),
            "--output-model-sets",
            str(paths["model_sets"]),
            "--output-models",
            str(paths["models"]),
            "--output-complex-loci",
            str(paths["complex_loci"]),
            "--output-summary",
            str(paths["summary"]),
            "--output-jsonl",
            str(paths["full_cards"]),
        ]
    )
    return cmd


def build_ai_cards_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    active_evidence = active_run_level_evidence_types(config, root)
    return [
        sys.executable,
        str(common_script("build_ai_cards")),
        "--input-jsonl",
        str(paths["full_cards"]),
        "--output-jsonl",
        str(paths["ai_cards"]),
        "--output-summary",
        str(paths["ai_summary"]),
        "--output-set-trace",
        str(paths["set_trace"]),
        "--output-auto-decisions",
        str(paths["auto_decisions"]),
        "--output-locus-review-router",
        str(paths["review_router"]),
        "--source-priority",
        ",".join(source_priority(config, root)),
        "--scoring-normalization",
        scoring_normalization(config),
        "--active-evidence-types",
        ",".join(active_evidence),
    ]


def run_api_command(config: Mapping[str, Any], paths: Mapping[str, Path], resume: bool) -> List[str]:
    api = api_section(config)
    cmd = [
        sys.executable,
        str(common_script("run_ai")),
        "--input-jsonl",
        str(paths["ai_cards"]),
        "--out-dir",
        str(paths["api_run_dir"]),
        "--profile",
        str(api.get("profile") or "deepseek_flash"),
        "--workers",
        str(api.get("workers", 1)),
        "--max-attempts",
        str(api.get("max_attempts", 2)),
        "--timeout",
        str(api.get("timeout", 180)),
        "--temperature",
        str(api.get("temperature", 0.0)),
        "--max-tokens",
        str(api.get("max_tokens", 8192)),
    ]
    if api.get("base_url"):
        cmd.extend(["--base-url", str(api["base_url"])])
    if api.get("dry_run_prompts"):
        cmd.append("--dry-run-prompts")
    if resume:
        cmd.append("--resume")
    return cmd


def final_calls_command(paths: Mapping[str, Path]) -> List[str]:
    return [
        sys.executable,
        str(common_script("final_calls")),
        "--set-trace",
        str(paths["set_trace"]),
        "--auto-decisions",
        str(paths["auto_decisions"]),
        "--api-decisions",
        str(paths["api_decisions"]),
        "--locus-review-router",
        str(paths["review_router"]),
        "--out-dir",
        str(paths["final_calls_dir"]),
    ]


def export_gff_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    export = export_section(config)
    cmd = [
        sys.executable,
        str(common_script("export_gff")),
        "--current-gff",
        current_gff_path(config, root),
        "--full-cards",
        str(paths["full_cards"]),
        "--final-calls",
        str(paths["final_calls"]),
        "--out-dir",
        str(paths["proposal_dir"]),
        "--source-label",
        str(export.get("source_label") or "GeneArbiter"),
        "--fail-on-truth-path",
    ]
    if export.get("sort_records"):
        cmd.append("--sort-records")
    return cmd


def reconcile_ids_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    export = export_section(config)
    cmd = [
        sys.executable,
        str(common_script("reconcile_ids")),
        "--proposal-gff",
        str(paths["proposal_gff"]),
        "--proposal-records",
        str(paths["proposal_records"]),
        "--out-dir",
        str(paths["final_naming_dir"]),
        "--new-gene-prefix",
        str(export.get("new_gene_prefix") or "AIFINAL_"),
        "--one-to-one-policy",
        str(export.get("one_to_one_policy") or "reuse_current_gene_id"),
        "--id-strategy",
        str(export.get("id_strategy") or "sequential"),
        "--transcript-id-strategy",
        str(export.get("transcript_id_strategy") or "template"),
        "--fail-on-truth-path",
    ]
    if export.get("id_strategy") in {"current_style_positional", "current_style_or_sequential"}:
        cmd.extend(["--current-gff", current_gff_path(config, root)])
    if export.get("transcript_template"):
        cmd.extend(["--transcript-template", str(export["transcript_template"])])
    if export.get("new_gene_start"):
        cmd.extend(["--new-gene-start", str(export["new_gene_start"])])
    if export.get("new_gene_width"):
        cmd.extend(["--new-gene-width", str(export["new_gene_width"])])
    if export.get("force_reconcile", True):
        cmd.append("--force")
    return cmd


def policy_export_command(config: Mapping[str, Any], paths: Mapping[str, Path]) -> List[str]:
    return [
        sys.executable,
        str(common_script("policy_export")),
        "--proposal-gff",
        str(paths["proposal_gff"]),
        "--proposal-records",
        str(paths["proposal_records"]),
        "--final-calls",
        str(paths["final_calls"]),
        "--models",
        str(paths["models"]),
        "--out-dir",
        str(paths["policy_dir"]),
        "--policy",
        policy_name(config),
        "--fail-on-truth-path",
    ]


def policy_reconcile_ids_command(config: Mapping[str, Any], root: Path, paths: Mapping[str, Path]) -> List[str]:
    export = export_section(config)
    policy = policy_export_section(config)
    cmd = [
        sys.executable,
        str(common_script("reconcile_ids")),
        "--proposal-gff",
        str(paths["policy_gff"]),
        "--proposal-records",
        str(paths["policy_records"]),
        "--out-dir",
        str(paths["policy_final_naming_dir"]),
        "--output-gff-name",
        str(policy.get("output_gff_name") or "policy_named_annotation.gff3"),
        "--clean-gff-name",
        str(policy.get("clean_gff_name") or "policy_named_annotation.auto_pass.clean.gff3"),
        "--gtf-name",
        str(policy.get("gtf_name") or "policy_named_annotation.auto_pass.hisat.gtf"),
        "--new-gene-prefix",
        str(policy.get("new_gene_prefix") or export.get("new_gene_prefix") or "AIFINAL_"),
        "--one-to-one-policy",
        str(policy.get("one_to_one_policy") or export.get("one_to_one_policy") or "reuse_current_gene_id"),
        "--id-strategy",
        str(policy.get("id_strategy") or export.get("id_strategy") or "sequential"),
        "--transcript-id-strategy",
        str(policy.get("transcript_id_strategy") or export.get("transcript_id_strategy") or "template"),
        "--fail-on-truth-path",
    ]
    if (policy.get("id_strategy") or export.get("id_strategy")) in {"current_style_positional", "current_style_or_sequential"}:
        cmd.extend(["--current-gff", current_gff_path(config, root)])
    transcript_template = policy.get("transcript_template") or export.get("transcript_template")
    proposal_source_label = policy.get("source_label") or export.get("source_label")
    if proposal_source_label:
        cmd.extend(["--proposal-source-label", str(proposal_source_label), "--clean-source-label", str(proposal_source_label)])
    if transcript_template:
        cmd.extend(["--transcript-template", str(transcript_template)])
    if policy.get("new_gene_start") or export.get("new_gene_start"):
        cmd.extend(["--new-gene-start", str(policy.get("new_gene_start") or export.get("new_gene_start"))])
    if policy.get("new_gene_width") or export.get("new_gene_width"):
        cmd.extend(["--new-gene-width", str(policy.get("new_gene_width") or export.get("new_gene_width"))])
    if policy.get("force_reconcile", export.get("force_reconcile", True)):
        cmd.append("--force")
    return cmd


def release_gff_qc_command(config: Mapping[str, Any], paths: Mapping[str, Path]) -> List[str]:
    section = release_gff_qc_section(config)
    cmd = [
        sys.executable,
        str(common_script("release_gff_qc")),
        "--input-gff",
        str(paths["release_source_gff"]),
        "--out-dir",
        str(paths["release_gff_qc_dir"]),
        "--output-gff-name",
        str(section.get("output_gff_name") or "annotation.trace.gff3"),
        "--clean-gff-name",
        str(section.get("clean_gff_name") or "annotation.clean.gff3"),
        "--summary-name",
        str(section.get("summary_name") or "release_gff_qc_summary.tsv"),
        "--mapping-name",
        str(section.get("mapping_name") or "release_gff_id_mapping.tsv"),
        "--report-name",
        str(section.get("report_name") or "release_gff_report.md"),
        "--transcript-template",
        str(section.get("transcript_template") or "{gene_id}.t{index:02d}"),
        "--exon-template",
        str(section.get("exon_template") or "{transcript_id}.exon{index:03d}"),
        "--fail-on-truth-path",
    ]
    if section.get("source_label"):
        cmd.extend(["--source-label", str(section["source_label"])])
    if section.get("force", True):
        cmd.append("--force")
    return cmd


def id_summary_command(config: Mapping[str, Any], paths: Mapping[str, Path]) -> List[str]:
    return [
        sys.executable,
        str(Path(__file__).with_name("id_summary.py")),
        "--workflow-mode",
        workflow_mode(config),
        "--full-cards",
        str(paths["full_cards"]),
        "--final-calls",
        str(paths["final_calls"]),
        "--proposal-records",
        str(paths["proposal_records"]),
        "--gene-id-mapping",
        str(paths["gene_id_mapping"]),
        "--transcript-id-mapping",
        str(paths["transcript_id_mapping"]),
        "--out-tsv",
        str(paths["id_summary"]),
    ]


def review_bundle_command(paths: Mapping[str, Path]) -> List[str]:
    return [
        sys.executable,
        str(Path(__file__).with_name("review_bundle.py")),
        "--final-calls",
        str(paths["final_calls"]),
        "--full-cards",
        str(paths["full_cards"]),
        "--out-dir",
        str(paths["review_bundle_dir"]),
    ]


def command_plan(config: Mapping[str, Any], root: Path, selected_steps: Sequence[str], resume_api: bool) -> List[Dict[str, Any]]:
    paths = output_paths(config, root)
    builders = {
        "normalize_gffs": lambda: normalize_gffs_command(config, root, paths),
        "validate_candidates": lambda: validate_candidates_command(config, root, paths),
        "audit_eligibility": lambda: audit_eligibility_command(config, root, paths),
        "cards": lambda: build_cards_command(config, root, paths),
        "ai_cards": lambda: build_ai_cards_command(config, root, paths),
        "api": lambda: run_api_command(config, paths, resume_api),
        "final_calls": lambda: final_calls_command(paths),
        "review_bundle": lambda: review_bundle_command(paths),
        "export_gff": lambda: export_gff_command(config, root, paths),
        "reconcile_ids": lambda: reconcile_ids_command(config, root, paths),
        "id_summary": lambda: id_summary_command(config, paths),
        "policy_export": lambda: policy_export_command(config, paths),
        "policy_reconcile_ids": lambda: policy_reconcile_ids_command(config, root, paths),
        "release_gff_qc": lambda: release_gff_qc_command(config, paths),
    }
    return [{"step": step, "command": builders[step]()} for step in selected_steps]


def expand_steps(
    raw_steps: str,
    skip_api: bool,
    reconcile: bool,
    mode: str,
    review: str,
    normalize_gffs: bool,
    candidate_qc: bool,
    policy_export: bool,
    policy_reconcile: bool,
    release_gff_qc: bool,
) -> List[str]:
    explicit = raw_steps != "all"
    if raw_steps == "all":
        if review == "review":
            steps = ["normalize_gffs", "cards", "ai_cards", "api", "final_calls", "review_bundle"]
        else:
            steps = [step for step in STEP_ORDER if step != "review_bundle"]
        if not normalize_gffs:
            steps = [step for step in steps if step != "normalize_gffs"]
        if not candidate_qc:
            steps = [step for step in steps if step not in {"validate_candidates", "audit_eligibility"}]
        if not policy_export:
            steps = [step for step in steps if step not in {"policy_export", "policy_reconcile_ids"}]
        elif not policy_reconcile:
            steps = [step for step in steps if step != "policy_reconcile_ids"]
        if not release_gff_qc:
            steps = [step for step in steps if step != "release_gff_qc"]
    else:
        steps = [item.strip() for item in raw_steps.split(",") if item.strip()]
        if "cards" in steps:
            card_index = steps.index("cards")
            if normalize_gffs and "normalize_gffs" not in steps:
                steps.insert(card_index, "normalize_gffs")
                card_index += 1
            if candidate_qc:
                for pre_step in ["validate_candidates", "audit_eligibility"]:
                    if pre_step not in steps:
                        steps.insert(card_index, pre_step)
                        card_index += 1
    unknown = sorted(set(steps) - set(STEP_ORDER))
    if unknown:
        raise SystemExit("Unknown step(s): " + ",".join(unknown))
    if skip_api and "api" in steps:
        steps.remove("api")
    if review == "review" and not explicit:
        downstream = {"export_gff", "reconcile_ids", "id_summary", "policy_export", "policy_reconcile_ids", "release_gff_qc"}
        steps = [step for step in steps if step not in downstream]
    correction_only = {"export_gff", "reconcile_ids", "policy_export", "policy_reconcile_ids", "release_gff_qc"}
    if mode != "correction":
        requested = correction_only & set(steps)
        if explicit and requested:
            raise SystemExit("Steps require correction mode with one current/base GFF: " + ",".join(sorted(requested)))
        steps = [step for step in steps if step not in correction_only]
    if not reconcile and "reconcile_ids" in steps:
        steps.remove("reconcile_ids")
    if not policy_reconcile and "policy_reconcile_ids" in steps:
        steps.remove("policy_reconcile_ids")
    if not release_gff_qc and "release_gff_qc" in steps:
        steps.remove("release_gff_qc")
    return steps


def ensure_parent_dirs(paths: Mapping[str, Path], mode: str, steps: Sequence[str]) -> None:
    keys = {"run_dir", "logs"}
    if "normalize_gffs" in steps:
        keys.update({"normalized_gff_dir", "normalized_gff_manifest", "normalized_gff_report"})
    if "validate_candidates" in steps or "audit_eligibility" in steps:
        keys.update({"candidate_qc_dir", "validation_summary", "eligibility_summary", "eligibility_reason_summary"})
    if "cards" in steps:
        keys.update({"models", "complex_loci", "summary", "units", "unit_edges", "model_sets", "full_cards"})
    if "ai_cards" in steps:
        keys.update({"ai_cards", "ai_summary", "set_trace", "auto_decisions", "review_router"})
    if "api" in steps:
        keys.add("api_run_dir")
    if "final_calls" in steps:
        keys.update({"final_calls_dir", "final_calls"})
    if "review_bundle" in steps:
        keys.update({"review_bundle_dir", "review_queue", "manual_review_template", "review_cards_dir"})
    if mode == "correction" and "export_gff" in steps:
        keys.update({"proposal_dir", "proposal_gff", "proposal_records"})
    if mode == "correction" and "reconcile_ids" in steps:
        keys.update({"final_naming_dir", "gene_id_mapping", "transcript_id_mapping", "clean_final_named_gff", "hisat_gtf"})
    if "id_summary" in steps:
        keys.add("id_summary")
    if mode == "correction" and "policy_export" in steps:
        keys.update({"policy_dir", "policy_gff", "policy_records", "policy_summary"})
    if mode == "correction" and "policy_reconcile_ids" in steps:
        keys.update({"policy_final_naming_dir", "policy_gene_id_mapping", "policy_transcript_id_mapping", "policy_clean_final_named_gff", "policy_hisat_gtf"})
    if mode == "correction" and "release_gff_qc" in steps:
        keys.update({"release_gff_qc_dir", "release_gff", "release_clean_gff", "release_gff_summary", "release_gff_mapping", "release_gff_report"})

    for key in keys:
        path = paths[key]
        if key in {"run_dir", "logs", "normalized_gff_dir", "candidate_qc_dir", "final_calls_dir", "review_bundle_dir", "review_cards_dir", "proposal_dir", "final_naming_dir", "api_run_dir", "policy_dir", "policy_final_naming_dir", "release_gff_qc_dir"}:
            path.mkdir(parents=True, exist_ok=True)
        elif path.suffix:
            path.parent.mkdir(parents=True, exist_ok=True)


def validate_paths(config: Mapping[str, Any], root: Path) -> List[str]:
    warnings: List[str] = []
    review_mode(config)
    scoring_normalization(config)
    if policy_name(config) not in {
        "supported_novel_only",
        "supported_novel_plus_risk_replacement",
        "supported_novel_multitool_only",
        "supported_novel_singletool_junction_supported",
        "strict_supported_novel_multitool_only",
        "strict_supported_novel_singletool_junction_supported",
    }:
        warnings.append("invalid_policy_export_policy:{0}".format(policy_name(config)))
    mode = workflow_mode(config)
    gffs = gff_entries(config, root)
    current_entries = [entry for entry in gffs if entry["role"] == "current"]
    candidate_entries = [entry for entry in gffs if entry["role"] == "candidate"]

    if len(gffs) < 2:
        warnings.append("too_few_gffs:provide at least two GFF files.")
    if config.get("gffs") is not None and mode == "correction" and len(gffs) > 4:
        warnings.append(f"too_many_gffs:ordered GFF workflow accepts at most four GFF files; found {len(gffs)}.")
    if mode == "correction":
        if len(current_entries) != 1:
            warnings.append(f"invalid_current_gff_count:correction mode requires exactly one current/base GFF; found {len(current_entries)}.")
        if len(candidate_entries) < 1:
            warnings.append("too_few_candidate_gffs:correction mode requires at least one candidate GFF in addition to current/base.")
    else:
        if current_entries:
            warnings.append("arbitration_mode_ignores_current_role:use correction mode if replacement mapping to current/base is required.")
        if len(candidate_entries) < 2:
            warnings.append("too_few_candidate_gffs:arbitration mode requires at least two candidate GFFs.")

    seen_names: set[str] = set()
    for entry in gffs:
        name = entry["name"]
        if name in seen_names:
            warnings.append(f"duplicate_gff_name:{name}")
        seen_names.add(name)
        if not Path(entry["path"]).exists():
            label = "current" if entry["role"] == "current" else "candidate"
            warnings.append(f"missing_{label}_gff:{entry['name']}:{entry['path']}")

    if candidate_qc_enabled(config, root):
        genome_fasta = genome_fasta_path(config, root)
        if not genome_fasta:
            warnings.append("missing_genome_fasta:candidate_qc requires genome_fasta.")
        elif not Path(genome_fasta).exists():
            warnings.append(f"missing_genome_fasta:{genome_fasta}")
    for entry in evidence_entries(config, root):
        evidence_path = Path(entry["path"])
        if not evidence_path.exists():
            warnings.append(f"missing_evidence:{entry['name']}:{entry['path']}")
            continue
        kind = entry["type"]
        if kind == "long_read_support_dir" and not evidence_path.is_dir():
            warnings.append(f"evidence_not_directory:{entry['name']}:{entry['path']}")
        elif kind in {"splice_junctions", "model_junction_support"} and not evidence_path.is_file():
            warnings.append(f"evidence_not_file:{entry['name']}:{entry['path']}")
        elif kind in {"splice_junctions", "model_junction_support"}:
            inferred = infer_evidence_type(str(evidence_path))
            if inferred != kind:
                warnings.append(
                    f"evidence_schema_mismatch:{entry['name']}:declared={kind}:detected={inferred}:{entry['path']}"
                )
    return warnings

def write_manifest(path: Path, config_path: Path, mode: str, steps: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at_utc": utc_now(),
        "config": str(config_path),
        "steps": list(steps),
        "workflow_mode": mode,
        "boundary": "No truth/manual evaluation is run by GeneArbiter. review.mode=review generates review materials and stops before final GFF export.",
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def format_command(command: Sequence[str]) -> str:
    return " ".join(json.dumps(item) if any(ch.isspace() for ch in item) else item for item in command)


class TeeWriter:
    """Write step output to both its log and the worker console."""

    def __init__(self, *handles: Any) -> None:
        self.handles = handles

    def write(self, text: str) -> int:
        for handle in self.handles:
            handle.write(text)
            handle.flush()
        return len(text)

    def flush(self) -> None:
        for handle in self.handles:
            handle.flush()


def run_python_script_in_process(command: Sequence[str], stdout: Any, stderr: Any, tee: bool = False) -> int:
    """Run a Python script command without spawning another frozen executable.

    Desktop builds launch the same executable in ``--worker`` mode.  In that
    environment ``sys.executable script.py`` would relaunch the GUI binary, so
    packaged Python step scripts are executed with runpy instead.  Normal CLI
    runs continue to use subprocesses unless GENEARBITER_IN_PROCESS_STEPS is set.
    """

    if len(command) < 2 or Path(command[0]).resolve() != Path(sys.executable).resolve():
        raise ValueError("Only sys.executable Python-script commands can run in process.")
    script = Path(command[1]).resolve()
    if not script.is_file():
        raise FileNotFoundError(f"Packaged workflow script is missing: {script}")

    old_argv = sys.argv[:]
    old_path = sys.path[:]
    display_out = TeeWriter(stdout, sys.stdout) if tee else stdout
    display_err = TeeWriter(stderr, sys.stderr) if tee else stderr
    try:
        sys.argv = [str(script), *command[2:]]
        sys.path.insert(0, str(script.parent))
        with contextlib.redirect_stdout(display_out), contextlib.redirect_stderr(display_err):
            try:
                runpy.run_path(str(script), run_name="__main__")
            except SystemExit as exc:
                if exc.code in {None, 0}:
                    return 0
                if isinstance(exc.code, int):
                    return exc.code
                print(str(exc.code), file=display_err)
                return 1
        return 0
    finally:
        sys.argv = old_argv
        sys.path[:] = old_path


def run_command(step: str, command: Sequence[str], log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = log_dir / f"{step}.stdout.log"
    stderr_path = log_dir / f"{step}.stderr.log"
    in_process = os.environ.get("GENEARBITER_IN_PROCESS_STEPS", "").strip().lower() in {"1", "true", "yes"}
    if in_process:
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            result_code = run_python_script_in_process(command, stdout, stderr, tee=step == "api")
        if result_code != 0:
            raise SystemExit(f"Step failed: {step}; returncode={result_code}; see {stderr_path}")
        return
    if step != "api":
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise SystemExit(f"Step failed: {step}; returncode={result.returncode}; see {stderr_path}")
        return

    def tee_stream(stream: Any, log_handle: Any, display_handle: Any) -> None:
        for line in iter(stream.readline, ""):
            log_handle.write(line)
            log_handle.flush()
            display_handle.write(line)
            display_handle.flush()
        stream.close()

    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        assert process.stderr is not None
        stdout_thread = threading.Thread(target=tee_stream, args=(process.stdout, stdout, sys.stdout), daemon=True)
        stderr_thread = threading.Thread(target=tee_stream, args=(process.stderr, stderr, sys.stderr), daemon=True)
        stdout_thread.start()
        stderr_thread.start()
        result_code = process.wait()
        stdout_thread.join()
        stderr_thread.join()
    if result_code != 0:
        raise SystemExit(f"Step failed: {step}; returncode={result_code}; see {stderr_path}")


def cmd_validate(args: argparse.Namespace) -> int:
    config_path = Path(args.config).expanduser().resolve()
    config = load_config(config_path)
    root = config_root(config, config_path)
    warnings = validate_paths(config, root)
    paths = output_paths(config, root)
    print(f"config\t{config_path}")
    print(f"root\t{root}")
    print(f"run_dir\t{paths['run_dir']}")
    if warnings:
        for warning in warnings:
            print(f"warning\t{warning}")
        return 1 if args.strict else 0
    print("status\tok")
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    config_path = Path(args.config).expanduser().resolve()
    config = load_config(config_path)
    root = config_root(config, config_path)
    reconcile = bool(export_section(config).get("reconcile_ids", True))
    review = review_mode(config, getattr(args, "review_mode", ""))
    steps = expand_steps(
        args.steps,
        args.skip_api,
        reconcile,
        workflow_mode(config),
        review,
        normalize_gffs_enabled(config),
        candidate_qc_enabled(config, root),
        policy_export_enabled(config),
        policy_reconcile_enabled(config),
        release_gff_qc_enabled(config),
    )
    plan = command_plan(config, root, steps, args.resume_api)
    for item in plan:
        print(f"## {item['step']}")
        print(format_command(item["command"]))
    return 0


def direct_named_path(value: str, label: str) -> tuple[str, str]:
    """Parse NAME=PATH while also accepting a bare path with a derived name."""
    if "=" in value:
        name, raw_path = value.split("=", 1)
        name = name.strip()
        raw_path = raw_path.strip()
        if not name or not raw_path:
            raise SystemExit(f"Invalid {label} value {value!r}; expected NAME=PATH.")
    else:
        raw_path = value.strip()
        name = Path(raw_path).stem
    if not raw_path:
        raise SystemExit(f"Invalid empty {label} path.")
    return safe_label(name), str(Path(raw_path).expanduser().resolve())


def direct_evidence_entry(value: str) -> Dict[str, str]:
    """Parse [TYPE=]PATH for the concise file-oriented CLI."""
    if "=" in value:
        kind, raw_path = value.split("=", 1)
        kind = kind.strip().lower()
        raw_path = raw_path.strip()
        aliases = {
            "junctions": "splice_junctions",
            "rna_junctions": "splice_junctions",
            "junction_support": "model_junction_support",
            "protein": "protein_gff",
            "long_read": "long_read_support_dir",
        }
        kind = aliases.get(kind, kind)
        if kind not in SUPPORTED_EVIDENCE_TYPES:
            allowed = ",".join(sorted(SUPPORTED_EVIDENCE_TYPES))
            raise SystemExit(f"Unsupported evidence type {kind!r}. Supported types: {allowed}")
    else:
        raw_path = value.strip()
        kind = "auto"
    if not raw_path:
        raise SystemExit("Evidence path cannot be empty.")
    path = str(Path(raw_path).expanduser().resolve())
    return {"name": safe_label(Path(raw_path).stem), "type": kind, "path": path}


def export_concise_run_results(config: Mapping[str, Any], config_path: Path) -> Dict[str, Path]:
    """Collect the three user-facing artifacts while retaining the full audit work tree."""
    try:
        from genearbiter.release_export import (
            PUBLIC_MAPPING_FIELDS,
            build_public_mapping_rows,
            read_tsv_rows,
            write_gzip_tsv,
        )
    except ModuleNotFoundError:
        from release_export import (  # type: ignore[no-redef]
            PUBLIC_MAPPING_FIELDS,
            build_public_mapping_rows,
            read_tsv_rows,
            write_gzip_tsv,
        )

    root = config_root(config, config_path)
    paths = output_paths(config, root)
    required = [
        paths["release_clean_gff"],
        paths["release_gff"],
        paths["gene_id_mapping"],
        paths["transcript_id_mapping"],
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Cannot collect final results; missing: " + ",".join(missing))

    public_dir = config_path.parent / "results"
    public_dir.mkdir()
    clean_path = public_dir / "annotation.clean.gff3"
    trace_path = public_dir / "annotation.trace.gff3"
    mapping_path = public_dir / "id_mapping.tsv.gz"
    summary_path = public_dir / "run_summary.txt"
    shutil.copy2(paths["release_clean_gff"], clean_path)
    shutil.copy2(paths["release_gff"], trace_path)

    _gene_fields, gene_rows = read_tsv_rows(paths["gene_id_mapping"])
    _tx_fields, tx_rows = read_tsv_rows(paths["transcript_id_mapping"])
    final_gene_ids = {row.get("final_gene_id", "") for row in gene_rows if row.get("final_gene_id")}
    mapping_rows = build_public_mapping_rows(gene_rows, tx_rows, final_gene_ids)
    for row in mapping_rows:
        row["included_in"] = "final"
    write_gzip_tsv(mapping_path, PUBLIC_MAPPING_FIELDS, mapping_rows)

    profile = str(api_section(config).get("profile") or "")
    candidate_qc = "enabled" if candidate_qc_enabled(config, root) else "disabled_no_genome_fasta"
    summary_path.write_text(
        "\n".join(
            [
                "status\tsuccess",
                f"decision_profile\t{profile}",
                f"candidate_coding_qc\t{candidate_qc}",
                f"clean_gff\t{clean_path.name}",
                f"trace_gff\t{trace_path.name}",
                f"id_mapping\t{mapping_path.name}",
                f"mapping_rows\t{len(mapping_rows)}",
                "mapping_scope\tchanged current lineages and selected source origins; unchanged retained IDs are omitted",
                f"audit_work_dir\t{paths['run_dir']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {"clean_gff": clean_path, "trace_gff": trace_path, "id_mapping": mapping_path, "summary": summary_path}


def cmd_run_files(args: argparse.Namespace) -> int:
    """Create a frozen run config from file arguments, then execute the normal workflow."""
    out_dir = Path(args.out_dir).expanduser().resolve()
    if out_dir.exists():
        raise SystemExit(f"Output directory already exists; choose a new directory: {out_dir}")

    current_name, current_path = direct_named_path(args.current, "current GFF")
    del current_name
    candidates = []
    for value in args.candidate:
        name, path = direct_named_path(value, "candidate GFF")
        candidates.append({"name": name, "path": path})
    if len({row["name"] for row in candidates}) != len(candidates):
        raise SystemExit("Candidate names must be unique.")

    config: Dict[str, Any] = {
        "project": {"name": args.project_name, "root": str(Path.cwd()), "run_id": out_dir.name},
        "strict_input_check": True,
        "gffs": [{"name": "current", "path": current_path}, *candidates],
        "evidence": [direct_evidence_entry(value) for value in args.evidence],
        "params": {"candidate_set_mode": "structure_consensus"},
        "api": {
            "profile": args.profile,
            "run_label": args.profile,
            "workers": args.workers,
            "max_attempts": args.max_attempts,
            "timeout": args.timeout,
            "temperature": 0.0,
            "max_tokens": args.max_tokens,
            "api_key_env": "" if args.profile == "local_rule" else args.api_key_env,
        },
        "outputs": {"run_dir": str(out_dir / "work"), "strategy": "genearbiter_default"},
        "export": {
            "source_label": "GeneArbiter",
            "sort_records": True,
            "reconcile_ids": True,
            "new_gene_prefix": args.new_gene_prefix,
            "one_to_one_policy": "reuse_current_gene_id",
            "transcript_template": "{gene_id}.t{index:02d}",
            "force_reconcile": True,
        },
        "review": {"mode": args.review_mode},
        "preprocess": {"normalize_gffs": True},
        "scoring": {"normalization": "run_level"},
        "release_gff_qc": {
            "enabled": True,
            "transcript_template": "{gene_id}.t{index:02d}",
            "exon_template": "{transcript_id}.exon{index:03d}",
            "output_gff_name": "annotation.trace.gff3",
            "clean_gff_name": "annotation.clean.gff3",
            "force": True,
        },
    }
    if args.genome:
        config["genome_fasta"] = str(Path(args.genome).expanduser().resolve())
    if args.base_url:
        config["api"]["base_url"] = args.base_url
    if args.model:
        config["api"]["model"] = args.model

    preflight_warnings = validate_paths(config, Path.cwd())
    if preflight_warnings:
        raise SystemExit("Input validation failed:\n" + "\n".join(preflight_warnings))
    if not args.genome:
        print(
            "warning\tgenome FASTA not supplied; CDS sequence validation and candidate coding QC are disabled",
            file=sys.stderr,
        )

    out_dir.mkdir(parents=True)
    config_path = out_dir / "genearbiter.config.json"
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"config\t{config_path}")
    run_args = argparse.Namespace(
        config=str(config_path),
        steps="all",
        skip_api=False,
        resume_api=True,
        review_mode=args.review_mode,
        dry_run=args.dry_run,
    )
    result = cmd_run(run_args)
    if result == 0 and not args.dry_run and args.review_mode == "auto":
        exported = export_concise_run_results(config, config_path)
        for name, path in exported.items():
            print(f"result_{name}\t{path}")
    return result


def cmd_run(args: argparse.Namespace) -> int:
    config_path = Path(args.config).expanduser().resolve()
    config = load_config(config_path)
    root = config_root(config, config_path)
    warnings = validate_paths(config, root)
    strict = bool(config.get("strict_input_check", True))
    if warnings and strict:
        raise SystemExit("Config validation failed:\n" + "\n".join(warnings))
    paths = output_paths(config, root)
    mode = workflow_mode(config)
    reconcile = bool(export_section(config).get("reconcile_ids", True))
    review = review_mode(config, args.review_mode)
    steps = expand_steps(
        args.steps,
        args.skip_api,
        reconcile,
        mode,
        review,
        normalize_gffs_enabled(config),
        candidate_qc_enabled(config, root),
        policy_export_enabled(config),
        policy_reconcile_enabled(config),
        release_gff_qc_enabled(config),
    )
    ensure_parent_dirs(paths, mode, steps)
    plan = command_plan(config, root, steps, args.resume_api)
    write_manifest(paths["run_dir"] / "genearbiter_command_manifest.json", config_path, mode, plan)
    if args.dry_run:
        for item in plan:
            print(f"## {item['step']}")
            print(format_command(item["command"]))
        return 0
    if "api" in steps and api_section(config).get("profile", "deepseek_flash") != "local_rule":
        key_env = str(api_section(config).get("api_key_env") or "DEEPSEEK_API_KEY")
        if not os.environ.get(key_env, "").strip() and not api_section(config).get("dry_run_prompts"):
            raise SystemExit(f"Missing API key environment variable: {key_env}")
    for item in plan:
        print(f"[{utc_now()}] run {item['step']}", file=sys.stderr)
        run_command(str(item["step"]), item["command"], paths["logs"])
    mode = workflow_mode(config)
    print(f"run_dir\t{paths['run_dir']}")
    print(f"final_calls\t{paths['final_calls']}")
    if "review_bundle" in steps:
        print(f"review_queue\t{paths['review_queue']}")
        print(f"manual_review_template\t{paths['manual_review_template']}")
    if mode == "correction" and "export_gff" in steps:
        print(f"proposal_gff\t{paths['proposal_gff']}")
    if mode == "correction" and "reconcile_ids" in steps:
        print(f"final_named_gff\t{paths['final_naming_dir'] / 'final_named_annotation.gff3'}")
        print(f"clean_final_named_gff\t{paths['clean_final_named_gff']}")
        print(f"hisat_gtf\t{paths['hisat_gtf']}")
    if "id_summary" in steps:
        print(f"id_mapping_summary\t{paths['id_summary']}")
    if mode == "correction" and "policy_export" in steps:
        print(f"policy_gff\t{paths['policy_gff']}")
        print(f"policy_records\t{paths['policy_records']}")
    if mode == "correction" and "policy_reconcile_ids" in steps:
        print(f"policy_final_named_gff\t{paths['policy_named_gff']}")
        print(f"policy_clean_final_named_gff\t{paths['policy_clean_final_named_gff']}")
        print(f"policy_hisat_gtf\t{paths['policy_hisat_gtf']}")
    if mode == "correction" and "release_gff_qc" in steps:
        print(f"release_gff\t{paths['release_gff']}")
        print(f"release_clean_gff\t{paths['release_clean_gff']}")
        print(f"release_gff_summary\t{paths['release_gff_summary']}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genearbiter",
        description="Evidence-constrained gene annotation arbitration and release workflow.",
    )
    parser.add_argument("--version", action="version", version=f"GeneArbiter {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate-config", help="Check config shape and input paths.")
    validate.add_argument("--config", required=True)
    validate.add_argument("--strict", action="store_true", help="Return non-zero when warnings are found.")
    validate.set_defaults(func=cmd_validate)

    add_extract_short_read_parser(sub)
    add_validate_evidence_parser(sub)

    run_files = sub.add_parser(
        "run-files",
        help="Run from GFF/evidence file arguments and save the generated config.",
    )
    run_files.add_argument("--current", required=True, help="Current/backbone GFF3 path.")
    run_files.add_argument(
        "--candidate",
        action="append",
        required=True,
        help="Candidate GFF3 as NAME=PATH; repeat for additional sources.",
    )
    run_files.add_argument(
        "--evidence",
        action="append",
        default=[],
        help="Evidence as [TYPE=]PATH; repeatable. Bare paths are type-detected.",
    )
    run_files.add_argument("--genome", default="", help="Matching genome FASTA; recommended for coding QC.")
    run_files.add_argument("--out-dir", required=True, help="New output directory; existing directories are refused.")
    run_files.add_argument("--project-name", default="genearbiter_run")
    run_files.add_argument("--new-gene-prefix", default="GeneArbiterG")
    run_files.add_argument(
        "--profile",
        choices=["local_rule", "deepseek_flash", "deepseek_v4_flash", "deepseek_flash_thinking", "deepseek_pro"],
        default="local_rule",
        help="Decision profile. local_rule is offline and deterministic.",
    )
    run_files.add_argument("--api-key-env", default="DEEPSEEK_API_KEY")
    run_files.add_argument("--base-url", default="")
    run_files.add_argument("--model", default="")
    run_files.add_argument("--workers", type=int, default=1)
    run_files.add_argument("--max-attempts", type=int, default=2)
    run_files.add_argument("--timeout", type=int, default=180)
    run_files.add_argument("--max-tokens", type=int, default=2048)
    run_files.add_argument("--review-mode", choices=["auto", "review"], default="auto")
    run_files.add_argument("--dry-run", action="store_true")
    run_files.set_defaults(func=cmd_run_files)

    plan = sub.add_parser("plan", help="Print commands that would be run.")
    plan.add_argument("--config", required=True)
    plan.add_argument("--steps", default="all", help="all or comma list: " + ",".join(STEP_ORDER))
    plan.add_argument("--skip-api", action="store_true")
    plan.add_argument("--resume-api", action="store_true", default=True)
    plan.add_argument("--review-mode", choices=["auto", "review"], default="", help="Override config review.mode.")
    plan.set_defaults(func=cmd_plan)

    run = sub.add_parser("run", help="Run configured workflow steps.")
    run.add_argument("--config", required=True)
    run.add_argument("--steps", default="all", help="all or comma list: " + ",".join(STEP_ORDER))
    run.add_argument("--skip-api", action="store_true", help="Skip API step; downstream final-calls require an existing decisions.jsonl.")
    run.add_argument("--resume-api", action="store_true", default=True)
    run.add_argument("--review-mode", choices=["auto", "review"], default="", help="Override config review.mode.")
    run.add_argument("--dry-run", action="store_true", help="Write manifest and print commands without executing them.")
    run.set_defaults(func=cmd_run)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
