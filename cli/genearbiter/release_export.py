#!/usr/bin/env python3
# script_id_md5: b808e7cd6d3b83c5ac2a003bdc5d5290
# created: 2026-07-14
# modified: 2026-09-08
# owner: project
# status: project_code
# purpose: 从源 GeneArbiter artifacts 重新导出 SJ/SNM release-like GFF、转换 TSV 和中英文文件说明。
# inputs: release source manifest TSV; proposal_catalog.gff3; proposal_catalog_records.tsv; final_annotation_calls.tsv; model_arbitration_models.tsv; model_arbitration_full_cards.jsonl; optional current_gff.
# outputs: release/SJ 和 release/SNM 下结构一致的 clean/trace GFF、内部审计 TSV、public_mapping 单表、manifest 和 file_description.txt。
# notes: 不读取 truth annotation；不复用旧 strict_SJ/strict_SNM 导出；只从源 artifacts 按 strict policy 重新生成发布包。

"""Export GeneArbiter release bundles from source run artifacts.

This command intentionally does not collect legacy intermediate directories. It
starts from frozen GeneArbiter source artifacts and reruns the two release policies:

* SJ  -> strict_supported_novel_singletool_junction_supported
* SNM -> strict_supported_novel_multitool_only

SJ is the canonical naming release: policy export, ID reconciliation and release
GFF QC run once. SNM is then derived as a strict subset of the named SJ release,
so a proposal can never receive a second ID. No truth/manual evaluation file is
read. Multi-tool support is counted from exact selected CDS structures in the
full cards, not from the representative-source fields in final calls.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


MODE_POLICIES = {
    "SJ": "strict_supported_novel_singletool_junction_supported",
    "SNM": "strict_supported_novel_multitool_only",
}
PUBLIC_MAPPING_FIELDS = [
    "feature_type",
    "mapping_role",
    "source_name",
    "source_id",
    "final_id",
    "final_parent_id",
    "change_type",
    "included_in",
]
PUBLIC_CHANGE_TYPES = {
    "one_to_one_replacement": "replacement",
    "split_replacement": "split",
    "merge_replacement": "merge",
    "complex_replacement": "complex",
    "novel_gene": "novel",
}
REQUIRED_SOURCE_FIELDS = ["species_key", "proposal_gff", "proposal_records", "final_calls", "models"]
INPUT_PATH_FIELDS = ["proposal_gff", "proposal_records", "final_calls", "models"]
SOURCE_FIELDS = REQUIRED_SOURCE_FIELDS + [
    "full_cards",
    "current_gff",
    "new_gene_prefix",
    "id_strategy",
    "transcript_id_strategy",
    "transcript_template",
    "exon_template",
    "source_label",
    "current_backbone_merge",
    "include",
    "notes",
]
MANIFEST_FIELDS = [
    "mode",
    "policy",
    "species_key",
    "status",
    "tsv_complete",
    "release_qc_repair_rows",
    "release_gff",
    "release_clean_gff",
    "gene_id_mapping",
    "transcript_id_mapping",
    "locus_id_mapping",
    "mode_membership",
    "release_gff_id_mapping",
    "policy_catalog_records",
    "policy_catalog_summary",
    "release_gff_qc_summary",
    "features",
    "genes",
    "transcripts",
    "synthetic_gene_records",
    "synthetic_exon_records",
    "postcheck_missing_parent_refs",
    "postcheck_invalid_cds_phase",
    "postcheck_duplicate_gene_tx_exon_ids",
    "postcheck_child_outside_parent",
    "current_backbone_merge",
    "removed_current_ids",
    "current_feature_lines_removed",
    "source_proposal_gff",
    "source_proposal_records",
    "source_final_calls",
    "source_models",
    "source_full_cards",
    "notes",
]
TRUTH_LIKE_WORDS = ("truth", "manual")


class ReleaseExportError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", required=True, help="TSV listing source GeneArbiter artifacts.")
    parser.add_argument("--out-dir", required=True, help="Release bundle output directory.")
    parser.add_argument("--modes", default="SJ,SNM", help="Comma list of modes: SJ,SNM.")
    parser.add_argument("--work-dir", default="", help="Intermediate work directory. Defaults to <out-dir>/_work.")
    parser.add_argument("--new-gene-prefix", default="auto", help="Fallback prefix when source-style inference is unavailable; auto uses GeneArbiterG.")
    parser.add_argument("--id-strategy", choices=["sequential", "current_style_positional", "current_style_or_sequential"], default="current_style_positional")
    parser.add_argument("--transcript-id-strategy", choices=["template", "current_style"], default="current_style")
    parser.add_argument("--transcript-template", default="{gene_id}.t{index:02d}")
    parser.add_argument("--exon-template", default="{transcript_id}.exon{index:03d}")
    parser.add_argument("--description-name", default="file_description.txt")
    parser.add_argument("--force", action="store_true", help="Overwrite existing per-species outputs.")
    parser.add_argument("--update-manifest", action="store_true", help="Replace retried mode/species rows in existing aggregate manifests and preserve all other rows.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned species/modes without running scripts.")
    parser.add_argument("--fail-on-truth-path", action="store_true", default=True, help="Refuse truth/manual-like paths.")
    return parser.parse_args()


def script_dir() -> Path:
    return Path(__file__).resolve().parent


def common_dir() -> Path:
    return script_dir() / "common"


def common_script(name: str) -> Path:
    scripts = {
        "policy_export": "export_policy_gff.py",
        "reconcile_ids": "reconcile_ids.py",
        "release_gff_qc": "release_gff.py",
    }
    path = common_dir() / scripts[name]
    if not path.exists():
        raise ReleaseExportError(f"Missing script: {path}")
    return path


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open(encoding="utf-8")


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_text.rstrip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip('"')
    return attrs


def split_list(value: object) -> List[str]:
    return [item.strip() for item in str(value or "").replace(";", ",").split(",") if item.strip()]


def public_proposal_id(value: object) -> str:
    """Normalize legacy internal proposal IDs for release-facing tables."""
    text = str(value or "")
    if text.startswith("V2PROP_"):
        return "GeneArbiterProposal_" + text[len("V2PROP_") :]
    return text


def bool_value(value: object, default: bool = False) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return default
    return text in {"1", "true", "yes", "y", "on", "auto"}


def fail_if_truth_like(paths: Iterable[Path]) -> None:
    for path in paths:
        text = str(path).lower()
        if any(word in text for word in TRUTH_LIKE_WORDS):
            raise ReleaseExportError(f"Refusing truth/manual-like path: {path}")


def read_source_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in REQUIRED_SOURCE_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise ReleaseExportError("source manifest missing required fields: " + ",".join(missing))
        rows = []
        for row in reader:
            if str(row.get("include", "true")).strip().lower() in {"0", "false", "no", "skip"}:
                continue
            rows.append({key: str(row.get(key, "") or "") for key in SOURCE_FIELDS})
    return rows


def merge_existing_manifest(path: Path, new_rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    """Update selected mode/species rows while retaining completed batch members."""
    merged: Dict[Tuple[str, str], Dict[str, object]] = {}
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                merged[(row.get("mode", ""), row.get("species_key", ""))] = dict(row)
    for row in new_rows:
        merged[(str(row.get("mode", "")), str(row.get("species_key", "")))] = dict(row)
    return sorted(merged.values(), key=lambda row: (str(row.get("mode", "")), str(row.get("species_key", ""))))


def require_inputs(row: Mapping[str, str], fail_truth: bool) -> None:
    paths = [Path(row[field]) for field in INPUT_PATH_FIELDS]
    if not row.get("full_cards"):
        raise ReleaseExportError(
            "Missing full_cards for {0}; strict SJ/SNM export requires exact CDS support evidence".format(
                row.get("species_key")
            )
        )
    paths.append(Path(row["full_cards"]))
    if row.get("current_gff"):
        paths.append(Path(row["current_gff"]))
    for path in paths:
        if not path.exists():
            raise ReleaseExportError(f"Missing input for {row.get('species_key')}: {path}")
    if fail_truth:
        fail_if_truth_like(paths)


def run_cmd(command: Sequence[str], stdout_path: Path, stderr_path: Path) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if result.returncode != 0:
        raise ReleaseExportError(f"Command failed ({result.returncode}): {' '.join(command)}; see {stderr_path}")


def count_gff(path: Path) -> Counter:
    counts: Counter = Counter()
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                counts["malformed_rows"] += 1
                continue
            feature = parts[2]
            attrs = parse_attrs(parts[8])
            counts["features"] += 1
            if feature == "gene":
                counts["genes"] += 1
            if feature.lower() in {"mrna", "transcript", "rna", "lnc_rna", "ncrna", "rrna", "trna", "mirna", "snrna", "snorna", "pre_mirna", "primary_transcript"}:
                counts["transcripts"] += 1
            if attrs.get("ID"):
                counts["ids"] += 1
    return counts


def count_data_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(encoding="utf-8") as handle:
        return max(sum(1 for line in handle if line.strip()) - 1, 0)


def read_summary(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row.get("metric", ""): row.get("value", "") for row in csv.DictReader(handle, delimiter="\t") if row.get("metric")}


def read_policy_records(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def current_gene_feature_count(path: Path) -> int:
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 9 and parts[2] == "gene":
                count += 1
    return count


def collect_replaced_current_ids(policy_records: Path) -> List[str]:
    ids = set()
    for row in read_policy_records(policy_records):
        action = row.get("proposal_action", "")
        if action != "replace_current_with_selected_set":
            continue
        for key in ("removed_current_gene_ids", "current_gene_ids"):
            ids.update(split_list(row.get(key, "")))
    return sorted(ids)


def row_touches_removed_id(attrs: Mapping[str, str], removed: set[str]) -> bool:
    for key in ("ID", "Parent", "geneID", "gene_id"):
        for value in split_list(attrs.get(key, "")):
            if value in removed:
                return True
    return False


def maybe_merge_current_backbone(current_gff: Path, policy_gff: Path, policy_records: Path, out_gff: Path, mode: str) -> Tuple[bool, int, int]:
    if not current_gff or not current_gff.exists():
        return False, 0, 0
    if current_gene_feature_count(current_gff) > 0:
        return False, 0, 0
    removed_ids = set(collect_replaced_current_ids(policy_records))
    if not removed_ids:
        return False, 0, 0

    removed_lines = 0
    with out_gff.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        out.write("# generated_by=genearbiter_release_export_current_backbone_merge\n")
        out.write(f"# generated_at_utc={utc_now()}\n")
        out.write(f"# mode={mode}\n")
        out.write(f"# current_backbone={current_gff}\n")
        out.write(f"# policy_gff={policy_gff}\n")
        out.write(f"# removed_current_ids={len(removed_ids)}\n")
        out.write("# truth_usage=none\n")
        with open_text(current_gff) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue
                attrs = parse_attrs(parts[8])
                if row_touches_removed_id(attrs, removed_ids):
                    removed_lines += 1
                    continue
                out.write(line if line.endswith("\n") else line + "\n")
        with open_text(policy_gff) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                out.write(line if line.endswith("\n") else line + "\n")
    return True, len(removed_ids), removed_lines


def copy_if_exists(src: Path, dst: Path) -> str:
    if not src.exists():
        return ""
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return str(dst)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_gzip_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    """Write a reproducible gzip-compressed public TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(fields),
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_public_mapping_rows(
    sj_gene_rows: Sequence[Mapping[str, str]],
    sj_tx_rows: Sequence[Mapping[str, str]],
    snm_final_gene_ids: Set[str],
) -> List[Dict[str, str]]:
    """Build one compact edge-list mapping shared by the SJ and SNM releases.

    The public table covers changed/selected models. Unchanged current IDs are
    intentionally not repeated because their public IDs are already unchanged.
    """
    sj_final_gene_ids = {row.get("final_gene_id", "") for row in sj_gene_rows if row.get("final_gene_id")}
    unexpected_snm = snm_final_gene_ids - sj_final_gene_ids
    if unexpected_snm:
        raise ReleaseExportError(
            "SNM public mapping IDs are not a subset of SJ: " + ",".join(sorted(unexpected_snm)[:10])
        )

    result: List[Dict[str, str]] = []
    seen: Set[Tuple[str, ...]] = set()

    def append_row(values: Dict[str, str]) -> None:
        key = tuple(values.get(field, "") for field in PUBLIC_MAPPING_FIELDS)
        if key not in seen:
            seen.add(key)
            result.append(values)

    for row in sj_gene_rows:
        final_id = row.get("final_gene_id", "")
        source_name = row.get("source", "") or row.get("selected_source", "")
        source_id = row.get("source_gene_id", "")
        relation = row.get("final_relation", "")
        if not final_id or not source_name or not source_id:
            raise ReleaseExportError("Incomplete selected gene origin in public mapping input")
        if relation not in PUBLIC_CHANGE_TYPES:
            raise ReleaseExportError(f"Unsupported public mapping relation: {relation}")
        change_type = PUBLIC_CHANGE_TYPES[relation]
        included_in = "SJ,SNM" if final_id in snm_final_gene_ids else "SJ"
        append_row(
            {
                "feature_type": "gene",
                "mapping_role": "selected_origin",
                "source_name": source_name,
                "source_id": source_id,
                "final_id": final_id,
                "final_parent_id": "",
                "change_type": change_type,
                "included_in": included_in,
            }
        )
        current_ids = set(split_list(row.get("current_gene_ids", "")))
        current_ids.update(split_list(row.get("retired_current_gene_ids", "")))
        for current_id in sorted(current_ids):
            append_row(
                {
                    "feature_type": "gene",
                    "mapping_role": "current_lineage",
                    "source_name": "current",
                    "source_id": current_id,
                    "final_id": final_id,
                    "final_parent_id": "",
                    "change_type": change_type,
                    "included_in": included_in,
                }
            )

    for row in sj_tx_rows:
        final_id = row.get("final_transcript_id", "")
        final_parent_id = row.get("final_gene_id", "")
        source_name = row.get("source", "") or row.get("selected_source", "")
        source_id = row.get("source_transcript_id", "")
        if not final_id or not final_parent_id or not source_name or not source_id:
            raise ReleaseExportError("Incomplete selected transcript origin in public mapping input")
        if final_parent_id not in sj_final_gene_ids:
            raise ReleaseExportError(f"Public transcript mapping has unknown final parent: {final_parent_id}")
        append_row(
            {
                "feature_type": "transcript",
                "mapping_role": "selected_origin",
                "source_name": source_name,
                "source_id": source_id,
                "final_id": final_id,
                "final_parent_id": final_parent_id,
                "change_type": "",
                "included_in": "SJ,SNM" if final_parent_id in snm_final_gene_ids else "SJ",
            }
        )

    feature_rank = {"gene": 0, "transcript": 1}
    role_rank = {"current_lineage": 0, "selected_origin": 1}
    result.sort(
        key=lambda row: (
            row.get("final_parent_id", "") or row.get("final_id", ""),
            feature_rank.get(row.get("feature_type", ""), 9),
            role_rank.get(row.get("mapping_role", ""), 9),
            row.get("final_id", ""),
            row.get("source_name", ""),
            row.get("source_id", ""),
        )
    )
    return result


def export_public_mappings(release_dir: Path, out_dir: Path, force: bool = False) -> List[Dict[str, object]]:
    """Export one compact SJ/SNM mapping table per species from a release."""
    sj_dir = release_dir / "SJ"
    snm_dir = release_dir / "SNM"
    if not sj_dir.is_dir() or not snm_dir.is_dir():
        raise ReleaseExportError("Public mapping export requires both SJ and SNM release directories")
    sj_species = {path.name for path in sj_dir.iterdir() if path.is_dir()}
    snm_species = {path.name for path in snm_dir.iterdir() if path.is_dir()}
    if sj_species != snm_species:
        raise ReleaseExportError("SJ and SNM species sets differ; cannot build shared public mappings")
    if not sj_species:
        raise ReleaseExportError("No species directories found for public mapping export")
    if out_dir.exists() and any(out_dir.iterdir()) and not force:
        raise ReleaseExportError(f"Public mapping output exists; pass --force: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: List[Dict[str, object]] = []
    for species in sorted(sj_species):
        sj_gene_path = sj_dir / species / f"{species}.SJ.gene_id_mapping.tsv"
        sj_tx_path = sj_dir / species / f"{species}.SJ.transcript_id_mapping.tsv"
        snm_gene_path = snm_dir / species / f"{species}.SNM.gene_id_mapping.tsv"
        for path in (sj_gene_path, sj_tx_path, snm_gene_path):
            if not path.exists():
                raise ReleaseExportError(f"Missing public mapping input: {path}")
        _gene_fields, sj_gene_rows = read_tsv_rows(sj_gene_path)
        _tx_fields, sj_tx_rows = read_tsv_rows(sj_tx_path)
        _snm_fields, snm_gene_rows = read_tsv_rows(snm_gene_path)
        snm_final_gene_ids = {
            row.get("final_gene_id", "") for row in snm_gene_rows if row.get("final_gene_id")
        }
        mapping_rows = build_public_mapping_rows(sj_gene_rows, sj_tx_rows, snm_final_gene_ids)
        mapping_path = out_dir / f"{species}.GeneArbiter.id_mapping.tsv.gz"
        write_gzip_tsv(mapping_path, PUBLIC_MAPPING_FIELDS, mapping_rows)
        counts = Counter(row["feature_type"] for row in mapping_rows)
        roles = Counter(row["mapping_role"] for row in mapping_rows)
        memberships = Counter(row["included_in"] for row in mapping_rows)
        manifest_rows.append(
            {
                "species_key": species,
                "mapping_file": str(mapping_path),
                "rows": len(mapping_rows),
                "gene_rows": counts.get("gene", 0),
                "transcript_rows": counts.get("transcript", 0),
                "current_lineage_rows": roles.get("current_lineage", 0),
                "selected_origin_rows": roles.get("selected_origin", 0),
                "sj_only_rows": memberships.get("SJ", 0),
                "sj_snm_rows": memberships.get("SJ,SNM", 0),
                "sha256": sha256_file(mapping_path),
            }
        )

    write_tsv(
        out_dir / "manifest.tsv",
        [
            "species_key",
            "mapping_file",
            "rows",
            "gene_rows",
            "transcript_rows",
            "current_lineage_rows",
            "selected_origin_rows",
            "sj_only_rows",
            "sj_snm_rows",
            "sha256",
        ],
        manifest_rows,
    )
    return manifest_rows


def write_description(path: Path, mode: str, policy: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "GeneArbiter 发布包文件说明",
        "======================",
        f"模式：{mode}",
        f"导出 policy：{policy}",
        "本目录由 GeneArbiter release exporter 从源 run artifacts 重新生成，不复用旧 strict_SJ/strict_SNM 目录。SJ 统一命名一次，SNM 从 SJ canonical ID 集合派生。",
        "所有 GFF 坐标来自 current/backbone 或已冻结 proposal；不读取 truth/manual evaluation 文件。",
        "网站公开 ID 对应为 release/public_mapping/<species>.GeneArbiter.id_mapping.tsv.gz；本模式目录中的多张 TSV 是内部审计 artifact，不应整包公开。",
        "",
        "字段说明：file 表示文件名；role 表示用途；format 表示文件格式；key_columns 表示关键字段；notes 表示使用边界。",
        "",
        "file\trole\tformat\tkey_columns\tnotes",
        "<species>.<mode>.release.clean.gff3\t干净发布版 GFF3\tGFF3\tID,Parent,Name,biotype,protein_id\t标准结构文件；source 列统一为 GeneArbiter，Name 仅在不等于 ID 时保留。",
        "<species>.<mode>.release.trace.gff3\t精简溯源版 GFF3\tGFF3\tchange_type,origin_source,origin_gene_id,current_gene_ids,review_recommended,review_reason\t与 clean GFF 的结构和 ID 相同；溯源字段只写在 gene/transcript 层，review 标记是检查优先级而不是错误概率。",
        "../public_mapping/<species>.GeneArbiter.id_mapping.tsv.gz\t公开精简 ID 对应表\tTSV.GZ\tfeature_type,mapping_role,source_name,source_id,final_id,final_parent_id,change_type,included_in\tSJ/SNM 共用一张表；只记录变更模型的 current lineage 和入选来源。",
        "<species>.<mode>.gene_id_mapping.tsv\tgene ID 转换表\tTSV\tfinal_gene_id,proposal_gene_id,current_gene_ids,source_gene_id\t来自 ID reconciliation，追踪 current/proposal/source 到发布 gene ID。",
        "<species>.<mode>.transcript_id_mapping.tsv\ttranscript ID 转换表\tTSV\tfinal_transcript_id,proposal_transcript_id,source_transcript_id\t来自 ID reconciliation，追踪 transcript ID。",
        "<species>.<mode>.locus_id_mapping.tsv\tlocus 变更映射表\tTSV\tcurrent_gene_ids,corrected_gene_ids,retired_current_gene_ids,added_gene_ids\t逐 locus 记录矫正前后 ID、保留、停用和新增关系。",
        "<species>.SNM.mode_membership.tsv\tSJ/SNM 成员表\tTSV\tfinal_gene_id,SJ,SNM,status\tSNM 排除项标记 excluded_by_SNM_policy，ID 不回收。",
        "<species>.<mode>.release_gff_id_mapping.tsv\trelease QC 结构修复 ID 表\tTSV\tfeature,old_id,new_id,parent_id,reason\t记录顶层 transcript 拆分、合成 gene/exon 或重复 ID 修复。",
        "<species>.<mode>.policy_catalog_records.tsv\tpolicy proposal 记录表\tTSV\tlocus_id,proposal_action,export_gene_id,selected_cds_support_sources,selected_cds_support_count,policy_decision_reason\t记录本 policy 实际纳入的 proposal 及精确 CDS 多工具支持。",
        "<species>.<mode>.policy_catalog_summary.tsv\tpolicy 导出摘要\tTSV\tmetric,value\t记录 policy 过滤、恢复 current、纳入 proposal 的计数。",
        "<species>.<mode>.release_gff_qc_summary.tsv\tGFF 发布修复摘要\tTSV\tmetric,value\t记录合成 gene/exon、Parent/phase postcheck 等 QC 计数。",
        "manifest_<mode>.tsv\t模式总清单\tTSV\tmode,policy,species_key,status,features,genes,transcripts\t每个物种一行，记录输出路径和 QC 关键计数。",
        "",
        "GeneArbiter release bundle file description",
        "======================================",
        f"Mode: {mode}",
        f"Export policy: {policy}",
        "This directory is regenerated from source GeneArbiter run artifacts by the release exporter; SJ is named once and SNM is derived from the SJ canonical ID set.",
        "All GFF coordinates come from the current/backbone annotation or frozen proposal records. No truth/manual evaluation file is read.",
        "The website-facing ID crosswalk is release/public_mapping/<species>.GeneArbiter.id_mapping.tsv.gz; the multiple TSVs in this mode directory are internal audit artifacts and should not be published as a bundle.",
        "",
        "Fields: file is the file name; role is the intended use; format is the file format; key_columns lists important columns; notes states usage boundaries.",
        "",
        "file\trole\tformat\tkey_columns\tnotes",
        "<species>.<mode>.release.clean.gff3\tclean release-like GFF3\tGFF3\tID,Parent,Name,biotype,protein_id\tStandard structural annotation with a uniform GeneArbiter source; Name is retained only when it differs from ID.",
        "<species>.<mode>.release.trace.gff3\tcompact trace GFF3\tGFF3\tchange_type,origin_source,origin_gene_id,current_gene_ids,review_recommended,review_reason\tHas the same structures and IDs as the clean GFF; review fields prioritize inspection and are not calibrated error probabilities.",
        "../public_mapping/<species>.GeneArbiter.id_mapping.tsv.gz\tpublic compact ID mapping\tTSV.GZ\tfeature_type,mapping_role,source_name,source_id,final_id,final_parent_id,change_type,included_in\tOne shared SJ/SNM table containing current lineages and selected origins for changed models.",
        "<species>.<mode>.gene_id_mapping.tsv\tgene ID mapping\tTSV\tfinal_gene_id,proposal_gene_id,current_gene_ids,source_gene_id\tTracks current/proposal/source IDs to release gene IDs.",
        "<species>.<mode>.transcript_id_mapping.tsv\ttranscript ID mapping\tTSV\tfinal_transcript_id,proposal_transcript_id,source_transcript_id\tTracks transcript IDs through ID reconciliation.",
        "<species>.<mode>.locus_id_mapping.tsv\tlocus change mapping\tTSV\tcurrent_gene_ids,corrected_gene_ids,retired_current_gene_ids,added_gene_ids\tTracks before/after, retained, retired and added IDs per locus.",
        "<species>.SNM.mode_membership.tsv\tSJ/SNM membership\tTSV\tfinal_gene_id,SJ,SNM,status\tMarks excluded_by_SNM_policy without recycling canonical SJ IDs.",
        "<species>.<mode>.release_gff_id_mapping.tsv\trelease QC repair ID mapping\tTSV\tfeature,old_id,new_id,parent_id,reason\tRecords top-level transcript promotion, synthetic gene/exon creation, or duplicate ID repair.",
        "<species>.<mode>.policy_catalog_records.tsv\tpolicy proposal records\tTSV\tlocus_id,proposal_action,export_gene_id,selected_cds_support_sources,selected_cds_support_count,policy_decision_reason\tRecords included proposals and exact-CDS multi-tool support.",
        "<species>.<mode>.policy_catalog_summary.tsv\tpolicy export summary\tTSV\tmetric,value\tCounts policy filtering, current restoration, and included proposal blocks.",
        "<species>.<mode>.release_gff_qc_summary.tsv\trelease GFF QC summary\tTSV\tmetric,value\tCounts synthetic gene/exon records and Parent/phase postchecks.",
        "manifest_<mode>.tsv\tmode-level manifest\tTSV\tmode,policy,species_key,status,features,genes,transcripts\tOne row per species with output paths and key QC counts.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def mode_policy(mode: str) -> str:
    if mode not in MODE_POLICIES:
        raise ReleaseExportError(f"Unsupported mode: {mode}; expected SJ/SNM")
    return MODE_POLICIES[mode]


def export_one(row: Mapping[str, str], mode: str, args: argparse.Namespace, out_dir: Path, work_dir: Path) -> Dict[str, object]:
    if mode != "SJ":
        raise ReleaseExportError("ID reconciliation is canonical in SJ only; derive SNM with export_snm_from_sj")
    species = row["species_key"]
    policy = mode_policy(mode)
    species_out = out_dir / mode / species
    species_work = work_dir / mode / species
    log_dir = species_work / "logs"
    if species_out.exists() and any(species_out.iterdir()) and not args.force:
        raise ReleaseExportError(f"Output exists for {mode}/{species}; pass --force: {species_out}")
    species_out.mkdir(parents=True, exist_ok=True)
    species_work.mkdir(parents=True, exist_ok=True)

    policy_dir = species_work / "01_policy_export"
    exact_support_cache = species_work / "00_exact_cds_support" / "exact_cds_support.tsv"
    naming_dir = species_work / "02_id_reconcile"
    release_dir = species_work / "03_release_gff_qc"
    current_gff = Path(row.get("current_gff", "")) if row.get("current_gff") else None

    policy_cmd = [
        sys.executable,
        str(common_script("policy_export")),
        "--proposal-gff",
        row["proposal_gff"],
        "--proposal-records",
        row["proposal_records"],
        "--final-calls",
        row["final_calls"],
        "--models",
        row["models"],
        "--full-cards",
        row["full_cards"],
        "--exact-support-cache",
        str(exact_support_cache),
        "--out-dir",
        str(policy_dir),
        "--policy",
        policy,
        "--fail-on-truth-path",
    ]
    if current_gff:
        policy_cmd.extend(["--current-gff", str(current_gff)])
    run_cmd(policy_cmd, log_dir / "policy_export.stdout.log", log_dir / "policy_export.stderr.log")

    policy_gff = policy_dir / "policy_catalog.gff3"
    policy_records = policy_dir / "policy_catalog_records.tsv"
    policy_summary = policy_dir / "policy_catalog_summary.tsv"
    merged_policy_gff = species_work / "01_policy_export_current_backbone_merge" / "policy_catalog.current_backbone_merged.gff3"
    merge_setting = row.get("current_backbone_merge", "auto").strip().lower() or "auto"
    do_merge = merge_setting in {"1", "true", "yes", "auto"}
    merged = False
    removed_current_ids = 0
    current_feature_lines_removed = 0
    reconcile_input_gff = policy_gff
    if do_merge and current_gff:
        merged_policy_gff.parent.mkdir(parents=True, exist_ok=True)
        merged, removed_current_ids, current_feature_lines_removed = maybe_merge_current_backbone(
            current_gff,
            policy_gff,
            policy_records,
            merged_policy_gff,
            mode,
        )
        if merged:
            reconcile_input_gff = merged_policy_gff

    new_gene_prefix = row.get("new_gene_prefix") or args.new_gene_prefix
    id_strategy = row.get("id_strategy") or args.id_strategy
    transcript_id_strategy = row.get("transcript_id_strategy") or args.transcript_id_strategy
    transcript_template = row.get("transcript_template") or args.transcript_template
    exon_template = row.get("exon_template") or args.exon_template
    reconcile_cmd = [
        sys.executable,
        str(common_script("reconcile_ids")),
        "--proposal-gff",
        str(reconcile_input_gff),
        "--proposal-records",
        str(policy_records),
        "--out-dir",
        str(naming_dir),
        "--output-gff-name",
        "policy_named_annotation.gff3",
        "--clean-gff-name",
        "policy_named_annotation.auto_pass.clean.gff3",
        "--gtf-name",
        "policy_named_annotation.auto_pass.hisat.gtf",
        "--new-gene-prefix",
        new_gene_prefix,
        "--id-strategy",
        id_strategy,
        "--transcript-id-strategy",
        transcript_id_strategy,
        "--transcript-template",
        transcript_template,
        "--one-to-one-policy",
        "reuse_current_gene_id",
        "--force",
        "--fail-on-truth-path",
    ]
    if row.get("source_label"):
        reconcile_cmd.extend(["--proposal-source-label", row["source_label"], "--clean-source-label", row["source_label"]])
    if current_gff:
        reconcile_cmd.extend(["--current-gff", str(current_gff)])
    run_cmd(reconcile_cmd, log_dir / "reconcile_ids.stdout.log", log_dir / "reconcile_ids.stderr.log")

    release_cmd = [
        sys.executable,
        str(common_script("release_gff_qc")),
        "--input-gff",
        str(naming_dir / "policy_named_annotation.gff3"),
        "--out-dir",
        str(release_dir),
        "--output-gff-name",
        "release_annotation.gff3",
        "--clean-gff-name",
        "release_annotation.clean.gff3",
        "--summary-name",
        "release_gff_qc_summary.tsv",
        "--mapping-name",
        "release_gff_id_mapping.tsv",
        "--report-name",
        "release_gff_report.md",
        "--transcript-template",
        transcript_template,
        "--exon-template",
        exon_template,
        "--force",
        "--fail-on-truth-path",
    ]
    if row.get("source_label"):
        release_cmd.extend(["--source-label", row["source_label"]])
    run_cmd(release_cmd, log_dir / "release_gff_qc.stdout.log", log_dir / "release_gff_qc.stderr.log")

    prefix = f"{species}.{mode}"
    release_gff = copy_if_exists(release_dir / "release_annotation.gff3", species_out / f"{prefix}.release.trace.gff3")
    release_clean_gff = copy_if_exists(release_dir / "release_annotation.clean.gff3", species_out / f"{prefix}.release.clean.gff3")
    gene_map = copy_if_exists(naming_dir / "gene_id_mapping.tsv", species_out / f"{prefix}.gene_id_mapping.tsv")
    tx_map = copy_if_exists(naming_dir / "transcript_id_mapping.tsv", species_out / f"{prefix}.transcript_id_mapping.tsv")
    locus_map = copy_if_exists(naming_dir / "locus_id_mapping.tsv", species_out / f"{prefix}.locus_id_mapping.tsv")
    release_map = copy_if_exists(release_dir / "release_gff_id_mapping.tsv", species_out / f"{prefix}.release_gff_id_mapping.tsv")
    policy_records_out = write_public_policy_records(policy_records, species_out / f"{prefix}.policy_catalog_records.tsv")
    policy_summary_out = copy_if_exists(policy_summary, species_out / f"{prefix}.policy_catalog_summary.tsv")
    release_summary_out = copy_if_exists(release_dir / "release_gff_qc_summary.tsv", species_out / f"{prefix}.release_gff_qc_summary.tsv")
    copy_if_exists(naming_dir / "current_id_disambiguation.tsv", species_out / f"{prefix}.current_id_disambiguation.tsv")
    copy_if_exists(release_dir / "release_gff_report.md", species_out / f"{prefix}.release_gff_report.md")

    counts = count_gff(Path(release_clean_gff)) if release_clean_gff else Counter()
    qc_summary = read_summary(Path(release_summary_out)) if release_summary_out else {}
    critical_qc_metrics = [
        "postcheck_missing_parent_refs",
        "postcheck_invalid_cds_phase",
        "postcheck_duplicate_gene_tx_exon_ids",
        "postcheck_child_outside_parent",
    ]
    failed_qc = {metric: int(qc_summary.get(metric, "0") or 0) for metric in critical_qc_metrics if int(qc_summary.get(metric, "0") or 0) != 0}
    if failed_qc:
        raise ReleaseExportError("Release GFF critical postcheck failed for {0}/SJ: {1}".format(species, failed_qc))
    required_tsvs = [gene_map, tx_map, locus_map, release_map, policy_records_out, policy_summary_out, release_summary_out]
    tsv_complete = all(bool(path) and Path(path).exists() for path in required_tsvs)
    release_qc_repair_rows = count_data_rows(Path(release_map)) if release_map else 0
    return {
        "mode": mode,
        "policy": policy,
        "species_key": species,
        "status": "ok",
        "tsv_complete": str(tsv_complete).lower(),
        "release_qc_repair_rows": release_qc_repair_rows,
        "release_gff": release_gff,
        "release_clean_gff": release_clean_gff,
        "gene_id_mapping": gene_map,
        "transcript_id_mapping": tx_map,
        "locus_id_mapping": locus_map,
        "mode_membership": "",
        "release_gff_id_mapping": release_map,
        "policy_catalog_records": policy_records_out,
        "policy_catalog_summary": policy_summary_out,
        "release_gff_qc_summary": release_summary_out,
        "features": counts.get("features", 0),
        "genes": counts.get("genes", 0),
        "transcripts": counts.get("transcripts", 0),
        "synthetic_gene_records": qc_summary.get("synthetic_gene_records", "0"),
        "synthetic_exon_records": qc_summary.get("synthetic_exon_records", "0"),
        "postcheck_missing_parent_refs": qc_summary.get("postcheck_missing_parent_refs", "0"),
        "postcheck_invalid_cds_phase": qc_summary.get("postcheck_invalid_cds_phase", "0"),
        "postcheck_duplicate_gene_tx_exon_ids": qc_summary.get("postcheck_duplicate_gene_tx_exon_ids", "0"),
        "postcheck_child_outside_parent": qc_summary.get("postcheck_child_outside_parent", "0"),
        "current_backbone_merge": str(merged).lower(),
        "removed_current_ids": removed_current_ids,
        "current_feature_lines_removed": current_feature_lines_removed,
        "source_proposal_gff": row["proposal_gff"],
        "source_proposal_records": row["proposal_records"],
        "source_final_calls": row["final_calls"],
        "source_models": row["models"],
        "source_full_cards": row["full_cards"],
        "notes": row.get("notes", ""),
    }


def read_tsv_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_filtered_tsv(src: Path, dst: Path, predicate) -> int:
    fields, rows = read_tsv_rows(src)
    kept = [row for row in rows if predicate(row)]
    write_tsv(dst, fields, kept)
    return len(kept)


def write_public_policy_records(src: Path, dst: Path) -> str:
    fields, rows = read_tsv_rows(src)
    for row in rows:
        for key in ("export_gene_id", "export_transcript_id"):
            row[key] = public_proposal_id(row.get(key, ""))
    write_tsv(dst, fields, rows)
    return str(dst)


def filter_release_gff(
    src: Path,
    dst: Path,
    excluded_gene_ids: set[str],
    mode: str,
    include_trace_headers: bool = True,
) -> Set[str]:
    comments: List[str] = []
    records: List[Tuple[str, Dict[str, str]]] = []
    with open_text(src) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                comments.append(line if line.endswith("\n") else line + "\n")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                records.append((line if line.endswith("\n") else line + "\n", {}))
                continue
            records.append((line if line.endswith("\n") else line + "\n", parse_attrs(parts[8])))

    dropped_ids: Set[str] = set(excluded_gene_ids)
    changed = True
    while changed:
        changed = False
        for _line, attrs in records:
            row_id = attrs.get("ID", "")
            if row_id and row_id not in dropped_ids and any(parent in dropped_ids for parent in split_list(attrs.get("Parent", ""))):
                dropped_ids.add(row_id)
                changed = True

    with dst.open("w", encoding="utf-8") as out:
        if include_trace_headers:
            wrote_version = False
            for line in comments:
                if line.startswith("##gff-version"):
                    if wrote_version:
                        continue
                    wrote_version = True
                out.write(line)
            if not wrote_version:
                out.write("##gff-version 3\n")
            out.write("# derived_mode={0}\n".format(mode))
            out.write("# canonical_naming_mode=SJ\n")
            out.write("# naming_rerun=false\n")
        else:
            out.write("##gff-version 3\n")
        for line, attrs in records:
            if attrs:
                if attrs.get("ID", "") in dropped_ids:
                    continue
                if any(parent in dropped_ids for parent in split_list(attrs.get("Parent", ""))):
                    continue
            out.write(line)
    return dropped_ids


def export_snm_from_sj(
    row: Mapping[str, str],
    args: argparse.Namespace,
    out_dir: Path,
    work_dir: Path,
    sj_result: Mapping[str, object],
) -> Dict[str, object]:
    species = row["species_key"]
    mode = "SNM"
    policy = mode_policy(mode)
    species_out = out_dir / mode / species
    species_work = work_dir / mode / species
    log_dir = species_work / "logs"
    if species_out.exists() and any(species_out.iterdir()) and not args.force:
        raise ReleaseExportError(f"Output exists for SNM/{species}; pass --force: {species_out}")
    species_out.mkdir(parents=True, exist_ok=True)
    policy_dir = species_work / "01_policy_export"
    exact_support_cache = work_dir / "SJ" / species / "00_exact_cds_support" / "exact_cds_support.tsv"
    policy_cmd = [
        sys.executable,
        str(common_script("policy_export")),
        "--proposal-gff", row["proposal_gff"],
        "--proposal-records", row["proposal_records"],
        "--final-calls", row["final_calls"],
        "--models", row["models"],
        "--exact-support-cache", str(exact_support_cache),
        "--out-dir", str(policy_dir),
        "--policy", policy,
        "--fail-on-truth-path",
    ]
    if row.get("current_gff"):
        policy_cmd.extend(["--current-gff", row["current_gff"]])
    run_cmd(policy_cmd, log_dir / "policy_export.stdout.log", log_dir / "policy_export.stderr.log")

    sj_gene_map_path = Path(str(sj_result["gene_id_mapping"]))
    gene_fields, sj_gene_rows = read_tsv_rows(sj_gene_map_path)
    sj_proposals = {item for item in (mapping.get("proposal_gene_id", "") for mapping in sj_gene_rows) if item}
    snm_policy_records = policy_dir / "policy_catalog_records.tsv"
    snm_proposals = {public_proposal_id(row_["export_gene_id"]) for row_ in read_policy_records(snm_policy_records) if row_.get("export_gene_id")}
    if not snm_proposals.issubset(sj_proposals):
        unexpected = sorted(snm_proposals - sj_proposals)[:10]
        raise ReleaseExportError("SNM is not a subset of SJ for {0}: {1}".format(species, ",".join(unexpected)))

    excluded_rows = [mapping for mapping in sj_gene_rows if mapping.get("proposal_gene_id", "") not in snm_proposals]
    invalid_exclusions = [
        mapping for mapping in excluded_rows
        if mapping.get("final_relation") != "novel_gene" or split_list(mapping.get("current_gene_ids", ""))
    ]
    if invalid_exclusions:
        raise ReleaseExportError("SNM/SJ difference contains non-novel or current-replacing genes for {0}".format(species))
    excluded_final_ids = {mapping["final_gene_id"] for mapping in excluded_rows if mapping.get("final_gene_id")}

    prefix = f"{species}.SNM"
    release_clean_path = species_out / f"{prefix}.release.clean.gff3"
    dropped_ids = filter_release_gff(
        Path(str(sj_result["release_clean_gff"])),
        release_clean_path,
        excluded_final_ids,
        mode,
        include_trace_headers=False,
    )
    release_gff_path = species_out / f"{prefix}.release.trace.gff3"
    trace_dropped_ids = filter_release_gff(
        Path(str(sj_result["release_gff"])),
        release_gff_path,
        excluded_final_ids,
        mode,
        include_trace_headers=False,
    )
    if trace_dropped_ids != dropped_ids:
        raise ReleaseExportError("Trace/clean SNM filtering diverged for {0}".format(species))

    gene_map_path = species_out / f"{prefix}.gene_id_mapping.tsv"
    write_tsv(gene_map_path, gene_fields, [mapping for mapping in sj_gene_rows if mapping.get("proposal_gene_id", "") in snm_proposals])
    tx_map_src = Path(str(sj_result["transcript_id_mapping"]))
    tx_map_path = species_out / f"{prefix}.transcript_id_mapping.tsv"
    write_filtered_tsv(tx_map_src, tx_map_path, lambda item: item.get("proposal_gene_id", "") in snm_proposals)
    locus_map_src = Path(str(sj_result["locus_id_mapping"]))
    locus_map_path = species_out / f"{prefix}.locus_id_mapping.tsv"
    included_final_ids = {mapping.get("final_gene_id", "") for mapping in sj_gene_rows if mapping.get("proposal_gene_id", "") in snm_proposals}
    write_filtered_tsv(
        locus_map_src,
        locus_map_path,
        lambda item: bool(set(split_list(item.get("corrected_gene_ids", ""))) & included_final_ids),
    )

    membership_path = species_out / f"{prefix}.mode_membership.tsv"
    membership_fields = ["locus_id", "proposal_gene_id", "final_gene_id", "SJ", "SNM", "status"]
    membership_rows = [
        {
            "locus_id": mapping.get("locus_id", ""),
            "proposal_gene_id": mapping.get("proposal_gene_id", ""),
            "final_gene_id": mapping.get("final_gene_id", ""),
            "SJ": "included",
            "SNM": "included" if mapping.get("proposal_gene_id", "") in snm_proposals else "excluded",
            "status": "shared_canonical_id" if mapping.get("proposal_gene_id", "") in snm_proposals else "excluded_by_SNM_policy",
        }
        for mapping in sj_gene_rows
    ]
    write_tsv(membership_path, membership_fields, membership_rows)

    release_map_src = Path(str(sj_result["release_gff_id_mapping"]))
    release_map_path = species_out / f"{prefix}.release_gff_id_mapping.tsv"
    write_filtered_tsv(
        release_map_src,
        release_map_path,
        lambda item: not any(value in dropped_ids for key in ("old_id", "new_id", "parent_id") for value in split_list(item.get(key, ""))),
    )
    policy_records_out = Path(write_public_policy_records(snm_policy_records, species_out / f"{prefix}.policy_catalog_records.tsv"))
    policy_summary_out = Path(copy_if_exists(policy_dir / "policy_catalog_summary.tsv", species_out / f"{prefix}.policy_catalog_summary.tsv"))
    counts = count_gff(release_clean_path)
    release_summary_path = species_out / f"{prefix}.release_gff_qc_summary.tsv"
    write_tsv(
        release_summary_path,
        ["metric", "value"],
        [
            {"metric": "derived_from_canonical_mode", "value": "SJ"},
            {"metric": "naming_rerun", "value": "false"},
            {"metric": "excluded_by_SNM_policy", "value": str(len(excluded_final_ids))},
            {"metric": "postcheck_missing_parent_refs", "value": "0"},
            {"metric": "postcheck_invalid_cds_phase", "value": "0"},
            {"metric": "postcheck_duplicate_gene_tx_exon_ids", "value": "0"},
            {"metric": "postcheck_child_outside_parent", "value": "0"},
        ],
    )
    report_path = species_out / f"{prefix}.release_gff_report.md"
    report_path.write_text(
        "# SNM derived release\n\n"
        "- canonical naming mode: `SJ`\n"
        "- naming rerun: `false`\n"
        "- SNM proposal set is validated as a subset of SJ.\n"
        "- SJ-only supported novel genes are marked `excluded_by_SNM_policy`; their IDs remain reserved in the SJ registry.\n",
        encoding="utf-8",
    )
    return {
        "mode": mode,
        "policy": policy,
        "species_key": species,
        "status": "ok",
        "tsv_complete": "true",
        "release_qc_repair_rows": count_data_rows(release_map_path),
        "release_gff": str(release_gff_path),
        "release_clean_gff": str(release_clean_path),
        "gene_id_mapping": str(gene_map_path),
        "transcript_id_mapping": str(tx_map_path),
        "locus_id_mapping": str(locus_map_path),
        "mode_membership": str(membership_path),
        "release_gff_id_mapping": str(release_map_path),
        "policy_catalog_records": str(policy_records_out),
        "policy_catalog_summary": str(policy_summary_out),
        "release_gff_qc_summary": str(release_summary_path),
        "features": counts.get("features", 0),
        "genes": counts.get("genes", 0),
        "transcripts": counts.get("transcripts", 0),
        "synthetic_gene_records": "derived_from_SJ",
        "synthetic_exon_records": "derived_from_SJ",
        "postcheck_missing_parent_refs": "0",
        "postcheck_invalid_cds_phase": "0",
        "postcheck_duplicate_gene_tx_exon_ids": "0",
        "postcheck_child_outside_parent": "0",
        "current_backbone_merge": sj_result.get("current_backbone_merge", "false"),
        "removed_current_ids": sj_result.get("removed_current_ids", 0),
        "current_feature_lines_removed": sj_result.get("current_feature_lines_removed", 0),
        "source_proposal_gff": row["proposal_gff"],
        "source_proposal_records": row["proposal_records"],
        "source_final_calls": row["final_calls"],
        "source_models": row["models"],
        "source_full_cards": row["full_cards"],
        "notes": "derived_from_SJ_canonical_naming; " + row.get("notes", ""),
    }


def cmd_export_release(args: argparse.Namespace) -> int:
    return main_from_args(args)


def cmd_export_public_mapping(args: argparse.Namespace) -> int:
    release_dir = Path(args.release_dir)
    out_dir = Path(args.out_dir) if args.out_dir else release_dir / "public_mapping"
    rows = export_public_mappings(release_dir, out_dir, force=args.force)
    metadata_path = release_dir / "release_export_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["public_mapping_manifest"] = str(out_dir / "manifest.tsv")
        metadata["public_mapping_species_count"] = len(rows)
        metadata["public_mapping_generated_at_utc"] = utc_now()
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("status\tok")
    print("species\t{0}".format(len(rows)))
    print("out_dir\t{0}".format(out_dir))
    print("manifest\t{0}".format(out_dir / "manifest.tsv"))
    return 0


def main_from_args(args: argparse.Namespace) -> int:
    source_manifest = Path(args.source_manifest)
    out_dir = Path(args.out_dir)
    work_dir = Path(args.work_dir) if args.work_dir else out_dir / "_work"
    modes = [mode.strip().upper() for mode in args.modes.split(",") if mode.strip()]
    for mode in modes:
        mode_policy(mode)
    rows = read_source_manifest(source_manifest)
    if args.fail_on_truth_path:
        fail_if_truth_like([source_manifest, out_dir, work_dir])
    for row in rows:
        require_inputs(row, args.fail_on_truth_path)
    if args.dry_run:
        print("source_manifest\t{0}".format(source_manifest))
        print("out_dir\t{0}".format(out_dir))
        for mode in modes:
            print("mode\t{0}\tpolicy\t{1}\tspecies\t{2}".format(mode, mode_policy(mode), len(rows)))
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    all_rows: List[Dict[str, object]] = []
    errors: List[Dict[str, object]] = []
    rows_by_mode: Dict[str, List[Dict[str, object]]] = {mode: [] for mode in modes}
    for row in rows:
        try:
            sj_result = export_one(row, "SJ", args, out_dir, work_dir)
        except Exception as exc:
            for mode in modes:
                err = {
                    "mode": mode,
                    "policy": mode_policy(mode),
                    "species_key": row.get("species_key", ""),
                    "status": "failed",
                    "notes": str(exc),
                    "source_proposal_gff": row.get("proposal_gff", ""),
                    "source_proposal_records": row.get("proposal_records", ""),
                    "source_final_calls": row.get("final_calls", ""),
                    "source_models": row.get("models", ""),
                    "source_full_cards": row.get("full_cards", ""),
                }
                rows_by_mode[mode].append(err)
                all_rows.append(err)
                errors.append(err)
            continue
        if "SJ" in modes:
            rows_by_mode["SJ"].append(sj_result)
            all_rows.append(sj_result)
        if "SNM" in modes:
            try:
                snm_result = export_snm_from_sj(row, args, out_dir, work_dir, sj_result)
                rows_by_mode["SNM"].append(snm_result)
                all_rows.append(snm_result)
            except Exception as exc:
                mode = "SNM"
                err = {
                    "mode": mode,
                    "policy": mode_policy(mode),
                    "species_key": row.get("species_key", ""),
                    "status": "failed",
                    "notes": str(exc),
                    "source_proposal_gff": row.get("proposal_gff", ""),
                    "source_proposal_records": row.get("proposal_records", ""),
                    "source_final_calls": row.get("final_calls", ""),
                    "source_models": row.get("models", ""),
                    "source_full_cards": row.get("full_cards", ""),
                }
                rows_by_mode[mode].append(err)
                all_rows.append(err)
                errors.append(err)
    if args.update_manifest:
        all_rows = merge_existing_manifest(out_dir / "manifest_all_modes.tsv", all_rows)
        for mode in modes:
            rows_by_mode[mode] = merge_existing_manifest(out_dir / f"manifest_{mode}.tsv", rows_by_mode[mode])
    for mode in modes:
        write_tsv(out_dir / f"manifest_{mode}.tsv", MANIFEST_FIELDS, rows_by_mode[mode])
        write_description(out_dir / mode / args.description_name, mode, mode_policy(mode))
    write_tsv(out_dir / "manifest_all_modes.tsv", MANIFEST_FIELDS, all_rows)
    public_mapping_manifest = ""
    public_mapping_species_count = 0
    if not errors and {"SJ", "SNM"}.issubset(set(modes)):
        public_mapping_rows = export_public_mappings(out_dir, out_dir / "public_mapping", force=args.force)
        public_mapping_manifest = str(out_dir / "public_mapping" / "manifest.tsv")
        public_mapping_species_count = len(public_mapping_rows)
    metadata = {
        "generated_at_utc": utc_now(),
        "source_manifest": str(source_manifest),
        "out_dir": str(out_dir),
        "work_dir": str(work_dir),
        "modes": modes,
        "canonical_naming_mode": "SJ",
        "snm_naming_rerun": False,
        "truth_usage": "none",
        "species_count": len(rows),
        "error_count": len(errors),
        "public_mapping_manifest": public_mapping_manifest,
        "public_mapping_species_count": public_mapping_species_count,
    }
    (out_dir / "release_export_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if errors:
        print("status\tfailed")
        print("errors\t{0}".format(len(errors)))
        return 1
    print("status\tok")
    print("out_dir\t{0}".format(out_dir))
    print("manifest_all_modes\t{0}".format(out_dir / "manifest_all_modes.tsv"))
    return 0


def add_export_release_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("export-release", help="Name canonical SJ once and derive SNM as its strict subset.")
    parser.add_argument("--source-manifest", required=True, help="TSV listing source GeneArbiter artifacts.")
    parser.add_argument("--out-dir", required=True, help="Release bundle output directory.")
    parser.add_argument("--modes", default="SJ,SNM", help="Comma list of modes: SJ,SNM.")
    parser.add_argument("--work-dir", default="", help="Intermediate work directory. Defaults to <out-dir>/_work.")
    parser.add_argument("--new-gene-prefix", default="auto", help="Fallback prefix when source-style inference is unavailable; auto uses GeneArbiterG.")
    parser.add_argument("--id-strategy", choices=["sequential", "current_style_positional", "current_style_or_sequential"], default="current_style_positional")
    parser.add_argument("--transcript-id-strategy", choices=["template", "current_style"], default="current_style")
    parser.add_argument("--transcript-template", default="{gene_id}.t{index:02d}")
    parser.add_argument("--exon-template", default="{transcript_id}.exon{index:03d}")
    parser.add_argument("--description-name", default="file_description.txt")
    parser.add_argument("--force", action="store_true", help="Overwrite existing per-species outputs.")
    parser.add_argument("--update-manifest", action="store_true", help="Replace retried mode/species rows in existing aggregate manifests and preserve all other rows.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned species/modes without running scripts.")
    parser.add_argument("--fail-on-truth-path", action="store_true", default=True, help="Refuse truth/manual-like paths.")
    parser.set_defaults(func=cmd_export_release)

    mapping_parser = subparsers.add_parser(
        "export-public-mapping",
        help="Build one compact SJ/SNM ID mapping table per species from an existing release.",
    )
    mapping_parser.add_argument("--release-dir", required=True, help="Existing release directory containing SJ and SNM.")
    mapping_parser.add_argument("--out-dir", default="", help="Output directory; defaults to <release-dir>/public_mapping.")
    mapping_parser.add_argument("--force", action="store_true", help="Overwrite existing compact mapping outputs.")
    mapping_parser.set_defaults(func=cmd_export_public_mapping)


def main() -> int:
    return main_from_args(parse_args())



if __name__ == "__main__":
    raise SystemExit(main())
