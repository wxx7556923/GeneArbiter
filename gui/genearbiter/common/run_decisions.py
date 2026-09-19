#!/usr/bin/env python3
# script_id_md5: 3b9b37e8a13e678aa585d2ab90382092
# created: 2026-06-24
# modified: 2026-07-10
# owner: project
# status: project_code
# purpose: 对 AI annotation card 运行受约束的 annotation-set 裁决。
# inputs: 命令行参数指定的 AI annotation card JSONL。
# outputs: 命令行参数指定的 decisions JSONL、summary TSV、prompt/raw response 文件。
# notes: 每个 locus 输出一个 listed annotation set；risk_flagged 只是检查优先级，不放弃选择。

"""Run constrained AI annotation-set decisions.

The input card is already the AI-facing biological summary. The model must
select one listed set_id. It may mark the selection as risk_flagged when
mutually exclusive candidate sets cannot be discriminated by the card evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple


PROFILES = {
    "local_rule": {
        "kind": "local",
        "label": "genearbiter_local_rule_v1",
        "model_env": "",
        "default_model": "genearbiter-local-rule-v1",
        "api_key_env": "",
    },
    "deepseek_flash": {
        "kind": "deepseek",
        "label": "deepseek_flash_non_thinking",
        "model_env": "DEEPSEEK_FLASH_MODEL",
        "default_model": "deepseek-chat",
        "api_key_env": "DEEPSEEK_API_KEY",
    },
    "deepseek_v4_flash": {
        "kind": "deepseek",
        "label": "deepseek_v4_flash_non_thinking",
        "model_env": "DEEPSEEK_V4_FLASH_MODEL",
        "default_model": "deepseek-v4-flash",
        "api_key_env": "DEEPSEEK_API_KEY",
    },
    "deepseek_flash_thinking": {
        "kind": "deepseek",
        "label": "deepseek_flash_thinking",
        "model_env": "DEEPSEEK_FLASH_THINKING_MODEL",
        "default_model": "deepseek-reasoner",
        "api_key_env": "DEEPSEEK_API_KEY",
    },
    "deepseek_pro": {
        "kind": "deepseek",
        "label": "deepseek_pro",
        "model_env": "DEEPSEEK_PRO_MODEL",
        "default_model": "deepseek-chat",
        "api_key_env": "DEEPSEEK_API_KEY",
    },
}

DECISIONS = {"select_annotation_set", "risk_flagged"}
CONFIDENCE = {"high", "medium", "low"}
REQUIRED_FIELDS = [
    "card_id",
    "decision",
    "selected_set_id",
    "confidence",
    "reason_tags",
    "risk_tags",
    "evidence_summary",
    "requires_validation",
]
OPTIONAL_FIELDS = ["normalization_notes"]
FORBIDDEN_EXTRA_FIELDS = {
    "coordinate",
    "coordinates",
    "new_coordinate",
    "new_coordinates",
    "exon",
    "exons",
    "cds",
    "CDS",
    "intron",
    "introns",
    "transcript",
    "transcripts",
    "gene",
    "genes",
    "gff",
    "gff3",
    "start",
    "end",
    "source",
    "tool",
    "model_id",
    "model_ids",
}

SUMMARY_FIELDS = [
    "profile",
    "model",
    "card_id",
    "api_status",
    "valid",
    "validation_errors",
    "decision",
    "selected_set_id",
    "confidence",
    "reason_tags",
    "risk_tags",
    "set_count",
    "component_count",
    "topology_class",
    "main_conflict",
    "discriminator_status",
    "recommended_default",
    "recommended_set_id",
]

SYSTEM_PROMPT = """你是基因注释 evidence card 裁决器。你只根据用户提供的 AI annotation card 做 annotation-set 选择。

总规则：
- 只返回一个 JSON object，不要 markdown。
- 只能选择 card.sets 中已经列出的 set_id，不能混合多个 set。
- 不能生成坐标、exon、CDS、intron、transcript、gene 或 GFF 记录。
- 默认必须做 select_annotation_set；只有 card 明确显示候选 set 互相排斥，且生物学证据无法区分时，才允许 risk_flagged。
- 即使输出 risk_flagged，也必须同时给出一个 listed selected_set_id，表示 provisional best annotation set；risk_flagged 只是人工复核/导出风险标签，不是放弃裁决。

判断顺序：
1. 先读 decision_focus、set_comparison 和 set_differences，确定本 locus 的主要冲突。
2. 再读 component_matrix、components 和每个 set 的 component_decision_summary，判断差异来自 split/merge、component absence、extra component、isoform/topology 差异还是 coding risk。
3. 对每个 set，综合比较 coding QC、phase consistency、splice/junction evidence、short-read support、long-read support、homolog protein support、component support counts 和 explicit risk tags。
4. 优先选择同时具备较好 coding 完整性、转录本/剪接支持、long-read 支持、protein 支持，并覆盖强证据 component 且没有 unsupported extra component 的 set。
5. coding disruption、frameshift、premature stop、phase inconsistency、unsupported junction、conflicting junction、unsupported extra component、missing strong component、unsupported split/merge risk 是重要负面证据。
6. candidate_priority 不是禁止选择；但 low_priority_* 表示该 set 的额外 component、缺失 component 或 split/merge topology 缺少直接支持，只有其它 set 证据更差或该结构有直接 splice/long-read/protein 支持时才优先选择。
7. RNA、protein 或 long-read 缺失不是 absence 的证明，可能只是未表达、数据不覆盖或库不完整，不能单独作为否定证据。
8. 多个候选生成来源一致只代表候选一致性，不等于生物学正确性。
9. 如果没有足够证据区分，但规则不允许 risk_flagged，也必须选择一个 listed set，并使用 low confidence 和 requires_validation=true。

risk_flagged 使用限制：
- 只有 set_comparison.recommended_default 为 risk_if_unresolved_conflict，且 discriminator_status 显示证据冲突或无法区分时，才允许 risk_flagged。
- 如果 recommended_default 是 select_best_listed_set，必须选择一个 listed set_id。
- risk_flagged 时 selected_set_id 也必须是 listed set_id；不要留空。
"""

USER_TEMPLATE = """请为这个 AI annotation card 返回受约束决策。

必须返回 JSON object，格式如下：
{
  "card_id": "必须逐字复制 card_id",
  "decision": "select_annotation_set | risk_flagged",
  "selected_set_id": "必须是 listed set_id；risk_flagged 时也必须给 provisional best set",
  "confidence": "high | medium | low",
  "reason_tags": ["英文短标签"],
  "risk_tags": ["英文短标签"],
  "evidence_summary": "一句话说明依据，只能引用 card 内信息，不要给出新坐标",
  "requires_validation": true
}

AI annotation card:
__CARD_JSON__
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-jsonl", required=True, type=Path, help="AI annotation card JSONL")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory for AI decisions")
    parser.add_argument("--profile", action="append", choices=sorted(PROFILES), help="AI profile to run. Repeatable.")
    parser.add_argument("--locus-id", action="append", default=[], help="Restrict to selected card_id; repeatable")
    parser.add_argument("--max-cards", type=int, default=0, help="0 means all selected cards")
    parser.add_argument("--force", action="store_true", help="Allow overwriting existing profile output files")
    parser.add_argument("--resume", action="store_true", help="Skip cards already present in decisions.jsonl")
    parser.add_argument("--dry-run-prompts", action="store_true", help="Write prompts and metadata without API calls")
    parser.add_argument("--base-url", default=os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com/chat/completions"))
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--max-tokens", type=int, default=2048)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--sleep-sec", type=float, default=0.2)
    parser.add_argument("--max-attempts", type=int, default=2)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def write_jsonl(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def append_jsonl(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")


def write_tsv(path: Path, rows: Sequence[Dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_endpoint(base_url: str) -> str:
    value = base_url.rstrip("/")
    if value.endswith("/chat/completions"):
        return value
    if value.endswith("/v1"):
        return value + "/chat/completions"
    return value + "/chat/completions"


def profile_model(profile: str) -> str:
    cfg = PROFILES[profile]
    env_name = cfg["model_env"]
    if env_name and os.environ.get(env_name, "").strip():
        return os.environ[env_name].strip()
    return cfg["default_model"]


def card_id(card: Dict[str, Any]) -> str:
    return str(card.get("card_id", ""))


def candidate_sets(card: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [row for row in card.get("sets") or [] if row.get("set_id")]


def set_ids(card: Dict[str, Any]) -> List[str]:
    return [str(row.get("set_id", "")) for row in candidate_sets(card)]


def select_cards(cards: Sequence[Dict[str, Any]], locus_ids: Sequence[str], max_cards: int) -> List[Dict[str, Any]]:
    wanted = set(locus_ids)
    out = [card for card in cards if not wanted or card_id(card) in wanted]
    if max_cards > 0:
        out = out[:max_cards]
    return out


def prompt_for_card(card: Dict[str, Any]) -> str:
    return USER_TEMPLATE.replace("__CARD_JSON__", json.dumps(card, ensure_ascii=False, sort_keys=True))


def response_to_json(text: str) -> Dict[str, Any]:
    stripped = text.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```(?:json)?\s*", "", stripped)
        stripped = re.sub(r"\s*```$", "", stripped)
    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", stripped, re.S)
        if not match:
            raise
        return json.loads(match.group(0))


def call_chat_completion(endpoint: str, api_key: str, model: str, system_prompt: str, user_prompt: str, temperature: float, max_tokens: int, timeout: int) -> Dict[str, Any]:
    payload = {
        "model": model,
        "messages": [{"role": "system", "content": system_prompt}, {"role": "user", "content": user_prompt}],
        "max_tokens": max_tokens,
        "response_format": {"type": "json_object"},
    }
    if model != "deepseek-reasoner":
        payload["temperature"] = temperature
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(endpoint, data=data, headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def comparison(card: Dict[str, Any]) -> Dict[str, Any]:
    return card.get("set_comparison") or {}


def risk_allowed_by_card(card: Dict[str, Any]) -> bool:
    cmp = comparison(card)
    return (
        str(cmp.get("recommended_default", "")) == "risk_if_unresolved_conflict"
        and str(cmp.get("discriminator_status", "")) == "conflicting"
        and str(cmp.get("main_conflict", "")) != "none"
    )


def confidence_from_comparison(card: Dict[str, Any], decision: str) -> str:
    if decision == "risk_flagged":
        return "low"
    status = str(comparison(card).get("discriminator_status", ""))
    if status == "clear":
        return "high"
    if status == "weak":
        return "medium"
    return "low"


def listed_or_first_set(card: Dict[str, Any], proposed: str) -> str:
    ids = set_ids(card)
    if proposed in ids:
        return proposed
    return ids[0] if ids else ""


def selected_set(card: Dict[str, Any], selected_set_id: str) -> Dict[str, Any]:
    for row in candidate_sets(card):
        if str(row.get("set_id", "")) == selected_set_id:
            return row
    return {}


def local_rule_decision(card: Dict[str, Any]) -> Dict[str, Any]:
    cmp = comparison(card)
    if risk_allowed_by_card(card):
        selected_id = listed_or_first_set(card, str(cmp.get("recommended_set_id", "")))
        return {
            "card_id": card_id(card),
            "decision": "risk_flagged",
            "selected_set_id": selected_id,
            "confidence": "low",
            "reason_tags": ["unresolved_mutually_exclusive_candidate_sets"],
            "risk_tags": ["risk_flagged_by_card_comparison"],
            "evidence_summary": "Card marks mutually exclusive candidate sets without a fully discriminating biological evidence signal; selected_set_id is the provisional best listed set and requires manual review.",
            "requires_validation": True,
        }

    selected_id = listed_or_first_set(card, str(cmp.get("recommended_set_id", "")))
    selected = selected_set(card, selected_id)
    reason_tags = [str(selected.get("interpretation", "selected_listed_annotation_set"))]
    if selected.get("component_fit"):
        reason_tags.append(str(selected.get("component_fit")))
    return {
        "card_id": card_id(card),
        "decision": "select_annotation_set",
        "selected_set_id": selected_id,
        "confidence": confidence_from_comparison(card, "select_annotation_set"),
        "reason_tags": sorted(set(reason_tags)),
        "risk_tags": [str(tag) for tag in (selected.get("risks") or [])],
        "evidence_summary": "Selected one listed annotation set using the card-level topology, interpreted evidence states, and risk tags.",
        "requires_validation": True,
    }


def normalize_decision(obj: Dict[str, Any], card: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    normalized = dict(obj)
    notes: List[str] = []
    if "requires_validation" not in normalized:
        normalized["requires_validation"] = True
        notes.append("missing_requires_validation_set_true")
    elif normalized.get("requires_validation") is not True:
        normalized["requires_validation"] = True
        notes.append("requires_validation_forced_true")

    selected = str(normalized.get("selected_set_id") or "")
    if str(normalized.get("decision", "")) == "risk_flagged" and not risk_allowed_by_card(card) and selected in set(set_ids(card)):
        normalized["decision"] = "select_annotation_set"
        normalized["confidence"] = "low"
        reason_tags = normalized.get("reason_tags")
        if not isinstance(reason_tags, list):
            reason_tags = []
        reason_tags.append("risk_flagged_normalized_to_select_annotation_set")
        normalized["reason_tags"] = sorted(set(str(tag) for tag in reason_tags if str(tag)))
        risk_tags = normalized.get("risk_tags")
        if not isinstance(risk_tags, list):
            risk_tags = []
        risk_tags.append("model_returned_disallowed_risk_flagged")
        normalized["risk_tags"] = sorted(set(str(tag) for tag in risk_tags if str(tag)))
        notes.append("risk_flagged_not_allowed_normalized_to_select_annotation_set")
    return normalized, notes


def validate_decision(obj: Dict[str, Any], card: Dict[str, Any]) -> List[str]:
    errors: List[str] = []
    allowed_top = set(REQUIRED_FIELDS) | set(OPTIONAL_FIELDS)
    for field in REQUIRED_FIELDS:
        if field not in obj:
            errors.append(f"missing_field:{field}")
    extra = sorted(set(obj) - allowed_top)
    if extra:
        errors.append("extra_top_level_fields:" + ",".join(extra))
    forbidden = sorted(set(obj) & FORBIDDEN_EXTRA_FIELDS)
    if forbidden:
        errors.append("forbidden_fields:" + ",".join(forbidden))

    if str(obj.get("card_id", "")) != card_id(card):
        errors.append("card_id_mismatch")
    decision = str(obj.get("decision", ""))
    if decision not in DECISIONS:
        errors.append("invalid_decision")
    if obj.get("confidence") not in CONFIDENCE:
        errors.append("invalid_confidence")
    if not isinstance(obj.get("reason_tags"), list):
        errors.append("reason_tags_not_list")
    if not isinstance(obj.get("risk_tags"), list):
        errors.append("risk_tags_not_list")
    if obj.get("requires_validation") is not True:
        errors.append("requires_validation_must_be_true")

    selected = str(obj.get("selected_set_id") or "")
    listed = set(set_ids(card))
    if decision == "select_annotation_set":
        if not selected:
            errors.append("select_annotation_set_requires_selected_set_id")
        elif selected not in listed:
            errors.append("selected_set_id_not_listed")
    elif decision == "risk_flagged":
        if not selected:
            errors.append("risk_flagged_requires_selected_set_id")
        elif selected not in listed:
            errors.append("selected_set_id_not_listed")
        if not risk_allowed_by_card(card):
            errors.append("risk_flagged_not_allowed_by_card_comparison")
    return errors


def decision_row(profile: str, model: str, obj: Dict[str, Any], card: Dict[str, Any], status: str, errors: Sequence[str]) -> Dict[str, Any]:
    cmp = comparison(card)
    topology = card.get("topology") or {}
    return {
        "profile": profile,
        "model": model,
        "card_id": card_id(card),
        "api_status": status,
        "valid": str(not errors and status in {"ok", "local_ok"}).lower(),
        "validation_errors": ";".join(errors) if errors else "ok",
        "decision": obj.get("decision", ""),
        "selected_set_id": obj.get("selected_set_id", ""),
        "confidence": obj.get("confidence", ""),
        "reason_tags": ";".join(str(tag) for tag in obj.get("reason_tags", [])) if isinstance(obj.get("reason_tags"), list) else "",
        "risk_tags": ";".join(str(tag) for tag in obj.get("risk_tags", [])) if isinstance(obj.get("risk_tags"), list) else "",
        "set_count": len(candidate_sets(card)),
        "component_count": len(card.get("components") or []),
        "topology_class": topology.get("class", ""),
        "main_conflict": cmp.get("main_conflict", ""),
        "discriminator_status": cmp.get("discriminator_status", ""),
        "recommended_default": cmp.get("recommended_default", ""),
        "recommended_set_id": cmp.get("recommended_set_id", ""),
    }


def decision_record_to_summary_row(profile: str, model: str, rec: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "profile": profile,
        "model": model,
        "card_id": rec.get("card_id", ""),
        "api_status": rec.get("api_status", ""),
        "valid": str(bool(rec.get("valid"))).lower(),
        "validation_errors": ";".join(rec.get("validation_errors") or []) if isinstance(rec.get("validation_errors"), list) else rec.get("validation_errors", ""),
        "decision": rec.get("decision", ""),
        "selected_set_id": rec.get("selected_set_id", ""),
        "confidence": rec.get("confidence", ""),
        "reason_tags": ";".join(rec.get("reason_tags", [])) if isinstance(rec.get("reason_tags"), list) else "",
        "risk_tags": ";".join(rec.get("risk_tags", [])) if isinstance(rec.get("risk_tags"), list) else "",
    }


def existing_card_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    out = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if rec.get("valid") is True and str(rec.get("api_status", "")) in {"ok", "local_ok"}:
                    out.add(str(rec.get("card_id", "")))
    return out


def compact_decisions_for_resume(path: Path) -> None:
    if not path.exists():
        return
    best: Dict[str, Dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            cid = str(rec.get("card_id", ""))
            if not cid:
                continue
            if rec.get("valid") is True and str(rec.get("api_status", "")) in {"ok", "local_ok"}:
                best[cid] = rec
            elif cid not in best:
                best[cid] = rec
    kept = [rec for rec in best.values() if rec.get("valid") is True and str(rec.get("api_status", "")) in {"ok", "local_ok"}]
    write_jsonl(path, sorted(kept, key=lambda row: str(row.get("card_id", ""))))


def run_one_job(job: Dict[str, Any]) -> Dict[str, Any]:
    profile = job["profile"]
    card = job["card"]
    cfg = PROFILES[profile]
    model = job["model"]
    raw_dir: Path = job["raw_dir"]
    prompt_dir: Path = job["prompt_dir"]
    cid = card_id(card)
    user_prompt = prompt_for_card(card)
    if job["write_prompts"]:
        (prompt_dir / f"{cid}.prompt.txt").write_text(user_prompt, encoding="utf-8")

    status = "ok"
    notes: List[str] = []
    obj: Dict[str, Any] = {}
    raw: Dict[str, Any] = {}

    if cfg["kind"] == "local":
        obj = local_rule_decision(card)
        status = "local_ok"
    elif job["dry_run_prompts"]:
        obj = {"card_id": cid}
        status = "dry_run"
        notes.append("api_not_called")
    else:
        for attempt in range(1, max(1, job["max_attempts"]) + 1):
            try:
                raw = call_chat_completion(job["endpoint"], job["api_key"], model, SYSTEM_PROMPT, user_prompt, job["temperature"], job["max_tokens"], job["timeout"])
                (raw_dir / f"{cid}.attempt{attempt}.raw_response.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                obj = response_to_json(raw["choices"][0]["message"]["content"])
                status = "ok"
                break
            except urllib.error.HTTPError as exc:
                status = "fail"
                body = exc.read().decode("utf-8", errors="replace")
                notes.append(f"http_error:{exc.code}:{body[:240]}")
            except Exception as exc:  # noqa: BLE001
                status = "fail"
                notes.append(f"exception:{type(exc).__name__}:{str(exc)[:240]}")
            if attempt < max(1, job["max_attempts"]):
                time.sleep(0.5)
        if status != "ok":
            notes.append(f"failed_after_attempts:{job['max_attempts']}")

    errors: List[str]
    if status in {"ok", "local_ok"}:
        obj, normalization_notes = normalize_decision(obj, card)
        if normalization_notes:
            notes.extend(["normalized:" + note for note in normalization_notes])
            obj["normalization_notes"] = normalization_notes
        errors = validate_decision(obj, card)
    else:
        errors = notes

    model_output_card_id = str(obj.get("card_id", "")) if isinstance(obj, dict) else ""
    record = {
        "profile": profile,
        "profile_label": cfg["label"],
        "model": model,
        "generated_at_utc": utc_now(),
        "api_status": status,
        "api_notes": notes,
        "valid": not errors and status in {"ok", "local_ok"},
        "validation_errors": errors,
        **obj,
    }
    record["card_id"] = cid
    if model_output_card_id and model_output_card_id != cid:
        record["model_output_card_id"] = model_output_card_id
    return {"record": record, "row": decision_row(profile, model, obj, card, status, errors)}


def run_profile(profile: str, cards: Sequence[Dict[str, Any]], args: argparse.Namespace) -> Tuple[Path, List[Dict[str, Any]]]:
    cfg = PROFILES[profile]
    model = profile_model(profile)
    out_dir = args.out_dir / profile
    decisions_path = out_dir / "decisions.jsonl"
    summary_path = out_dir / "decision_summary.tsv"
    raw_dir = out_dir / "raw_api_responses"
    prompt_dir = out_dir / "prompts"

    if out_dir.exists() and decisions_path.exists() and not args.force and not args.resume:
        raise SystemExit(f"Refusing to overwrite existing profile output: {out_dir}\nUse --force or --resume.")
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(exist_ok=True)
    prompt_dir.mkdir(exist_ok=True)
    if args.force and decisions_path.exists() and not args.resume:
        decisions_path.unlink()
    if args.resume:
        compact_decisions_for_resume(decisions_path)

    cards_to_run = list(cards)
    if args.resume:
        done = existing_card_ids(decisions_path)
        cards_to_run = [card for card in cards_to_run if card_id(card) not in done]

    api_key = ""
    endpoint = normalize_endpoint(args.base_url)
    if cfg["kind"] != "local" and not args.dry_run_prompts:
        key_env = cfg["api_key_env"]
        api_key = os.environ.get(key_env, "").strip()
        if not api_key:
            raise SystemExit(f"Missing API key environment variable for {profile}: {key_env}")

    metadata = {
        "generated_at_utc": utc_now(),
        "script": str(Path(__file__).name),
        "profile": profile,
        "profile_label": cfg["label"],
        "model": model,
        "input_jsonl": str(args.input_jsonl),
        "selected_card_count": len(cards),
        "cards_run_this_invocation": len(cards_to_run),
        "dry_run_prompts": args.dry_run_prompts,
        "api_key_written_to_disk": False,
        "incremental_decisions_jsonl": True,
        "endpoint": endpoint if cfg["kind"] != "local" else "",
        "card_schema": "ai_annotation_card_v1",
        "scientific_boundary": "AI always selects one listed annotation set; risk_flagged is only a review/export-risk label for unresolved mutually exclusive candidates. No coordinate generation or deletion decision.",
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out_dir / "system_prompt.txt").write_text(SYSTEM_PROMPT, encoding="utf-8")

    jobs = [
        {
            "profile": profile,
            "model": model,
            "card": card,
            "raw_dir": raw_dir,
            "prompt_dir": prompt_dir,
            "write_prompts": args.dry_run_prompts,
            "dry_run_prompts": args.dry_run_prompts,
            "endpoint": endpoint,
            "api_key": api_key,
            "temperature": args.temperature,
            "max_tokens": args.max_tokens,
            "timeout": args.timeout,
            "max_attempts": args.max_attempts,
        }
        for card in cards_to_run
    ]

    progress_counts: Counter[str] = Counter()

    def print_progress(idx: int, total: int, result: Dict[str, Any]) -> None:
        record = result["record"]
        status = str(record.get("api_status", ""))
        valid = bool(record.get("valid"))
        progress_counts[f"status:{status}"] += 1
        if status in {"ok", "local_ok"}:
            progress_counts["success"] += 1
        if not valid:
            progress_counts["invalid"] += 1
        percent = (idx / total * 100.0) if total else 100.0
        print(
            f"[{profile} {idx}/{total} {percent:5.1f}%] "
            f"card={record.get('card_id')} status={status} valid={str(valid).lower()} "
            f"ok={progress_counts.get('success', 0)} "
            f"fail={progress_counts.get('status:fail', 0)} "
            f"invalid={progress_counts.get('invalid', 0)}",
            file=sys.stderr,
            flush=True,
        )

    def persist_result(result: Dict[str, Any]) -> None:
        append_jsonl(decisions_path, [result["record"]])

    if not jobs:
        print(f"[{profile} 0/0 100.0%] no cards to run; existing decisions will be reused", file=sys.stderr, flush=True)
    if args.workers <= 1:
        for idx, job in enumerate(jobs, 1):
            result = run_one_job(job)
            persist_result(result)
            print_progress(idx, len(jobs), result)
            if cfg["kind"] != "local" and args.sleep_sec > 0:
                time.sleep(args.sleep_sec)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_to_job = {executor.submit(run_one_job, job): job for job in jobs}
            for idx, future in enumerate(as_completed(future_to_job), 1):
                result = future.result()
                persist_result(result)
                print_progress(idx, len(jobs), result)

    all_records = sorted(read_jsonl(decisions_path), key=lambda row: str(row.get("card_id", ""))) if decisions_path.exists() else []
    if not decisions_path.exists():
        decisions_path.write_text("", encoding="utf-8")
    summary_rows = [decision_record_to_summary_row(profile, model, rec) for rec in all_records]

    write_tsv(summary_path, summary_rows, SUMMARY_FIELDS)
    write_profile_report(out_dir / "decision_report.md", profile, model, summary_rows)
    return decisions_path, summary_rows


def write_profile_report(path: Path, profile: str, model: str, rows: Sequence[Dict[str, Any]]) -> None:
    valid_counts = Counter(str(row.get("valid")) for row in rows)
    decision_counts = Counter(str(row.get("decision")) for row in rows)
    confidence_counts = Counter(str(row.get("confidence")) for row in rows)
    api_counts = Counter(str(row.get("api_status")) for row in rows)
    lines = [
        f"# AI annotation-set decision report: {profile}",
        "",
        f"- generated_at_utc: {utc_now()}",
        f"- model: {model}",
        f"- rows: {len(rows)}",
        f"- api_status_counts: {dict(api_counts)}",
        f"- valid_counts: {dict(valid_counts)}",
        f"- decision_counts: {dict(decision_counts)}",
        f"- confidence_counts: {dict(confidence_counts)}",
        "",
        "Boundary: decisions always select one listed annotation set; risk_flagged only marks unresolved mutually exclusive candidates for review/export caution. No coordinates are generated.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def selected_key(row: Dict[str, Any]) -> str:
    decision = str(row.get("decision", ""))
    if decision == "select_annotation_set":
        return "set:" + str(row.get("selected_set_id", ""))
    if decision == "risk_flagged":
        return "risk_flagged:set:" + str(row.get("selected_set_id", ""))
    return "invalid"


def write_cross_profile_comparison(out_dir: Path, profile_rows: Dict[str, List[Dict[str, Any]]]) -> None:
    by_card: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
    provider_summary: List[Dict[str, Any]] = []
    for profile, rows in profile_rows.items():
        valid_rows = [row for row in rows if str(row.get("valid")) == "true"]
        decision_counts = Counter(str(row.get("decision")) for row in valid_rows)
        provider_summary.append(
            {
                "profile": profile,
                "row_count": len(rows),
                "valid_count": len(valid_rows),
                "invalid_count": len(rows) - len(valid_rows),
                "select_annotation_set": decision_counts.get("select_annotation_set", 0),
                "risk_flagged": decision_counts.get("risk_flagged", 0),
            }
        )
        for row in rows:
            by_card[str(row.get("card_id", ""))][profile] = row

    comparison_rows: List[Dict[str, Any]] = []
    profiles = sorted(profile_rows)
    for cid in sorted(by_card):
        selections = []
        decisions = []
        fields: Dict[str, Any] = {"card_id": cid}
        for profile in profiles:
            row = by_card[cid].get(profile, {})
            valid = str(row.get("valid", "false")) == "true" if row else False
            key = selected_key(row) if row else ""
            fields[f"{profile}_valid"] = str(valid).lower()
            fields[f"{profile}_decision"] = row.get("decision", "")
            fields[f"{profile}_selected_key"] = key
            fields[f"{profile}_confidence"] = row.get("confidence", "")
            if valid:
                selections.append(key)
                decisions.append(str(row.get("decision", "")))
        fields["valid_profile_count"] = sum(1 for profile in profiles if fields.get(f"{profile}_valid") == "true")
        fields["all_valid_profiles_same_selected_key"] = str(len(set(selections)) <= 1 if selections else False).lower()
        fields["all_valid_profiles_same_decision"] = str(len(set(decisions)) <= 1 if decisions else False).lower()
        fields["selected_key_values"] = ";".join(selections)
        fields["decision_values"] = ";".join(decisions)
        comparison_rows.append(fields)

    comparison_fields = ["card_id", "valid_profile_count", "all_valid_profiles_same_selected_key", "all_valid_profiles_same_decision", "selected_key_values", "decision_values"]
    for profile in profiles:
        comparison_fields.extend([f"{profile}_valid", f"{profile}_decision", f"{profile}_selected_key", f"{profile}_confidence"])
    write_tsv(out_dir / "ai_decision_cross_profile_comparison.tsv", comparison_rows, comparison_fields)
    write_tsv(out_dir / "ai_decision_provider_summary.tsv", provider_summary, ["profile", "row_count", "valid_count", "invalid_count", "select_annotation_set", "risk_flagged"])


def main() -> int:
    args = parse_args()
    profiles = args.profile or ["local_rule", "deepseek_flash", "deepseek_flash_thinking", "deepseek_pro"]
    cards = select_cards(read_jsonl(args.input_jsonl), args.locus_id, args.max_cards)
    if not cards:
        print("No AI annotation cards selected; writing empty API decisions and continuing.", file=sys.stderr)
    no_sets = [card_id(card) for card in cards if not candidate_sets(card)]
    if no_sets:
        raise SystemExit("Input cards do not contain listed sets: " + ",".join(no_sets[:10]))
    args.out_dir.mkdir(parents=True, exist_ok=True)

    profile_rows: Dict[str, List[Dict[str, Any]]] = {}
    for profile in profiles:
        _, rows = run_profile(profile, cards, args)
        profile_rows[profile] = rows
    if len(profile_rows) > 1:
        write_cross_profile_comparison(args.out_dir, profile_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
