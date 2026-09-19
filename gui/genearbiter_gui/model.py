"""Pure-Python input model shared by the GUI and its unit tests.

The GUI only constructs the public ``genearbiter run-files`` invocation.  All
annotation arbitration, evidence interpretation, ID reconciliation and GFF
release logic remains in the bundled GeneArbiter core package.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ProfileChoice:
    label: str
    profile: str
    model_name: str


# Only profiles already implemented by the bundled DeepSeek client are public.
PROFILE_CHOICES = (
    ProfileChoice("DeepSeek Chat（非思考）", "deepseek_flash", "deepseek-chat"),
    ProfileChoice("DeepSeek Reasoner（思考）", "deepseek_flash_thinking", "deepseek-reasoner"),
)


@dataclass
class CandidateInput:
    name: str
    path: str


@dataclass
class RunRequest:
    current_gff: str
    candidates: list[CandidateInput]
    output_dir: str
    api_key: str
    profile: str = "deepseek_flash"
    workers: int = 4
    genome_fasta: str = ""
    splice_junctions: str = ""
    protein_gff: str = ""
    long_read_support_dir: str = ""


@dataclass
class ValidationResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def source_name_from_path(path: str, fallback_index: int) -> str:
    """Create a stable CLI source label from a selected filename."""

    name = Path(path).name
    for suffix in (".gff3.gz", ".gff.gz", ".gff3", ".gff"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    safe = re.sub(r"[^0-9A-Za-z_.-]+", "_", name).strip("_.-")
    return safe or f"candidate_{fallback_index}"


def _check_file(value: str, label: str, errors: list[str]) -> None:
    if value and not Path(value).expanduser().is_file():
        errors.append(f"{label}不存在或不是文件：{value}")


def validate_request(request: RunRequest) -> ValidationResult:
    """Validate GUI-level invariants before starting the core workflow."""

    result = ValidationResult()
    if not request.current_gff.strip():
        result.errors.append("请选择当前注释 GFF3。")
    else:
        _check_file(request.current_gff, "当前注释 GFF3", result.errors)

    if not 1 <= len(request.candidates) <= 3:
        result.errors.append("候选 GFF3 必须为 1–3 个。")
    names: list[str] = []
    for index, candidate in enumerate(request.candidates, start=1):
        name = candidate.name.strip()
        if not name:
            result.errors.append(f"候选 {index} 缺少来源名。")
        elif "=" in name:
            result.errors.append(f"候选来源名不能包含 '='：{name}")
        names.append(name)
        if not candidate.path.strip():
            result.errors.append(f"候选 {index} 缺少 GFF3 路径。")
        else:
            _check_file(candidate.path, f"候选 {index} GFF3", result.errors)
    if len({name.casefold() for name in names}) != len(names):
        result.errors.append("候选来源名必须唯一。")

    if not request.output_dir.strip():
        result.errors.append("请设置新的输出目录。")
    elif Path(request.output_dir).expanduser().exists():
        result.errors.append(f"输出目录已存在，请换一个新目录：{request.output_dir}")

    if not request.api_key.strip():
        result.errors.append("请输入 DeepSeek API Key。")
    if request.profile not in {choice.profile for choice in PROFILE_CHOICES}:
        result.errors.append(f"不支持的模型配置：{request.profile}")
    if request.workers < 1 or request.workers > 64:
        result.errors.append("并发请求数必须在 1–64 之间。")

    _check_file(request.genome_fasta, "参考基因组 FASTA", result.errors)
    _check_file(request.splice_junctions, "RNA 剪接证据", result.errors)
    _check_file(request.protein_gff, "同源蛋白证据", result.errors)
    if request.long_read_support_dir and not Path(request.long_read_support_dir).expanduser().is_dir():
        result.errors.append(f"长读长证据不存在或不是目录：{request.long_read_support_dir}")

    if not request.genome_fasta:
        result.warnings.append("未提供 FASTA：将跳过候选 CDS 序列和编码质量检查。")
    if not any((request.splice_junctions, request.protein_gff, request.long_read_support_dir)):
        result.warnings.append("未提供 RNA、蛋白或长读长证据；运行仍可继续，但生物学证据有限。")
    return result


def build_worker_arguments(request: RunRequest) -> list[str]:
    """Translate one validated request to public GeneArbiter CLI arguments."""

    def absolute(value: str) -> str:
        return str(Path(value).expanduser().resolve())

    args = ["run-files", "--current", absolute(request.current_gff)]
    for candidate in request.candidates:
        args.extend(["--candidate", f"{candidate.name}={absolute(candidate.path)}"])
    evidence: Iterable[tuple[str, str]] = (
        ("splice_junctions", request.splice_junctions),
        ("protein_gff", request.protein_gff),
        ("long_read_support_dir", request.long_read_support_dir),
    )
    for kind, path in evidence:
        if path:
            args.extend(["--evidence", f"{kind}={absolute(path)}"])
    if request.genome_fasta:
        args.extend(["--genome", absolute(request.genome_fasta)])
    args.extend(
        [
            "--out-dir",
            absolute(request.output_dir),
            "--profile",
            request.profile,
            "--api-key-env",
            "DEEPSEEK_API_KEY",
            "--workers",
            str(request.workers),
        ]
    )
    return args
