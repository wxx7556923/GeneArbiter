from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from genearbiter_gui.model import (
    PROFILE_CHOICES,
    CandidateInput,
    RunRequest,
    build_worker_arguments,
    source_name_from_path,
    validate_request,
)


class GuiModelTests(unittest.TestCase):
    def test_public_profiles_only_expose_deepseek(self) -> None:
        profiles = {choice.profile for choice in PROFILE_CHOICES}
        self.assertEqual(profiles, {"deepseek_flash", "deepseek_flash_thinking"})
        self.assertNotIn("local_rule", profiles)

    def test_source_name_is_safe(self) -> None:
        self.assertEqual(source_name_from_path("C:/data/My Tool.gff3", 1), "My_Tool")

    def test_request_builds_public_cli_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            current = root / "current.gff3"
            candidate = root / "candidate.gff3"
            rna = root / "junctions.tsv"
            for path in (current, candidate, rna):
                path.write_text("test\n", encoding="utf-8")
            request = RunRequest(
                current_gff=str(current),
                candidates=[CandidateInput("tool_a", str(candidate))],
                splice_junctions=str(rna),
                output_dir=str(root / "new_run"),
                api_key="secret",
                profile="deepseek_flash",
                workers=3,
            )
            validation = validate_request(request)
            self.assertTrue(validation.ok, validation.errors)
            args = build_worker_arguments(request)
            self.assertEqual(args[0], "run-files")
            self.assertIn("tool_a=" + str(candidate), args)
            self.assertIn("splice_junctions=" + str(rna), args)
            self.assertNotIn("secret", args)
            self.assertEqual(args[args.index("--workers") + 1], "3")

    def test_existing_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            current = root / "current.gff3"
            candidate = root / "candidate.gff3"
            current.write_text("test\n", encoding="utf-8")
            candidate.write_text("test\n", encoding="utf-8")
            request = RunRequest(
                current_gff=str(current),
                candidates=[CandidateInput("tool_a", str(candidate))],
                output_dir=str(root),
                api_key="secret",
            )
            result = validate_request(request)
            self.assertFalse(result.ok)
            self.assertTrue(any("输出目录已存在" in value for value in result.errors))

    def test_release_guide_and_build_metadata_are_public_safe(self) -> None:
        root = Path(__file__).resolve().parents[1]
        guide = (root / "WINDOWS_USER_GUIDE.md").read_text(encoding="utf-8")
        build = (root / "build_windows.ps1").read_text(encoding="utf-8")
        self.assertIn("1–3", guide)
        self.assertNotIn("1–13", guide)
        self.assertIn("WINDOWS_USER_GUIDE.md", build)
        self.assertNotIn("pip freeze", build)


if __name__ == "__main__":
    unittest.main()
