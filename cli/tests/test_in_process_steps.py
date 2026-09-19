from __future__ import annotations

import io
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from genearbiter.cli import run_python_script_in_process
from genearbiter.common.run_decisions import run_profile


class InProcessStepTests(unittest.TestCase):
    def test_run_python_script(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            script = Path(temp_dir) / "step.py"
            script.write_text("import sys\nprint('out:' + sys.argv[1])\nprint('err', file=sys.stderr)\n", encoding="utf-8")
            stdout = io.StringIO()
            stderr = io.StringIO()
            code = run_python_script_in_process([sys.executable, str(script), "ok"], stdout, stderr)
            self.assertEqual(code, 0)
            self.assertIn("out:ok", stdout.getvalue())
            self.assertIn("err", stderr.getvalue())

    def test_system_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            script = Path(temp_dir) / "step.py"
            script.write_text("raise SystemExit(7)\n", encoding="utf-8")
            code = run_python_script_in_process([sys.executable, str(script)], io.StringIO(), io.StringIO())
            self.assertEqual(code, 7)

    def test_empty_api_card_set_writes_empty_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            args = SimpleNamespace(
                out_dir=Path(temp_dir),
                force=False,
                resume=True,
                base_url="",
                dry_run_prompts=False,
                input_jsonl=Path(temp_dir) / "empty_cards.jsonl",
                workers=1,
                sleep_sec=0,
                max_attempts=1,
                timeout=1,
                temperature=0.0,
                max_tokens=256,
            )
            decisions, rows = run_profile("local_rule", [], args)
            self.assertTrue(decisions.is_file())
            self.assertEqual(decisions.read_text(encoding="utf-8"), "")
            self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
