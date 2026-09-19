from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

from genearbiter.cli import main


FIXTURES = Path(__file__).resolve().parent / "fixtures"


class RunFilesCliTests(unittest.TestCase):
    def test_file_oriented_cli_exports_concise_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            out_dir = Path(temp_dir) / "run"
            code = main(
                [
                    "run-files",
                    "--current",
                    str(FIXTURES / "tiny_current.gff3"),
                    "--candidate",
                    f"ToolA={FIXTURES / 'tiny_tool_a.gff3'}",
                    "--candidate",
                    f"ToolB={FIXTURES / 'tiny_tool_b.gff3'}",
                    "--out-dir",
                    str(out_dir),
                ]
            )
            self.assertEqual(code, 0)
            public_dir = out_dir / "results"
            clean = public_dir / "annotation.clean.gff3"
            trace = public_dir / "annotation.trace.gff3"
            mapping = public_dir / "id_mapping.tsv.gz"
            self.assertTrue(clean.is_file())
            self.assertTrue(trace.is_file())
            self.assertTrue(mapping.is_file())
            self.assertNotIn("risk", clean.read_text(encoding="utf-8").lower())
            self.assertNotIn("review", clean.read_text(encoding="utf-8").lower())
            with gzip.open(mapping, "rt", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.readline().rstrip("\n").split("\t"),
                    [
                        "feature_type",
                        "mapping_role",
                        "source_name",
                        "source_id",
                        "final_id",
                        "final_parent_id",
                        "change_type",
                        "included_in",
                    ],
                )
            config = json.loads((out_dir / "genearbiter.config.json").read_text(encoding="utf-8"))
            self.assertEqual(config["api"]["profile"], "local_rule")
            self.assertEqual(Path(config["outputs"]["run_dir"]), out_dir / "work")

    def test_existing_output_directory_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            with self.assertRaises(SystemExit):
                main(
                    [
                        "run-files",
                        "--current",
                        str(FIXTURES / "tiny_current.gff3"),
                        "--candidate",
                        f"ToolA={FIXTURES / 'tiny_tool_a.gff3'}",
                        "--out-dir",
                        temp_dir,
                    ]
                )

    def test_unrecognized_evidence_is_not_silently_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            parent = Path(temp_dir)
            unknown = parent / "expression.tsv"
            unknown.write_text("gene_id\ttpm\nG1\t1.0\n", encoding="utf-8")
            out_dir = parent / "run"
            with self.assertRaises(SystemExit):
                main(
                    [
                        "run-files",
                        "--current",
                        str(FIXTURES / "tiny_current.gff3"),
                        "--candidate",
                        f"ToolA={FIXTURES / 'tiny_tool_a.gff3'}",
                        "--evidence",
                        str(unknown),
                        "--out-dir",
                        str(out_dir),
                    ]
                )
            self.assertFalse(out_dir.exists())


if __name__ == "__main__":
    unittest.main()
