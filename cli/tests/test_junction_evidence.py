#!/usr/bin/env python3

"""Regression tests for the GeneArbiter short-read junction evidence contract."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from genearbiter.evidence_validation import check_junction_coordinates
from genearbiter.short_read import merge_junctions, supported_junctions, write_merged_junctions


class JunctionEvidenceTests(unittest.TestCase):
    def write_star(self, path: Path, unique: int, multi: int = 0) -> None:
        path.write_text(f"chr1\t201\t299\t1\t1\t0\t{unique}\t{multi}\t40\n", encoding="utf-8")

    def test_multisample_output_uses_canonical_comma_list(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first.SJ.out.tab"
            second = root / "second.SJ.out.tab"
            output = root / "merged_splice_junctions.tsv"
            self.write_star(first, 5)
            self.write_star(second, 6, 2)

            merged = merge_junctions([f"sample02={second}", f"sample01={first}"])
            write_merged_junctions(output, merged)

            with output.open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["junction_read_count"], "13")
            self.assertEqual(row["sample_support_count"], "2")
            self.assertEqual(row["supporting_samples"], "sample01,sample02")
            self.assertEqual(check_junction_coordinates(output)[0].split(":", 1)[0], "ok_junction_coordinates")

    def test_support_thresholds_are_both_required(self) -> None:
        enough_both = ("chr1", 201, 299, "+")
        enough_reads_only = ("chr1", 401, 499, "+")
        enough_samples_only = ("chr1", 601, 699, "+")
        merged = {
            enough_both: {"read_count": 8, "samples": {"s1", "s2"}},
            enough_reads_only: {"read_count": 8, "samples": {"s1"}},
            enough_samples_only: {"read_count": 2, "samples": {"s1", "s2"}},
        }

        supported = supported_junctions(merged, min_read_count=3, min_sample_count=2)

        self.assertEqual(set(supported), {enough_both})

    def test_validator_rejects_noncanonical_sample_list_and_bed_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "merged_splice_junctions.tsv"
            path.write_text(
                "chrom\tintron_start_1based\tintron_end_1based\tbed_start_0based\tbed_end_0based\tstrand\t"
                "junction_read_count\tsample_support_count\tsupporting_samples\n"
                "chr1\t201\t299\t200\t300\t+\t5\t2\tsample02;sample01\n",
                encoding="utf-8",
            )

            self.assertTrue(check_junction_coordinates(path)[0].startswith("bad_junction_coordinates:"))

if __name__ == "__main__":
    unittest.main()
