from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


COMMON = Path(__file__).resolve().parents[1] / "genearbiter" / "common"
if str(COMMON) not in sys.path:
    sys.path.insert(0, str(COMMON))

from reconcile_ids import (
    build_gene_mapping,
    build_locus_mapping,
    build_tx_mapping,
    collect_current_genes,
    parse_current_style_gene_id,
    rewrite_gff,
)
from export_policy_gff import (
    GeneBlock as PolicyGeneBlock,
    Record as PolicyRecord,
    block_cds_signatures,
    cds_signature,
    include_strict_novel_block,
    split_list,
)
from release_gff import Record, build_release_records, validate_release_records, write_gff
from genearbiter.release_export import (
    PUBLIC_MAPPING_FIELDS,
    build_public_mapping_rows,
    filter_release_gff,
    write_gzip_tsv,
    write_public_policy_records,
)


def record(gene_id: str, locus: str, start: int, end: int, action: str, current: str = "") -> dict[str, str]:
    return {
        "export_gene_id": gene_id,
        "export_transcript_id": gene_id + ".tmp1",
        "locus_id": locus,
        "proposal_action": action,
        "removed_current_gene_ids": current,
        "current_gene_ids": current,
        "seqid": "Chr02",
        "start": str(start),
        "end": str(end),
        "strand": "+",
        "source": "tool",
    }


class ReleaseNamingTests(unittest.TestCase):
    def test_gff3_encoded_multigene_replacement_list_is_decoded(self) -> None:
        self.assertEqual(
            split_list("OsMH_01G0000900%3BOsMH_01G0001000"),
            ["OsMH_01G0000900", "OsMH_01G0001000"],
        )

    @staticmethod
    def novel_policy_block() -> PolicyGeneBlock:
        gene = PolicyRecord(
            "Chr01\tGeneArbiter\tgene\t100\t300\t.\t+\t.\tID=p1;locus_id=L1;proposal_action=add_novel_selected_set",
            ["Chr01", "GeneArbiter", "gene", "100", "300", ".", "+", ".", "ID=p1;locus_id=L1;proposal_action=add_novel_selected_set"],
            {"ID": "p1", "locus_id": "L1", "proposal_action": "add_novel_selected_set"},
            0,
        )
        transcript = PolicyRecord(
            "Chr01\tGeneArbiter\tmRNA\t100\t300\t.\t+\t.\tID=t1;Parent=p1",
            ["Chr01", "GeneArbiter", "mRNA", "100", "300", ".", "+", ".", "ID=t1;Parent=p1"],
            {"ID": "t1", "Parent": "p1"},
            1,
        )
        cds1 = PolicyRecord(
            "Chr01\tGeneArbiter\tCDS\t100\t150\t.\t+\t0\tParent=t1",
            ["Chr01", "GeneArbiter", "CDS", "100", "150", ".", "+", "0", "Parent=t1"],
            {"Parent": "t1"},
            2,
        )
        cds2 = PolicyRecord(
            "Chr01\tGeneArbiter\tCDS\t250\t300\t.\t+\t0\tParent=t1",
            ["Chr01", "GeneArbiter", "CDS", "250", "300", ".", "+", "0", "Parent=t1"],
            {"Parent": "t1"},
            3,
        )
        return PolicyGeneBlock("p1", gene, [gene, transcript, cds1, cds2])

    def test_snm_uses_exact_selected_cds_multitool_support(self) -> None:
        block = self.novel_policy_block()
        signature = block_cds_signatures(block)[0]
        call = {
            "best_evidence_level": "strong",
            "manual_review_unsupported_locus": "false",
            "member_real_sources": "ANNEVO",
        }
        include, reason = include_strict_novel_block(
            block,
            call,
            {},
            {},
            {"L1": {signature: {"ANNEVO", "Tiberius"}}},
            "strict_supported_novel_multitool_only",
        )
        self.assertTrue(include)
        self.assertEqual(reason, "strict_keep_supported_exact_cds_multitool_novel")

    def test_snm_does_not_trust_final_call_member_source_count(self) -> None:
        block = self.novel_policy_block()
        signature = block_cds_signatures(block)[0]
        call = {
            "best_evidence_level": "strong",
            "manual_review_unsupported_locus": "false",
            "member_real_sources": "ANNEVO,Tiberius",
        }
        include, reason = include_strict_novel_block(
            block,
            call,
            {},
            {},
            {"L1": {signature: {"ANNEVO"}}},
            "strict_supported_novel_multitool_only",
        )
        self.assertFalse(include)
        self.assertEqual(reason, "strict_drop_supported_singletool_or_unknown_novel")

    def test_exact_cds_support_distinguishes_phase(self) -> None:
        phase_zero = cds_signature("Chr01", "+", [(100, 150, "0"), (250, 300, "0")])
        phase_one = cds_signature("Chr01", "+", [(100, 150, "1"), (250, 300, "0")])
        self.assertNotEqual(phase_zero, phase_one)

    def test_split_inherits_one_id_and_inserts_position_id(self) -> None:
        records = [
            record("V2PROP_a", "locus_split", 100, 145, "replace_current_with_selected_set", "OsMH_02G0390000"),
            record("V2PROP_b", "locus_split", 155, 199, "replace_current_with_selected_set", "OsMH_02G0390000"),
        ]
        current = {
            "OsMH_02G0389900": {"seqid": "Chr02", "start": 1, "end": 90},
            "OsMH_02G0390000": {"seqid": "Chr02", "start": 100, "end": 199},
            "OsMH_02G0390100": {"seqid": "Chr02", "start": 300, "end": 390},
        }
        mapping, rows, _next = build_gene_mapping(
            records,
            set(),
            "reuse_current_gene_id",
            "NEW_",
            1,
            6,
            "current_style_positional",
            current,
        )
        self.assertEqual(mapping["V2PROP_a"]["final_gene_id"], "OsMH_02G0390000")
        self.assertEqual(mapping["V2PROP_b"]["final_gene_id"], "OsMH_02G0390050")
        self.assertEqual(rows[0]["proposal_gene_id"], "GeneArbiterProposal_a")
        locus = build_locus_mapping(rows)[0]
        self.assertEqual(locus["retained_current_gene_ids"], "OsMH_02G0390000")
        self.assertEqual(locus["added_gene_ids"], "OsMH_02G0390050")

    def test_split_anchor_uses_cds_contribution_not_gene_span(self) -> None:
        records = [
            record("V2PROP_left", "locus_split", 100, 400, "replace_current_with_selected_set", "OsMH_02G0390000"),
            record("V2PROP_right", "locus_split", 401, 900, "replace_current_with_selected_set", "OsMH_02G0390000"),
        ]
        current = {
            "OsMH_02G0389900": {"seqid": "Chr02", "start": 1, "end": 90, "strand": "+", "cds": ()},
            "OsMH_02G0390000": {
                "seqid": "Chr02",
                "start": 100,
                "end": 900,
                "strand": "+",
                "cds": ((120, 200), (500, 800)),
            },
            "OsMH_02G0390100": {"seqid": "Chr02", "start": 1000, "end": 1090, "strand": "+", "cds": ()},
        }
        proposal_cds = {
            "V2PROP_left": ((120, 200),),
            "V2PROP_right": ((500, 800),),
        }
        mapping, _rows, _next = build_gene_mapping(
            records,
            set(),
            "reuse_current_gene_id",
            "NEW_",
            1,
            6,
            "current_style_positional",
            current,
            proposal_cds,
        )
        self.assertEqual(mapping["V2PROP_right"]["final_gene_id"], "OsMH_02G0390000")
        self.assertEqual(mapping["V2PROP_right"]["id_assignment_status"], "deterministic_cds_anchor")

    def test_public_id_style_keeps_version_suffix(self) -> None:
        self.assertEqual(parse_current_style_gene_id("Prupe.1G000100_v2.0.a1"), ("Prupe.1G", 100, 6, "_v2.0.a1"))
        self.assertEqual(parse_current_style_gene_id("GH_scaffold49195_objG0001"), ("GH_scaffold49195_objG", 1, 4, ""))

    def test_parentless_mrna_is_current_gene_anchor(self) -> None:
        content = """##gff-version 3
A01\tcurrent\tmRNA\t100\t300\t.\t+\t.\tID=GH_A01G0001
A01\tcurrent\tCDS\t120\t180\t.\t+\t0\tParent=GH_A01G0001
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "cotton.gff3"
            source.write_text(content, encoding="utf-8")
            genes = collect_current_genes(str(source))
        self.assertIn("GH_A01G0001", genes)
        self.assertEqual(genes["GH_A01G0001"]["cds"], ((120, 180),))

    def test_auto_style_learns_seqid_embedding_without_species_config(self) -> None:
        records = [record("V2PROP_new", "locus_new", 100, 200, "add_novel_selected_set")]
        records[0]["seqid"] = "scaffold0099"
        current = {}
        for scaffold in ("scaffold0001", "scaffold0002"):
            for index in range(1, 7):
                gene_id = "Bna{0}G{1:07d}ZS".format(scaffold, index * 100)
                current[gene_id] = {"seqid": scaffold, "start": index * 1000, "end": index * 1000 + 100, "strand": "+", "cds": ()}
        mapping, _rows, _next = build_gene_mapping(
            records, set(), "reuse_current_gene_id", "auto", 1, 6, "current_style_or_sequential", current
        )
        self.assertEqual(mapping["V2PROP_new"]["final_gene_id"], "Bnascaffold0099G0000100ZS")
        self.assertEqual(mapping["V2PROP_new"]["id_assignment_status"], "deterministic_new_source_style_id")

    def test_auto_style_continues_dominant_static_css_pattern(self) -> None:
        records = [record("V2PROP_new", "locus_new", 500, 600, "add_novel_selected_set")]
        records[0]["seqid"] = "Chr1"
        numbers = [6089, 43312, 36687, 3312]
        current = {
            "CSS{0:07d}".format(number): {"seqid": "Chr1", "start": index * 100, "end": index * 100 + 50, "strand": "+", "cds": ()}
            for index, number in enumerate(numbers, start=1)
        }
        mapping, _rows, _next = build_gene_mapping(
            records, set(), "reuse_current_gene_id", "auto", 1, 6, "current_style_or_sequential", current
        )
        final_id = mapping["V2PROP_new"]["final_gene_id"]
        self.assertRegex(final_id, r"^CSS\d{7}$")
        self.assertNotIn(final_id, current)

    def test_auto_style_does_not_mint_ncbi_loc_ids(self) -> None:
        records = [record("V2PROP_new", "locus_new", 500, 600, "add_novel_selected_set")]
        records[0]["seqid"] = "NC_000001.1"
        current = {
            "gene-LOC101204932": {"seqid": "NC_000001.1", "start": 100, "end": 150, "strand": "+", "cds": ()},
            "gene-LOC101205177": {"seqid": "NC_000001.1", "start": 200, "end": 250, "strand": "+", "cds": ()},
        }
        mapping, _rows, _next = build_gene_mapping(
            records, set(), "reuse_current_gene_id", "auto", 1, 6, "current_style_or_sequential", current
        )
        self.assertEqual(mapping["V2PROP_new"]["final_gene_id"], "GeneArbiterG000001")

    def test_novel_gene_uses_free_number_between_neighbors(self) -> None:
        records = [record("V2PROP_novel", "locus_novel", 200, 240, "add_novel_selected_set")]
        current = {
            "OsMH_02G0068800": {"seqid": "Chr02", "start": 100, "end": 150},
            "OsMH_02G0069000": {"seqid": "Chr02", "start": 300, "end": 350},
        }
        mapping, _rows, _next = build_gene_mapping(
            records, set(), "reuse_current_gene_id", "NEW_", 1, 6, "current_style_positional", current
        )
        self.assertEqual(mapping["V2PROP_novel"]["final_gene_id"], "OsMH_02G0068900")

    def test_current_style_transcript_id(self) -> None:
        records = [record("V2PROP_a", "locus", 100, 145, "replace_current_with_selected_set", "OsMH_02G0390000")]
        gene_mapping = {
            "V2PROP_a": {
                "final_gene_id": "OsMH_02G0390000",
                "final_relation": "one_to_one_replacement",
                "current_gene_ids": "OsMH_02G0390000",
            }
        }
        mapping, _rows = build_tx_mapping(records, gene_mapping, set(), "{gene_id}.t{index:02d}", "current_style")
        self.assertEqual(mapping["V2PROP_a.tmp1"]["final_transcript_id"], "OsMH_02T0390000.1")

    def test_named_proposal_source_is_genearbiter(self) -> None:
        content = "Chr02\ttemporary\tgene\t100\t145\t.\t+\t.\tID=V2PROP_a;proposal_catalog=true\n"
        gene_mapping = {
            "V2PROP_a": {
                "final_gene_id": "OsMH_02G0390000",
                "final_relation": "one_to_one_replacement",
                "reuse_current_gene_id": "true",
                "current_gene_ids": "OsMH_02G0390000",
                "retired_current_gene_ids": "",
            }
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "proposal.gff3"
            target = Path(temp_dir) / "named.gff3"
            source.write_text(content, encoding="utf-8")
            rewrite_gff(str(source), str(target), gene_mapping, {}, "GeneArbiter")
            result = target.read_text(encoding="utf-8")
        self.assertIn("\tGeneArbiter\tgene\t", result)
        self.assertIn("proposal_gene_id=GeneArbiterProposal_a", result)
        self.assertNotIn("V2PROP", result)

    def test_snm_filter_removes_complete_sj_only_subtree(self) -> None:
        content = """##gff-version 3
Chr01\tGeneArbiter\tgene\t1\t100\t.\t+\t.\tID=OsMH_01G0000100
Chr01\tGeneArbiter\tmRNA\t1\t100\t.\t+\t.\tID=OsMH_01T0000100.1;Parent=OsMH_01G0000100
Chr01\tGeneArbiter\texon\t1\t100\t.\t+\t.\tID=e1;Parent=OsMH_01T0000100.1
Chr01\tcurrent\tgene\t200\t300\t.\t+\t.\tID=OsMH_01G0000200
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "sj.gff3"
            target = Path(temp_dir) / "snm.gff3"
            source.write_text(content, encoding="utf-8")
            dropped = filter_release_gff(source, target, {"OsMH_01G0000100"}, "SNM")
            result = target.read_text(encoding="utf-8")
        self.assertIn("OsMH_01T0000100.1", dropped)
        self.assertNotIn("OsMH_01G0000100", result)
        self.assertIn("OsMH_01G0000200", result)
        self.assertIn("naming_rerun=false", result)

    def test_snm_clean_filter_keeps_only_standard_header(self) -> None:
        content = """##gff-version 3
# generated_by=GeneArbiter release_gff
Chr01\tGeneArbiter\tgene\t1\t100\t.\t+\t.\tID=g1
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "sj.clean.gff3"
            target = Path(temp_dir) / "snm.clean.gff3"
            source.write_text(content, encoding="utf-8")
            filter_release_gff(source, target, set(), "SNM", include_trace_headers=False)
            result = target.read_text(encoding="utf-8").splitlines()

        self.assertEqual(result, ["##gff-version 3", "Chr01\tGeneArbiter\tgene\t1\t100\t.\t+\t.\tID=g1"])

    def test_public_policy_records_hide_legacy_version_label(self) -> None:
        content = "locus_id\texport_gene_id\texport_transcript_id\nL1\tV2PROP_L1_g01\tV2PROP_L1_g01.t01\n"
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "internal.tsv"
            target = Path(temp_dir) / "public.tsv"
            source.write_text(content, encoding="utf-8")
            write_public_policy_records(source, target)
            result = target.read_text(encoding="utf-8")
        self.assertIn("GeneArbiterProposal_L1_g01", result)
        self.assertNotIn("V2PROP", result)

    def test_public_mapping_is_one_compact_edge_list_for_both_modes(self) -> None:
        gene_rows = [
            {
                "final_relation": "split_replacement",
                "final_gene_id": "G1",
                "current_gene_ids": "OLD1;OLD2",
                "retired_current_gene_ids": "OLD2",
                "source": "ANNEVO",
                "source_gene_id": "source_g1",
            },
            {
                "final_relation": "novel_gene",
                "final_gene_id": "G2",
                "current_gene_ids": "",
                "retired_current_gene_ids": "",
                "source": "Helixer",
                "source_gene_id": "source_g2",
            },
        ]
        transcript_rows = [
            {
                "final_gene_id": "G1",
                "final_transcript_id": "T1",
                "source": "ANNEVO",
                "source_transcript_id": "source_t1",
            },
            {
                "final_gene_id": "G2",
                "final_transcript_id": "T2",
                "source": "Helixer",
                "source_transcript_id": "source_t2",
            },
        ]
        rows = build_public_mapping_rows(gene_rows, transcript_rows, {"G1"})

        self.assertEqual(len(rows), 6)
        self.assertEqual(
            {row["mapping_role"] for row in rows},
            {"current_lineage", "selected_origin"},
        )
        self.assertEqual(
            {row["source_id"] for row in rows if row["mapping_role"] == "current_lineage"},
            {"OLD1", "OLD2"},
        )
        self.assertEqual(
            {row["included_in"] for row in rows if row["final_id"] in {"G1", "T1"}},
            {"SJ,SNM"},
        )
        self.assertEqual(
            {row["included_in"] for row in rows if row["final_id"] in {"G2", "T2"}},
            {"SJ"},
        )
        self.assertTrue(all(row["change_type"] == "" for row in rows if row["feature_type"] == "transcript"))

    def test_public_mapping_gzip_has_only_public_columns(self) -> None:
        rows = [
            {
                "feature_type": "gene",
                "mapping_role": "selected_origin",
                "source_name": "ANNEVO",
                "source_id": "source_g1",
                "final_id": "G1",
                "final_parent_id": "",
                "change_type": "replacement",
                "included_in": "SJ,SNM",
            }
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            first = Path(temp_dir) / "first.tsv.gz"
            second = Path(temp_dir) / "second.tsv.gz"
            write_gzip_tsv(first, PUBLIC_MAPPING_FIELDS, rows)
            write_gzip_tsv(second, PUBLIC_MAPPING_FIELDS, rows)

            import gzip

            content = gzip.open(first, "rt", encoding="utf-8").read().splitlines()
            self.assertEqual(first.read_bytes(), second.read_bytes())

        self.assertEqual(content[0].split("\t"), PUBLIC_MAPPING_FIELDS)
        self.assertEqual(len(content), 2)
        self.assertNotIn("proposal", content[1])

    def test_release_qc_emits_multi_parent_exon_once(self) -> None:
        records = [
            Record("Chr01", "current", "gene", 1, 100, ".", "+", ".", {"ID": "g1"}, 1),
            Record("Chr01", "current", "mRNA", 1, 100, ".", "+", ".", {"ID": "t1", "Parent": "g1"}, 2),
            Record("Chr01", "current", "mRNA", 1, 100, ".", "+", ".", {"ID": "t2", "Parent": "g1"}, 3),
            Record("Chr01", "current", "exon", 1, 100, ".", "+", ".", {"ID": "e1", "Parent": "t1,t2"}, 4),
        ]
        released, _mapping, counts = build_release_records(
            records, "{gene_id}.t{index:02d}", "{transcript_id}.exon{index:03d}", "GeneArbiter"
        )
        self.assertEqual(sum(record.attrs.get("ID") == "e1" for record in released), 1)
        self.assertEqual([record.attrs.get("ID") for record in released], ["g1", "t1", "t2", "e1"])
        self.assertEqual(counts.get("postcheck_duplicate_gene_tx_exon_ids", 0), 0)

    def test_release_qc_repairs_unambiguous_child_strand_typo(self) -> None:
        records = [
            Record("Chr01", "current", "gene", 1, 100, ".", "+", ".", {"ID": "g1"}, 1),
            Record("Chr01", "current", "mRNA", 1, 100, ".", "+", ".", {"ID": "t1", "Parent": "g1"}, 2),
            Record("Chr01", "current", "five_prime_UTR", 1, 10, ".", "-", ".", {"ID": "u1", "Parent": "t1"}, 3),
        ]
        released, _mapping, counts = build_release_records(
            records, "{gene_id}.t{index:02d}", "{transcript_id}.exon{index:03d}", "GeneArbiter"
        )
        utr = next(record for record in released if record.attrs.get("ID") == "u1")
        self.assertEqual(utr.strand, "+")
        self.assertEqual(utr.attrs["original_release_strand"], "-")
        self.assertEqual(counts.get("postcheck_child_outside_parent", 0), 0)

    def test_release_qc_separates_inherited_source_hierarchy_warning(self) -> None:
        records = [
            Record("chr1", "Titan", "gene", 1, 100, ".", "+", ".", {"ID": "g1"}, 1),
            Record("chr1", "Titan", "lncRNA", 1, 100, ".", "-", ".", {"ID": "t1", "Parent": "g1"}, 2),
            Record("chr1", "GeneArbiter", "gene", 200, 300, ".", "+", ".", {"ID": "g2"}, 3),
            Record("chr1", "GeneArbiter", "mRNA", 200, 300, ".", "-", ".", {"ID": "t2", "Parent": "g2"}, 4),
        ]
        counts = validate_release_records(records, "GeneArbiter")
        self.assertEqual(counts["postcheck_inherited_source_child_outside_parent"], 1)
        self.assertEqual(counts["postcheck_child_outside_parent"], 1)

    def test_release_clean_gff_uses_feature_attribute_whitelist(self) -> None:
        records = [
            Record(
                "Chr01",
                "GeneArbiter",
                "gene",
                1,
                100,
                ".",
                "+",
                ".",
                {
                    "ID": "g1",
                    "Name": "gene one",
                    "biotype": "protein_coding",
                    "locus_id": "locus_1",
                    "selected_source": "Helixer",
                    "original_release_span": "2-99",
                },
                1,
            ),
            Record(
                "Chr01",
                "GeneArbiter",
                "mRNA",
                1,
                100,
                ".",
                "+",
                ".",
                {"ID": "t1", "Parent": "g1", "Name": "tx one", "transcript_biotype": "protein_coding", "source_model_id": "tool:t1"},
                2,
            ),
            Record(
                "Chr01",
                "GeneArbiter",
                "CDS",
                1,
                100,
                ".",
                "+",
                "0",
                {"ID": "cds1", "Parent": "t1", "protein_id": "p1", "original_release_feature_id": "old_cds"},
                3,
            ),
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "release.clean.gff3"
            write_gff(output, records, Path("input.gff3"), clean=True)
            lines = output.read_text(encoding="utf-8").splitlines()

        self.assertEqual(lines[0], "##gff-version 3")
        self.assertEqual(sum(line.startswith("#") for line in lines), 1)
        self.assertEqual(lines[1].split("\t")[8], "ID=g1;Name=gene%20one;biotype=protein_coding")
        self.assertEqual(lines[2].split("\t")[8], "ID=t1;Parent=g1;Name=tx%20one;biotype=protein_coding")
        self.assertEqual(lines[3].split("\t")[8], "ID=cds1;Parent=t1;protein_id=p1")

    def test_trace_release_gff_keeps_only_compact_trace_attributes(self) -> None:
        record = Record(
            "Chr01",
            "Helixer",
            "gene",
            1,
            100,
            ".",
            "+",
            ".",
            {
                "ID": "g1",
                "Name": "g1",
                "locus_id": "locus_1",
                "selected_source": "Helixer",
                "source_gene_id": "helixer_g1",
                "final_id_relation": "novel_gene",
                "risk_tags": "set_contains_not_assessed_models",
            },
            1,
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "release.gff3"
            write_gff(output, [record], Path("input.gff3"), clean=False)
            content = output.read_text(encoding="utf-8")

        self.assertEqual(sum(line.startswith("#") for line in content.splitlines()), 1)
        self.assertIn("\tGeneArbiter\tgene\t", content)
        self.assertIn("change_type=novel", content)
        self.assertIn("origin_source=Helixer", content)
        self.assertIn("origin_gene_id=helixer_g1", content)
        self.assertIn("review_recommended=true", content)
        self.assertIn("review_reason=junction_not_assessed", content)
        self.assertNotIn("Name=g1", content)
        self.assertNotIn("locus_id", content)
        self.assertNotIn("selected_source", content)
        self.assertNotIn("risk_tags", content)


if __name__ == "__main__":
    unittest.main()
