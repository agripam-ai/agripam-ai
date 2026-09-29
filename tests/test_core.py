import unittest
import tempfile
from pathlib import Path

from agripam.core import (analyze_ttc_genome, exact_target_copy_counts, genome_qc,
                          genome_summary, parse_fasta, reverse_complement,
                          scan_pam_targets_with_ambiguity, scan_ttc_targets, score_pam)
from agripam.workflow import _keyword_hits, _pam_gene_context, _write_fasta, run_genome_workflow, tool_inventory, zip_results
from agripam.pam_discovery import infer_pams_from_mobile_references, supported_pams
from agripam.mobile_references import build_tiered_mobile_reference
from agripam.system_assessment import assess_editing_systems
from agripam.editor_targeting import approximate_offtarget_counts
from agripam.system_assessment import assess_introduced_editors


class CoreTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AAGT"), "ACTT")

    def test_scan_forward_ttc(self):
        sequence = "GGGTTC" + "A" * 35 + "CCC"
        hits = scan_ttc_targets(sequence, 35)
        forward = [hit for hit in hits if hit["strand"] == "+"]
        self.assertEqual(len(forward), 1)
        self.assertEqual(forward[0]["start_1based"], 7)
        self.assertEqual(forward[0]["protospacer"], "A" * 35)

    def test_score_lookup(self):
        rows = [{
            "pam3": "TTC",
            "mobile_nearby": "1",
            "predicted_probability": "0.97",
            "predicted_percent": "97",
            "rank_within_context": "1",
        }]
        result = score_pam("ttc", True, rows)
        self.assertEqual(result["rank"], 1)
        self.assertAlmostEqual(result["probability"], 0.97)

    def test_invalid_sequence(self):
        with self.assertRaises(ValueError):
            scan_ttc_targets("TTCNNNN")

    def test_parse_fasta_and_summary(self):
        records = parse_fasta(">contig1\nACGTNN\n>contig2\nGGCC\n")
        summary = genome_summary(records)
        self.assertEqual(summary["contigs"], 2)
        self.assertEqual(summary["total_length_bp"], 10)
        self.assertEqual(summary["ambiguous_bases"], 2)
        self.assertEqual(summary["n50_bp"], 6)
        self.assertEqual(summary["l50_contigs"], 1)

    def test_genome_qc_flags_fragmentation_and_ambiguity(self):
        records = [(f"c{i}", "ACGTN" * 30) for i in range(21)]
        checks = {row["check"]: row for row in genome_qc(records)}
        self.assertEqual(checks["sequence_content"]["status"], "review")
        self.assertEqual(checks["contiguity"]["status"], "review")
        self.assertEqual(checks["N50"]["status"], "review")

    def test_exact_target_copy_counts_both_strands(self):
        records = [("c1", "AAAACGTTTTACGTAAAA")]
        counts = exact_target_copy_counts(records, ["ACGT", "AAAA"])
        self.assertEqual(counts["ACGT"], 2)
        self.assertGreaterEqual(counts["AAAA"], 2)

    def test_bounded_approximate_offtarget_screen(self):
        targets = [
            {"editor": "SpCas9", "protospacer": "A" * 20},
            {"editor": "SpCas9", "protospacer": "A" * 19 + "C"},
            {"editor": "SpCas9", "protospacer": "A" * 18 + "CC"},
            {"editor": "SpCas9", "protospacer": "A" * 17 + "CCC"},
        ]
        result = approximate_offtarget_counts(targets)
        primary = result["SpCas9\t" + "A" * 20]
        self.assertEqual(primary["count"], 2)
        self.assertEqual(primary["nearest_mismatches"], 1)

    def test_effector_compatibility_retains_uncertainty(self):
        targets = [{"editor": "dCas9/CRISPRi", "exact_uniqueness": "unique",
                    "off_target_screen_status": "complete_substitutions_only", "approx_offtargets_le2": 0}]
        rows = assess_introduced_editors(targets, 0, "unavailable", False, True)
        crispri = next(row for row in rows if row["editor"] == "dCas9/CRISPRi")
        self.assertIn("conditional", crispri["compatibility_tier"])
        self.assertIn("not an editing-efficiency", crispri["uncertainty"])

    def test_genome_analysis_ranks_candidates(self):
        records = parse_fasta(">c1\nGGGTTC" + "ACGT" * 9 + "AAA\n")
        rows = [{"pam3": "TTC", "mobile_nearby": "0", "predicted_probability": "0.8127", "predicted_percent": "81.27", "rank_within_context": "1"}]
        summary, candidates = analyze_ttc_genome(records, 35, rows)
        self.assertEqual(summary["ttc_candidates"], 1)
        self.assertEqual(candidates[0]["candidate_id"], "AGP000001")

    def test_workflow_normalizes_long_fasta_contig_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genome.fna"
            records = _write_fasta(">CM181814.1 [taxon withheld] strain NBUS12 chromosome\nACGTACGT\n", path)
            self.assertEqual(records[0][0], "CM181814.1")
            self.assertEqual(path.read_text().splitlines()[0], ">CM181814.1")

    def test_reproducible_workflow_package(self):
        root = Path(__file__).resolve().parents[1]
        sequence = "ATG" + "GCT" * 800 + "TAA"
        work, manifest = run_genome_workflow(">contig\n" + sequence, "test", root, include_heavy=False)
        self.assertEqual(manifest["product"], "RhizoForge-Select")
        self.assertTrue((work / "pipeline_status.tsv").exists())
        if not (work / "proteins.faa").exists():
            status = (work / "pipeline_status.tsv").read_text()
            self.assertIn("coding_sequence_annotation\tunavailable", status)
            self.assertIn("Install Bakta, Prokka or Prodigal", status)
        self.assertGreater(len(zip_results(work)), 0)

    def test_tool_inventory_has_required_components(self):
        rows = tool_inventory(Path(__file__).resolve().parents[1])
        self.assertIn("CRISPRCasTyper", {row["component"] for row in rows})

    def test_pam_target_is_linked_to_overlapping_gene(self):
        sequence = "GGGTTC" + "A" * 35 + "CCC"
        features = [{"contig": "c1", "type": "CDS", "start": 1, "end": 50,
                     "strand": "+", "id": "cds1", "gene": "ftsZ",
                     "product": "cell division protein FtsZ", "attributes": ""}]
        rows = _pam_gene_context([("c1", sequence)], features)
        self.assertEqual(rows[0]["gene_relation"], "within CDS")
        self.assertEqual(rows[0]["gene"], "ftsZ")
        self.assertNotEqual(rows[0]["growth_relevance_screen"], "not flagged")

    def test_strain_specific_pam_inference_ranks_all_triplets(self):
        spacer = "ACGTACGTACGTACGTACGTACGTACGTAC"
        mobile = ">phage_1\nGGGTTT" + spacer + "AACCGG\n"
        spacers = [{"spacer_id": "s1", "array_id": "a1", "spacer": spacer,
                    "trusted": True, "subtype": "I-C"}]
        summary, observations, ranking, logo = infer_pams_from_mobile_references(spacers, mobile)
        self.assertEqual(summary["exact_protospacer_matches"], 1)
        self.assertEqual(observations[0]["pam3_normalized"], "TTT")
        self.assertEqual(len(ranking), 64)
        self.assertEqual(len(logo), 12)
        self.assertEqual(ranking[0]["pam3"], "TTT")
        self.assertEqual(supported_pams(ranking), [])

    def test_generalized_pam_scanner(self):
        sequence = "GGGAAA" + "C" * 35 + "TTT"
        hits = scan_pam_targets_with_ambiguity(sequence, ["AAA"], 35)
        self.assertEqual(hits[0]["pam"], "AAA")
        self.assertEqual(hits[0]["protospacer"], "C" * 35)

    def test_host_aware_mobile_reference_keeps_taxonomic_tiers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "metadata" / "caudoviricetes_refseq_sequence_manifest.tsv"
            phages = root / "data" / "phages" / "processed" / "caudoviricetes_refseq_exact_deduplicated.fna"
            plasmids = root / "data" / "plasmids" / "ncbi_raw_chunks" / "part.fasta"
            manifest.parent.mkdir(parents=True)
            phages.parent.mkdir(parents=True)
            plasmids.parent.mkdir(parents=True)
            manifest.write_text(
                "accession\tvirus_name\thost_name\n"
                "P1\tBacillus velezensis phage\tBacillus velezensis\n"
                "P2\tBacillus subtilis phage\tBacillus subtilis\n"
                "P3\tBacillus cereus phage\tBacillus cereus\n"
            )
            phages.write_text(">P1\nACGT\n>P2\nACGT\n>P3\nACGT\n")
            plasmids.write_text(">PL1 Bacillus velezensis plasmid\nACGT\n")
            fasta, metadata = build_tiered_mobile_reference(root, "Bacillus velezensis")
            self.assertEqual(metadata["tier_2_same_species"], 2)
            self.assertEqual(metadata["tier_3_species_complex"], 1)
            self.assertEqual(metadata["tier_4_same_genus"], 1)
            self.assertIn("evidence_tier=4", fasta)

    def test_system_assessment_reports_completeness_and_unknown_activity(self):
        features = [
            {"contig": "c1", "start": 100, "end": 1000, "strand": "+", "id": "recA",
             "gene": "recA", "product": "DNA recombination protein RecA", "attributes": "gene=recA"},
            {"contig": "c1", "start": 1500, "end": 2200, "strand": "+", "id": "recF",
             "gene": "recF", "product": "DNA repair protein RecF", "attributes": "gene=recF"},
        ]
        rows = assess_editing_systems(features, {"c1": 5000}, 0, "complete", [], "complete")
        hdr = next(row for row in rows if row["system_model"] == "Homology-directed repair route")
        self.assertEqual(hdr["completeness"], "minimum component model satisfied")
        self.assertEqual(hdr["expression_status"], "not evaluated from genome sequence")
        self.assertIn("co-localized", hdr["locus_organization"])

    def test_short_gene_symbol_does_not_match_inside_unrelated_word(self):
        feature = {"contig": "c1", "start": 1, "end": 100, "strand": "+", "id": "x",
                   "product": "tyrosine recombinase XerC", "attributes": "product=tyrosine recombinase XerC"}
        hits = _keyword_hits([feature], {"repair": ("recO",)})
        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
