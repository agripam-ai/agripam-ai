import unittest
from pathlib import Path

import pandas as pd

from agripam import biobricks as bb
from agripam import parts_io

ROOT = Path(__file__).resolve().parents[1]

XML = """<?xml version="1.0"?><rsbpml><part_list><part><part_id>1</part_id><part_name>BBa_TEST1</part_name>
<part_short_desc>test promoter</part_short_desc><part_type>Regulatory</part_type><release_status>Released</release_status>
<part_results>Works</part_results><part_author>A. Author</part_author>
<sequences><seq_data>ttgacg
gctagc</seq_data></sequences></part></part_list></rsbpml>"""


class FakeResponse:
    def __init__(self, text): self.text = text
    def read(self): return self.text.encode()
    def __enter__(self): return self
    def __exit__(self, *a): return False


class BioBrickTests(unittest.TestCase):
    def test_forbidden_sites_are_found_on_both_strands(self):
        f = bb.check_sequence("p", "AAAGAATTCAAATCTAGAAA", ["RFC10 (BioBrick)"])
        detail = next(x for x in f if x["level"] == "error")["detail"]
        self.assertIn("EcoRI at 4", detail)
        self.assertIn("XbaI", detail)
        ok = bb.check_sequence("p", "AAAAAAAAAT", ["RFC10 (BioBrick)"])
        self.assertEqual(ok[0]["level"], "ok")

    def test_iupac_and_reverse_strand_motifs(self):
        self.assertEqual(len(bb.find_motif("AAGGTCTCAA", "GGTCTC")), 1)          # BsaI forward
        self.assertEqual(bb.find_motif("AAGAGACCAA", "GGTCTC")[0]["strand"], "-")  # reverse complement
        self.assertEqual(len(bb.find_motif("CCAGG", "CCWGG")), 1)                 # degenerate motif
        with self.assertRaises(ValueError):
            bb.find_motif("ACGT", "XX1")

    def test_assembly_scars_and_annotation(self):
        parts = [{"part_id": "a", "sequence": "AAAA", "role": "promoter"}, {"part_id": "b", "sequence": "CCCC", "role": "cds"}]
        built = bb.assemble(parts, scar="TACTAGAG")
        self.assertEqual(built["sequence"], "AAAATACTAGAGCCCC")
        self.assertEqual([(n["start"], n["end"]) for n in built["annotations"]], [(1, 4), (5, 12), (13, 16)])
        self.assertEqual(bb.assemble(parts, scar="")["sequence"], "AAAACCCC")
        self.assertTrue(bb.assemble(parts, add_biobrick_ends=True)["sequence"].startswith(bb.PREFIX))

    def test_homology_arms_and_integration(self):
        genome = [("c1", "".join("ACGT"[(i * 7) % 4] for i in range(400)))]
        arms = bb.homology_arms(genome, "c1", 200, 50)
        self.assertEqual(len(arms["left"]), 50)
        self.assertEqual(arms["left"], genome[0][1][150:200])
        self.assertEqual(arms["right"], genome[0][1][200:250])
        cass = bb.assemble([{"part_id": "a", "sequence": "GATTACAGATTACA", "role": "cds"}])
        fwd = bb.integration_construct(cass, arms, "+")
        rev = bb.integration_construct(cass, arms, "-")
        self.assertEqual(len(fwd["sequence"]), 50 + 14 + 50)
        self.assertEqual(rev["sequence"][50:64], bb.revcomp("GATTACAGATTACA"))
        self.assertEqual(rev["annotations"][1]["strand"], "-")
        with self.assertRaises(ValueError):
            bb.homology_arms(genome, "c1", 10, 50)
        pcr = bb.diagnostic_pcr(arms, 14, outside=100)
        self.assertEqual(pcr["edited_bp"] - pcr["wild_type_bp"], 14)

    def test_sgrna_cassette_validates_spacer(self):
        pr = {"part_id": "pr", "sequence": "TTGACA", "role": "promoter"}
        sc = {"part_id": "sc", "sequence": "GTTTTAGAGC", "role": "scaffold"}
        built = bb.sgrna_cassette(pr, "ACGTACGTACGTACGTACGT", sc, None)
        self.assertEqual(built["sequence"], "TTGACA" + "ACGTACGTACGTACGTACGT" + "GTTTTAGAGC")
        with self.assertRaises(ValueError):
            bb.sgrna_cassette(pr, "ACGT", sc, None)
        # a scar must never sit inside the guide RNA: it is used only before the terminator
        tm = {"part_id": "tm", "sequence": "TTTTTTT", "role": "terminator"}
        with_term = bb.sgrna_cassette(pr, "ACGTACGTACGTACGTACGT", sc, tm, scar="TACTAGAG")
        self.assertEqual(with_term["sequence"], "TTGACA" + "ACGTACGTACGTACGTACGT" + "GTTTTAGAGC" + "TACTAGAG" + "TTTTTTT")
        self.assertEqual(bb.assemble([pr, sc], junction_scars=["AA"])["sequence"], "TTGACAAAGTTTTAGAGC")

    def test_host_range_and_chassis_motifs(self):
        seq = "ACGT" * 30
        f = bb.check_sequence("p", seq, [], chassis="Bacillus subtilis", hosts="E. coli")
        self.assertTrue(any(x["check"] == "host range" and x["level"] == "warning" for x in f))
        f = bb.check_sequence("p", seq, [], chassis="Bacillus subtilis", hosts="E. coli; B. subtilis")
        self.assertTrue(any(x["check"] == "host range" and x["level"] == "ok" for x in f))
        f = bb.check_sequence("p", "AAGCTAAA", [], chassis_motifs={"HindIII-like": "AAGCTT"})
        self.assertFalse(any(x["check"] == "chassis restriction motif" for x in f))
        f = bb.check_sequence("p", "AAAGCTTAAA", [], chassis_motifs={"HindIII": "AAGCTT"})
        self.assertTrue(any(x["check"] == "chassis restriction motif" for x in f))

    def test_registry_parse_and_fetch_with_fake_opener(self):
        part = bb.parse_registry_xml(XML)
        self.assertEqual((part["part_id"], part["role"], part["sequence"]), ("BBa_TEST1", "promoter", "TTGACGGCTAGC"))
        got = bb.fetch_registry_part("BBa_TEST1", opener=lambda req, timeout: FakeResponse(XML))
        self.assertEqual(got["sequence"], "TTGACGGCTAGC")
        with self.assertRaises(ValueError):
            bb.fetch_registry_part("not an id; rm -rf")
        def offline(req, timeout): raise OSError("offline")
        with self.assertRaises(RuntimeError):
            bb.fetch_registry_part("BBa_TEST1", opener=offline)

    def test_genbank_and_fasta(self):
        built = bb.assemble([{"part_id": "a", "name": "A", "sequence": "ACGT" * 20, "role": "promoter"}])
        gb = bb.to_genbank("c1", built["sequence"], built["annotations"])
        self.assertIn("1..80", gb)
        self.assertTrue(gb.rstrip().endswith("//"))
        self.assertTrue(bb.to_fasta("c1", "ACGT" * 40).startswith(">c1 160 bp"))

    def test_tracker_funnel_controls_and_correlation(self):
        cons = [{"construct_id": f"K{i}", "edit_type": "CRISPRi", "target": "g", "parts": "p", "length_bp": 100,
                 "predicted_score": 10 + i} for i in range(6)]
        rows = bb.tracker_rows(cons, "exp1")
        self.assertEqual(sum(r["role_in_experiment"] == "test" for r in rows), 6)
        self.assertEqual(sorted(r["role_in_experiment"] for r in rows if r["role_in_experiment"] != "test"),
                         ["low-ranked", "negative", "parental"])
        for i, r in enumerate(r for r in rows if r["role_in_experiment"] == "test"):
            r.update(stage_assembled="2026-01-01", measured_value=5 + i, normalized_effect=1 + i * 0.5, replicates=3)
        result = bb.analyze_tracker(rows)
        self.assertEqual(result["funnel"]["designed"], 6)
        self.assertEqual(result["funnel"]["assembled"], 6)
        self.assertAlmostEqual(result["spearman_rho"], 1.0)
        self.assertTrue(any("parental control" in w for w in result["warnings"]))

    def test_starter_library_and_workbook_roundtrip(self):
        parts = parts_io.read_parts_tsv(ROOT / "data/biobricks/starter_library.tsv")
        ids = {p["part_id"] for p in parts}
        self.assertIn("BBa_J23100", ids)
        self.assertIn("SCAFFOLD_SpCas9_unverified", ids)
        frame = pd.DataFrame([{"part_id": "X1", "role": "wrongrole", "sequence": "ACGTN"},
                              {"part_id": "EXAMPLE_x", "role": "cds", "sequence": "ACGT"}])
        parts, issues = parts_io.parts_from_frame(frame)
        self.assertEqual(parts, [])
        self.assertTrue(any("other than A, C, G, T" in i["message"] for i in issues))
        data = parts_io.tracker_to_xlsx(bb.tracker_rows([{"construct_id": "K1", "edit_type": "reporter", "target": "t",
                                                          "parts": "p", "length_bp": 5}], "g"))
        rows = parts_io.read_tracker(__import__("io").BytesIO(data))
        self.assertEqual(rows[0]["construct_id"], "K1")


if __name__ == "__main__":
    unittest.main()
