import unittest
from pathlib import Path

import pandas as pd

try:
    import openpyxl  # noqa: F401
    HAVE_XLSX = True
except ImportError:  # the app degrades gracefully without it
    HAVE_XLSX = False

from agripam import syncom_io

ROOT = Path(__file__).resolve().parents[1]


def sheets_from(strains, fungi=None, pairs=None):
    out = {"Strains": pd.DataFrame(strains)}
    if fungi:
        out["Fungi_compatibility"] = pd.DataFrame(fungi)
    if pairs:
        out["Bacteria_compatibility"] = pd.DataFrame(pairs)
    return out


class SynComIoTests(unittest.TestCase):
    def test_validation_catches_bad_values(self):
        sheets = sheets_from(
            [{"strain": "A", "biosafety_hold": "no", "chitinase": 9},
             {"strain": "B", "biosafety_hold": "no", "chitinase": "abc"}],
            [{"strain": "A", "fungus": "F1", "symbol": "?"}])
        bank, issues = syncom_io.build_bank(sheets)
        self.assertIsNone(bank)
        text = " ".join(i["message"] for i in issues)
        self.assertIn("outside 0-5", text)
        self.assertIn("not a number", text)
        self.assertIn("must be one of", text)

    def test_blank_biosafety_warns_but_hold_blocks_ranking(self):
        sheets = sheets_from(
            [{"strain": "A", "biosafety_hold": "", "chitinase": 4, "glucanase": 1},
             {"strain": "B", "biosafety_hold": "yes", "chitinase": 1, "glucanase": 4},
             {"strain": "C", "biosafety_hold": "no", "chitinase": 3, "glucanase": 3}],
            [{"strain": s, "fungus": "F1", "symbol": "++"} for s in "ABC"])
        bank, issues = syncom_io.build_bank(sheets)
        self.assertTrue(any(i["level"] == "warning" and "blank" in i["message"] for i in issues))
        result = syncom_io.analyze(bank, ["A", "B", "C"], ["F1"])
        b = next(r for r in result["ranking"] if r["strain"] == "B")
        self.assertFalse(b["safety_eligible"])
        self.assertEqual(b["rank"], "")

    def test_custom_function_columns_and_threshold(self):
        sheets = sheets_from([{"strain": "A", "biosafety_hold": "no", "my_trait": 3},
                              {"strain": "B", "biosafety_hold": "no", "my_trait": 1}])
        bank, _ = syncom_io.build_bank(sheets)
        low = syncom_io.analyze(bank, ["A", "B"], [], threshold=2)
        high = syncom_io.analyze(bank, ["A", "B"], [], threshold=4)
        self.assertEqual(low["missing"], [])
        self.assertEqual(high["missing"], ["my trait"])

    def test_pairwise_inhibition_lowers_fit(self):
        sheets = sheets_from(
            [{"strain": "A", "biosafety_hold": "no", "x": 3}, {"strain": "B", "biosafety_hold": "no", "x": 3}],
            pairs=[{"strain_a": "A", "strain_b": "B", "outcome": "inhibits"}])
        bank, _ = syncom_io.build_bank(sheets)
        ranking = syncom_io.analyze(bank, ["A", "B"], [])["ranking"]
        self.assertEqual(ranking[0]["fit_score"], 0.0)

    def _fungal_sheets(self):
        strains = [
            {"strain": "B1", "kind": "bacterium", "biosafety_hold": "no", "growth": 4, "biocontrol": 0, "editing_precedent": "species"},
            {"strain": "B2", "kind": "bacterium", "biosafety_hold": "no", "growth": 0, "biocontrol": 3, "editing_precedent": "none"},
            {"strain": "F1", "kind": "fungus", "biosafety_hold": "no", "growth": 0, "biocontrol": 5, "editing_precedent": "genus"},
            {"strain": "F2", "kind": "fungus", "biosafety_hold": "no", "growth": 0, "biocontrol": 4, "editing_precedent": "none"}]
        sheets = sheets_from(strains,
                             fungi=[{"strain": "B1", "fungus": "F1", "symbol": "++"}, {"strain": "B2", "fungus": "F1", "symbol": "-"},
                                    {"strain": "B1", "fungus": "F2", "symbol": "++"}])
        sheets["Fungus_compatibility"] = pd.DataFrame([{"fungus_a": "F1", "fungus_b": "F2", "symbol": "-"}])
        return sheets

    def test_fungi_are_scored_members(self):
        bank, issues = syncom_io.build_bank(self._fungal_sheets())
        self.assertEqual(issues, [])
        self.assertEqual(bank["kind"]["F1"], "fungus")
        result = syncom_io.analyze(bank, fungi=["F1"])
        self.assertEqual(result["members"][0], "F1")            # the fungus is a member
        self.assertIn("B1", result["members"])                  # the compatible bacterium is added
        self.assertNotIn("B2", result["members"])               # B2 is inhibited by F1
        self.assertEqual(result["kinds"]["F1"], "fungus")
        self.assertNotIn("biocontrol", result["missing"])      # the fungus supplies biocontrol
        ranking = {r["strain"]: r for r in result["ranking"]}
        self.assertEqual(ranking["F1"]["kind"], "fungus")

    def test_fungal_functions_count_in_dispensability(self):
        bank, _ = syncom_io.build_bank(self._fungal_sheets())
        result = syncom_io.analyze(bank, ["B1", "B2"], ["F1"])
        ranking = {r["strain"]: r for r in result["ranking"]}
        # B2's biocontrol (3) is covered by the fungus (5), so nothing is lost if B2's activity is lost
        self.assertEqual(ranking["B2"]["dispensability_score"], 100.0)

    def test_incompatible_fungi_raise_a_warning(self):
        bank, _ = syncom_io.build_bank(self._fungal_sheets())
        result = syncom_io.analyze(bank, ["B1"], ["F1", "F2"])
        self.assertTrue(any("F1 and F2 inhibit each other" in w for w in result["warnings"]))
        self.assertTrue(any(r["fit_notes"] for r in result["ranking"] if r["strain"] == "F1"))

    def test_legacy_partner_fungi_are_not_members(self):
        bank, _ = syncom_io.build_bank(sheets_from(
            [{"strain": "A", "biosafety_hold": "no", "x": 3}], [{"strain": "A", "fungus": "EXT", "symbol": "++"}]))
        result = syncom_io.analyze(bank, ["A"], ["EXT"])
        self.assertEqual(result["members"], ["A"])
        self.assertEqual(result["kinds"], {"A": "bacterium"})

    @unittest.skipUnless(HAVE_XLSX, "openpyxl not installed")
    def test_template_and_example_workbooks(self):
        template = syncom_io.read_workbook(ROOT / "data/syncom/SynCom_input_template.xlsx")
        bank, issues = syncom_io.build_bank(template)
        self.assertIsNone(bank)  # only EXAMPLE rows: must be replaced first
        example = syncom_io.read_workbook(ROOT / "data/syncom/SynCom_example_filled.xlsx")
        bank, issues = syncom_io.build_bank(example)
        self.assertEqual(issues, [])
        result = syncom_io.analyze(bank, fungi=["F1"])
        self.assertIn("F1", result["members"])      # the fungal anchor is a scored member
        self.assertNotIn("S4", result["members"])   # inhibited by the fungus
        self.assertNotIn("S5", result["members"])   # biosafety hold
        self.assertEqual(result["missing"], [])
        self.assertEqual(syncom_io.analyze(bank, fungi=[])["ranking"][0]["strain"], "S6")   # bacteria only
        self.assertGreater(len(syncom_io.results_to_xlsx(result, issues)), 1000)


if __name__ == "__main__":
    unittest.main()
