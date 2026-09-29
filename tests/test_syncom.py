import unittest
from pathlib import Path

from agripam import syncom

ROOT = Path(__file__).resolve().parents[1]


class SynComTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = syncom.load_bank(ROOT / "data/syncom")

    def test_biosafety_hold_never_ranked_as_chassis(self):
        ranking = syncom.rank_chassis(["B1", "B39", "B56"], ["FI20"], self.bank)
        b1 = next(r for r in ranking if r["strain"] == "B1")
        self.assertFalse(b1["safety_eligible"])
        self.assertEqual(b1["rank"], "")
        self.assertEqual(b1["composite"], 0.0)

    def test_redundant_member_is_preferred_over_sole_provider(self):
        ranking = syncom.rank_chassis(["B1", "B39", "B56"], ["FI20"], self.bank)
        self.assertEqual(ranking[0]["strain"], "B56")
        b39 = next(r for r in ranking if r["strain"] == "B39")
        self.assertIn("chitinase", b39["functions_lost_if_removed"])

    def test_blank_cells_are_not_zero(self):
        self.assertIsNone(self.bank["traits"]["54"]["catalase"])

    def test_incompatible_strains_are_not_assembled(self):
        community = syncom.assemble(self.bank, ["FI20"], size=4)
        self.assertNotIn("B34", community)
        self.assertNotIn("54", community)

    def test_weight_sensitivity_is_deterministic(self):
        a = syncom.weight_sensitivity(["B39", "B56"], ["FI20"], self.bank)
        b = syncom.weight_sensitivity(["B39", "B56"], ["FI20"], self.bank)
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
