import csv
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_hawkins2020_cross_species_validation import build


class HawkinsTransferTests(unittest.TestCase):
    def test_reconstructs_33_guides_with_outcome_separation(self):
        root = Path(__file__).resolve().parents[1]
        raw = root / "data" / "external_validation" / "raw" / "HAWKINS2020_GFP_CROSS_SPECIES"
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            metrics = build(raw, output, raw / "gfp.origmap.tsv", resamples=40)
            self.assertEqual(metrics["n_guides"], 33)
            self.assertAlmostEqual(metrics["spearman_rho"], 0.45254010695187163)
            with (output / "frozen_ecoli_transfer_scores.tsv").open() as handle:
                prediction_fields = csv.DictReader(handle, delimiter="\t").fieldnames
            with (output / "heldout_bsubtilis_outcomes.tsv").open() as handle:
                outcome_fields = csv.DictReader(handle, delimiter="\t").fieldnames
            self.assertNotIn("outcome", prediction_fields)
            self.assertNotIn("prediction", outcome_fields)


if __name__ == "__main__":
    unittest.main()
