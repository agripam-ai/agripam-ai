import json
from pathlib import Path
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/external_validation/processed/CRISPRHAL_TEV_RELEASED_HOLDOUT"


class CrisprHALReleasedHoldoutTests(unittest.TestCase):
    def test_outcomes_are_absent_from_blinded_inputs_and_frozen_predictions(self):
        blinded = pd.read_csv(DATA / "prediction_inputs_blinded.tsv", sep="\t")
        frozen = pd.read_csv(DATA / "frozen_predictions.tsv", sep="\t")
        outcomes = pd.read_csv(DATA / "heldout_outcomes.tsv", sep="\t")
        self.assertNotIn("outcome", blinded.columns)
        self.assertNotIn("outcome", frozen.columns)
        self.assertIn("outcome", outcomes.columns)
        self.assertEqual(set(blinded.guide_id), set(frozen.guide_id))
        self.assertEqual(set(blinded.guide_id), set(outcomes.guide_id))

    def test_reported_holdout_metrics_and_scope(self):
        metrics = json.loads((DATA / "validation_metrics.json").read_text())
        self.assertEqual(metrics["exact_sequence_overlap_train_test"], 0)
        self.assertGreater(metrics["primary"]["spearman_rho"], metrics["gc_baseline"]["spearman_rho"])
        self.assertIn("not external-study", metrics["benchmark_role"])
        self.assertEqual(metrics["uncertainty_and_null"]["resamples"], 2000)
        similarity = metrics["sequence_similarity_audit"]
        self.assertEqual(similarity["distant_subset_n"], 5031)
        self.assertGreater(similarity["distant_subset_metrics"]["spearman_rho"], 0.6)


if __name__ == "__main__":
    unittest.main()
