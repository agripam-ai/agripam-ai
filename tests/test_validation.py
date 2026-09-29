import math
import unittest

import pandas as pd

from agripam.validation import (ValidationConfig, cluster_bootstrap_intervals, evaluate_ablations,
                                evaluate_predictions, grouped_holdout_manifest,
                                leave_one_group_out_influence, permutation_null_test,
                                permuted_outcome_control, random_effects_spearman_meta,
                                top_k_permutation_test)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame({
            "dataset_id": ["TEST_FIXTURE"] * 6, "organism": ["synthetic"] * 6,
            "genome_id": ["train"] * 3 + ["held"] * 3,
            "guide_id": [f"g{i}" for i in range(6)],
            "prediction": [0.05, 0.15, 0.25, 0.70, 0.80, 0.95],
            "outcome": [0.0, 0.1, 0.2, 0.75, 0.85, 1.0],
        })

    def test_metrics_and_calibration(self):
        self.data["gene_id"] = ["a"] * 3 + ["b"] * 3
        result = evaluate_predictions(self.data, ValidationConfig(positive_threshold=0.5))
        self.assertAlmostEqual(result["spearman_rho"], 1.0)
        self.assertAlmostEqual(result["roc_auc"], 1.0)
        self.assertAlmostEqual(result["average_precision"], 1.0)
        self.assertIn("brier_score", result)
        self.assertAlmostEqual(result["median_within_gene_spearman"], 1.0)

    def test_grouped_holdout_is_disjoint(self):
        split = grouped_holdout_manifest(self.data, ["held"])
        self.assertEqual(split["train_genomes"], ["train"])
        self.assertEqual(split["test_genomes"], ["held"])
        self.assertTrue(split["outcomes_hidden_until_predictions_frozen"])

    def test_no_threshold_does_not_invent_classification_metrics(self):
        result = evaluate_predictions(self.data)
        self.assertNotIn("roc_auc", result)
        self.assertIn("classification_note", result)
        self.assertEqual(len(result["continuous_top_k"]), 3)
        self.assertGreater(result["continuous_top_k"][0]["mean_lift_over_all"], 0)

    def test_permutation_control_is_deterministic(self):
        config = ValidationConfig(positive_threshold=0.5, seed=7)
        self.assertEqual(permuted_outcome_control(self.data, config), permuted_outcome_control(self.data, config))

    def test_cluster_bootstrap_uses_whole_genes(self):
        self.data["gene_id"] = ["a"] * 3 + ["b"] * 3
        result = cluster_bootstrap_intervals(self.data, config=ValidationConfig(resamples=50, seed=3))
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["cluster_unit"], "gene_id")
        self.assertAlmostEqual(result["median_within_cluster_spearman"]["median"], 1.0)

    def test_within_gene_permutation_null_is_deterministic(self):
        self.data["gene_id"] = ["a"] * 3 + ["b"] * 3
        config = ValidationConfig(resamples=50, seed=4)
        first = permutation_null_test(self.data, config=config)
        second = permutation_null_test(self.data, config=config)
        self.assertEqual(first, second)
        self.assertLessEqual(first["pooled_spearman"]["empirical_two_sided_p"], 1.0)

    def test_sensitivity_reports_delta_from_full_model(self):
        frame = self.data.copy()
        frame["prediction_sensitivity_permuted_x"] = list(reversed(frame["prediction"]))
        rows = evaluate_ablations(frame, ["prediction", "prediction_sensitivity_permuted_x"])
        self.assertEqual(rows[0]["delta_spearman_vs_full"], 0.0)
        self.assertLess(rows[1]["delta_spearman_vs_full"], 0.0)

    def test_influence_meta_analysis_and_top_k_null(self):
        frame = pd.DataFrame({
            "dataset_id": ["TEST"] * 12, "organism": ["synthetic"] * 12,
            "genome_id": ["g"] * 12, "guide_id": [f"x{i}" for i in range(12)],
            "gene_id": ["a"] * 4 + ["b"] * 4 + ["c"] * 4,
            "prediction": list(range(12)), "outcome": list(range(12)),
        })
        influence = leave_one_group_out_influence(frame)
        self.assertEqual(influence["status"], "complete")
        self.assertEqual(len(influence["rows"]), 3)
        meta = random_effects_spearman_meta(frame)
        self.assertEqual(meta["status"], "complete")
        self.assertGreater(meta["pooled_spearman"], 0.99)
        top_k = top_k_permutation_test(frame, config=ValidationConfig(resamples=50, seed=5))
        self.assertEqual(top_k["status"], "complete")
        self.assertEqual(len(top_k["rows"]), 3)
        self.assertGreater(top_k["rows"][0]["observed_mean_lift"], 0)


if __name__ == "__main__":
    unittest.main()
