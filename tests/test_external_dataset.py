import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.prepare_gse196911_validation import normalize
from scripts.attach_gse196911_sequences import attach
from scripts.generate_yu2024_frozen_predictions import build_model_matrix


class ExternalDatasetTests(unittest.TestCase):
    def test_gse196911_normalization_declares_sign_and_withholds_predictions(self):
        raw = pd.DataFrame({
            "guides": ["purA_1", "purB_2"], "type": ["targeting", "targeting"],
            "OD02.logFC": [-1.0, -2.0], "OD02.pvalue": [0.1, 0.01],
            "OD06.logFC": [-2.0, -3.0], "OD06.pvalue": [0.01, 0.001],
            "OD1.logFC": [-4.0, -5.0], "OD1.pvalue": [0.001, 0.0001],
        })
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.tsv"
            raw.to_csv(path, sep="\t", index=False)
            normalized = normalize(path)
        self.assertEqual(normalized["outcome"].tolist(), [4.0, 5.0])
        self.assertNotIn("prediction", normalized.columns)
        self.assertTrue(normalized["prediction_status"].str.startswith("withheld").all())

    def test_supplement_sequence_matching_uses_all_outcomes(self):
        geo = pd.DataFrame({"guides": ["purA_1"], "OD02.logFC": [-1.0], "OD06.logFC": [-2.0], "OD1.logFC": [-3.0]})
        supplement = pd.DataFrame({
            "gene_name": ["purA"], "OD02.logFC": [-1.0], "OD06.logFC": [-2.0], "OD1.logFC": [-3.0],
            "sequence": ["A" * 20], "sequence_30nt": ["C" * 30], "seq_60nt": ["G" * 60], "PAM": ["AGG"],
            "genome_pos": [10], "guide_strand": ["+"], "guide_GC_content": [0],
            "distance_start_codon": [10], "distance_start_codon_perc": [1.0],
        })
        joined = attach(geo, supplement)
        self.assertEqual(joined.loc[0, "sequence"], "A" * 20)

    def test_frozen_prediction_matrix_excludes_outcomes(self):
        blinded = pd.DataFrame({"guides": ["g1"], "sequence": ["A" * 20], "OD1.logFC": [-9.0]})
        feature_values = {
            "sequence": ["A" * 20], "sequence_30nt": ["ATCG" * 7 + "AT"],
            "guide_GC_content": [0.0], "distance_start_codon": [1],
            "distance_start_codon_perc": [1.0], "homopolymers": [1],
            "MFE_hybrid_full": [-1.0], "MFE_hybrid_seed": [-1.0],
            "MFE_homodimer_guide": [-1.0], "MFE_monomer_guide": [-1.0],
            "if_promoter": [0], "OD1_edgeR.batch": [-99.0],
        }
        headers = ["guide_GC_content", "sequence_1_A"]
        _, matrix = build_model_matrix(blinded, pd.DataFrame(feature_values), headers)
        self.assertEqual(matrix.columns.tolist(), headers)
        self.assertNotIn("OD1.logFC", matrix.columns)


if __name__ == "__main__":
    unittest.main()
