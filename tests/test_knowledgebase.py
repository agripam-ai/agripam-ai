import unittest

import pandas as pd

from agripam.knowledgebase import delete_records, merge_records, with_record_ids


COLUMNS = ["organism", "strain", "editing_route", "source_id", "notes"]


class KnowledgebaseTests(unittest.TestCase):
    def test_first_upload_initializes_empty_database(self):
        incoming = pd.DataFrame([["Bacillus", "B1", "Cas9", "PMID1", "first"]], columns=COLUMNS)
        merged, summary = merge_records(pd.DataFrame(columns=COLUMNS), incoming, COLUMNS)
        self.assertEqual(len(merged), 1)
        self.assertEqual(summary, {"added": 1, "updated": 0, "unchanged": 0})

    def test_later_upload_supplements_and_updates_without_duplicates(self):
        current = pd.DataFrame([["Bacillus", "B1", "Cas9", "PMID1", "old"]], columns=COLUMNS)
        incoming = pd.DataFrame([
            ["Bacillus", "B1", "Cas9", "PMID1", "reviewed"],
            ["Pseudomonas", "P2", "CRISPRi", "PMID2", "new"],
        ], columns=COLUMNS)
        merged, summary = merge_records(current, incoming, COLUMNS)
        self.assertEqual(len(merged), 2)
        self.assertEqual(summary, {"added": 1, "updated": 1, "unchanged": 0})
        self.assertEqual(merged.loc[merged["source_id"] == "PMID1", "notes"].item(), "reviewed")

    def test_delete_selected_records(self):
        current = pd.DataFrame([
            ["Bacillus", "B1", "Cas9", "PMID1", "one"],
            ["Pseudomonas", "P2", "CRISPRi", "PMID2", "two"],
        ], columns=COLUMNS)
        selected = with_record_ids(current).iloc[0]["record_id"]
        remaining = delete_records(current, [selected])
        self.assertEqual(remaining["source_id"].tolist(), ["PMID2"])

    def test_identity_matching_ignores_case_and_surrounding_spaces(self):
        current = pd.DataFrame([["Bacillus", "B1", "Cas9", "PMID1", "old"]], columns=COLUMNS)
        incoming = pd.DataFrame([[" bacillus ", "b1", "cas9", "pmid1", "reviewed"]], columns=COLUMNS)
        merged, summary = merge_records(current, incoming, COLUMNS)
        self.assertEqual(len(merged), 1)
        self.assertEqual(summary["updated"], 1)
        self.assertEqual(merged.iloc[0]["notes"], "reviewed")


if __name__ == "__main__":
    unittest.main()
