import unittest
from agripam.community_ranking import rank_communities


class CommunityRankingTests(unittest.TestCase):
    def bank(self):
        return {'traits': {'A': {'P': 5, 'K': 0}, 'B': {'P': 0, 'K': 5},
                           'C': {'P': None, 'K': None}, 'D': {'P': 5, 'K': 5}},
                'meta': {s: {'biosafety_hold': 'yes' if s == 'D' else 'no', 'biosafety_reason': ''}
                         for s in 'ABCD'}, 'functions': {'P': 'phosphorus', 'K': 'potassium'},
                'fungi': {}, 'pairs': {frozenset('AB'): 'mixes'}, 'threshold': 2}

    def test_complementary_pair_wins_and_hold_excluded(self):
        result = rank_communities(self.bank(), ['P', 'K'], max_size=3)
        self.assertEqual(result['rows'][0]['members'], ['A', 'B'])
        self.assertTrue(all('D' not in r['members'] for r in result['rows']))
        self.assertTrue(all(len(r['members']) >= 2 for r in result['rows']))

    def test_inhibition_and_unknown_are_distinct(self):
        bank = self.bank()
        bank['pairs'][frozenset('AB')] = 'inhibits'
        result = rank_communities(bank, ['P', 'K'], max_size=3)
        self.assertTrue(all(not {'A', 'B'}.issubset(r['members']) for r in result['rows']))
        self.assertTrue(all(r['untested_pairs'] > 0 for r in result['rows']))

    def test_unknown_function_not_called_absent(self):
        bank = self.bank()
        bank['traits']['A']['K'] = None
        row = next(r for r in rank_communities(bank, ['K'])['rows'] if r['members'] == ['A', 'C'])
        self.assertEqual(row['not_measured'], 'potassium')
        self.assertEqual(row['below_threshold'], '')

    def test_held_anchor_and_bounded_search(self):
        self.assertFalse(rank_communities(self.bank(), ['P'], anchors=['D'])['rows'])
        result = rank_communities(self.bank(), ['P', 'K'], limit=3)
        self.assertFalse(result['exhaustive'])
        self.assertLessEqual(result['evaluated'], 3)
