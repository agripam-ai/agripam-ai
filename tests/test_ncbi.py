import unittest
from agripam.ncbi import validate_assembly_accession

class NcbiTests(unittest.TestCase):
    def test_valid_accession(self):
        self.assertEqual(validate_assembly_accession("gcf_000597985.1"), "GCF_000597985.1")

    def test_invalid_accession(self):
        with self.assertRaises(ValueError):
            validate_assembly_accession("GCF_123; rm -rf")

if __name__ == "__main__":
    unittest.main()
