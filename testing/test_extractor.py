# coding: utf-8
import unittest as ut

import tRNAsorter.extractor as extractor

class ExtractorTestCase(ut.TestCase):
    def setUp(self):
        self.extractor = extractor.Extractor()

    def test_pass_pair_check(self):
        result = self.extractor.pair_check("CACCCCTGATAAGGGTG")
        self.assertTrue(result)

    def test_fail_pair_check(self):
        result = self.extractor.pair_check("CAGGCCTGATAAGGGTG")
        self.assertFalse(result)

    def test_pass_get_anticodon(self):
        result = self.extractor.get_anticodon("CACCCCTGATAAGGGTG")
        self.assertEqual(result, "GAT")

    def test_fail_no_matching_precursor_get_anticodon(self):
        result = self.extractor.get_anticodon("CACCCCAGATAAGGGTG")
        self.assertEqual(result, "")

    def test_fail_no_matching_postcursor_get_anticodon(self):
        result = self.extractor.get_anticodon("CACCCCTGATTAGGGTG")
        self.assertEqual(result, "")

    def test_pass_one_extract_anticodon(self):
        result = self.extractor.extract_anticodon("AGGCTTGTAGCTCAGGTGGTtAGAGCGCACCCCTGATAAGGGTGAGGtCGGTGGTTCAAGTCCACTCAGGCCTACCA")
        self.assertEqual(result, ["GAT"])

    def test_pass_two_extract_anticodon(self):
        result = self.extractor.extract_anticodon("GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGtCGGCGGTTCGATCCCGTCATCACCCACCA")
        self.assertEqual(result, ["ACA", "TAC"])

    def test_fail_extract_anticodon(self):
        result = self.extractor.extract_anticodon("CGGGATGTAGCACAGTTGGCTAGCTCACCACGTTGGGACATGGAGGTCGGAAATTCGAGTCTTCTCATCCTGACCA")
        self.assertEqual(result, [])

if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(ExtractorTestCase)
    ut.TextTestRunner(verbosity=2).run(suite)
