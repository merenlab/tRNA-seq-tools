import unittest as ut
import sys

sys.path.append("../")

import sort_tRNA as sort


class SortTestCase(ut.TestCase):
    def setUp(self):
        self.sorter = sort.Sorter()

    def test_pass_is_tRNA(self):
        result = self.sorter.is_tRNA("AACCGTTGAACTGAAAGGTTCCTGGGGTTCGAATCCCCATCTCTCCGCCA")
        self.assertTrue(result[0])
    
    def test_fail_is_tRNA(self):
        result = self.sorter.is_tRNA("GAGTACCAAGATCGGAAGAGCACACGTCTAGTTCTACAGTCCGACGATCATCCTTTGG")
        self.assertFalse(result[0])

    def test_stats_file(self):
        with open(")

if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(SortTestCase)
    ut.TextTestRunner(verbosity=2).run(suite)
