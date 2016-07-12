# coding: utf-8
import unittest as ut

import tRNAsorter.sorter as sort

class SortTestCase(ut.TestCase):
    def setUp(self):
        self.sorter = sort.Sorter()
        self.maxDiff = None
        self.cur_specs = sort.SeqSpecs()

    def test_pass_is_tRNA(self):
        result = self.sorter.is_tRNA("AACCGTTGAACTGAAAGGTTCCTGGGGTTCGAATCCCCATCTCTCCGCCA")
        self.assertTrue(result[0])
    
    def test_fail_is_tRNA(self):
        result = self.sorter.is_tRNA("GAGTACCAAGATCGGAAGAGCACACGTCTAGTTCTACAGTCCGACGATCATCCTTTGG")
        self.assertFalse(result[0])

    def test_pass_check_full_length(self):
        self.cur_specs.seq = "GCGCGAGTGGCTCAGGGGTGGAGCACCACCTTGCCAAGGTGGGGGCCGCGGGTTCGAATCCCGTCTCGCGCTCCA"
        self.cur_specs.length = len(self.cur_specs.seq)
        cur_specs_res = self.sorter.check_full_length(self.cur_specs)
        self.assertTrue(cur_specs_res.full_length)

    def test_fail_too_short_check_full_length(self):
        self.cur_specs.seq = "CTAGAGTTCAAATCTCCCTTCCGCTACCA"
        self.cur_specs.length = len(self.cur_specs.seq)
        cur_specs_res = self.sorter.check_full_length(self.cur_specs)
        self.assertFalse(cur_specs_res.full_length)

    def test_fail_no_match_check_full_length(self):
        self.cur_specs.seq = "TGGTGGGTATAGCGCAGTTGGTTAGTGCGCCAGATTGTGGCTCTGGAGGCCGAGGGTTCGAATCCCTTTACCCACCCCA"
        self.cur_specs.length = len(self.cur_specs.seq)
        cur_specs_res = self.sorter.check_full_length(self.cur_specs)
        self.assertFalse(cur_specs_res.full_length)

    def test_split_3_trailer_with_trailer(self):
        self.cur_specs.seq = "CTCCAGGTTCGAGTCCTGGTAGAACAACCAA"
        self.cur_specs.length = len(self.cur_specs.seq)
        cur_specs_res = self.sorter.split_3_trailer(self.cur_specs, 1)
        self.assertEqual(cur_specs_res.seq, "CTCCAGGTTCGAGTCCTGGTAGAACAACCA")
        self.assertEqual(cur_specs_res.three_trailer, "A")
        self.assertEqual(cur_specs_res.length, 30)
        self.assertEqual(cur_specs_res.trailer_length, 1)
        
    def test_split_3_trailer_without_trailer(self):
        self.cur_specs.seq = "CTCCAGGTTCGAGTCCTGGTAGAACAACCA"
        self.cur_specs.length = len(self.cur_specs.seq)
        cur_specs_res = self.sorter.split_3_trailer(self.cur_specs, 0)
        self.assertEqual(cur_specs_res.seq, "CTCCAGGTTCGAGTCCTGGTAGAACAACCA")
        self.assertEqual(cur_specs_res.three_trailer, "")
        self.assertEqual(cur_specs_res.length, 30)
        self.assertEqual(cur_specs_res.trailer_length, 0)
    
    def test_stats_file(self):
        with open("sort_stats_gold") as gold_file:
            with open("../sort_stats") as test_file:
                gold_lines = gold_file.readlines()
                test_lines = test_file.readlines()
                self.assertEqual(gold_lines, test_lines)
                    

if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(SortTestCase)
    ut.TextTestRunner(verbosity=2).run(suite)
