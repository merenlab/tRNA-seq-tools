# coding: utf-8
import unittest as ut

import tempfile
import os
import subprocess
import csv
import shutil

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
        self.assertEqual(result, ["TAC", "ACA"])

    def test_fail_extract_anticodon(self):
        result = self.extractor.extract_anticodon("CGGGATGTAGCACAGTTGGCTAGCTCACCACGTTGGGACATGGAGGTCGGAAATTCGAGTCTTCTCATCCTGACCA")
        self.assertEqual(result, [])
    
    def check_output_file(self, file_path):
        tempdir = tempfile.mkdtemp()
        try:
            os.chdir(tempdir)
            subprocess.call(["trna-sort", "-n", "testing-run", file_path])
            gold_file = open("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file", "r")
            with open("testing-run_TAB_NO_TRAILER", "r") as test_file:
                test_file_reader = csv.DictReader(test_file, delimiter="\t")
                gold_file_reader = csv.DictReader(gold_file, delimiter="\t")
                for row in test_file_reader:
                    testing_row = {}
                    id_string_list = row["ID"].split(" ")
                    testing_row["ID"] = id_string_list[0]
                    testing_row["Anticodon"] = row["Anticodon"]
                    self.assertIn(testing_row, gold_file_reader)
        finally:
            shutil.rmtree(tempdir)
           

    def test_check_output_bactBact(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/bactBact_CF-tRNAs.fa")

    def test_check_output_bactFrag(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/bactFrag_638R-tRNAs.fa")

    def test_check_output_e_coli_1(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/eschColi_042-tRNAs.fa")

    def test_check_output_e_coli_2(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/eschColi_K_12_MG1655-tRNAs.fa")

    def test_check_output_mycoTube(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/mycoTube_0B076XDR-tRNAs.fa")

    def test_check_output_pseuPuti(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/pseuPuti_H8234-tRNAs.fa")

    def test_check_output_rhodPalu(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/rhodPalu_HAA2-tRNAs.fa")

    def test_check_output_saccDegr(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/saccDegr_2_40-tRNAs.fa")

    def test_check_output_salmEnte(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/salmEnte_ENTERICA_ENTER_111-tRNAs.fa")

    def test_check_output_streSuis(self):
        self.check_output_file("/groups/merenlab/scui/meren-lab-tools/testing/extractor_test_file_creator/streSuis_GZ1-tRNAs.fa")



if __name__ == "__main__":
    suite = ut.TestLoader().loadTestsFromTestCase(ExtractorTestCase)
    ut.TextTestRunner(verbosity=2).run(suite)
