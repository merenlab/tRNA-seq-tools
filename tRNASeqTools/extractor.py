# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to extract info from tRNA sequences."""

import csv
import Levenshtein as lev

import tRNASeqTools
import tRNASeqTools.filters as filters

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = tRNASeqTools.__version__
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"

class ExtractorStats:
    """This class handles keeping track of extraction statistics."""
    
    def __init__(self, guidelines):
        """Initializes statistics."""
        self.total_seqs = 0
        self.type_I_seqs = 0
        self.type_II_seqs = 0
        self.type_I_match_dict = {}
        self.type_II_match_dict = {}
        self.subseq_match = 0
        for x in guidelines[0]:
            self.type_I_match_dict[x] = 0
        for x in guidelines[1]:
            self.type_II_match_dict[x] = 0

    def format_line(self, label, value, level, padding = 55):
        """Handles indenting/formatting lines for statistics."""
        levels_dict = {1:"%s%s\t%s\n" % 
            (label, " " + " " * (padding - len(label)), value),
            2:"\t%s%s\t%s\n" % 
                (label, " " + " " * (padding - (4 + len(label))), value),
            3:"\t\t%s%s\t%s\n" % 
                (label, " " + " " * (padding - (12 + len(label))), value)}
        return levels_dict[level]


    def write_stats(self, out_file_path):
        """Writes statistics to an output file.""" 
        with open(out_file_path, "w") as outfile:
            outfile.write(self.format_line("Total seqs", "%d" % 
                self.total_seqs,1))
            
            outfile.write(self.format_line("Type I Seqs", "%d" % 
                self.type_I_seqs,2))
            for key in self.type_I_match_dict:
                outfile.write(self.format_line("Match at dist " + str(key), 
                    "%d" % self.type_I_match_dict[key],3))
            
            outfile.write(self.format_line("Type II Seqs", "%d" % 
                self.type_II_seqs,2))
            for key in self.type_II_match_dict:
                outfile.write(self.format_line("Match at dist " + str(key), 
                    "%d" % self.type_II_match_dict[key],3))
            
            outfile.write(self.format_line("Sub Seq Matches", "%d" % 
                self.subseq_match,2))


class Extractor:
    """This class handles the extraction of info from tRNAs"""
    
    def __init__(self):
        """Initializes variables for the extractor"""
        self.extractor_stats_file = ""
        self.loop_guidelines = filters.IsTRNA("").getAnticodonGuidelines()
        self.extractor_stats = ExtractorStats([self.loop_guidelines[2], self.loop_guidelines[3]])

        self.allowed_pairings = {"G":["C", "T"], "T":["A", "G"], "C":["G"], "A":["T"]}


    def pair_check(self, a_arm):
        """Checks a given anticodon arm for valid pairing"""
        pair_seg_length = 5
        total_mismatch = 0
        
        for x in range(pair_seg_length):
            if a_arm[-(x + 1)] not in self.allowed_pairings[a_arm[x]]:
                total_mismatch += 1
        return total_mismatch < self.loop_guidelines[0] + 1


    def get_anticodon(self, a_arm):
        """Takes a given anticodon and returns the anticodon"""
        a_loop = a_arm[5:11]
        anticodon = ""
        if a_loop[1] in self.loop_guidelines[1][0] and a_loop[5] in self.loop_guidelines[1][1]:
            anticodon = a_loop[2:5]
        return anticodon


    def extract_anticodon(self, seq, fullLength):
        """Takes a given sequence and checks rules to find the anticodon, and
        returns the anticodon if there is one.
        """
        self.extractor_stats.total_seqs += 1
        length = len(seq)
        anticodon_arm_start = 24 + 17
        anticodon_arm_end = 24
        anticodon_list = []

        # handles type I full-length seqs
        if fullLength and length < 78 or not fullLength and length > 50:
            for x in self.loop_guidelines[2]:
                if self.pair_check(seq[-(anticodon_arm_start + x):-(anticodon_arm_end + x)]):
                    anticodon = self.get_anticodon(seq[-(anticodon_arm_start + x):-(anticodon_arm_end + x)])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_I_match_dict[x] += 1
                        self.extractor_stats.type_I_seqs += 1
                        anticodon_list.append(anticodon)

        # handles type II full-length seqs
        if fullLength and length > 81 or not fullLength and length > 67:
            for x in self.loop_guidelines[3]:
                if self.pair_check(seq[-(anticodon_arm_start + x):-(x + anticodon_arm_end)]):
                    anticodon = self.get_anticodon(seq[-(anticodon_arm_start + x): -(anticodon_arm_end + x)])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_II_match_dict[x] += 1
                        self.extractor_stats.type_II_seqs += 1
                        anticodon_list.append(anticodon)

        return anticodon_list


    def match_unassigned_sequences(self, file_name, max_seq_width, fieldnames):
        """takes unmatched sequences and attempts to match them to already
        matched sequences"""
        
        match_dict = {}
        unassigned_rows = []
        assigned_rows = []

        with open(file_name, "r") as readfile:
            readfile_reader = csv.DictReader(readfile, delimiter="\t")
            for row in readfile_reader:
                if row["Anticodon"]:
                    match_dict[row["Seq"].strip("-")] = row["Anticodon"]
                    assigned_rows.append(row)
                else:
                    unassigned_rows.append(row)

        for row in unassigned_rows:
            seq = row["Seq"].strip("-")
            for key in list(match_dict.keys()):
                if seq in key:
                    self.extractor_stats.subseq_match += 1
                    row["Anticodon"] = match_dict[key]
                    
        with open(file_name, "w") as writefile:
            writefile_writer = csv.DictWriter(writefile, fieldnames=fieldnames,
                delimiter="\t")
            writefile_writer.writeheader()
            writerows = assigned_rows + unassigned_rows
            for row in writerows:
                writefile_writer.writerow(row)


