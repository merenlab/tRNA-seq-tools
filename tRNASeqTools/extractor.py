# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to extract info from tRNA sequences."""

import csv
import Levenshtein as lev


__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class ExtractorStats:
    """This class handles keeping track of extraction statistics."""

    def __init__(self):
        """Initializes statistics."""
        self.total_seqs = 0
        self.type_I_seqs = 0
        self.type_II_seqs = 0
        self.type_I_match_dict = {8:0, 9:0}
        self.type_II_match_dict = {}
        self.subseq_match = 0
        for x in range (16, 27):
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

        self.extractor_stats = ExtractorStats()


    def pair_check(self, a_arm):
        """Checks a given anticodon arm for valid pairing"""
        pair_seg_length = 5
        total_mismatch = 0
        allowed_pairings = {"G":["C", "T"],
            "T":["A", "G"],
            "C":["G"],
            "A":["T"]}
        
        for x in range(pair_seg_length):
            if a_arm[-(x + 1)] not in allowed_pairings[a_arm[x]]:
                total_mismatch += 1
        return total_mismatch < 2


    def get_anticodon(self, a_arm):
        """Takes a given anticodon and returns the anticodon"""
        a_loop = a_arm[5:11]
        anticodon = ""
        if a_loop[1] == "T" and (a_loop[5] == "A" or a_loop[5] == "G"):
            anticodon = a_loop[2:5]
        return anticodon


    def extract_anticodon(self, seq):
        """Takes a given sequence and checks rules to find the anticodon, and
        returns the anticodon if there is one.
        """
        self.extractor_stats.total_seqs += 1
        length = len(seq)
        anticodon_arm_start = 24 + 8 + 17
        anticodon_arm_end = 24 + 8
        anticodon_list = []

        if lev.distance("GTTC", seq[-24:-20]) < 2:
            # handles type I full-length seqs
            if length < 78:
                # checks if there is an anticodon at a distance of 9
                if self.pair_check(seq[-(anticodon_arm_start +
                        1):-(anticodon_arm_end + 1)]):
                    anticodon = self.get_anticodon(seq[-(anticodon_arm_start + 1):
                        -(anticodon_arm_end + 1)])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_I_match_dict[9] += 1
                        self.extractor_stats.type_I_seqs += 1
                        anticodon_list.append(anticodon)
                # checks if there is an anticodon at a distance of 8
                if self.pair_check(seq[-anticodon_arm_start: -anticodon_arm_end]):
                    anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_I_match_dict[8] += 1
                        self.extractor_stats.type_I_seqs += 1
                        anticodon_list.append(anticodon)
            # handles type II full-length seqs
            elif length > 81:
                for x in range(16, 27):
                    anticodon_arm_start = 24 + x + 17
                    anticodon_arm_end = 24 + x
                 
                    if self.pair_check(seq[-anticodon_arm_start:
                        -anticodon_arm_end]):
                        anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                        if len(anticodon) > 0:
                            self.extractor_stats.type_II_match_dict[x] += 1
                            self.extractor_stats.type_II_seqs += 1
                            anticodon_list.append(anticodon)
                        else:
                            continue
        else:
            print("error: GTTC didn't match")
            print(seq[-24:-20])

        return anticodon_list

    
    def extract_anticodon_not_full_length(self, seq):
        """Takes a given sequence (not full length) and checks rules to find
        the anticodon, and returns the anticodon if there is one.
        """
        self.extractor_stats.total_seqs += 1
        length = len(seq)
        anticodon_arm_start = 24 + 8 + 17
        anticodon_arm_end = 24 + 8
        anticodon_list = []

        if lev.distance("GTTC", seq[-24:-20]) < 2:
            if length > 50:
                # checks if there is an anticodon at a distance of 9
                if self.pair_check(seq[-(anticodon_arm_start +
                        1):-(anticodon_arm_end + 1)]):
                    anticodon = self.get_anticodon(seq[-(anticodon_arm_start + 1):
                        -(anticodon_arm_end + 1)])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_I_match_dict[9] += 1
                        self.extractor_stats.type_I_seqs += 1
                        anticodon_list.append(anticodon)
                # checks if there is an anticodon at a distance of 8
                if self.pair_check(seq[-anticodon_arm_start: -anticodon_arm_end]):
                    anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                    if len(anticodon) > 0:
                        self.extractor_stats.type_I_match_dict[8] += 1
                        self.extractor_stats.type_I_seqs += 1
                        anticodon_list.append(anticodon)
            if length > 67:
                for x in range(16, 27):
                    anticodon_arm_start = 24 + x + 17
                    anticodon_arm_end = 24 + x
                    
                    if self.pair_check(seq[-anticodon_arm_start:
                        -anticodon_arm_end]):
                        anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                        if len(anticodon) > 0:
                            self.extractor_stats.type_II_match_dict[x] += 1
                            self.extractor_stats.type_II_seqs += 1
                            anticodon_list.append(anticodon)
                        else:
                            continue
        else:
            print("error: GTTC didn't match")
            print(seq[-24:-20])

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


