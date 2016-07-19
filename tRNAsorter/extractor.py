# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to extract info from tRNA sequences."""

import os 
import sys
import csv
import Levenshtein as lev

import Oligotyping.lib.fastalib as u

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

class Extractor:
    """This class handles the extraction of info from tRNAs"""

    def __init__(self):
        """Initializes variables for the extractor"""
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
        length = len(seq)
        anticodon_arm_start = 24 + 8 + 17
        anticodon_arm_end = 24 + 8

        if lev.distance("GTTC", seq[-24:-20]) < 2:
            if length < 78:
                if self.pair_check(seq[-anticodon_arm_start: -anticodon_arm_end]):
                    anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                    #print "8 passed: " + anticodon
                    return anticodon
                elif self.pair_check(seq[-(anticodon_arm_start +
                        1):-(anticodon_arm_end + 1)]):
                    anticodon = self.get_anticodon(seq[-(anticodon_arm_start + 1):
                        -(anticodon_arm_end + 1)])
                    #print "9 passed: " + anticodon
                    return anticodon
            elif length > 81:
                for x in range(17, 27):
                    anticodon_arm_start = 24 + x + 17
                    anticodon_arm_end = 24 + x
                    if self.pair_check(seq[-anticodon_arm_start:
                        -anticodon_arm_end]):
                        anticodon = self.get_anticodon(seq[-anticodon_arm_start: -anticodon_arm_end])
                        if len(anticodon) > 0:
                           #print str(x) + " passed: " + anticodon
                            return anticodon
                        else:
                            continue
            else:
                print "length fell between type I and type II"
                    
        else:
            print "error: GTTC didn't match"
            print seq[-24:-20]


