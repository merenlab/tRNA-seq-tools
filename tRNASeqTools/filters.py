#Filters:

import os
import re

import tRNASeqTools.extractor as extractor

class IsTRNA:

    def __init__(self, output_path):
        self.ANTICODON_LOOP_GUIDELINES = [0, (('T'), ('A', 'G')), [], [], 0]
        self.allowed_pairings = {"G":("C", "T"), "T":("A", "G"), "C":("G"), "A":("T")}
        self.sub_size = 24
        self.output_path = output_path + "/filteredSequences/"
        self.T_LOOP_AND_ACCEPTOR_GUIDELINES = [[], 0, 0]
        self.SET_UP_FILTERS = {"Allow_one_mismatch_in_the_anticodon_pairs": self.change_anticodon_loop_guidelines(0, 1), #Canonical 
                               "Positions_34_and_37": self.change_anticodon_loop_guidelines(1, (('T'), ('A', 'G'))), #Canonical 
##                               "Anticodon_arm_starting_pos_at_11": self.change_anticodon_loop_guidelines(4, 11),
##                               "Anticodon_arm_starting_pos_at_24": self.change_anticodon_loop_guidelines(4, 24), #Canonical 
##                               "Type_I_length_between_8_and_9": self.change_anticodon_loop_guidelines(2, [8, 9]), #Canonical 
##                               "Type_II_length_between_16_and_27": self.change_anticodon_loop_guidelines(3, range(16, 27)), #Canonical 
                               "T_region_length_between_4_and_20": self.change_anticodon_loop_guidelines(2, range(4, 21)),
##                               "require_T_Loop_G_at_0": self.change_T_and_acc_guidelines(0, (0, "G")), #Canonical 
##                               "require_T_Loop_T_at_1": self.change_T_and_acc_guidelines(0, (1, "T")), #Canonical 
##                               "require_T_Loop_T_at_2": self.change_T_and_acc_guidelines(0, (2, "T")), #Canonical 
##                               "require_T_Loop_C_at_3": self.change_T_and_acc_guidelines(0, (3, "C")), #Canonical 
##                               "require_T_Loop_C_at_8": self.change_T_and_acc_guidelines(0, (8, "C")), #Canonical 
                               "require_acceptor_C_at_-3": self.change_T_and_acc_guidelines(0, (-3, "C")), #Canonical 
                               "require_acceptor_C_at_-2": self.change_T_and_acc_guidelines(0, (-2, "C")), #Canonical 
                               "require_acceptor_A_at_-1": self.change_T_and_acc_guidelines(0, (-1, "A")), #Canonical 
##                               "Allow_one_mismatch_in_T-loop_and_acceptor": self.change_T_and_acc_guidelines(1, 1), #Canonical 
                               "Allow_no_mismatches_in_T-loop_and_acceptor": self.change_T_and_acc_guidelines(1, 0),
                               "Require_Acceptor_Stem_Matching_with_one_mismatch":  self.change_T_and_acc_guidelines(2, (True, 2))
                               }
        self.FILTERS = {"Longer_than_30": lambda seq: len(seq) > 30, #Canonical: 24
                        "Shorter_than_200": lambda seq: len(seq) < 61, #Canonical
                        "Anticodon_is_known": self.isAnticodonKnown,
                        "T_loop_and_acceptor_is_acceptable":  self.t_loop_and_acceptor, #Canonical
                        "D_region_and_T_region_acceptable_length": self.check_D_and_T_region_lengths
                        }

        for elem in self.SET_UP_FILTERS:
            self.SET_UP_FILTERS[elem]
        if self.T_LOOP_AND_ACCEPTOR_GUIDELINES[2] == 0:
            self.T_LOOP_AND_ACCEPTOR_GUIDELINES = [False]

        self.ANTICODON_LOOP_GUIDELINES[4] = 11

        self.D_region_range = []
        for i in range(4, 17):
            self.D_region_range.append(i)
        for i in range(len(self.D_region_range)):
            self.D_region_range[i] += 7 + 6
        self.T_region_range = []
        for i in range(4, 21):
            self.T_region_range.append(i)
        for j in range(len(self.T_region_range)):
            self.T_region_range[j] = -(self.T_region_range[j]  + 10 + 5)
        
        FILTER_DESCRIPTIONS = {"Allow_one_mismatch_in_the_anticodon_pairs": "This filter allows a single mismatch when paring the anticodon stem",
                               "Positions_34_and_37": "This filter requires that position 34 (right before the anticodon) is a T, and position 37 (right after the anticodon) is an A or G",
                               "Type_I_length_between_8_and_9": "For a sequence to be counted as a Type I sequence, the distance between the G at the end of the T-Stem loop and the anticodon must be 8 or 9 (this is the 4-5 of the V-loop + 4)",
                               "Type_II_length_between_16_and_27": "For a sequence to be counted as a Type II sequence, the distance between the G at the end of the T-Stem loop and the anticodon must be 16-27 (this is the 12-23 of the V-loop + 4)",
                               "require_T_Loop_G_at_0": "This filter requires that the 5' T-stem ends with a G",
                               "require_T_Loop_T_at_1": "This filter requires that the T-loop starts with a T",
                               "require_T_Loop_T_at_2": "This filter requires that the T-loop's second position is a T",
                               "require_T_Loop_C_at_3": "This filter requires that the T-loop's third position is a C",
                               "require_T_Loop_C_at_8": "This filter requires that the 3' T-stem stars with a C",
                               "require_acceptor_C_at_minus_3": "This filter requires that the acceptor ends with CNN",
                               "require_acceptor_C_at_minus_2": "This filter requires that the acceptor ends with NCN",
                               "require_acceptor_A_at_minus_1": "This filter requires that the acceptor ends with NNA",
                               "Allow_one_mismatch_in_T-loop_and_acceptor": "This filter allows a single mismatch within the requirements for the T-loop and the acceptor, as listed above",
                                "Longer_than_24": "This filter requires that the sequence be longer than 24 nucleotides",
                               "Shorter_than_200": "This filter requires that the sequence be shorter than 200 nucleotides",
                                "Anticodon_is_known": "This filter discards all sequences where the anticodon was not found with the current parameters",
                                "T-Loop_and_acceptor are acceptable": "This filter discards all sequences where the T-loop and Acceptor are not found using the current parameters"}

    def getFilters(self):
        return self.FILTERS

    def getSetUpFilters(self):
        return self.SET_UP_FILTERS
    
    def istRNA(self, seq, name):
        problems = []
        self.anticodon = []
        self.name = name
        for filt in self.FILTERS:
            if not self.FILTERS[filt](seq):
                open(self.output_path + filt, "a").write(self.name + "\n" + seq + "\n")
                return filt
        return ""

    def change_anticodon_loop_guidelines(self, i, changeItTo):
        self.ANTICODON_LOOP_GUIDELINES[i] = changeItTo

    def getAnticodonGuidelines(self):
        return self.ANTICODON_LOOP_GUIDELINES
    
    def isAnticodonKnown(self, seq):
        if self.anticodon == []:
            self.extractor = extractor.Extractor()
            fullLength = len(seq) > 70 and len(seq) < 100 and seq[7] == "T" and seq[13] == "A"
            self.anticodon = self.extractor.extract_anticodon(seq, fullLength)
        return not self.anticodon == []
            
    def check_D_and_T_region_lengths(self, seq):
        if self.isAnticodonKnown(seq):
            for anti in self.anticodon:
                for pos34 in self.ANTICODON_LOOP_GUIDELINES[1][0]:
                    for elem in [m.start() for m in re.finditer(pos34 + anti, seq)]:
                        if elem in self.D_region_range and elem - len(seq) in self.T_region_range:
                            return True
        return False     

    def get_t_loop_and_acceptor_guidelines(self):
        return self.T_LOOP_AND_ACCEPTOR_GUIDELINES

    def change_T_and_acc_guidelines(self, i, changeItTo):
        if type(self.T_LOOP_AND_ACCEPTOR_GUIDELINES[i]) == int:
            self.T_LOOP_AND_ACCEPTOR_GUIDELINES[i] = changeItTo
        else:
            self.T_LOOP_AND_ACCEPTOR_GUIDELINES[i].append(changeItTo)

    def t_loop_and_acceptor(self, seq):
        length = len(seq)
        shortestMissed = [0]
        for elem in self.T_LOOP_AND_ACCEPTOR_GUIDELINES:
            shortestMissed.append(0)
        for i in range(length - self.sub_size + 1):
            sub_str = seq[-(i + self.sub_size):(length - i)]
            missed = []
            for position_tuple in self.T_LOOP_AND_ACCEPTOR_GUIDELINES[0]:
                if sub_str[position_tuple[0]] != position_tuple[1]:
                    missed.append(position_tuple)
            if len(missed) < self.T_LOOP_AND_ACCEPTOR_GUIDELINES[1] + 1:
                if self.T_LOOP_AND_ACCEPTOR_GUIDELINES[2][0]:
                    misses = 0 
                    for j in range(7):
                        if seq[j] not in self.allowed_pairings[seq[-5 - j]]:
                            misses += 1
                    if misses < self.T_LOOP_AND_ACCEPTOR_GUIDELINES[2][1]:
                        return True
                    else:
                        open(self.output_path + "Require_Acceptor_Stem_Matching_with_one_mismatch", "a").write(self.name + "\n" + seq + "\n")
        
                else:
                    return True
            if len(missed) < len(shortestMissed):
                shortestMissed = missed
        for elem in shortestMissed:
            fileName = "require_"
            if elem[0] > 0:
                fileName += "T_Loop_"
            else:
                fileName += "acceptor_"
            fileName += elem[1] + "_at_" +  str(elem[0])
            open(self.output_path + fileName, "a").write(self.name + "\n" + seq + "\n")
        open(self.output_path + "Allow_one_mismatch_in_T-loop_and_acceptor", "a").write(self.name + "\n" + seq + "\n")
        return False
