#Filters:

import os

class IsTRNA:

    def __init__(self, output_path):
        self.ANTICODON_LOOP_GUIDELINES = [0, 0, 0, 0]
        self.allowed_pairings = {"G":("C", "T"), "T":("A", "G"), "C":("G"), "A":("T")}
        self.sub_size = 24
        self.output_path = output_path + "/filteredSequences/"
        self.T_LOOP_AND_ACCEPTOR_GUIDELINES = [[], 0]
        self.SET_UP_FILTERS = {"Allow_one_mismatch_in_the_anticodon_pairs": self.change_anticodon_loop_guidelines(0, 1),
                               "Positions_34_and_37": self.change_anticodon_loop_guidelines(1, (('T'), ('A', 'G'))),
                               "Type_I_length_between_8_and_9": self.change_anticodon_loop_guidelines(2, [8, 9]),
                               "Type_II_length_between_16_and_27": self.change_anticodon_loop_guidelines(3, range(16, 27)),
                               "require_T_Loop_G_at_0": self.change_T_and_acc_guidelines(0, (0, "G")),
                               "require_T_Loop_T_at_1": self.change_T_and_acc_guidelines(0, (1, "T")),
                               "require_T_Loop_T_at_2": self.change_T_and_acc_guidelines(0, (2, "T")),
                               "require_T_Loop_C_at_3": self.change_T_and_acc_guidelines(0, (3, "C")),
                               "require_T_Loop_C_at_8": self.change_T_and_acc_guidelines(0, (8, "C")),
                               "require_acceptor_C_at_-3": self.change_T_and_acc_guidelines(0, (-3, "C")),
                               "require_acceptor_C_at_-2": self.change_T_and_acc_guidelines(0, (-2, "C")),
                               "require_acceptor_A_at_-1": self.change_T_and_acc_guidelines(0, (-1, "A")),
                               "Allow_one_mismatch_in_T-loop_and_acceptor": self.change_T_and_acc_guidelines(1, 1),
                               }
        self.FILTERS = {"Longer_than_24": lambda seq: len(seq) > 24,
                        "Shorter_than_200": lambda seq: len(seq) < 200,
             ##           "Anticodon_is_known": self.isAnticodonKnown,
                        "T-Loop_and_acceptor are acceptable":  self.t_loop_and_acceptor}

        for elem in self.SET_UP_FILTERS:
            self.SET_UP_FILTERS[elem]
        
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
        anticodon = ""
        length = len(seq)
        anticodon_arm_start = 24 + 17
        anticodon_arm_end = 24
        fullLength = len(seq) > 70 and len(seq) < 100 and seq[7] == "T" and seq[13] == "A"
        if fullLength and length < 78 or not fullLength and length > 50:
            for x in self.ANTICODON_LOOP_GUIDELINES[2]:
                if pair_check(seq[-(anticodon_arm_start + x):-(anticodon_arm_end + x)]):
                    if get_anticodon(seq[-(anticodon_arm_start + x):-(anticodon_arm_end + x)]):
                        return True
        if fullLength and length > 81 or not fullLength and length > 67:
            for x in self.ANTICODON_LOOP_GUIDELINES[3]:
                if pair_check(seq[-(anticodon_arm_start + x):-(anticodon_arm_end + x)]):
                    if get_anticodon(seq[-(anticodon_arm_start + x): -(anticodon_arm_end + x)]):
                        return True
        return False

    def pair_check(self, a_arm):
        total_mismatch = 0
        
        for x in range(5):
            if a_arm[-(x + 1)] not in self.allowed_pairings[a_arm[x]]:
                total_mismatch += 1
        return total_mismatch < self.ANTICODON_LOOP_GUIDELINES[0] + 1


    def get_anticodon(self, a_arm):
        a_loop = a_arm[5:11]
        anticodon = ""
        if a_loop[1] in self.ANTICODON_LOOP_GUIDELINES[1][0] and a_loop[5] in self.ANTICODON_LOOP_GUIDELINES[1][1]:
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
        shortestMissed = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(length - self.sub_size + 1):
            sub_str = seq[-(i + self.sub_size):(length - i)]
            missed = []
            for position_tuple in self.T_LOOP_AND_ACCEPTOR_GUIDELINES[0]:
                if sub_str[position_tuple[0]] != position_tuple[1]:
                    missed.append(position_tuple)
            if len(missed) < self.T_LOOP_AND_ACCEPTOR_GUIDELINES[1] + 1:
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
