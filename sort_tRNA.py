import os 
import sys
import Oligotyping.lib.fastalib as u 
import argparse
import csv
import Levenshtein as lev


class SorterStats:
    """This class handles keeping track of sort statistics."""

    def __init__(self):
        """Initializes statistics."""
        self.total_seqs = 0
        self.total_rejected = 0
        self.total_passed = 0
        self.num_trailer = 0
        self.total_full_length = 0

        self.no_divergence = 0
        self.t_loop_divergence = 0
        self.div_at_0 = 0
        self.div_at_1 = 0
        self.div_at_2 = 0
        self.div_at_3 = 0
        self.div_at_8 = 0
        self.acceptor_divergence = 0
        self.div_at_neg_1 = 0
        self.div_at_neg_2 = 0
        self.div_at_neg_3 = 0

        self.t_loop_seq_rejected = 0
        self.acceptor_seq_rejected = 0
        self.both_rejected = 0
        self.short_rejected = 0
       

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
            outfile.write(self.format_line("Total full-length", "%d" %
                self.total_full_length, 1))
            outfile.write(self.format_line("With trailer", "%d" %
                self.num_trailer, 1))
            outfile.write(self.format_line("Total passed", "%d" % 
                self.total_passed, 1))

            outfile.write(self.format_line("No divergence", "%d" % 
                self.no_divergence, 2))
            outfile.write(self.format_line("T-loop divergence", "%d" % 
                self.t_loop_divergence, 2))
            outfile.write(self.format_line("Divergence at pos 0", "%d" % 
                self.div_at_0, 3))
            outfile.write(self.format_line("Divergence at pos 1", "%d" % 
                self.div_at_1, 3))
            outfile.write(self.format_line("Divergence at pos 2", "%d" % 
                self.div_at_2, 3))
            outfile.write(self.format_line("Divergence at pos 3", "%d" % 
                self.div_at_3, 3))
            outfile.write(self.format_line("Divergence at pos 8", "%d" %
                self.div_at_8, 3))
            outfile.write(self.format_line("Acceptor divergence", "%d" % 
                self.acceptor_divergence, 2))
            outfile.write(self.format_line("Divergence at pos -3", "%d" % 
                self.div_at_neg_3, 3))
            outfile.write(self.format_line("Divergence at pos -2", "%d" % 
                self.div_at_neg_2, 3))
            outfile.write(self.format_line("Divergence at pos -1", "%d" % 
                self.div_at_neg_1, 3))

            outfile.write(self.format_line("Total failed", "%d" % 
                self.total_rejected, 1))
            outfile.write(self.format_line("T-loop seq rejected", "%d" % 
                self.t_loop_seq_rejected, 2))
            outfile.write(self.format_line("Acceptor seq rejected", "%d" % 
                self.acceptor_seq_rejected, 2))
            outfile.write(self.format_line("Both rejected", "%d" % 
                self.both_rejected, 2))
            outfile.write(self.format_line("Short rejected", "%d" % 
                self.short_rejected, 2))
 
       
class SeqSpecs:
    """Class to store seq information during sort."""

    def __init__(self):
        """Initializes seq information categories."""
        self.length = 0
        self.mis_count = 100
        self.t_loop_error = True
        self.acceptor_error = True
        self.seq = ""
        self.seq_sub = ""
        self.full_length = False
        self.t_loop_seq = ""
        self.acceptor_seq = ""
        self.three_trailer = ""
        self.trailer_length = 0


    def gen_id_string(self, id):
        """Generates ID string for _PASSED and _FAILED files"""
        mod_id_list = []
        mod_id_list.append(id)
        mod_id_list.append("mismatches:" + str(self.mis_count))
        mod_id_list.append("t_loop_error:" + str(self.t_loop_error))
        mod_id_list.append("acceptor_error:" + str(self.acceptor_error))
        mod_id_list.append("seq_sub:" + self.seq_sub)
        if self.full_length:
            mod_id_list.append("full length: Yes")
        else:
            mod_id_list.append("full length: Maybe")
        mod_id_list.append("t_loop:" + self.t_loop_seq)
        mod_id_list.append("acceptor:" + self.acceptor_seq)
        mod_id_list.append("5_trailer:" + self.three_trailer)
        return "|".join(mod_id_list)


    def write_specs(self, writer, id):
        """Uses the passed in DictWriter to write to a csv output file
        (_TAB_NO_TRAILER or _TAB_TRAILER).
        """
        temp_dict = {"ID" : id, "Seq" : self.seq, "3-trailer" :
            self.three_trailer, "t-loop" : self.t_loop_seq, "acceptor" :
            self.acceptor_seq, "full-length" : str(self.full_length), 
            "Seq_length" : str(self.length), "Trailer_length" :
            str(self.trailer_length)}
        writer.writerow(temp_dict) 


class Sorter:
    """Class that handles the sorting of the seqs."""
    
    def __init__(self):
        """Initializes variables for input/output files and statistics."""
        self.passed_seqs_write_fasta = ""
        self.rejected_seqs_write_fasta = ""
        self.stats_write_file = ""
        self.read_fasta = ""
        self.no_trailer_tabfile = ""
        self.trailer_tabfile = ""

        self.fieldnames = ["ID", "Seq", "3-trailer", "t-loop", "acceptor",
            "full-length", "Seq_length", "Trailer_length"]
        self.max_seq_width = len(self.fieldnames[1])

        self.stats = SorterStats()


    def set_file_names(self, args):
        """Takes command line arguments from args and assigns input/output file
        variables.
        """
        self.passed_seqs_write_fasta = u.FastaOutput(args.sample_name +
            "_PASSED")
        self.rejected_seqs_write_fasta = u.FastaOutput(args.sample_name +
            "_FAILED")
        self.stats_write_file = args.sample_name + "_STATS"
        self.read_fasta = u.SequenceSource(args.readfile)
        self.no_trailer_tabfile = args.sample_name + "_TAB_NO_TRAILER"
        self.trailer_tabfile = args.sample_name + "_TAB_TRAILER"


    def check_divergence_pos(self, cur_seq_specs):
        """Takes a SeqSpecs class and updates statistics on divergence
        position.
        """
        if cur_seq_specs.t_loop_error:
            self.stats.t_loop_divergence += 1
            if cur_seq_specs.seq_sub[0] != "G":
                self.stats.div_at_0 += 1
            elif cur_seq_specs.seq_sub[1] != "T":
                self.stats.div_at_1 += 1
            elif cur_seq_specs.seq_sub[2] != "T":
                self.stats.div_at_2 += 1
            elif cur_seq_specs.seq_sub[3] != "C":
                self.stats.div_at_3 += 1
            elif cur_seq_specs.seq_sub[8] != "C":
                self.stats.div_at_8 += 1
        elif cur_seq_specs.acceptor_error:
            self.stats.acceptor_divergence += 1
            if cur_seq_specs.seq_sub[-3] != "C":
                self.stats.div_at_neg_3 += 1
            elif cur_seq_specs.seq_sub[-2] != "C":
                self.stats.div_at_neg_2 += 1
            elif cur_seq_specs.seq_sub[-1] != "A":
                self.stats.div_at_neg_1 += 1
        else:
            self.stats.no_divergence += 1


    def check_full_length(self, cur_seq_specs):
        """Takes a SeqSpecs class and checks whether or not the sequence is a
        full-length sequence, returns an updated SeqSpecs class.
        """
        if cur_seq_specs.length > 70 and cur_seq_specs.length < 100:
            if cur_seq_specs.seq[7] == "T" and cur_seq_specs.seq[13] == "A":
                self.stats.total_full_length += 1
                cur_seq_specs.full_length = True
        return cur_seq_specs


    def split_3_trailer(self, cur_seq_specs, i):
        """Takes a SeqSpecs class and splits the 3-trailer from the seqience,
        returns an updated SeqSpecs class.
        """
        full_seq = cur_seq_specs.seq
        cur_seq_specs.seq = full_seq[:(cur_seq_specs.length - i)]
        cur_seq_specs.three_trailer = full_seq[(cur_seq_specs.length - i):]
        cur_seq_specs.length = len(cur_seq_specs.seq)
        cur_seq_specs.trailer_length = len(cur_seq_specs.three_trailer)
        return cur_seq_specs


    def handle_pass_seq(self, cur_seq_specs, i):
        """Consdolidates all the methods run specifically for a confirmed passed
        sequence.
        """
        self.stats.total_passed += 1
        self.check_divergence_pos(cur_seq_specs)
        cur_seq_specs = self.check_full_length(cur_seq_specs)
        cur_seq_specs = self.split_3_trailer(cur_seq_specs, i)
        return cur_seq_specs


    def is_tRNA(self, seq):
        """Takes a sequence and determines whether or not it matches the
        criterion for being a tRNA
        """
        length = len(seq)
        sub_size = 24
        t_loop_error = True
        acceptor_error = True
        cur_seq_specs = SeqSpecs()

        # Start the sliding window at the last 24 bases, and move to the left
        # one at a time
        for i in xrange(length - sub_size + 1):
            sub_str = seq[-(i + sub_size):(length - i)]
            t_loop_seq = sub_str[0:9]
            acceptor_seq = sub_str[-3:]
            t_loop_dist = (lev.distance("GTTC", sub_str[0:4]) 
                + lev.distance("C", sub_str[8]))
            acceptor_dist = lev.distance("CCA", sub_str[-3:])
            mis_count = t_loop_dist + acceptor_dist
           
            if t_loop_dist < 1:    
                t_loop_error = False
            else:
                t_loop_error = True
            if acceptor_dist < 1:
                acceptor_error = False
            else:
                acceptor_error = True
            if mis_count < cur_seq_specs.mis_count:
                cur_seq_specs.length = length
                cur_seq_specs.mis_count = mis_count
                cur_seq_specs.t_loop_error = t_loop_error
                cur_seq_specs.acceptor_error = acceptor_error
                cur_seq_specs.seq = seq
                cur_seq_specs.seq_sub = sub_str
                cur_seq_specs.t_loop_seq = t_loop_seq
                cur_seq_specs.acceptor_seq = acceptor_seq
            if mis_count < 2:
                cur_seq_specs = self.handle_pass_seq(cur_seq_specs, i)
                res_tup = (True, cur_seq_specs)
                return res_tup
        
        # Handles a failed sequence
        if cur_seq_specs.t_loop_error and cur_seq_specs.acceptor_error:
            if length < 24:
                self.stats.short_rejected += 1
            else:
                self.stats.both_rejected += 1
        elif cur_seq_specs.acceptor_error and not cur_seq_specs.t_loop_error:
            self.stats.acceptor_seq_rejected += 1
        elif cur_seq_specs.t_loop_error and not cur_seq_specs.acceptor_error:
            self.stats.t_loop_seq_rejected += 1
        self.stats.total_rejected += 1
        res_tup = (False, cur_seq_specs)
        return res_tup


    def fix_spacing_csv(self):
        """Fixes spacing and indentation on the csv files (_TAB_NO_TRAILER and
        _TAB_TRAILER)
        """
        # Declare and assign variables for output files
        tabfile = open(self.no_trailer_tabfile, "w")
        trailer_tabfile = open(self.trailer_tabfile, "w")
       
        with open("tab_passed", "r") as temp_tabfile:
            temp_tabfile_reader = csv.DictReader(temp_tabfile, delimiter="\t")
            tabfile_writer = csv.DictWriter(tabfile, fieldnames=self.fieldnames,
                delimiter="\t")
            trailer_tabfile_writer = csv.DictWriter(trailer_tabfile,
                fieldnames=self.fieldnames, delimiter="\t")
            tabfile_writer.writeheader()
            trailer_tabfile_writer.writeheader()
            trailer_count = 0

            for row in temp_tabfile_reader:
                row["Seq"] = ("-" * (self.max_seq_width - int(row["Seq_length"]))) + row["Seq"]
                if row["Trailer_length"] == "0":
                    row["3-trailer"] = "-"
                    tabfile_writer.writerow(row)   
                else:
                    trailer_tabfile_writer.writerow(row)
                    trailer_count += 1

        tabfile.close()
        trailer_tabfile.close()
       
        self.stats.num_trailer = trailer_count
        os.remove("tab_passed")


    def write_to_outputs(self, spec_writer):
        """Writes the sort results to output files."""
        self.stats.total_seqs += 1
        is_tRNA_result = self.is_tRNA(self.read_fasta.seq.upper()) 
        mod_id = is_tRNA_result[1].gen_id_string(
            self.read_fasta.id.split('|')[0]) 

        if is_tRNA_result[0]:
            self.passed_seqs_write_fasta.write_id(mod_id)
            self.passed_seqs_write_fasta.write_seq(is_tRNA_result[1].seq)
            is_tRNA_result[1].write_specs(spec_writer, 
                self.read_fasta.id.split('|')[0])

            if len(is_tRNA_result[1].seq) > self.max_seq_width:
                self.max_seq_width = len(is_tRNA_result[1].seq)
        else:
            self.rejected_seqs_write_fasta.write_id(mod_id)
            self.rejected_seqs_write_fasta.write_seq(is_tRNA_result[1].seq)
    

    def write_sorted(self, readfile):
        """Reads in csv files and rewrites them, sorted by seq length."""
        sort_list = []
        with open(readfile, "r") as temp_tabfile:  
            temp_tabfile_reader = csv.DictReader(temp_tabfile, delimiter="\t")
            for row in temp_tabfile_reader:
                sort_list.append(row)

        sort_list.sort(key=lambda dict: dict["Seq_length"], reverse=True)
       
        with open(readfile, "w") as tabfile:
            tabfile_writer = csv.DictWriter(tabfile, fieldnames=self.fieldnames, delimiter="\t")
            tabfile_writer.writeheader()
            tabfile_writer.writerows(sort_list)
    

    def run(self, args):
        """Run the sorter."""
        print "sort started"
        self.set_file_names(args)

        with open("tab_passed", "w") as temp_tabfile:
            spec_writer = csv.DictWriter(temp_tabfile, fieldnames=self.fieldnames, 
                delimiter="\t")
            spec_writer.writeheader()
        
            while self.read_fasta.next():
                self.write_to_outputs(spec_writer)

        self.fix_spacing_csv()

        if args.length_sort:
            self.write_sorted(self.no_trailer_tabfile)
            self.write_sorted(self.trailer_tabfile)
        self.stats.write_stats(self.stats_write_file)
        print "sort finished"

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sort tRNAs")
    parser.add_argument("readfile", help="name of read file")
    parser.add_argument("-n", "--sample_name", help="sample name (to be used" +
        " for naming output files")
    parser.add_argument("-s", "--length_sort", help="sort sequences based on length (excluding trailer)", 
        action="store_true")

    args = parser.parse_args()
    
    if args.length_sort:
        print "length sorting turned on"

    try:
        sorter = Sorter()
        sorter.run(args)
    except:
        exit(1)

