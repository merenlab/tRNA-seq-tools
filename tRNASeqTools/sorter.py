# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to deal with tRNA sequences."""

import os
import sys
import Levenshtein as lev

import tRNASeqTools.fastalib as u
import tRNASeqTools.dbops as dbops
import tRNASeqTools.utils as utils
import tRNASeqTools.terminal as terminal
import tRNASeqTools.extractor as extractor
import tRNASeqTools.filesnpaths as filesnpaths

from tRNASeqTools.errors import ConfigError


__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


pp = terminal.pretty_print


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


    def gen_sql_query_info_tuple(self):
        info_string_list = []
        info_string_list.append(self.total_seqs)
        info_string_list.append(self.total_rejected)
        info_string_list.append(self.total_passed)
        info_string_list.append(self.num_trailer)
        info_string_list.append(self.total_full_length)
        info_string_list.append(self.no_divergence)
        info_string_list.append(self.t_loop_divergence)
        info_string_list.append(self.div_at_0)
        info_string_list.append(self.div_at_1)
        info_string_list.append(self.div_at_2)
        info_string_list.append(self.div_at_3)
        info_string_list.append(self.div_at_8)
        info_string_list.append(self.acceptor_divergence)
        info_string_list.append(self.div_at_neg_1)
        info_string_list.append(self.div_at_neg_2)
        info_string_list.append(self.div_at_neg_3)
        info_string_list.append(self.t_loop_seq_rejected)
        info_string_list.append(self.acceptor_seq_rejected)
        info_string_list.append(self.both_rejected)
        info_string_list.append(self.short_rejected)
        return tuple(info_string_list)


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
        self.anticodon = ""


    def gen_sql_query_info_tuple(self, id):
        """Generates tuple of values to add into database to use in SQL query."""
        info_string_list = []
        info_string_list.append(id)
        info_string_list.append(self.seq)

        if self.trailer_length == 0:
            info_string_list.append(None)
        else:
            info_string_list.append(self.three_trailer)

        info_string_list.append(self.t_loop_seq)
        info_string_list.append(self.acceptor_seq)
        info_string_list.append(str(self.full_length))
        info_string_list.append(str(self.length))
        info_string_list.append(str(self.trailer_length))

        if len(self.anticodon) == 0:
            info_string_list.append(None)
        else:
            info_string_list.append(self.anticodon)

        return tuple(info_string_list)


class Sorter:
    def __init__(self, args):
        """Class that handles the sorting of the seqs."""

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.sample_name = A('sample_name')
        self.input_fasta_path = A('input_fasta')
        self.output_db_path = A('output_db_path')

        self.run = terminal.Run()
        self.progress = terminal.Progress()

        self.sort_stats = SorterStats()
        self.extractor = extractor.Extractor()
        self.db = None
        self.seq_count_dict = {}


    def sanity_check(self):
        """Takes command line arguments from args and assigns input/output file
        variables.
        """

        if not self.sample_name:
            raise ConfigError('You must provide a sample name (which should uniquely describe this profile).')

        if not self.output_db_path:
            raise ConfigError('You must provide an output databaes file path (which should end with exteniosn ".db").')

        if not self.output_db_path.endswith('.db'):
            raise ConfigError('The output database file name must end with ".db".')

        if filesnpaths.is_file_exists(self.output_db_path, dont_raise=True):
            raise ConfigError("The output file already exists. We don't like overwriting stuff here :/")

        if '-' in self.sample_name:
            self.sample_name = self.sample_name.replace('-', '_')
            self.run.warning('I just replaced all "-" characters with "_" characters in your sample name. This program\
                              does not like "-" characters in sample names.')

        utils.check_sample_id(self.sample_name)
        filesnpaths.is_output_file_writable(self.output_db_path)
        filesnpaths.is_file_fasta_formatted(self.input_fasta_path)

        self.input_fasta_path = os.path.abspath(self.input_fasta_path)


    def check_divergence_pos(self, cur_seq_specs):
        """Takes a SeqSpecs class and updates statistics on divergence
        position.
        """
        if cur_seq_specs.t_loop_error:
            self.sort_stats.t_loop_divergence += 1
            if cur_seq_specs.seq_sub[0] != "G":
                self.sort_stats.div_at_0 += 1
            elif cur_seq_specs.seq_sub[1] != "T":
                self.sort_stats.div_at_1 += 1
            elif cur_seq_specs.seq_sub[2] != "T":
                self.sort_stats.div_at_2 += 1
            elif cur_seq_specs.seq_sub[3] != "C":
                self.sort_stats.div_at_3 += 1
            elif cur_seq_specs.seq_sub[8] != "C":
                self.sort_stats.div_at_8 += 1
        elif cur_seq_specs.acceptor_error:
            self.sort_stats.acceptor_divergence += 1
            if cur_seq_specs.seq_sub[-3] != "C":
                self.sort_stats.div_at_neg_3 += 1
            elif cur_seq_specs.seq_sub[-2] != "C":
                self.sort_stats.div_at_neg_2 += 1
            elif cur_seq_specs.seq_sub[-1] != "A":
                self.sort_stats.div_at_neg_1 += 1
        else:
            self.sort_stats.no_divergence += 1


    def check_full_length(self, cur_seq_specs):
        """Takes a SeqSpecs class and checks whether or not the sequence is a
        full-length sequence, returns an updated SeqSpecs class. Also determines
        the anticodon if it is a full-length seq
        """
        if cur_seq_specs.length > 70 and cur_seq_specs.length < 100:
            if cur_seq_specs.seq[7] == "T" and cur_seq_specs.seq[13] == "A":
                self.sort_stats.total_full_length += 1
                cur_seq_specs.full_length = True
        return cur_seq_specs


    def assign_anticodons(self, cur_seq_specs):
        """Takes a SeqSpecs class and tries to assign an anticodon to it"""
        if cur_seq_specs.full_length:
            anticodon = ",".join(self.extractor.extract_anticodon(cur_seq_specs.seq))
            cur_seq_specs.anticodon = anticodon
        else:
            anticodon = ",".join(self.extractor.extract_anticodon_not_full_length(cur_seq_specs.seq))
            cur_seq_specs.anticodon = anticodon
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


    def check_seq_count(self, cur_seq_specs):
        """Checks the sequence in current SeqSpecs to see if the seq has
        already been encountered already, increments count if it has,
        adds it to the dict.
        """
        if cur_seq_specs.seq in self.seq_count_dict:
            self.seq_count_dict[cur_seq_specs.seq] += 1
        else:
            self.seq_count_dict[cur_seq_specs.seq] = 0

        # put the call to a hash function here so we can save the hashed seq in
        # cur_seq_specs
        return cur_seq_specs


    def handle_pass_seq(self, cur_seq_specs, i):
        """Consdolidates all the methods run specifically for a confirmed passed
        sequence.
        """
        self.sort_stats.total_passed += 1
        self.check_divergence_pos(cur_seq_specs)
        cur_seq_specs = self.split_3_trailer(cur_seq_specs, i)
        cur_seq_specs = self.check_full_length(cur_seq_specs)
        cur_seq_specs = self.assign_anticodons(cur_seq_specs)
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
        for i in range(length - sub_size + 1):
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
                self.sort_stats.short_rejected += 1
            else:
                self.sort_stats.both_rejected += 1
        elif cur_seq_specs.acceptor_error and not cur_seq_specs.t_loop_error:
            self.sort_stats.acceptor_seq_rejected += 1
        elif cur_seq_specs.t_loop_error and not cur_seq_specs.acceptor_error:
            self.sort_stats.t_loop_seq_rejected += 1
        self.sort_stats.total_rejected += 1
        res_tup = (False, cur_seq_specs)

        return res_tup


    def gen_sql_query_info_tuple(self):
        info_string_list = []
        info_string_list.append(self.total_seqs)
        info_string_list.append(self.total_rejected)

        return tuple(info_string_list)


    def process(self):
        """Run the sorter."""

        self.sanity_check()

        # creating an empty profile databsae
        profile_db = dbops.tRNADatabase(self.output_db_path)
        profile_db.create(meta_values={'sample_name': self.sample_name})

        # a list buffer to keep results
        results_buffer = []

        # an arbitrary max size to store and reset the buffer
        memory_max = 2000000

        input_fasta = u.SequenceSource(self.input_fasta_path)

        self.run.info('Hi', terminal.get_date(), mc='green')
        self.run.info('Sample name', self.sample_name)
        self.run.info('Input FASTA', self.input_fasta_path)

        table_for_tRNA_seqs = dbops.TableFortRNASequences(self.output_db_path)

        self.progress.new('Profiling tRNAs')
        self.progress.update('...')
        while next(input_fasta):
            self.sort_stats.total_seqs += 1
            is_tRNA_result = self.is_tRNA(input_fasta.seq.upper())

            if is_tRNA_result[0]:
                results_buffer.append(('%s_%d' % (self.sample_name, input_fasta.pos), is_tRNA_result[1]))

            if sys.getsizeof(results_buffer) > memory_max:
                self.progress.update('Writing %d items in the buffer to the DB ...' % len(results_buffer))
                table_for_tRNA_seqs.append_sequences(results_buffer)
                results_buffer = []

            if self.sort_stats.total_seqs % 1000 == 0:
                t, p = self.sort_stats.total_seqs, self.sort_stats.total_passed
                self.progress.update('%s :: %s (num tRNAs :: num raw reads so far): %.2f%% ...' %\
                                        (pp(p), pp(t), p * 100 / t))

        self.progress.update('Writing %d items in the buffer to the DB ...' % len(results_buffer))
        table_for_tRNA_seqs.append_sequences(results_buffer)
        results_buffer = []

        self.progress.update('Writing stats ...')
        dbops.TableFortRNAProfilingStats(self.output_db_path).insert(self.sort_stats)

        self.progress.end()

        self.run.info('Total raw seqs processed', self.sort_stats.total_seqs)
        self.run.info('Total tRNA seqs recovered', self.sort_stats.total_passed)
        self.run.info('Total full length tRNA seqs', self.sort_stats.total_full_length)
        self.run.info('Output DB path', self.output_db_path)
        self.run.info('Bye', terminal.get_date(), mc='green')

        self.run.quit()
