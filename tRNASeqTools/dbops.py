# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import time

import tRNASeqTools.db as db
import tRNASeqTools.tables as t
import tRNASeqTools.terminal as terminal


__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class tRNADatabase:
    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            self.db = db.DB(self.db_path, t.__profile_db__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            self.run.info('tRNA Profile database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Sample name', self.meta['sample_name'], quiet=self.quiet)
        else:
            self.db = None


    def create(self, meta_values={}):
        self.db = db.DB(self.db_path, t.__profile_db__version__, new_database=True)

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # creating empty default tables
        self.db.create_table(t.profile_table_name, t.profile_table_structure, t.profile_table_types)
        self.db.create_table(t.stats_table_name, t.stats_table_structure, t.stats_table_types)

        self.disconnect()

        self.run.info('Profile database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


    def gen_anticodon_profile(self, only_full_length, min_seq_length, 
        max_seq_length, anticodons): 
        """returns an anticodon profile from the database"""

        anticodon_count_dict = {}

        where_clause = "Three_trailer IS NULL AND Anticodon IS NOT NULL"
        if only_full_length:
            where_clause += " AND Full_length = 'True'"
        if min_seq_length:
            where_clause += (" AND Seq_length >= " + str(min_seq_length))
        if max_seq_length:
            where_clause += (" AND Seq_length <= " + str(max_seq_length))
        if anticodons:
            spec_anticodons_list = anticodons.split(",")
            or_string = " OR ".join(["Anticodon IN ('" + spec_anticodon + "')"
                for spec_anticodon in spec_anticodons_list])
            where_clause += " AND (" + or_string + ")"

        profile_dict_no_trailer = self.db.get_some_rows_from_table_as_dict(t.profile_table_name,
            where_clause)


        for key in list(profile_dict_no_trailer.keys()):
            anticodon_list = profile_dict_no_trailer[key]["Anticodon"].split(",")
                
            for anticodon in anticodon_list:
                if anticodon in anticodon_count_dict:
                    anticodon_count_dict[anticodon] += 1
                else:
                    anticodon_count_dict[anticodon] = 1
        return anticodon_count_dict


class TableFortRNASequences:
    """A class to populate the profile databse with tRNA results"""
    
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run


    def append_sequences(self, tRNA_profile_list):
        """Insert a seq and its info into the profile table in a tRNA datbase."""

        values = []
        for sequence_id, sequence_object in tRNA_profile_list:
            values.append(sequence_object.gen_sql_query_info_tuple(sequence_id))

        profile_db = tRNADatabase(self.db_path)
        profile_db.db._exec_many("""INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?)""" % (t.profile_table_name), values)
        profile_db.disconnect()


class TableFortRNAProfilingStats:
    """A class to populate and process stats generated during the tRNA profiling"""
    
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress


    def insert(self, sorter_stats):
        """Insert a seq and its info into the profile table in a tRNA datbase."""

        profile_db = tRNADatabase(self.db_path)
        profile_db.db._exec("""INSERT INTO %s VALUES (%s)""" % \
                        (t.stats_table_name, (", ".join(['?'] *
                         len(t.stats_table_structure)))),
                         sorter_stats.gen_sql_query_info_tuple())

        profile_db.disconnect()
