# -*- coding: utf-8
# pylint: disable=line-too-long

import tRNASeqTools.tables as t
import tRNASeqTools.db as db

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class tRNADatabase:
    """To access or create a tRNA Database."""
    
    def __init__(self, db_path, skip_init=False):
        """Initializes variables."""
        self.db = db.DB(db_path)
        self.db_path = db_path

        if not skip_init:
            self.create()

    def create(self):
        """Creates a table in a tRNA database."""
        self.db.create_self()
        self.db.create_table(t.profile_table_name, t.profile_table_structure, t.profile_table_types)
        self.db.create_table(t.stats_table_name, t.stats_table_structure, t.stats_table_types)


    def insert_seq(self, seq_data, id):
        """Insert a seq and its info into the profile table in a tRNA datbase."""
        self.db._exec("""INSERT INTO %s VALUES (%s)""" %
            (t.profile_table_name, (", ".join(['?'] * len(t.profile_table_structure)))),
            seq_data.gen_sql_query_info_tuple(id))


    def insert_stats(self, sorter_stats):
        """Insert stats into stats table in a tRNA database."""
        self.db._exec("""INSERT INTO %s VALUES (%s)""" % \
                        (t.stats_table_name, (", ".join(['?'] *
                         len(t.stats_table_structure)))),
                         sorter_stats.gen_sql_query_info_tuple())


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
