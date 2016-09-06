# -*- coding: utf-8
# pylint: disable=line-too-long


import os
import sys
import sqlite3
import tables as t
import db

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class tRNADatabase:
    """To access or create a tRNA Database."""
    
    def __init__(self, db_path):
        """Initializes variables."""
        self.db = db.DB(db_path)
        self.db_path = db_path
        self.table_name = t.tRNA_profiling_table_name
        self.table_structure = t.tRNA_profiling_table_structure
        self.table_types = t.tRNA_profiling_table_types
        
        self.create()


    def create(self):
        """Creates a table in a tRNA database."""
        self.db.create_table(self.table_name, self.table_structure, 
            self.table_types)

    def insert_seq(self, seq_data, id):
        """Insert a seq and its info into the table in a tRNA datbase."""
        self.db._exec("""INSERT INTO %s VALUES (%s)""" %
            (self.table_name, (", ".join(['?'] * len(self.table_structure)))),
            seq_data.gen_sql_query_info_tuple(id))
