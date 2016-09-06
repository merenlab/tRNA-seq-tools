# -*- coding: utf-8
# pylint: disable=line-too-long

import tables as t
import db

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class DB_Accessor:
    """This class handles accessing the tRNA database."""

    def __init__(self, db_path):
        """Initializes access to tRNA database."""
        self.db = db.DB(db_path)
        self.db_path = db_path

    #def write_to_file(self, dict, 

    def run(self):
        #db_dict = self.db.get_table_as_dict(t.tRNA_profiling_table_name)
        #print len(db_dict)
        #print db_dict["test_657"]

        db_dict_no_trailer = self.db.get_some_rows_from_table_as_dict(t.tRNA_profiling_table_name,
            "three_trailer IS NULL")
        print len(db_dict_no_trailer)
        print db_dict_no_trailer["test_657"]


        db_dict_trailer = self.db.get_some_rows_from_table_as_dict(t.tRNA_profiling_table_name,
            "three_trailer IS NOT NULL")
        print len(db_dict_trailer)
        print db_dict_trailer["test_1823"]
            
