# -*- coding: utf-8
# pylint: disable=line-too-long
"""table schemas for database"""


__author__ = "Steven Cui"
__copyright__ = "Copyright 2017, Meren Lab"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.erengmail.com"


profile_db_version = 1


###############################################################################
#
#   TABLE DESCRIPTIONS
#
###############################################################################


# table schemas for tRNA profiling table
profile_table_name      = "profile"
profile_table_structure = ["ID",   "Seq" , "Three_trailer", "T_loop", "Acceptor", "Full_length", "Seq_length", "Trailer_length", "Anticodon"]
profile_table_types     = ["text", "text",     "text"     ,  "text" ,   "text"  ,     "text"   ,    "int"    ,      "int"      ,   "text"   ]
