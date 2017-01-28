# -*- coding: utf-8
# pylint: disable=line-too-long
"""table schemas for database"""

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


###############################################################################
#
#   TABLE DESCRIPTIONS
#
###############################################################################


# table schemas for tRNA profiling table
profile_table_name      = "profile"
profile_table_structure = ["ID",   "Seq" , "Three_trailer", "T_loop", "Acceptor", "Full_length", "Seq_length", "Trailer_length", "Anticodon"]
profile_table_types     = ["text", "text",     "text"     ,  "text" ,   "text"  ,     "text"   ,    "int"    ,      "int"      ,   "text"   ]

# table schemas for stats table
stats_table_name = "stats"
stats_table_structure = ["Total_seqs", "Total_full_length", "With_trailer",
    "Total_passed", "No_divergence", "T_loop_divergence", "Pos_0_divergence",
    "Pos_1_divergence", "Pos_2_divergence", "Pos_3_divergence", 
    "Pos_8_divergence", "Acceptor_divergence", "Pos_neg_3_divergence", 
    "Pos_neg_2_divergence", "Pos_neg_1_divergence", "Total_failed", "T_loop_reject",
    "Acceptor_reject", "Both_reject", "Short_reject"]
stats_table_types = ["int"] * len(stats_table_structure)
