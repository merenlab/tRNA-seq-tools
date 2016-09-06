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
tRNA_profiling_table_name = "tRNA_profiling"
tRNA_profiling_table_structure = ["ID", "Seq", "three_trailer", "t_loop",
    "acceptor", "full_length", "Seq_length", "Trailer_length", "Anticodon"]
tRNA_profiling_table_types = ["text", "text", "text", "text", "text", "text",
    "int", "int", "text"]
