# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes to extract info from tRNA sequences."""

import os 
import csv
import Levenshtein as lev


__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class ExtractorStats:
    """This class handles keeping track of extraction statistics."""

    def __init__(self):
        """Initializes statistics."""
        
