# -*- coding: utf-8
# pylint: disable=line-too-long

"""Lonely, helper functions that are broadly used and don't fit anywhere"""


import string

from tRNASeqTools.errors import ConfigError

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, Meren Lab"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


digits = string.digits
allowed_chars = string.ascii_letters + string.digits + '_' + '-' + '.'

def check_sample_id(sample_id):
    if sample_id:
        if sample_id[0] in digits:
            raise ConfigError("Sample names can't start with digits. Long story. Please specify a sample name\
                                that starts with an ASCII letter (you may want to check '-s' parameter to set\
                                a sample name if your client permits (otherwise you are going to have to edit\
                                your input files)).")

        allowed_chars_for_samples = allowed_chars.replace('-', '').replace('.', '')
        if len([c for c in sample_id if c not in allowed_chars_for_samples]):
            raise ConfigError("Sample name ('%s') contains characters that anvio does not like. Please\
                                limit the characters that make up the project name to ASCII letters,\
                                digits, and the underscore character ('_')." % sample_id)


