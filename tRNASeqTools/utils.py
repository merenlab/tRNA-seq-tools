# -*- coding: utf-8
# pylint: disable=line-too-long

"""Lonely, helper functions that are broadly used and don't fit anywhere"""

import os
import string
import textwrap

import tRNASeqTools.filesnpaths as filesnpaths

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


def store_dict_as_TAB_delimited_file(d, output_path, headers=None, file_obj=None):
    if not file_obj and not os.access(os.path.dirname(os.path.abspath(output_path)), os.W_OK):
        raise ConfigError("The output file path '%s' is not writable..." % output_path)

    if not file_obj:
        f = open(output_path, 'w')
    else:
        f = file_obj

    if not headers:
        headers = ['key'] + sorted(list(d.values())[0].keys())

    f.write('%s\n' % '\t'.join(headers))

    for k in sorted(d.keys()):
        line = [str(k)]
        for header in headers[1:]:
            try:
                val = d[k][header]
            except KeyError:
                raise ConfigError("Header ('%s') is not found in the dict :/" % (header))
            except TypeError:
                raise ConfigError("Your dictionary is not properly formatted to be exported\
                                    as a TAB-delimited file :/ You ask for '%s', but it is not\
                                    even a key in the dictionary" % (header))

            line.append(str(val) if not isinstance(val, type(None)) else '')

        f.write('%s\n' % '\t'.join(line))

    f.close()
    return output_path


def store_dict_as_FASTA_file(d, output_file_path, report_unique_sequences=False, wrap_from=200):
    filesnpaths.is_output_file_writable(output_file_path)
    output = open(output_file_path, 'w')

    props_h = sorted(list(list(d.values())[0]['props'].keys()))

    if report_unique_sequences:
        seqs_sorted_by_frequency = [x[1] for x in sorted([(len(d[seq]['ids']), seq) for seq in d], reverse=True)]

        for seq in seqs_sorted_by_frequency:
            frequency = len(d[seq]['ids'])
            seq_id = d[seq]['ids'].pop()
            output.write('>%s %s|frequency:%d\n' % (seq_id, '|'.join(['%s:%s' % (t[0], t[1]) for t in d[seq]['props'].items()]), frequency))
            output.write('%s\n' % textwrap.fill(seq, wrap_from, break_on_hyphens=False))
    else:
        for seq in d:
            props = '|'.join(['%s:%s' % (t[0], t[1]) for t in d[seq]['props'].items()])
            for seq_id in d[seq]['ids']:
                output.write('>%s %s\n' % (seq_id, props))
                output.write('%s\n' % textwrap.fill(seq, wrap_from, break_on_hyphens=False))

    output.close()
    return True


