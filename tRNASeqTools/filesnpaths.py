# -*- coding: utf-8
# pylint: disable=line-too-long
"""File/Path operations"""

import os
import mmap
import time
import shutil
import tempfile

import tRNASeqTools.fastalib as u
from tRNASeqTools.terminal import Run
from tRNASeqTools.terminal import Progress
from tRNASeqTools.errors import FilesNPathsError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, Meren Lab"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def is_file_exists(file_path, dont_raise=False):
    if not file_path:
        raise FilesNPathsError("No input file is declared...")
    if not os.path.exists(file_path):
        if dont_raise:
            return False
        else:
            raise FilesNPathsError("No such file: '%s' :/" % file_path)
    return True


def is_output_file_writable(file_path):
    if not file_path:
        raise FilesNPathsError("No output file is declared...")
    if not os.access(os.path.dirname(os.path.abspath(file_path)), os.W_OK):
        raise FilesNPathsError("You do not have permission to generate the output file '%s'" % file_path)
    return True


def is_output_dir_writable(dir_path):
    if not dir_path:
        raise FilesNPathsError("No output directory path is declared...")
    if not os.path.isdir(dir_path):
        raise FilesNPathsError("'%s' is not a directory..." % dir_path)
    if not os.access(os.path.abspath(dir_path), os.W_OK):
        raise FilesNPathsError("You do not have permission to generate files in '%s'" % dir_path)
    return True


def is_dir_empty(dir_path):
    if not dir_path:
        raise FilesNPathsError("is_dir_empty: No directory path is declared...")
    if not os.path.isdir(dir_path):
        raise FilesNPathsError("is_dir_empty: '%s' is not a directory..." % dir_path)

    return False if len(os.listdir(dir_path)) else True


def is_file_tab_delimited(file_path, separator='\t', expected_number_of_fields=None):
    is_file_exists(file_path)
    f = open(file_path, 'rU')

    while True:
        line = f.readline().strip(' ')
        if line.startswith('#'):
            continue
        else:
            break

    if len(line.split(separator)) == 1 and expected_number_of_fields != 1:
        raise FilesNPathsError("File '%s' does not seem to have TAB characters.\
                            Did you export this file on MAC using EXCEL? :(" % file_path)

    f.seek(0)
    num_fields_set = set([len(line.split(separator)) for line in f.readlines()])
    if len(num_fields_set) != 1:
        raise FilesNPathsError("Not all lines in the file '%s' have equal number of fields..." % file_path)

    if expected_number_of_fields:
        num_fields_in_file = list(num_fields_set)[0]
        if num_fields_in_file != expected_number_of_fields:
            raise FilesNPathsError("The expected number of fileds for '%s' is %d. Yet, it has %d\
                                     of them :/" % (file_path, expected_number_of_fields, num_fields_in_file))

    f.close()
    return True


def is_file_fasta_formatted(file_path):
    is_file_exists(file_path)

    try:
        f = u.SequenceSource(file_path)
    except u.FastaLibError as e:
        raise FilesNPathsError("Someone is not happy with your FASTA file '%s' (this is\
                            what the lib says: '%s'." % (file_path, e))

    f.close()

    return True


def is_program_exists(program):
    """adapted from http://stackoverflow.com/a/377028"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    raise FilesNPathsError("'%s' is not found" % program)


def get_temp_directory_path():
    return tempfile.mkdtemp()


def get_temp_file_path():
    f = tempfile.NamedTemporaryFile(delete=False)
    temp_file_name = f.name
    f.close()
    return temp_file_name


def get_num_lines_in_file(file_path):
    if os.stat(file_path).st_size == 0:
        return 0

    f = open(file_path, "rU+")
    buf = mmap.mmap(f.fileno(), 0)
    num_lines = 0
    readline = buf.readline
    while readline():
        num_lines += 1

    return num_lines


def check_output_directory(output_directory, ok_if_exists=False):
    if not output_directory:
        raise FilesNPathsError("Sorry. You must declare an output directory path.")

    output_directory = os.path.abspath(output_directory)

    if os.path.exists(output_directory) and not ok_if_exists:
        raise FilesNPathsError("The output directory already exists. We do not like overwriting stuff.")

    return output_directory


def gen_output_directory(output_directory, progress=Progress(verbose=False), run=Run(), delete_if_exists=False):
    if os.path.exists(output_directory) and delete_if_exists and not is_dir_empty(output_directory):
        try:
            run.warning('filesnpaths::gen_output_directory: the client asked the existing directory \
                         "%s" to be removed.. Just so you know :/ (You have 5 seconds to press\
                         CTRL + C).' % output_directory)
            time.sleep(5)
            shutil.rmtree(output_directory)
        except:
            progress.end()
            raise FilesNPathsError("I was instructed to remove this directory, but I failed: '%s' :/" % output_directory)

    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except:
            progress.end()
            raise FilesNPathsError("Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                            (output_directory))
    if not os.access(output_directory, os.W_OK):
        progress.end()
        raise FilesNPathsError("You do not have write permission for the output directory: '%s'" % output_directory)

    return output_directory


def get_name_from_file_path(file_path, postfix_separator="."):
    """Return a decent name for a given file at file at `file_path`"""

    splits = os.path.basename(os.path.abspath(file_path)).split(postfix_separator)

    return splits[0] if len(splits) in [1, 2] else '.'.join(splits[:-1])
