# -*- coding: utf-8
# pylint: disable=line-too-long

"""Exceptions"""

import textwrap

from tRNASeqTools.ttycolors import color_text

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, Meren Lab"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def remove_spaces(text):
    while True:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


class tRNASeqToolsError(Exception, object):
    def __init__(self, e=None):
        Exception.__init__(self)
        return

    def __str__(self):
        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        error_lines = ['%s%s' % (l, ' ' * (max_len - len(l))) for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s: %s' % (color_text(self.error_type, 'red'), error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(self.error_type) + 2), error_line))

        return '\n\n' + '\n'.join(error_message) + '\n\n'

    def clear_text(self):
        return '%s: %s' % (self.error_type, self.e)


class ConfigError(tRNASeqToolsError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Config Error'
        tRNASeqToolsError.__init__(self)


class TerminalError(tRNASeqToolsError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'Terminal Error'
        tRNASeqToolsError.__init__(self)


class FilesNPathsError(tRNASeqToolsError):
    def __init__(self, e=None):
        self.e = remove_spaces(e)
        self.error_type = 'File/Path Error'
        tRNASeqToolsError.__init__(self)
