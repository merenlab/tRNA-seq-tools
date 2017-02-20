# -*- coding: utf-8
# pylint: disable=line-too-long

import os
import sys
import copy
import pkg_resources

# Make sure the Python environment hasn't changed since the installation
try:
    if sys.version_info.major != 3:
        sys.stderr.write("tRNA-seq-tools expects to be run in a Python 3 environment, howeever, your active Python major version is %d.\n" % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(tRNA-seq-tools failed to learn about your Python version, it will pretend as if nothing happened, but this is so not normal)\n\n")


# a dict of communtiy parameters and flags
D = {
    'profile-db': (
            ['-p', '--profile-db'],
            {'metavar': "PROFILE_DB",
             'required': True,
             'help': "tRNA-seq-tools profile database"}
                ),
}

# two functions that works with the dictionary above.
def A(param_id):
    return D[param_id][0]

def K(param_id, params_dict={}):
    kwargs = copy.deepcopy(D[param_id][1])
    for key in params_dict:
        kwargs[key] = params_dict[key]

    return kwargs


import tRNASeqTools.tables as t
import tRNASeqTools.terminal as terminal

run = terminal.Run()

def set_version():
    tRNASeqTools_version = 'unknown'

    try:
        tRNASeqTools_version = pkg_resources.require("tRNASeqTools")[0].version
    except:
        # maybe anvi'o is not installed but it is being run from the codebase dir?
        # some hacky stuff to get version from the setup.py
        try:
            setup_py_path = os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/../setup.py'
            version_string = [l.strip() for l in open(setup_py_path).readlines() if l.strip().startswith('tRNASeqTools_version')][0]
            tRNASeqTools_version = version_string.split('=')[1].strip().strip("'").strip('"')
        except:
            pass

    return tRNASeqTools_version, \
           t.profile_db_version


def print_version():
    run.info("tRNASeqTools version", __version__, mc='green')
    run.info("Profile DB version", __profile__version__)


__version__, \
__profile__version__ = set_version()


if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()
