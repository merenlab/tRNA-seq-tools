# -*- coding: utf-8
# pylint: disable=line-too-long

import sys

# Make sure the Python environment hasn't changed since the installation
try:
    if sys.version_info.major != 3:
        sys.stderr.write("tRNA-seq-tools expects to be run in a Python 3 environment, howeever, your active Python major version is %d.\n" % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(tRNA-seq-tools failed to learn about your Python version, it will pretend as if nothing happened, but this is so not normal)\n\n")
