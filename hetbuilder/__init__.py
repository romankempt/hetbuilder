"""Heterostructure interface builder."""

import sys
from distutils.version import LooseVersion

if sys.version_info[0] == 2:
    raise ImportError("Requires Python3. This is Python2.")

__version__ = "0.2.0"


from hetbuilder.algorithm import *
