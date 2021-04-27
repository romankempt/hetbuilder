"""hetbuilder"""

import sys
from distutils.version import LooseVersion
from pathlib import Path

if sys.version_info[0] == 2:
    raise ImportError("Requires Python3. This is Python2.")

__version__ = "0.5.8"

PROJECT_ROOT_DIR = Path(__file__).absolute().parent

from hetbuilder.algorithm import Interface, CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

__all__ = ["__version__", "Interface", "CoincidenceAlgorithm", "InteractivePlot"]

