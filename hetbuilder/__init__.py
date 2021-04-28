"""hetbuilder"""

import sys
from distutils.version import LooseVersion
from pathlib import Path

if sys.version_info[0] == 2:
    raise ImportError("Requires Python3. This is Python2.")

__version__ = "0.6.0"

PROJECT_ROOT_DIR = Path(__file__).absolute().parent

from hetbuilder.algorithm import Interface, CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot
from hetbuilder_backend import (
    double2dVector,
    double1dVector,
    int1dVector,
    int2dVector,
    CppAtomsClass,
    CppCoincidenceAlgorithmClass,
    CppInterfaceClass,
    get_number_of_omp_threads,
)

__all__ = ["__version__", "Interface", "CoincidenceAlgorithm", "InteractivePlot"]

