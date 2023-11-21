import ase.io
from pathlib import Path

from hetbuilder import PROJECT_ROOT_DIR
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot
from hetbuilder.log import set_verbosity_level

import numpy as np

import time
import os

# # we set up the algorithm class
# alg = CoincidenceAlgorithm(bottom, top)
# # we run the algorithm for a choice of parameters
# results = alg.run(Nmax = 12, Nmin = 0, angles = [0, 10, 15, 20, 25, 30], tolerance = 0.2, weight= 0.5)

# iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=0.5)
# iplot.plot_results()

