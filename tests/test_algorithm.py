import ase.io
from pathlib import Path

from hetbuilder import *
from hetbuilder.algorithm import *
from hetbuilder.plotting import InteractivePlot

import numpy as np


def test_coincidence_algorithm():

    bottom = ase.io.read(PROJECT_ROOT_DIR.joinpath("../tests/MoS2_2H_1l.xyz"))

    top = ase.io.read(PROJECT_ROOT_DIR.joinpath("../tests/WS2_2H_1l.xyz"))

    alg = CoincidenceAlgorithm(bottom, top)
    results = alg.run(Nmin=-5, Nmax=5, angles=[0, 10, 20, 30], tolerance=0.1)
    assert results is not None, "Found no results at all."
    assert len(results) == 1, "For these settings there should only be one result."
    example = results[0]
    assert np.array(example.M).sum() == 3, "M matrix should be unity."
