import ase.io
from pathlib import Path

from hetbuilder import PROJECT_ROOT_DIR
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

import numpy as np

import time


def test_performance():
    bottom = ase.io.read(PROJECT_ROOT_DIR.joinpath("../tests/MoS2_2H_1l.xyz"))
    top = ase.io.read(PROJECT_ROOT_DIR.joinpath("../tests/WS2_2H_1l.xyz"))
    alg = CoincidenceAlgorithm(bottom, top)
    timings = []
    ncombs = []
    for i in [5, 10, 15]:
        subtimings = []
        ncombinations = ((2 * (i)) ** 4) * 1
        ncombs.append(ncombinations)
        for j in range(5):
            start = time.time()
            results = alg.run(Nmin=0, Nmax=i, angles=[0], tolerance=0.1)
            end = time.time()
            subtimings.append(end - start)
        subtime = sum(subtimings) / len(subtimings)
        timings.append(subtime)
    return (timings, ncombs)

