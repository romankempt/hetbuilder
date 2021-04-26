import ase.io
from hetbuilder.algorithm import *
from hetbuilder.plotting import InteractivePlot

bottom = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/MoS2_2H_1l.xyz"
)

top = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/WS2_2H_1l.xyz"
)

alg = CoincidenceAlgorithm(bottom, top)
results = alg.run(Nmin=-15, Nmax=15, angles=[0, 10, 11, 12, 13, 14, 15], tolerance=0.1)

if results is not None:
    for s in results:
        print(s.angle, s.M)
    ipl = InteractivePlot(bottom, top, results, weight=0.5)
    ipl.plot_results()
