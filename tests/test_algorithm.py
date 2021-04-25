import ase.io
from hetbuilder.algorithm import *

bottom = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/MoS2_2H_1l.xyz"
)

top = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/WS2_2H_1l.xyz"
)

alg = CoincidenceAlgorithm(bottom, top)
results = alg.run(Nmin=-10, Nmax=10, angles=[0, 5, 10, 15, 20, 25, 30], tolerance=0.2)

from ase.visualize import view

if results is not None:
    for s in results:
        print(s.M)
        print(s.angle)

    view(results[1].stack)
