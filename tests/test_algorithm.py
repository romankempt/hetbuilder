import ase.io
from hetbuilder.algorithm import *

atoms = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/graphene.xyz"
)

alg = CoincidenceAlgorithm(atoms, atoms)
results = alg.run(Nmax=3, Nmin=0, angles=[0])

from ase.visualize import view

for s in results:
    view(s.stack)

