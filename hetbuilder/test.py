import ase.io
from hetbuilder.algorithm import *

atoms = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/graphene.xyz"
)

alg = CoincidenceAlgorithm(atoms, atoms)
results = alg.run(Nmax=3, Nmin=0, angles=[0])

example = results[0]
one = example.bottom
two = example.top
stack = example.stack
print(results[0].angle)

import numpy as np

# N = np.array([[-5, -4, 0], [-1, -2, 0], [0, 0, 1]])

# M = np.array([[k for k in l] for l in example.M])


from ase.visualize import view
from ase.build.supercells import lattice_points_in_supercell, clean_matrix


# test1 = make_supercell(atoms, N, wrap=False)
# test2 = cpp_make_supercell(bottom, Nc)
# test2 = cpp_atoms_to_ase_atoms(test2)

# view(test1)
# view(test2)
# view(one)
# view(two)

