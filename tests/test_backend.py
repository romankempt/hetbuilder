import ase.io
from ase.atoms import Atoms

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

from hetbuilder.algorithm import *
from ase.build import make_supercell

atoms = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/graphene.xyz"
)

N = [[3, 1, 0], [-1, 2, 0], [0, 0, 1]]
atoms = make_supercell(atoms, N)
cppatoms = ase_atoms_to_cpp_atoms(atoms)

import spglib

cell = (atoms.cell, atoms.get_scaled_positions(), atoms.numbers)
spgcell = spglib.standardize_cell(
    cell, to_primitive=True, no_idealize=False, symprec=1e-5, angle_tolerance=5
)
atoms = Atoms(
    cell=spgcell[0],
    scaled_positions=spgcell[1],
    numbers=spgcell[2],
    pbc=[True, True, True],
)
cppatoms.standardize(1, 0, 1e-5, 5)
cppatoms = cpp_atoms_to_ase_atoms(cppatoms)

# newcell = atoms.cell * 2
# N = double2dVector([double1dVector(k) for k in newcell])
# cppatoms.scale_cell(N)


# atoms.set_cell(newcell, scale_atoms=True)
