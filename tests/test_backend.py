import ase.io
from ase.atoms import Atoms
from ase.build import make_supercell

from pathlib import Path

from hetbuilder import PROJECT_ROOT_DIR
from hetbuilder.algorithm import cpp_atoms_to_ase_atoms, ase_atoms_to_cpp_atoms

from hetbuilder_backend import (
    double2dVector,
    double1dVector,
    int1dVector,
    int2dVector,
    CppAtomsClass,
    CppCoincidenceAlgorithmClass,
    CppInterfaceClass,
    get_number_of_omp_threads,
    cpp_make_supercell,
)

import spglib

from ase.utils.structure_comparator import SymmetryEquivalenceCheck


def test_backend_supercell():
    for i in [
        "../tests/MoS2_2H_1l.xyz",
        "../tests/WS2_2H_1l.xyz",
        "../tests/graphene.xyz",
    ]:
        atoms = ase.io.read(PROJECT_ROOT_DIR.joinpath(i))
        cppatoms = ase_atoms_to_cpp_atoms(atoms)
        N = [[3, 1, 0], [-1, 2, 0], [0, 0, 1]]
        atoms = make_supercell(atoms, N)

        N = int2dVector([int1dVector(k) for k in N])
        cppatoms = cpp_make_supercell(cppatoms, N)
        cppatoms = cpp_atoms_to_ase_atoms(cppatoms)
        comp = SymmetryEquivalenceCheck()
        is_equal = comp.compare(atoms, cppatoms)
        assert (
            is_equal
        ), "Supercell generation in backend and from ASE do not yield same result."


def test_spglib_standardize():
    for i in [
        "../tests/MoS2_2H_1l.xyz",
        "../tests/WS2_2H_1l.xyz",
        "../tests/graphene.xyz",
    ]:
        atoms = ase.io.read(PROJECT_ROOT_DIR.joinpath(i))
        N = [[3, 1, 0], [-1, 2, 0], [0, 0, 1]]
        sc = make_supercell(atoms, N)
        cppatoms = ase_atoms_to_cpp_atoms(sc)
        cppatoms.standardize(1, 0, 1e-5, 5)

        cell = (sc.cell, sc.get_scaled_positions(), sc.numbers)
        spgcell = spglib.standardize_cell(
            cell, to_primitive=True, no_idealize=False, symprec=1e-5, angle_tolerance=5
        )
        atoms = Atoms(
            cell=spgcell[0],
            scaled_positions=spgcell[1],
            numbers=spgcell[2],
            pbc=[True, True, True],
        )
        cppatoms = cpp_atoms_to_ase_atoms(cppatoms)
        comp = SymmetryEquivalenceCheck()
        is_equal = comp.compare(atoms, cppatoms)
        assert (
            is_equal
        ), "Standardization in backend and from spglib do not yield same result."
