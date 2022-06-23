""" Test functions for local testing, by copying the shared library hetbuilder_backend.so to hetbuilder/ """

from dataclasses import dataclass
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

from pathlib import Path

import ase.io
from ase.build import mx2

import numpy as np

from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from ase.build import make_supercell
from ase.atoms import Atoms

from hetbuilder.algorithm import ase_atoms_to_cpp_atoms
from timeit import default_timer as time
from itertools import islice

from random import randint


def test_algorithm():

    graphene_cell = cellpar_to_cell(np.array([2.46, 2.46, 100.0, 90.0, 90.0, 120]))
    graphene_positions = np.array([[0.0, 1.42028166, -1.6775], [0.0, 0.0, -1.6775]])
    atoms1 = Atoms(
        symbols=["C", "C"], positions=graphene_positions, cell=graphene_cell, pbc=True
    )

    atoms2 = mx2("WS2")
    atoms2.pbc = True
    atoms2.cell[2, 2] = 100

    alg = CoincidenceAlgorithm(atoms1, atoms2)
    results = alg.run(
        tolerance=0.1, Nmax=20, angle_limits=(0, 30), angle_stepsize=1, verbosity=2
    )
    if results is not None:
        ip = InteractivePlot(atoms1, atoms2, results, 0.5)
        ip.plot_results()
    else:
        print("nope")


def test_scaling_ase(M=4, N=5):
    atoms1 = mx2("MoS2")
    atoms1.pbc = True
    atoms1.cell[2, 2] = 100
    atoms2 = mx2("MoS2")
    atoms2.pbc = True
    atoms2.cell[2, 2] = 100
    atoms2.set_cell(atoms2.cell * 1.1, scale_atoms=True)

    for m in range(1, M):
        atoms = make_supercell(atoms1, m * np.eye(3))

        atoms_list = [atoms] * N

        for j in atoms_list:
            if randint(0, 10) > 0:
                j.rotate(randint(0, 90), "z", rotate_cell=True)
                j.rattle()

        atoms_list[randint(0, len(atoms1) - 1)] = atoms_list[0]
        atoms_list[randint(0, len(atoms1) - 1)].translate([0, 0, 100])
        atoms_list[randint(0, len(atoms1) - 1)] = atoms2

        def compare(a1, other):
            # other can be list
            comp = SymmetryEquivalenceCheck()
            return comp.compare(a1, other)

        def del_dups_ase(lst):
            """O(n**2) algorithm, O(1) in memory"""
            pos = 0
            for item in lst:
                if not compare(item, islice(lst, pos)):
                    # we haven't seen `item` yet
                    lst[pos] = item
                    pos += 1
            del lst[pos:]

        ase_timings = []
        for k in range(10):
            t1 = time()
            del_dups_ase(atoms_list)
            t2 = time()
            ase_timings.append(t2 - t1)

        ase_time = np.average(ase_timings) / len(ase_timings) / 1000
        print(f"M={m} N={N} \t ASE {ase_time:e} ms")


def test_scaling_cpp(M=4, N=5):

    atoms1 = mx2("MoS2")
    atoms1.pbc = True
    atoms1.cell[2, 2] = 100
    atoms2 = mx2("MoS2")
    atoms2.pbc = True
    atoms2.cell[2, 2] = 100
    atoms2.set_cell(atoms2.cell * 1.1, scale_atoms=True)

    for m in range(1, M):
        atoms = make_supercell(atoms1, m * np.eye(3))

        atoms_list = [atoms] * N

        for j in atoms_list:
            if randint(0, 10) > 0:
                j.rotate(randint(0, 90), "z", rotate_cell=True)
                j.rattle()

        atoms_list[randint(0, len(atoms1) - 1)] = atoms_list[0]
        atoms_list[randint(0, len(atoms1) - 1)].translate([0, 0, 100])
        atoms_list[randint(0, len(atoms1) - 1)] = atoms2

        cpp1 = [ase_atoms_to_cpp_atoms(j) for j in atoms_list]

        def del_dups_cpp(lst):
            """O(n**2) algorithm, O(1) in memory"""
            pos = 0
            for item in lst:
                if not all([(item.compare(item2)) for item2 in islice(lst, pos)]):
                    # we haven't seen `item` yet
                    lst[pos] = item
                    pos += 1
            del lst[pos:]

        cpp_timings = []
        for k in range(10):
            t1 = time()
            del_dups_cpp(cpp1)
            t2 = time()
            cpp_timings.append(t2 - t1)

        cpp_time = np.average(cpp_timings) / len(cpp_timings) / 1000
        print(f"M={m} N={N} \t CPP {cpp_time:e} ms")


if __name__ == "__main__":
    test_algorithm()
