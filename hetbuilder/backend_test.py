from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

from pathlib import Path

import ase.io
from ase.build import mx2

import numpy as np

from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from hetbuilder.algorithm import ase_atoms_to_cpp_atoms
from timeit import default_timer as time
from itertools import islice
from ase.build import make_supercell
from random import randint


def test_algorithm():

    atoms1 = mx2("MoS2")
    atoms1.pbc = True
    atoms1.cell[2, 2] = 100

    atoms2 = mx2("MoS2")
    atoms2.pbc = True
    atoms2.cell[2, 2] = 100
    atoms2.set_cell(atoms2.cell * 1.5, scale_atoms=True)

    alg = CoincidenceAlgorithm(atoms1, atoms2)
    results = alg.run(tolerance=0.1, Nmax=15, angle_limits=(0, 60))
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

        def compare(a1, other):
            # other can be list
            comp = SymmetryEquivalenceCheck()
            return comp.compare(a1, other)

        def del_dups_cpp(lst):
            """O(n**2) algorithm, O(1) in memory"""
            pos = 0
            for item in lst:
                if not all([compare(item, item2) for item2 in islice(lst, pos)]):
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


# test_scaling_ase()
test_scaling_cpp(M=20, N=50)

# test_algorithm()

