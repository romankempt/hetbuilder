from ase.utils import rotate
from hetbuilder_backend import *


import logging
import pretty_errors
from rich.logging import RichHandler
import sys
from contextlib import redirect_stdout

pretty_errors.configure(
    separator_character="*",
    filename_display=pretty_errors.FILENAME_EXTENDED,
    line_number_first=True,
    display_link=True,
    lines_before=2,
    lines_after=1,
    line_color=pretty_errors.RED + "> " + pretty_errors.default_config.line_color,
    code_color="  " + pretty_errors.default_config.line_color,
    truncate_code=True,
    display_locals=True,
)
pretty_errors.blacklist("c:/python")


class DuplicateFilter(logging.Filter):
    def filter(self, record):
        # add other fields if you need more granular comparison, depends on your app
        current_log = (record.module, record.levelno, record.msg)
        if current_log != getattr(self, "last_log", None):
            self.last_log = current_log
            return True
        return False


def setup_custom_logger(name):
    formatter = logging.Formatter(fmt="{message:s}", style="{")
    handler = RichHandler(
        show_time=False, markup=True, rich_tracebacks=True, show_path=False
    )
    logging.basicConfig(stream=sys.stdout)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.addFilter(DuplicateFilter())
    return logger


logger = setup_custom_logger("root")


def set_verbosity_level(verbosity):
    logger = logging.getLogger("root")
    if verbosity == 0:
        level = "WARNING"
    elif verbosity == 1:
        level = "INFO"
    else:
        level = "DEBUG"
    logger.setLevel(level)


set_verbosity_level(2)

import ase.io
from ase.atoms import Atoms

atoms = ase.io.read(
    "/mnt/c/Users/rkempt/Repositories/heterostructure_builder/tests/graphene.xyz"
)


def ase_atoms_to_cpp_atoms(atoms):
    lattice = atoms.cell.copy()
    positions = atoms.positions.copy()
    atomic_numbers = int1dVector([k for k in atoms.numbers])

    lattice = double2dVector([double1dVector(k) for k in lattice])
    positions = double2dVector([double1dVector(k) for k in positions])
    return CppAtomsClass(lattice, positions, atomic_numbers)


def cpp_atoms_to_ase_atoms(cppatoms):
    lattice = [[j for j in k] for k in cppatoms.lattice]
    positions = [[j for j in k] for k in cppatoms.positions]
    numbers = [i for i in cppatoms.atomic_numbers]
    atoms = Atoms(
        numbers=numbers, positions=positions, cell=lattice, pbc=[True, True, True]
    )
    return atoms


bottom = ase_atoms_to_cpp_atoms(atoms)
top = ase_atoms_to_cpp_atoms(atoms)


def run_coincidence_algorithm(
    bottom,
    top,
    Nmax=10,
    Nmin=-10,
    angles=[0, 10, 20, 30],
    tolerance=0.01,
    weight=0.5,
    distance=4,
    no_idealize=False,
    symprec=1e-5,
    angle_tolerance=5,
):
    angles = double1dVector(angles)
    no_idealize = int(no_idealize)
    alg = CppCoincidenceAlgorithmClass(
        bottom,
        top,
        Nmax,
        Nmin,
        angles,
        tolerance,
        weight,
        distance,
        no_idealize,
        symprec,
        angle_tolerance,
    )
    results = alg.run()
    return results


nthreads = get_number_of_omp_threads()
print("Nthreads: ", nthreads)
results = run_coincidence_algorithm(bottom, top, Nmax=5, Nmin=-5)
print("Found results: ", len(results))
example = results[5]
one = cpp_atoms_to_ase_atoms(example.bottom)
two = cpp_atoms_to_ase_atoms(example.top)
stack = cpp_atoms_to_ase_atoms(example.stack)
print(results[5].angle)

import numpy as np

# N = np.array([[-5, -4, 0], [-1, -2, 0], [0, 0, 1]])

# M = np.array([[k for k in l] for l in example.M])
N = np.array([[k for k in l] for l in example.N])
Nc = int2dVector([int1dVector(k) for k in N])
print(N)

from ase.visualize import view
from ase.build.supercells import lattice_points_in_supercell, clean_matrix


def make_supercell(prim, P, wrap=True, tol=1e-5):
    r"""Generate a supercell by applying a general transformation (*P*) to
    the input configuration (*prim*).

    The transformation is described by a 3x3 integer matrix
    `\mathbf{P}`. Specifically, the new cell metric
    `\mathbf{h}` is given in terms of the metric of the input
    configuration `\mathbf{h}_p` by `\mathbf{P h}_p =
    \mathbf{h}`.

    Parameters:

    prim: ASE Atoms object
        Input configuration.
    P: 3x3 integer matrix
        Transformation matrix `\mathbf{P}`.
    wrap: bool
        wrap in the end
    tol: float
        tolerance for wrapping
    """

    supercell_matrix = P
    supercell = clean_matrix(supercell_matrix @ prim.cell)

    # cartesian lattice points
    lattice_points_frac = lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell)

    lattice_points = clean_matrix(lattice_points)

    superatoms = Atoms(cell=supercell, pbc=prim.pbc)

    for lp in lattice_points:
        shifted_atoms = prim.copy()
        shifted_atoms.positions += lp
        superatoms.extend(shifted_atoms)

    return superatoms


test1 = make_supercell(atoms, N, wrap=False)
test2 = cpp_make_supercell(bottom, Nc)
test2 = cpp_atoms_to_ase_atoms(test2)

# view(test1)
# view(test2)
view(one)
# view(two)

