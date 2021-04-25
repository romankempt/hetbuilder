import ase.io
from ase.atoms import Atoms
from ase.spacegroup import Spacegroup

from hetbuilder.log import *
from hetbuilder.atom_checks import check_atoms, recenter

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


class Interface:
    def __init__(self, interface: "CppInterfaceClass" = None) -> None:
        bottom = cpp_atoms_to_ase_atoms(interface.bottom)
        top = cpp_atoms_to_ase_atoms(interface.top)
        stack = cpp_atoms_to_ase_atoms(interface.stack)
        self.bottom = recenter(bottom)
        self.top = recenter(top)
        self.stack = recenter(stack)
        self.M = [[j for j in k] for k in interface.M]
        self.N = [[j for j in k] for k in interface.N]
        self.spacegroup = Spacegroup(interface.spacegroup)
        self.angle = interface.angle
        # if monoclinic by bravais lattice, turn b to c


class CoincidenceAlgorithm:
    def __init__(self, bottom: "ase.atoms.Atoms", top: "ase.atoms.Atoms") -> None:
        self.bottom = check_atoms(bottom)
        self.top = check_atoms(top)

    def run(
        self,
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
        bottom = ase_atoms_to_cpp_atoms(self.bottom)
        top = ase_atoms_to_cpp_atoms(self.top)

        if (self.bottom == self.top) and (0 in angles):
            logger.warning("The bottom and top structure seem to be identical.")
            logger.warning(
                "Removing the angle 0° from the search because all lattice points would match."
            )
            angles = [k for k in angles if abs(k) > 1e-4]

        angles = double1dVector(angles)
        no_idealize = int(no_idealize)

        ncombinations = ((Nmax - Nmin + 1) ** 4) * len(angles)

        nthreads = get_number_of_omp_threads()
        logger.info("Using {:d} OpenMP threads.".format(nthreads))
        logger.info("Running through {:d} grid points...".format(ncombinations))
        alg = CppCoincidenceAlgorithmClass(bottom, top)
        results = alg.run(
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
        if len(results) == 0:
            logger.error("Could not find any coincidence pairs for these parameters.")

        elif len(results) > 0:
            logger.info("Found {:d} results.".format(len(results)))
            interfaces = [Interface(k) for k in results]
            return interfaces
        else:
            return None
