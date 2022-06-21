import ase.io
from ase.atoms import Atoms
from ase.spacegroup import Spacegroup
from ase.geometry import permute_axes
from ase.geometry.analysis import Analysis

from itertools import combinations_with_replacement

from dataclasses import dataclass

import numpy as np

from scipy.linalg import polar

from hetbuilder.log import *
from hetbuilder.atom_checks import check_atoms, recenter

import sys

from hetbuilder.hetbuilder_backend import (
    double2dVector,
    double1dVector,
    int1dVector,
    int2dVector,
    CppAtomsClass,
    CppCoincidenceAlgorithmClass,
    CppInterfaceClass,
    get_number_of_omp_threads,
)


def ase_atoms_to_cpp_atoms(atoms: "ase.atoms.Atoms") -> "CppAtomsClass":
    """Converts :class:`~ase.atoms.Atoms` to the C++ CppAtomsClass."""
    lattice = atoms.cell.copy()
    positions = atoms.positions.copy()
    atomic_numbers = int1dVector([k for k in atoms.numbers])
    lattice = double2dVector([double1dVector(k) for k in lattice])
    positions = double2dVector([double1dVector(k) for k in positions])
    return CppAtomsClass(lattice, positions, atomic_numbers)


def cpp_atoms_to_ase_atoms(cppatoms: "CppAtomsClass") -> "ase.atoms.Atoms":
    """Converts the C++ CppAtomsClass to :class:`~ase.atoms.Atoms`"""
    lattice = [[j for j in k] for k in cppatoms.lattice]
    positions = [[j for j in k] for k in cppatoms.positions]
    numbers = [i for i in cppatoms.atomic_numbers]
    atoms = Atoms(
        numbers=numbers, positions=positions, cell=lattice, pbc=[True, True, True]
    )
    return atoms


def check_angles(
    angle_stepsize: float = 1, angle_limits: tuple = (0, 180), angles: list = []
) -> list:
    """ Helper function to assert correct input of angles."""
    if len(angles) == 0:
        a1 = angle_limits[0]
        a2 = angle_limits[1]
        assert a2 > a1, "Second angle must be larger than first one."
        assert angle_stepsize > 0, "Angle stepsize must be larger than zero."
        assert angle_stepsize < abs(
            a2 - a1
        ), "Angle stepsize must be larger then difference between angles."
        angles = list(np.arange(a1, a2, step=angle_stepsize)) + [a2]
        logger.info(
            "Searching {:d} angles between {:.1f} and {:.1f} degree with a stepsize of {:.1f} degree.".format(
                len(angles), a1, a2, angle_stepsize
            )
        )
        return angles
    elif angles != None:
        msg = ", ".join([str(k) for k in angles])
        logger.info("Calculating the following angles: {} in degree.".format(msg))
        return list(angles)
    else:
        logger.error("Angle specifications not recognized.")


def get_unique_bond_lengths(atoms: "ase.atoms.AToms") -> dict:
    """ Returns a dictionary containing average bond values of a structure. """
    symbs = set(atoms.get_chemical_symbols())
    pairs = list(combinations_with_replacement(symbs, 2))
    ana = Analysis(atoms)
    unique_bonds = {}
    for p in pairs:
        bonds = ana.get_bonds(*p)
        if bonds == [[]]:
            continue
        bond_values = ana.get_values(bonds)
        avg = np.average(bond_values)
        unique_bonds[p] = avg
    return unique_bonds


@dataclass
class Interface:
    """Exposes the C++ implementation of the CppInterfaceClass.
    
    Attributes:
        bottom (ase.atoms.Atoms): Lower layer as supercell.
        top (ase.atoms.Atoms): Upper layer as supercell.
        stack (ase.atoms.Atoms): Combined lower and upper layer as supercell.
        M (numpy.ndarray): Supercell matrix M.
        N (numpy.ndarray): Supercell matrix N.
        angle (float): Twist angle in degree.
        stress (float): Stress measure of the unit cell.
        strain (float): Strain measure of the bond lengths.

    """

    def __init__(
        self, interface: "CppInterfaceClass" = None, weight=0.5, **kwargs
    ) -> None:
        bottom = cpp_atoms_to_ase_atoms(interface.bottom)
        top = cpp_atoms_to_ase_atoms(interface.top)
        stack = cpp_atoms_to_ase_atoms(interface.stack)

        self.bottom = recenter(bottom)
        self.top = recenter(top)
        self.stack = recenter(stack)
        self.M = [[j for j in k] for k in interface.M]
        self.N = [[j for j in k] for k in interface.N]
        self.angle = interface.angle
        self._weight = weight
        self._stress = None
        self.bbl = kwargs.get("bottom_bond_lengths", None)
        self.tbl = kwargs.get("top_bond_lengths", None)

    def __repr__(self):
        return "{}(M={}, N={}, angle={:.1f}, stress={:.1f})".format(
            self.__class__.__name__, self.M, self.N, self.angle, self.stress,
        )

    @property
    def stress(self) -> float:
        """Returns the stress measure."""
        return self.measure_stress()

    @property
    def strain(self) -> float:
        """Returns the strain measure."""
        return self.measure_strain()

    def measure_stress(self) -> float:
        """Measures the stress on both unit cells."""
        A = self.bottom.cell.copy()[:2, :2]
        B = self.top.cell.copy()[:2, :2]
        C = A + self._weight * (B - A)
        T1 = C @ np.linalg.inv(A)
        T2 = C @ np.linalg.inv(B)

        def measure(P):
            eps = P - np.identity(2)
            meps = np.sqrt(
                (
                    eps[0, 0] ** 2
                    + eps[1, 1] ** 2
                    + eps[0, 0] * eps[1, 1]
                    + eps[1, 0] ** 2
                )
                / 4
            )
            return meps

        U1, P1 = polar(T1)  # this one goes counterclockwise
        U2, P2 = polar(T2)  # this one goes clockwise
        # u is rotation, p is strain
        meps1 = measure(P1)
        meps2 = measure(P2)
        stress = meps1 + meps2
        # return (stress, P1 - np.identity(2), P2 - np.identity(2))
        return stress

    def measure_strain(self) -> float:
        """Measures the average strain on bond lengths on both substructures."""
        bond_lengths = get_unique_bond_lengths(self.stack)
        bottom_strain = []
        top_strain = []
        for k, b in bond_lengths.items():
            for k2, b2 in self.bbl.items():
                if (k2 == k) or (k2[::-1] == k):
                    d = np.abs((b2 - b)) / b * 100
                    # print("bottom", k2, d)
                    bottom_strain.append(d)
            for k3, b3 in self.tbl.items():
                if (k3 == k) or (k3[::-1] == k):
                    d = np.abs((b3 - b)) / b * 100
                    # print("top", k3, d)
                    top_strain.append(d)
        strain = np.average(bottom_strain) + np.average(top_strain)
        return strain


class CoincidenceAlgorithm:
    """Exposes the C++ implementation of the CppCoincidenceAlgorithmClass.
    
    Args:
        bottom (ase.atoms.Atoms): Lower layer, needs to be two-dimensional.
        top (ase.atoms.Atoms): Upper layer, needs to be two-dimensional.   
    """

    def __init__(self, bottom: "ase.atoms.Atoms", top: "ase.atoms.Atoms") -> None:
        self.bottom = check_atoms(bottom)
        self.top = check_atoms(top)
        self.bdl = get_unique_bond_lengths(bottom)
        self.tbl = get_unique_bond_lengths(top)

    def __repr__(self):
        return "{}(bottom={}, top={})".format(
            self.__class__.__name__, self.bottom, self.top
        )

    def run(
        self,
        Nmax: int = 10,
        Nmin: int = 0,
        angles: list = [],
        angle_limits: tuple = (0, 90),
        angle_stepsize: float = 1.0,
        tolerance: float = 0.1,
        weight: float = 0.5,
        distance: float = 4,
        standardize: bool = False,
        no_idealize: bool = False,
        symprec: float = 1e-5,
        angle_tolerance: float = 5,
        verbosity: int = 0,
    ) -> list:
        """Executes the coincidence lattice algorithm.

        Args:
            Nmax (int): Maximum number of translations. Defaults to 10.
            Nmin (int): Minimum number of translations. Defaults to 0.
            angles (list): List of angles in degree to search. Takes precedence over angle_limits and angle_stepsize.
            angle_limits (tuple): Lower and upper bound of angles too look through with given step size by angle_stepsize. Defaults to (0, 90) degree.
            angle_stepsize (float): Increment of angles to look through. Defaults to 1.0 degree.
            tolerance (float): Tolerance criterion to accept lattice match. Corresponds to a distance in Angström. Defaults to 0.1.
            weight (float): The coincidence unit cell is C = A + weight * (B-A). Defaults to 0.5.
            distance (float): Interlayer distance of the stacks. Defaults to 4.0 Angström.
            standardize (bool): Perform spglib standardization. Defaults to true.
            no_idealize (bool): Does not idealize unit cell parameters in the spglib standardization routine. Defaults to False.
            symprec (float): Symmetry precision for spglib. Defaults to 1e-5 Angström.
            angle_tolerance (float): Angle tolerance fo the spglib `spgat` routines. Defaults to 5.
            verbosity (int): Debug level for printout of Coincidence Algorithm. Defaults to 0.

        Returns:
            list : A list of :class:`~hetbuilder.algorithm.Interface`.

        """
        bottom = ase_atoms_to_cpp_atoms(self.bottom)
        top = ase_atoms_to_cpp_atoms(self.top)
        angles = check_angles(
            angle_limits=angle_limits, angle_stepsize=angle_stepsize, angles=angles
        )
        if (self.bottom == self.top) and (0 in angles):
            logger.warning("The bottom and top structure seem to be identical.")
            logger.warning(
                "Removing the angle 0° from the search because all lattice points would match."
            )
            angles = [k for k in angles if abs(k) > 1e-4]
            assert len(angles) > 0, "List of angles contains no values."

        assert Nmin < Nmax, "Nmin must be smaller than Nmax."
        assert Nmin >= 0, "Nmin must be larger than or equal 0."
        assert Nmax > 0, "Nmax must be larger than 0."
        assert distance > 0, "Interlayer distance must be larger than zero."
        assert tolerance > 0, "Tolerance must be larger than zero."
        assert (
            angle_tolerance >= 0
        ), "Angle tolerance must be larger than or equal zero."
        assert (symprec) > 0, "Symmetry precision must be larger than zero."
        assert (weight >= 0) and (weight <= 1), "Weight factor must be between 0 and 1."
        assert verbosity in [0, 1, 2], "Verbose must be 0, 1, or 2."

        angles = double1dVector(angles)
        no_idealize = int(no_idealize)

        ncombinations = ((2 * (Nmax - Nmin)) ** 4) * len(angles)

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
            standardize,
            no_idealize,
            symprec,
            angle_tolerance,
            verbosity,
        )
        if len(results) == 0:
            logger.error("Could not find any coincidence pairs for these parameters.")
            return None
        elif len(results) > 0:
            if len(results) == 1:
                logger.info("Found 1 result.")
            else:
                logger.info("Found {:d} results.".format(len(results)))

            interfaces = [
                Interface(
                    k,
                    weight=weight,
                    bottom_bond_lengths=self.bdl,
                    top_bond_lengths=self.tbl,
                )
                for k in results
            ]
            return interfaces
