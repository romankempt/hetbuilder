import ase.io
from ase.atoms import Atoms
from ase.spacegroup import Spacegroup
from ase.geometry import permute_axes

import numpy as np

from scipy.linalg import polar

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


class Interface:
    """Exposes the C++ implementation of the CppInterfaceClass.
    
    Attributes:
        bottom (ase.atoms.Atoms): Lower layer as supercell.
        top (ase.atoms.Atoms): Upper layer as supercell.
        stack (ase.atoms.Atoms): Combined lower and upper layer as supercell.
        M (numpy.ndarray): Supercell matrix M.
        N (numpy.ndarray): Supercell matrix N.
        spacegroup (ase.spacegroup.Spacegroup): Space group of the stack.
        angle (float): Twist angle in degree.
        stress (float): Stress measure.

    """

    def __init__(self, interface: "CppInterfaceClass" = None, weight=0.5) -> None:
        bottom = cpp_atoms_to_ase_atoms(interface.bottom)
        top = cpp_atoms_to_ase_atoms(interface.top)
        stack = cpp_atoms_to_ase_atoms(interface.stack)

        self.bottom = recenter(bottom)
        self.top = recenter(top)
        self.stack = recenter(stack)
        self.stack = stack
        self.M = [[j for j in k] for k in interface.M]
        self.N = [[j for j in k] for k in interface.N]
        self.spacegroup = Spacegroup(interface.spacegroup)
        self.angle = interface.angle
        self._weight = weight
        self._stress = None

    def __repr__(self):
        return "{}(M={}, N={}, spacegroup='{}', angle={:.1f}, stress={:.1f})".format(
            self.__class__.__name__,
            self.M,
            self.N,
            self.spacegroup.symbol,
            self.angle,
            self.stress,
        )

    @property
    def stress(self) -> float:
        """Returns the stress measure."""
        return self.measure_stress()

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


class CoincidenceAlgorithm:
    """Exposes the C++ implementation of the CppCoincidenceAlgorithmClass.
    
    Args:
        bottom (ase.atoms.Atoms): Lower layer, needs to be two-dimensional.
        top (ase.atoms.Atoms): Upper layer, needs to be two-dimensional.   
        check (bool): Runs checks on the input structures.
    """

    def __init__(
        self, bottom: "ase.atoms.Atoms", top: "ase.atoms.Atoms", check=True
    ) -> None:
        if check:
            self.bottom = check_atoms(bottom)
            self.top = check_atoms(top)
        else:
            self.bottom = bottom
            self.top = top

    def __repr__(self):
        return "{}(bottom={}, top={})".format(
            self.__class__.__name__, self.bottom, self.top
        )

    def run(
        self,
        Nmax: int = 10,
        Nmin: int = 0,
        angles: list[float] = [],
        tolerance: float = 0.1,
        weight: float = 0.5,
        distance: float = 4,
        no_idealize: bool = False,
        symprec: float = 1e-5,
        angle_tolerance: float = 5,
    ) -> list[Interface]:
        """Executes the coincidence lattice algorithm.

        Args:
            Nmax (int): Maximum number of translations. Defaults to 10.
            Nmin (int): Minimum number of translations. Defaults to -10.
            angles (list): List of angles in degree to search.
            tolerance (float): Tolerance criterion to accept lattice match. Corresponds to a distance in Angström. Defaults to 0.01.
            weight (float): The coincidence unit cell is C = A + weight * (B-A). Defaults to 0.5.
            distance (float): Interlayer distance of the stacks. Defaults to 4.0 Angström.
            no_idealize (bool): Does not idealize unit cell parameters in the spglib standardization routine. Defaults to False.
            symprec (float): Symmetry precision for spglib. Defaults to 1e-5 Angström.
            angle_tolerance (float): Angle tolerance fo the spglib `spgat` routines. Defaults to 5.

        Returns:
            list : A list of :class:`~hetbuilder.algorithm.Interface`.

        """
        bottom = ase_atoms_to_cpp_atoms(self.bottom)
        top = ase_atoms_to_cpp_atoms(self.top)
        if angles == []:
            angles = np.arange(0, 180, 1, dtype=float).tolist() + [180.0]

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
            no_idealize,
            symprec,
            angle_tolerance,
        )
        if len(results) == 0:
            logger.error("Could not find any coincidence pairs for these parameters.")
            return None
        elif len(results) > 0:
            if len(results) == 1:
                logger.info("Found 1 result.")
            else:
                logger.info("Found {:d} results.".format(len(results)))
            interfaces = [Interface(k, weight=weight) for k in results]
            return interfaces
