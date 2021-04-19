import numpy as np

scm = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]])


def lattice_points_in_supercell(supercell_matrix):
    """Find all lattice points contained in a supercell.

    Adapted from pymatgen, which is available under MIT license:
    The MIT License (MIT) Copyright (c) 2011-2012 MIT & The Regents of the
    University of California, through Lawrence Berkeley National Laboratory
    """

    diagonals = np.array(
        [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1],
        ]
    )
    d_points = np.dot(diagonals, supercell_matrix)

    mins = np.min(d_points, axis=0)
    maxes = np.max(d_points, axis=0) + 1

    ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
    br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
    cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))

    tvects = frac_points[
        np.all(frac_points < 1 - 1e-10, axis=1) & np.all(frac_points >= -1e-10, axis=1)
    ]
    assert len(tvects) == round(abs(np.linalg.det(supercell_matrix)))
    return tvects


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
    from ase.atoms import Atoms

    supercell_matrix = P
    supercell = supercell_matrix @ prim.cell

    # cartesian lattice points
    lattice_points_frac = lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell)

    superatoms = Atoms(cell=supercell, pbc=prim.pbc)

    for lp in lattice_points:
        shifted_atoms = prim.copy()
        shifted_atoms.positions += lp
        superatoms.extend(shifted_atoms)

    return superatoms


def rotate(atoms, angle, v=(0, 0, 1), center=(0, 0, 0), rotate_cell=True):
    import numpy as np
    from numpy import cos, sin

    normv = np.linalg.norm(v)
    atoms = atoms.copy()
    otheratoms = atoms.copy()
    angle *= np.pi / 180
    v /= normv
    c = cos(angle)
    s = sin(angle)

    p = atoms.positions - center
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    otheratoms.positions[:] = (R @ otheratoms.positions.T).T

    atoms.positions[:] = (
        c * p - np.cross(p, s * v) + np.outer(np.dot(p, v), (1.0 - c) * v) + center
    )
    if rotate_cell:
        rotcell = atoms.get_cell()
        rotcell[:] = (
            c * rotcell
            - np.cross(rotcell, s * v)
            + np.outer(np.dot(rotcell, v), (1.0 - c) * v)
        )
        atoms.set_cell(rotcell)
        othercell = otheratoms.get_cell()
        othercell[:] = (R @ othercell.T).T
        otheratoms.set_cell(othercell)
    return atoms, otheratoms

