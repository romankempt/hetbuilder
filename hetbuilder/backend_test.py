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
        tolerance=0.1, Nmax=10, angle_limits=(0, 30), angle_stepsize=1, verbosity=2
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


from ase.gui.gui import GUI
import ase.gui.ui as ui
from ase.gui.images import Images

# test_scaling_ase()
# test_scaling_cpp(M=10, N=50)

# test_algorithm()


from ase.neighborlist import natural_cutoffs, NeighborList


def get_bonds(atoms):
    radii = np.array(natural_cutoffs(atoms))
    nl = NeighborList(radii * 1.5, skin=0, self_interaction=True, bothways=True)
    nl.update(atoms)
    nbonds = nl.nneighbors + nl.npbcneighbors

    # bonds = np.zeros((nbonds, 6))
    bonds = []
    for a in range(len(atoms)):
        indices, offsets = nl.get_neighbors(a)
        for i, offset in zip(indices, offsets):
            startvector = atoms.positions[a]
            endvector = atoms.positions[i] + offset @ atoms.get_cell()
            if np.sum((endvector - startvector) ** 2) > 0:
                # print(a, i, offset, startvector, endvector)
                bonds.append([startvector, endvector])

    return bonds


import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms


import numpy as np
from math import sqrt
from itertools import islice

from ase.io.formats import string2index
from ase.utils import rotate
from ase.data import covalent_radii, atomic_numbers
from ase.data.colors import jmol_colors


class PlottingVariables:
    # removed writer - self
    def __init__(
        self,
        atoms,
        rotation="",
        show_unit_cell=2,
        radii=None,
        bbox=None,
        colors=None,
        scale=20,
        maxwidth=500,
        extra_offset=(0.0, 0.0),
        do_not_shift=True,
    ):
        self.numbers = atoms.get_atomic_numbers()
        self.colors = colors
        if colors is None:
            ncolors = len(jmol_colors)
            self.colors = jmol_colors[self.numbers.clip(max=ncolors - 1)]

        if radii is None:
            radii = covalent_radii[self.numbers]
        elif isinstance(radii, float):
            radii = covalent_radii[self.numbers] * radii
        else:
            radii = np.array(radii)

        natoms = len(atoms)

        if isinstance(rotation, str):
            rotation = rotate(rotation)

        cell = atoms.get_cell()
        disp = atoms.get_celldisp().flatten()

        if show_unit_cell > 0:
            L, T, D = cell_to_lines(self, cell)
            cell_vertices = np.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        cell_vertices[c1, c2, c3] = np.dot([c1, c2, c3], cell) + disp
            cell_vertices.shape = (8, 3)
            cell_vertices = np.dot(cell_vertices, rotation)
        else:
            L = np.empty((0, 3))
            T = None
            D = None
            cell_vertices = None

        nlines = len(L)

        positions = np.empty((natoms + nlines, 3))
        R = atoms.get_positions()
        positions[:natoms] = R
        positions[natoms:] = L

        r2 = radii ** 2
        for n in range(nlines):
            d = D[T[n]]
            if (
                (((R - L[n] - d) ** 2).sum(1) < r2)
                & (((R - L[n] + d) ** 2).sum(1) < r2)
            ).any():
                T[n] = -1

        positions = np.dot(positions, rotation)
        R = positions[:natoms]

        if bbox is None:
            X1 = (R - radii[:, None]).min(0)
            X2 = (R + radii[:, None]).max(0)
            if show_unit_cell == 2:
                X1 = np.minimum(X1, cell_vertices.min(0))
                X2 = np.maximum(X2, cell_vertices.max(0))
            M = (X1 + X2) / 2
            S = 1.05 * (X2 - X1)
            w = scale * S[0]
            if w > maxwidth:
                w = maxwidth
                scale = w / S[0]
            h = scale * S[1]
            offset = np.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])
        else:
            w = (bbox[2] - bbox[0]) * scale
            h = (bbox[3] - bbox[1]) * scale
            offset = np.array([bbox[0], bbox[1], 0]) * scale

        offset[0] = offset[0] - extra_offset[0]
        offset[1] = offset[1] - extra_offset[1]
        self.w = w + extra_offset[0]
        self.h = h + extra_offset[1]

        positions *= scale
        if not do_not_shift:
            positions -= offset

        if nlines > 0:
            D = np.dot(D, rotation)[:, :2] * scale

        if cell_vertices is not None:
            cell_vertices *= scale
            cell_vertices -= offset

        cell = np.dot(cell, rotation)
        cell *= scale

        self.cell = cell
        self.positions = positions
        self.D = D
        self.T = T
        self.cell_vertices = cell_vertices
        self.natoms = natoms
        self.d = 2 * scale * radii
        self.constraints = atoms.constraints

        # extension for partial occupancies
        self.frac_occ = False
        self.tags = None
        self.occs = None

        try:
            self.occs = atoms.info["occupancy"]
            self.tags = atoms.get_tags()
            self.frac_occ = True
        except KeyError:
            pass


class Matplotlib(PlottingVariables):
    def __init__(
        self,
        atoms,
        ax,
        rotation="",
        radii=None,
        colors=None,
        scale=1,
        offset=(0, 0),
        **parameters,
    ):
        PlottingVariables.__init__(
            self,
            atoms,
            rotation=rotation,
            radii=radii,
            colors=colors,
            scale=scale,
            extra_offset=offset,
            **parameters,
        )

        self.ax = ax
        self.figure = ax.figure
        self.ax.set_aspect("equal")

    def write(self):
        self.write_body()
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(0, self.h)

    def write_body(self):
        patch_list = make_patch_list(self)
        for patch in patch_list:
            self.ax.add_patch(patch)


def cell_to_lines(writer, cell):
    # XXX this needs to be updated for cell vectors that are zero.
    # Cannot read the code though!  (What are T and D? nn?)
    nlines = 0
    nsegments = []
    for c in range(3):
        d = sqrt((cell[c] ** 2).sum())
        n = max(2, int(d / 0.3))
        nsegments.append(n)
        nlines += 4 * n

    positions = np.empty((nlines, 3))
    T = np.empty(nlines, int)
    D = np.zeros((3, 3))

    n1 = 0
    for c in range(3):
        n = nsegments[c]
        dd = cell[c] / (4 * n - 2)
        D[c] = dd
        P = np.arange(1, 4 * n + 1, 4)[:, None] * dd
        T[n1:] = c
        for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            n2 = n1 + n
            positions[n1:n2] = P + i * cell[c - 2] + j * cell[c - 1]
            n1 = n2

    return positions, T, D


def make_patch_list(writer):
    from matplotlib.path import Path
    from matplotlib.patches import Circle, PathPatch, Wedge

    indices = writer.positions[:, 2].argsort()
    patch_list = []
    for a in indices:
        xy = writer.positions[a, :2]
        if a < writer.natoms:
            r = writer.d[a] / 2
            if writer.frac_occ:
                site_occ = writer.occs[str(writer.tags[a])]
                # first an empty circle if a site is not fully occupied
                if (np.sum([v for v in site_occ.values()])) < 1.0:
                    # fill with white
                    fill = "#ffffff"
                    patch = Circle(xy, r, facecolor=fill, edgecolor="black")
                    patch_list.append(patch)

                start = 0
                # start with the dominant species
                for sym, occ in sorted(
                    site_occ.items(), key=lambda x: x[1], reverse=True
                ):
                    if np.round(occ, decimals=4) == 1.0:
                        patch = Circle(
                            xy, r, facecolor=writer.colors[a], edgecolor="black"
                        )
                        patch_list.append(patch)
                    else:
                        # jmol colors for the moment
                        extent = 360.0 * occ
                        patch = Wedge(
                            xy,
                            r,
                            start,
                            start + extent,
                            facecolor=jmol_colors[atomic_numbers[sym]],
                            edgecolor="black",
                        )
                        patch_list.append(patch)
                        start += extent

            else:
                if (
                    (xy[1] + r > 0)
                    and (xy[1] - r < writer.h)
                    and (xy[0] + r > 0)
                    and (xy[0] - r < writer.w)
                ):
                    patch = Circle(xy, r, facecolor=writer.colors[a], edgecolor="black")
                    patch_list.append(patch)
        else:
            a -= writer.natoms
            c = writer.T[a]
            if c != -1:
                hxy = writer.D[c]
                patch = PathPatch(Path((xy + hxy, xy - hxy)))
                patch_list.append(patch)
    return patch_list


def plot_atoms(atoms, ax=None, **parameters):
    """Plot an atoms object in a matplotlib subplot.

    Parameters
    ----------
    atoms : Atoms object
    ax : Matplotlib subplot object
    rotation : str, optional
        In degrees. In the form '10x,20y,30z'
    show_unit_cell : int, optional, default 2
        Draw the unit cell as dashed lines depending on value:
        0: Don't
        1: Do
        2: Do, making sure cell is visible
    radii : float, optional
        The radii of the atoms
    colors : list of strings, optional
        Color of the atoms, must be the same length as
        the number of atoms in the atoms object.
    scale : float, optional
        Scaling of the plotted atoms and lines.
    offset : tuple (float, float), optional
        Offset of the plotted atoms and lines.
    """
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]

    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()
    Matplotlib(atoms, ax, **parameters).write()
    return ax


def test_get_bonds():
    atoms = mx2("MoS2")
    atoms.pbc = True
    atoms.cell[2, 2] = 100
    bonds = get_bonds(atoms)

    fig = plt.figure()
    # ax = fig.add_subplot(projection="3d")
    ax = fig.add_subplot(111)
    for i, k in bonds:
        x, y, z = list(zip(i, k))
        print(x, y, z)
        ax.plot(x, y, color="black")
    plot_atoms(
        atoms, ax=ax, radii=0.3, offset=(0, 0), scale=1, do_not_shift=True,
    )
    ax.scatter(
        atoms.positions[:, 0],
        atoms.positions[:, 1],
        # atoms.positions[:, 2],
        c="crimson",
    )
    plt.show()

    # for i, k in zip(a, b):
    #    print(i, k)


test_get_bonds()
