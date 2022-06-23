from gettext import find
import matplotlib.pyplot as plt

import numpy as np
from math import sqrt
from itertools import islice

from ase.io.formats import string2index
from ase.utils import rotate
from ase.data import covalent_radii, atomic_numbers
from ase.data.colors import jmol_colors

from hetbuilder.algorithm import ase_atoms_to_cpp_atoms, cpp_atoms_to_ase_atoms

from hetbuilder.hetbuilder_backend import cpp_make_supercell, int1dVector, int2dVector
from ase.neighborlist import find_mic


def cpp_make_supercell(atoms, N):
    """C++ implementation of making a supercell.
    
    Also returns a list of indices holding equivalent atoms.
    """
    cppatoms = ase_atoms_to_cpp_atoms(atoms)
    N = int2dVector([int1dVector(k) for k in N])
    sc = cpp_make_supercell(cppatoms, N)
    index_map = sc.get_index_mapping()
    index_map = [int(k) for k in index_map]
    atoms = cpp_atoms_to_ase_atoms(sc)
    return atoms, index_map


class PlottingVariables:
    """ Modified from https://gitlab.com/ase/ase/-/blob/master/ase/visualize/plot.py to add bonds."""

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
        bonds=None,
    ):
        self.numbers = atoms.get_atomic_numbers()
        self.colors = colors
        self.offset = None
        show_bonds = True
        if bonds == None:
            show_bonds = False
        self.bonds = bonds

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
    """ Modified from https://gitlab.com/ase/ase/-/blob/master/ase/visualize/plot.py to add bonds."""

    def __init__(
        self,
        atoms,
        ax,
        rotation="",
        radii=None,
        colors=None,
        scale=1,
        offset=(0, 0),
        bonds=None,
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
            bonds=bonds,
            **parameters,
        )
        self.atoms = atoms
        self.ax = ax
        self.figure = ax.figure
        self.ax.set_aspect("equal")
        show_bonds = True if bonds is not None else False
        if show_bonds:
            self.show_bonds()

    def write(self):
        self.write_body()
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(0, self.h)

    def write_body(self):
        patch_list = make_patch_list(self)
        for patch in patch_list:
            self.ax.add_patch(patch)

    def show_bonds(self):
        from matplotlib import cm, colors

        cmap = colors.ListedColormap(["gray", "black"])

        pmin = np.min(self.atoms.positions[:, 2])
        pmax = np.max(self.atoms.positions[:, 2])
        for i, k in self.bonds:
            p1 = self.positions[i]
            p2 = self.positions[k]
            dr, l = find_mic(p2 - p1, self.cell)
            c = self.atoms.positions[k][2]
            c = (c - pmin) / (pmax - pmin)
            c = cmap(c)

            if l < 2.5:
                e1 = p1 + dr / 2
                e2 = p2 - dr / 2
                x1, y1, _ = list(zip(p1, e1))
                x2, y2, _ = list(zip(p2, e2))
                self.ax.plot(x1, y1, color=c, zorder=0)
                self.ax.plot(x2, y2, color=c, zorder=0)


def cell_to_lines(writer, cell):
    """ Taken from https://gitlab.com/ase/ase/-/blob/master/ase/visualize/plot.py."""
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
    """ Taken from https://gitlab.com/ase/ase/-/blob/master/ase/visualize/plot.py."""
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
    """ Taken from https://gitlab.com/ase/ase/-/blob/master/ase/visualize/plot.py.
    Plot an atoms object in a matplotlib subplot.

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
