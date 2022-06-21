from dataclasses import dataclass, astuple
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, Normalize, LinearSegmentedColormap
from matplotlib import path, patches

import numpy as np

from itertools import product

from ase.visualize.plot import plot_atoms

from hetbuilder.log import *
from hetbuilder.algorithm import Interface

from collections import namedtuple

from ase.neighborlist import natural_cutoffs, NeighborList

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import tkinter as tk

mutedblack = "#1a1a1a"


def plot_stack(stack: "ase.atoms.Atoms" = None, supercell_data: "namedtuple" = None):
    """Wrapper of the :func:~`ase.visualize.plot.plot_atoms` function."""
    fig = plt.gcf()
    axes = plt.gca()
    canvas = fig.canvas
    axes.clear()
    axes.set_yticks([])
    axes.set_xticks([])
    axes.set_xlabel("")
    axes.set_ylabel("")
    description = r"#{:d}, {:d} atoms, {:.1f} % deformation, $\theta=${:.2f}°".format(
        supercell_data.index,
        supercell_data.natoms,
        supercell_data.stress + supercell_data.strain,
        supercell_data.angle,
    )
    axes.set_title(description, fontsize=12)
    plot_atoms(stack, axes, radii=0.3)

    axes.set_frame_on(False)
    canvas.draw()


def plot_grid(
    basis: "numpy.ndarray" = None,
    supercell_matrix: "numpy.ndarray" = None,
    Nmax: int = None,
    **kwargs,
):
    """Plots lattice points of a unit cell."""
    axes = plt.gca()
    basis = basis[:2, :2].copy()
    a1 = basis[0, :].copy()
    a2 = basis[1, :].copy()
    # p = product(range(-Nmax, Nmax + 1), range(-Nmax, Nmax + 1))
    # points = np.array([n[0] * a1 + n[1] * a2 for n in p])
    # axes.scatter(points[:, 0], points[:, 1], **kwargs)
    axes.axline(0 * a1, 1 * a1, **kwargs)
    axes.axline(0 * a2, 1 * a2, **kwargs)
    for n in range(-Nmax, Nmax + 1):
        if n != 0:
            axes.axline(n * a1 + 0 * a2, n * a1 + n * a2, **kwargs)
            axes.axline(0 * a1 + n * a2, n * a1 + n * a2, **kwargs)


def plot_unit_cell_patch(cell: "numpy.ndarray", **kwargs):
    """Plots a face patch for a unit cell."""
    axes = plt.gca()
    cell = cell.copy()
    path1 = [
        (0, 0),
        (cell[0, :]),
        (cell[0, :] + cell[1, :]),
        (cell[1, :]),
        (0, 0),
    ]
    path1 = path.Path(path1)
    patch = patches.PathPatch(path1, **kwargs)
    axes.add_patch(patch)
    path2 = np.array(path1.vertices)
    xlim = (np.min(path2[:, 0]) - 4, np.max(path2[:, 0] + 4))
    ylim = (np.min(path2[:, 1]) - 4, np.max(path2[:, 1] + 4))
    axes.axis("equal")
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)


def plot_lattice_points(
    basis1: "numpy.ndarray" = None,
    basis2: "numpy.ndarray" = None,
    supercell_data: "namedtuple" = None,
    weight: float = None,
):
    """Plots lattice points of both bases on top of each other, as well as the coincidence supercell."""
    fig = plt.gcf()
    axes = plt.gca()
    canvas = fig.canvas
    axes.clear()
    axes.set_yticks([])
    axes.set_xticks([])
    axes.set_xlabel("")
    axes.set_ylabel("")
    (
        natoms,
        m1,
        m2,
        m3,
        m4,
        n1,
        n2,
        n3,
        n4,
        angle,
        stress,
        strain,
        index,
    ) = supercell_data
    sc1 = np.array([[m1, m2], [m3, m4]])
    sc2 = np.array([[n1, n2], [n3, n4]])
    Nmax = max(abs(j) for j in [m1, m2, m3, m4, n1, n2, n3, n4]) * 2

    # first cell
    plot_grid(
        basis=basis1,
        supercell_matrix=sc1,
        color="tab:red",
        # facecolor="tab:red",
        # edgecolor="tab:red",
        alpha=0.25,
        # s=2,
        lw=1,
        Nmax=Nmax,
    )
    A = sc1 @ basis1

    # second cell
    plot_grid(
        basis=basis2,
        supercell_matrix=sc2,
        color="tab:blue",
        # facecolor="tab:blue",
        # edgecolor="tab:blue",
        alpha=0.25,
        # s=2,
        lw=1,
        Nmax=Nmax,
    )
    B = sc2 @ basis2

    # supercell lattice points
    C = A + weight * (B - A)
    plot_unit_cell_patch(
        C,
        facecolor="tab:purple",
        alpha=0.5,
        edgecolor=mutedblack,
        linewidth=1,
        linestyle="--",
    )

    axes.set_frame_on(False)

    scdata = """M = ({: 2d}, {: 2d}, {: 2d}, {: 2d})\nN = ({: 2d}, {: 2d}, {: 2d}, {: 2d})""".format(
        m1, m2, m3, m4, n1, n2, n3, n4
    )
    axes.set_title(scdata, fontsize=12)
    canvas.draw()


def rand_jitter(arr, jitter):
    stdev = jitter * (max(arr) - min(arr)) + 0.01
    return arr + np.random.randn(len(arr)) * stdev


@dataclass
class SuperCellData:
    natoms: int
    m1: int
    m2: int
    m3: int
    m4: int
    n1: int
    n2: int
    n3: int
    n4: int
    angle: float
    stress: float
    strain: float
    index: int

    def __iter__(self):
        return iter(astuple(self))


class InteractivePlot:
    """ Interactive visualization of the results via matplotlib. 
    
    Args:
        bottom (ase.atoms.Atoms): Lower layer as primitive.
        top (ase.atoms.Atoms): Upper layer as primitive.
        results (list): List of :class:`~hetbuilder.algorithm.Interface` returned from the coincidence lattice search.
        weight (float, optional): Weight of the supercell.
    """

    def __init__(
        self,
        bottom: "ase.atoms.Atoms" = None,
        top: "ase.atoms.Atoms" = None,
        results: list = None,
        weight: float = 0.5,
    ) -> None:
        self.bottom = bottom
        self.top = top
        self.results = results
        self._weight = weight

    def __repr__(self):
        return "{}(nresults={})".format(self.__class__.__name__, len(self.results))

    def plot_results(self):
        """ Plots results interactively.

        Generates a matplotlib interface that allows to select the reconstructed stacks and save them to a file.    
        """
        results = self.results
        data = np.array(
            [[i.stress + i.strain, len(i.stack)] for i in results], dtype=float
        )
        color = [i.angle for i in results]
        norm = Normalize(vmin=0, vmax=90, clip=True)
        cmap = LinearSegmentedColormap.from_list(
            "",
            [
                "darkgreen",
                "tab:green",
                "lightgreen",
                "lightblue",
                "tab:blue",
                "royalblue",
            ],
        )
        mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        color = [mapper.to_rgba(v) for v in color]

        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=([6.4 * 3, 6.4]))
        ax1.scatter(
            data[:, 0], data[:, 1], color=color, alpha=0.75, picker=3.5, marker=".",
        )
        clb = plt.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1, ticks=[0, 30, 60, 90]
        )
        clb.set_label(
            r"Twist angle $\theta$ [°]", rotation=270, labelpad=12, fontsize=12
        )
        clb.ax.set_yticklabels(["0", "30", "60", "90"])

        ax1.set_xlim(0.00, np.max(data[:, 0]) * 1.1)
        ax1.set_ylim(np.min(data[:, 1]) * 0.95, np.max(data[:, 1]) * 1.05)
        ax1.set_ylim(0, ax1.get_ylim()[1] + 10)
        ax1.set_xlabel(r"deformation [%]", fontsize=12)
        ax1.set_ylabel("Number of atoms", fontsize=14)
        ax1.set_title("Click a point to select a structure.", fontsize=12)
        ax1.grid(axis="both", color="lightgray", linestyle="-", linewidth=1, alpha=0.2)

        ax2.set_yticks([])
        ax2.set_xticks([])
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        ax2.set_frame_on(False)

        ax3.set_yticks([])
        ax3.set_xticks([])
        ax3.set_xlabel("")
        ax3.set_ylabel("")
        ax3.set_frame_on(False)
        axbutton = plt.axes([0.8, 0.05, 0.1, 0.05])

        fig.canvas.mpl_connect("pick_event", self.__onpick)

        def __save(stack):
            try:
                sc = self.current_scdata
                form = self.current_stack.get_chemical_formula()
                info_string = []
                info_string.append(form)
                M = [sc.m1, sc.m2, sc.m3, sc.m4]
                N = [sc.n1, sc.n2, sc.n3, sc.n4]
                a = sc.angle
                s1 = sc.stress
                s2 = sc.strain
                info_string.append("M = [{} {} // {} {}]".format(*M))
                info_string.append("N = [{} {} // {} {}]".format(*N))
                info_string.append(f"angle = {a:.4f} [degree]")
                info_string.append(f"stress = {s1:.4f} [%]")
                info_string.append(f"strain = {s2:.4f} [%]")
                index = sc.index
                name = "{}_{:.2f}_degree_{}.in".format(form, a, index)
                stack.write(name, info_str=info_string, format="aims")
                logger.info("Saved structure to {}".format(name))
            except Exception as excpt:
                logger.error("You need to select a point first.")

        save = Button(axbutton, " Save this structure. ")
        save.on_clicked(lambda x: __save(self.current_stack))
        plt.show()

    def __onpick(self, event):
        point = event.artist
        mouseevent = event.mouseevent
        index = event.ind[0]
        fig = point.properties()["figure"]
        axes = fig.axes
        stack = self.results[index].stack.copy()
        M = np.array(self.results[index].M)
        N = np.array(self.results[index].N)
        angle = self.results[index].angle
        stress = self.results[index].stress
        strain = self.results[index].strain
        m1, m2, m3, m4 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
        n1, n2, n3, n4 = N[0, 0], N[0, 1], N[1, 0], N[1, 1]

        scdata = SuperCellData(
            natoms=int(len(stack)),
            m1=int(m1),
            m2=int(m2),
            m3=int(m3),
            m4=int(m4),
            n1=int(n1),
            n2=int(n2),
            n3=int(n3),
            n4=int(n4),
            angle=float(angle),
            stress=float(stress),
            strain=float(strain),
            index=int(index),
        )
        self.current_scdata = scdata
        self.current_stack = stack

        lower = self.bottom.copy()
        upper = self.top.copy()
        upper.rotate(angle, v="z", rotate_cell=True)
        basis1 = lower.cell.copy()[:2, :2]
        basis2 = upper.cell.copy()[:2, :2]

        plt.sca(fig.axes[1])
        plot_lattice_points(
            basis1=basis1, basis2=basis2, supercell_data=scdata, weight=self._weight,
        )
        plt.sca(fig.axes[2])
        plot_stack(stack=stack, supercell_data=scdata)

