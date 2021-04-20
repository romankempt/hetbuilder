#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <tuple>

#include "math_functions.h"
#include "logging_functions.h"
#include "atom_class.h"

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

/**
 * Find all lattice points contained in a supercell.

    Adapted from pymatgen, which is available under MIT license:
    The MIT License (MIT) Copyright (c) 2011-2012 MIT & The Regents of the
    University of California, through Lawrence Berkeley National Laboratory
 */
double2dvec_t lattice_points_in_supercell(int2dvec_t &SuperCellMatrix)
{
    int2dvec_t diagonals = {{0, 0, 0},
                            {0, 0, 1},
                            {0, 1, 0},
                            {0, 1, 1},
                            {1, 0, 0},
                            {1, 0, 1},
                            {1, 1, 0},
                            {1, 1, 1}};
    int2dvec_t dpoints;
    int1dvec_t dotproduct;
    for (int row = 0; row < diagonals.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<int, int>(diagonals[row], SuperCellMatrix);
        dpoints.push_back(dotproduct);
    }

    int k = dpoints.size() - 1;
    int1dvec_t mins = dpoints[0];
    int1dvec_t maxes = dpoints[k];

    int minrowsum = dpoints[0][0] + dpoints[0][1] + dpoints[0][2];
    int maxrowsum = dpoints[k][0] + dpoints[k][1] + dpoints[k][2];
    int rowsum;
    for (int row = 0; row < dpoints.size(); row++)
    {
        rowsum = dpoints[row][0] + dpoints[row][1] + dpoints[row][2];
        if (rowsum < minrowsum)
        {
            mins = dpoints[row];
        }
        if (rowsum > maxrowsum)
        {
            maxes = dpoints[row];
        }
    }

    int2dvec_t ar, br, cr;
    int1dvec_t subvec(3, 0);
    for (int a = mins[0]; a < maxes[0]; a++)
    {
        subvec = {a, 0, 0};
        ar.push_back(subvec);
    }
    for (int b = mins[1]; b < maxes[1]; b++)
    {
        subvec = {0, b, 0};
        br.push_back(subvec);
    }
    for (int c = mins[2]; c < maxes[2]; c++)
    {
        subvec = {0, 0, c};
        cr.push_back(subvec);
    }

    int2dvec_t allpoints;
    for (int i = 0; i < ar.size(); i++)
    {
        for (int j = 0; j < br.size(); j++)
        {
            for (int k = 0; k < cr.size(); k++)
            {
                subvec[0] = ar[i][0] + br[j][0] + cr[k][0];
                subvec[1] = ar[i][1] + br[j][1] + cr[k][1];
                subvec[2] = ar[i][2] + br[j][2] + cr[k][2];
                allpoints.push_back(subvec);
            }
        }
    }

    // convert integer matrix to doubles
    double2dvec_t allpoints_double;
    for (int row = 0; row < allpoints.size(); row++)
    {
        double1dvec_t doubleVec(allpoints[row].begin(), allpoints[row].end());
        allpoints_double.push_back(doubleVec);
    };

    double determinant = get_3x3_matrix_determinant<int>(SuperCellMatrix);
    double2dvec_t InvSuperCellMatrix = invert_3x3_matrix<int>(SuperCellMatrix);
    double2dvec_t fracpoints;
    std::vector<double> dp;
    for (int row = 0; row < allpoints.size(); row++)
    {
        dp = vec1x3_dot_3x3_matrix<double, double>(allpoints_double[row], InvSuperCellMatrix);
        fracpoints.push_back(dp);
    }

    double2dvec_t tvects;
    double fa, fb, fc;
    double1dvec_t fvec;
    for (int row = 0; row < fracpoints.size(); row++)
    {
        fa = fracpoints[row][0];
        fb = fracpoints[row][1];
        fc = fracpoints[row][2];
        if ((fa <= (1 - 1e-10) && (fa >= (-1e-10))) && (fb <= (1 - 1e-10) && (fb >= (-1e-10))) && (fc <= (1 - 1e-10) && (fc >= (-1e-10))))
        {
            fvec = {fa, fb, fc};
            tvects.push_back(fvec);
        }
    }
    try
    {
        int detsize = (int)determinant;
        if (detsize != determinant)
        {
            throw "Determinant of supercell does not match number of lattice points.";
        }
    }
    catch (const char *msg)
    {
        std::cout << msg << std::endl;
        tvects = {};
    }

    return tvects;
};

/**
 * Generate a supercell by applying a SuperCellMatrix to
    the input atomic configuration prim.
*/
Atoms make_supercell(Atoms &prim, int2dvec_t &SuperCellMatrix)
{
    double2dvec_t fracpoints = lattice_points_in_supercell(SuperCellMatrix);
    double2dvec_t cell = prim.lattice;
    double2dvec_t supercell = matrix3x3_dot_matrix3x3<int, double>(SuperCellMatrix, cell);

    double2dvec_t lattice_points;
    double1dvec_t dotproduct;
    for (int row = 0; row < fracpoints.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<double, double>(fracpoints[row], supercell);
        lattice_points.push_back(dotproduct);
    }

    double2dvec_t new_positions = prim.positions;
    int1dvec_t new_numbers = prim.atomic_numbers;
    int1dvec_t new_spin = prim.spins;
    int1dvec_t new_equiv = prim.equivalent_atoms;
    for (int i = 0; i < prim.num_atom; i++)
    {
        double1dvec_t atom_pos = prim.positions[i];
        int number = prim.atomic_numbers[i];
        int spin = prim.spins[i];
        int equivalent = prim.equivalent_atoms[i];
        for (int row = 0; row < lattice_points.size(); row++)
        {
            double1dvec_t lp = lattice_points[row];
            if (!((lp[0] == 0) && (lp[1] == 0) && (lp[2] == 0)))
            {
                for (int k = 0; k < 3; k++)
                {
                    lp[k] += atom_pos[k];
                }
                new_positions.push_back(lp);
                new_numbers.push_back(number);
                new_spin.push_back(spin);
                new_equiv.push_back(equivalent);
            }
        }
    }
    Atoms superatoms = {supercell, new_positions, new_numbers, new_spin, new_equiv};
    return superatoms;
};

Atoms rotate_atoms_around_z(Atoms &atoms, const double &theta)
{
    double t = M_PI * theta / 180.0;
    double c = std::cos(t);
    double s = std::sin(t);
    double2dvec_t R = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};

    double2dvec_t positions = atoms.positions;
    double2dvec_t rotpositions;
    for (int row = 0; row < positions.size(); row++)
    {
        std::vector<double> vec = positions[row];
        std::vector<double> rotvec = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                rotvec[i] += (R[i][j] * vec[j]);
            }
        }
        rotpositions.push_back(rotvec);
    }

    double2dvec_t tcell = transpose_matrix3x3<double>(atoms.lattice);
    double2dvec_t rotcelltranspose = matrix3x3_dot_matrix3x3(R, tcell);
    double2dvec_t rotcell = transpose_matrix3x3<double>(rotcelltranspose);

    Atoms rotatoms(rotcell, rotpositions, atoms.atomic_numbers, atoms.spins, atoms.equivalent_atoms);
    return rotatoms;
};

/* def stack_atoms(atom1, atom2, weight=0.5, distance=4):
    """ Stacks two layered structures on top of each other.
    
    Args:
        atom1 (atoms): Lower layer.
        atom2 (atoms): Upper layer.
        weight (float, optional): Value between 0 and 1, defaults to 0.5. The unit cell of the reconstructed stack is :math:`B + w \cdot (T - B)`.
        distance (int, optional): Interlayer distance in Angström. Defaults to 4.
    
    Returns:
        atoms: Reconstructed stack.
    """

    bottom = atom1.copy()
    top = atom2.copy()
    c1 = np.linalg.norm(bottom.cell[2])
    c2 = np.linalg.norm(top.cell[2])
    cell1 = bottom.cell.copy()
    cell2 = top.cell.copy()
    cell1[2] /= c1
    cell2[2] /= c2
    cell = cell1 + weight * (cell2 - cell1)
    cell[2] /= np.linalg.norm(cell[2])
    cell1 = cell.copy()
    cell2 = cell.copy()
    cell1[2] *= c1
    cell2[2] *= c2

    bottom.set_cell(cell1, scale_atoms=True)
    top.set_cell(cell2, scale_atoms=True)

    zeroshift = np.min(bottom.get_positions()[:, 2])
    bottom.translate([0, 0, -zeroshift])
    zeroshift = np.min(top.get_positions()[:, 2])
    top.translate([0, 0, -zeroshift])
    bottom_thickness = np.max(bottom.get_positions()[:, 2]) - np.min(
        bottom.get_positions()[:, 2]
    )
    top.translate([0, 0, bottom_thickness])
    top.translate([0, 0, distance])
    bottom.extend(top)
    stack = recenter(bottom)
    return stack */

void translate_atoms_z(Atoms &atoms, const double shift)
{
    double2dvec_t pos1 = atoms.positions;
    double2dvec_t new_pos = {};
    for (int row = 0; row < pos1.size(); row++)
    {
        double1dvec_t subvec = {pos1[row][0],
                                pos1[row][1],
                                pos1[row][2] + shift};
        new_pos.push_back(subvec);
    };
    atoms.positions = new_pos;
};

std::tuple<double, double> get_min_max_z(Atoms &atoms)
{
    double2dvec_t pos = atoms.positions;
    double min_z = pos[0][2];
    double max_z = pos[pos.size()][2];
    for (int row = 0; row < pos.size(); row++)
    {
        double z = pos[row][2];
        if (z < min_z)
        {
            min_z = z;
        }
        if (z > max_z)
        {
            max_z = z;
        }
    }
    return std::make_tuple(min_z, max_z);
};

Atoms stack_atoms(Atoms &bottom, Atoms &top, double weight, double distance)
{
    auto [min_z1, max_z1] = get_min_max_z(bottom);
    auto [min_z2, max_z2] = get_min_max_z(bottom);
    translate_atoms_z(bottom, -min_z1);
    double bottom_thickness = max_z1 - min_z1;
    double top_thickness = max_z2 - min_z2;
    translate_atoms_z(top, -min_z2 + bottom_thickness + distance);

    double2dvec_t latticeA = bottom.lattice;
    double2dvec_t latticeB = top.lattice;
    double2dvec_t newcell(3, std::vector<double>(3, 0));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            newcell[i][j] = latticeA[i][j] + weight * (latticeB[i][j] - latticeA[i][j]);
        }
    }
    newcell[2][2] = bottom_thickness + top_thickness + distance + 50.0;
    bottom.scale_cell(newcell);
    top.scale_cell(newcell);
    Atoms stack = bottom + top;
    return stack;
};