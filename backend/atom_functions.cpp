#include <iostream>
#include <vector>
#include <cmath>
#include <array>

#include "math_functions.h"
#include "logging_functions.h"

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

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
    std::vector<int> dotproduct;
    for (int row = 0; row < diagonals.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<int, int>(diagonals[row], SuperCellMatrix);
        dpoints.push_back(dotproduct);
    }

    int k = dpoints.size() - 1;
    std::vector<int> mins = dpoints[0];
    std::vector<int> maxes = dpoints[k];

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
    std::vector<int> subvec(3, 0);
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
        std::vector<double> doubleVec(allpoints[row].begin(), allpoints[row].end());
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
    std::vector<double> fvec;
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

struct Atoms
{
    double2dvec_t lattice;
    double2dvec_t positions;
    std::vector<int> atomic_numbers;
    int num_atom;
    std::vector<int> spins;
    std::vector<int> equivalent_atoms;
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
    std::vector<double> dotproduct;
    for (int row = 0; row < fracpoints.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<double, double>(fracpoints[row], supercell);
        lattice_points.push_back(dotproduct);
    }

    double2dvec_t new_positions = prim.positions;
    std::vector<int> new_numbers = prim.atomic_numbers;
    std::vector<int> new_spin = prim.spins;
    std::vector<int> new_equiv = prim.equivalent_atoms;
    for (int i = 0; i < prim.num_atom; i++)
    {
        std::vector<double> atom_pos = prim.positions[i];
        int number = prim.atomic_numbers[i];
        int spin = prim.spins[i];
        int equivalent = prim.equivalent_atoms[i];
        for (int row = 0; row < lattice_points.size(); row++)
        {
            std::vector<double> lp = lattice_points[row];
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
    Atoms superatoms = {};
    superatoms.lattice = supercell;
    superatoms.positions = new_positions;
    superatoms.atomic_numbers = new_numbers;
    superatoms.num_atom = new_numbers.size();
    superatoms.equivalent_atoms = new_equiv;
    superatoms.spins = new_spin;
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

    Atoms rotatoms = {};
    rotatoms.positions = rotpositions;
    rotatoms.lattice = rotcell;
    rotatoms.num_atom = atoms.num_atom;
    rotatoms.atomic_numbers = atoms.atomic_numbers;
    rotatoms.spins = atoms.spins;
    rotatoms.equivalent_atoms = atoms.equivalent_atoms;

    return rotatoms;
};