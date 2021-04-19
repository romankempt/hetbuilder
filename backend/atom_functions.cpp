#include <iostream>
#include <vector>
#include <cmath>
#include <array>

using std::cout, std::cin, std::endl;
using std::sin, std::cos, std::sqrt, std::pow, std::abs;
using std::vector;

typedef vector<vector<int>> int2dvec_t;

template <typename T>
std::vector<T> vec1x3_dot_3x3_matrix(std::vector<T> &a, std::vector<std::vector<T>> &matrix)
{
    std::vector<T> b(3, 0);
    for (int i = 0; i < a.size(); i++)
    {
        b[i] = a[0] * matrix[0][i] + a[1] * matrix[1][i] + a[2] * matrix[2][i];
    }
    return b;
};

template <typename T>
double get_3x3_matrix_determinant(vector<vector<T>> &mat)
{
    double determinant = 0;

    //finding determinant
    for (int i = 0; i < 3; i++)
        determinant = determinant + (mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] - mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]));

    return determinant;
};

template <typename T>
vector<vector<double>> invert_3x3_matrix(vector<vector<T>> &mat)
{
    double determinant = get_3x3_matrix_determinant(mat);
    vector<vector<double>> minv(3, vector<double>(3, 0)); // inverse of matrix m
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            minv[i][j] = ((mat[(j + 1) % 3][(i + 1) % 3] * mat[(j + 2) % 3][(i + 2) % 3]) - (mat[(j + 1) % 3][(i + 2) % 3] * mat[(j + 2) % 3][(i + 1) % 3])) / determinant;
    }
    return minv;
}

/**
 * Find all lattice points contained in a supercell.

    Adapted from pymatgen, which is available under MIT license:
    The MIT License (MIT) Copyright (c) 2011-2012 MIT & The Regents of the
    University of California, through Lawrence Berkeley National Laboratory
 */
vector<vector<double>> lattice_points_in_supercell(int2dvec_t &SuperCellMatrix)
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
    vector<int> dotproduct;
    for (int row = 0; row < diagonals.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix(diagonals[row], SuperCellMatrix);
        dpoints.push_back(dotproduct);
    }

    int k = dpoints.size() - 1;
    vector<int> mins = dpoints[0];
    vector<int> maxes = dpoints[k];

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
    vector<int> subvec(3, 0);
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
    vector<vector<double>> allpoints_double;
    for (int row = 0; row < allpoints.size(); row++)
    {
        vector<double> doubleVec(allpoints[row].begin(), allpoints[row].end());
        allpoints_double.push_back(doubleVec);
    };

    double determinant = get_3x3_matrix_determinant(SuperCellMatrix);
    vector<vector<double>> InvSuperCellMatrix = invert_3x3_matrix(SuperCellMatrix);
    vector<vector<double>> fracpoints;
    vector<double> dp;
    for (int row = 0; row < allpoints.size(); row++)
    {
        dp = vec1x3_dot_3x3_matrix(allpoints_double[row], InvSuperCellMatrix);
        fracpoints.push_back(dp);
    }

    vector<vector<double>> tvects;
    double fa, fb, fc;
    vector<double> fvec;
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
        cout << msg << endl;
        tvects = {};
    }

    return tvects;
};