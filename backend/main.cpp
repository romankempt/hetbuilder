#include <Python.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <map>

#include "spglib.h"
#include "symmetry.h"

#include "logging.h"
#include "math_functions.h"
#include "atom_functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cout, std::cin, std::endl;
using std::sin, std::cos, std::sqrt, std::pow, std::abs;
using std::vector;

typedef vector<vector<int>> int2dvec_t;

/**
 * Solves the equation |Am - R(theta)Bn| < tolerance for a given angle theta.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t findCoincidences(const double (&A)[2][2], const double (&B)[2][2], const double &theta, const int &Nmin, const int &Nmax, const double &tolerance)
{
    double Am[2] = {}, Bn[2] = {};
    int vecM[2] = {}, vecN[2] = {};
    double RBn[2] = {};
    double norm;
    int match;
    bool all_equal;
    int2dvec_t coincidences;

    const int nCombinations = (int)pow((Nmax - Nmin + 1), 4);
    cout << "Doing " << nCombinations << " combinations." << endl;

#pragma omp parallel for default(none) shared(A, B, theta, Nmin, Nmax, tolerance, Am, Bn, vecM, vecN, RBn, norm, match, all_equal, coincidences) schedule(static) ordered collapse(4)
    for (int i = Nmin; i < (Nmax + 1); i++)
    {
        for (int j = Nmin; j < (Nmax + 1); j++)
        {
            for (int k = Nmin; k < (Nmax + 1); k++)
            {
                for (int l = Nmin; l < (Nmax + 1); l++)
                {

                    vecM[0] = i;
                    vecM[1] = j;
                    vecN[0] = k;
                    vecN[1] = l;
                    BasisDotVector(A, vecM, Am);
                    BasisDotVector(B, vecN, Bn);
                    RotateVector(Bn, theta, RBn);
                    norm = getDistance(Am, RBn);
                    match = norm < tolerance;
                    all_equal = (i == j) && (j == k) && (k == l);
                    if (match && !all_equal)
                    {
                        vector<int> row{i, j, k, l};
#pragma omp ordered
                        coincidences.push_back(row);
                    };
                }
            }
        }
    }
    if (coincidences.size() > 0)
    {
        return coincidences;
    }
    else
    {
        return {};
    };
};

/**
 * Constructs the independent pairs (m1,m2,m3,m4) and (n1,n2,n3,n4).
 * 
 * Loop is OpenMP parallel.
 * 
 * All pairs with an absolute greatest common divisor different from 1 are removed,
 * because they correspond to scalar multiples of other smaller super cells.
 */
int2dvec_t findUniquePairs(int2dvec_t &coincidences)
{
    int2dvec_t uniquePairs;
    int m1, m2, m3, m4, n1, n2, n3, n4;
    int detM, detN;
    int _gcd;

#pragma omp parallel for shared(m1, m2, m3, m4, n1, n2, n3, n4, detM, detN, _gcd) schedule(static) ordered collapse(2)
    for (std::size_t i = 0; i < coincidences.size(); i++)
    {
        for (std::size_t j = 0; j < coincidences.size(); j++)
        {
            m1 = coincidences[i][0];
            m2 = coincidences[i][1];
            n1 = coincidences[i][2];
            n2 = coincidences[i][3];

            m3 = coincidences[j][0];
            m4 = coincidences[j][1];
            n3 = coincidences[j][2];
            n4 = coincidences[j][3];

            detM = m1 * m4 - m2 * m3;
            detN = n1 * n4 - n2 * n3;

            if ((detM > 0) && (detN > 0) && (j != i))
            {
                vector<int> subvec{m1, m2, m3, m4, n1, n2, n3, n4};
                _gcd = findGCD(subvec, 8);
                if (abs(_gcd) == 1)
                {
#pragma omp ordered
                    uniquePairs.push_back(subvec);
                };
            };
        };
    };

    // return {} if no unique pairs were found
    if (uniquePairs.size() > 0)
    {
        return uniquePairs;
    }
    else
    {
        return {};
    }
};

int main()
{
    double latticeA[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
    double latticeB[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};

    double position[][3] = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}};
    int types[] = {1, 1};
    int num_atom = 2, num_primitive_atom;
    double symprec = 1e-5;

    //test_spg_get_symmetry();
    int SuperCellMatrix[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
    vector<vector<double>> fracpoints;
    fracpoints = lattice_points_in_supercell(SuperCellMatrix);
    //print_2d_vector<int>(fracpoints);

    const double basisA[2][2] = {{latticeA[0][0], latticeA[0][1]}, {latticeA[1][0], latticeA[1][1]}};
    const double basisB[2][2] = {{latticeB[0][0], latticeB[0][1]}, {latticeB[1][0], latticeB[1][1]}};

    const int Nmax = 3;
    const int Nmin = -Nmax;
    const double tolerance = 0.01;

    double theta;
    double thetaMin, thetaMax;
    unsigned int numberOfAngles;
    int2dvec_t coincidences;
    int2dvec_t uniquePairs;

    std::map<double, int2dvec_t> AnglesMN;

#ifdef _OPENMP
    log_number_of_threads();
#endif

    numberOfAngles = 2;
    thetaMin = 0;
    thetaMax = 90;
    for (unsigned int i = 1; i <= numberOfAngles; i++)
    {
        theta = (thetaMax - thetaMin) * (i / (double)numberOfAngles);
        printf("i %d  angle %.2f of %d \n", i, theta, numberOfAngles);
        coincidences = findCoincidences(basisA, basisB, theta, Nmin, Nmax, tolerance);
        uniquePairs = findUniquePairs(coincidences);
        if (uniquePairs.size() > 0)
        {
            AnglesMN[theta] = uniquePairs;
        };
    };
    //print_map_key_2d_vector_of_ints(AnglesMN);

    //cout << "####################" << endl;
    //cout << "Coincidences: " << coincidences.size() << endl;
    //cout << "Unique Pairs " << uniquePairs.size() << endl;
    //print_2d_vector(results);
}