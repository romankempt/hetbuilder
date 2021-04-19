#include <Python.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <map>

#include "spglib.h"
#include "symmetry.h"

#include "logging_functions.h"
#include "math_functions.h"
#include "atom_functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

/**
 * Solves the equation |Am - R(theta)Bn| < tolerance for a given angle theta.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t find_coincidences(const double2dvec_t &A, const double2dvec_t &B, const double &theta, const int &Nmin, const int &Nmax, const double &tolerance)
{

    int2dvec_t coincidences;
    const int nCombinations = (int)pow((Nmax - Nmin + 1), 4);
    std::cout << "Doing " << nCombinations << " combinations." << std::endl;

#pragma omp parallel for default(none) shared(A, B, theta, Nmin, Nmax, tolerance, coincidences) schedule(static) ordered collapse(4)
    for (int i = Nmin; i < (Nmax + 1); i++)
    {
        for (int j = Nmin; j < (Nmax + 1); j++)
        {
            for (int k = Nmin; k < (Nmax + 1); k++)
            {
                for (int l = Nmin; l < (Nmax + 1); l++)
                {
                    std::vector<int> vecM = {i, j};
                    std::vector<int> vecN = {k, l};
                    std::vector<double> Am;
                    std::vector<double> Bn;
                    std::vector<double> RBn;
                    double norm;
                    int match;
                    bool all_equal;
                    Am = basis_2x2_dot_2d_vector<double, int>(A, vecM);
                    Bn = basis_2x2_dot_2d_vector<double, int>(B, vecN);
                    RBn = rotate_2d_vector<double>(Bn, theta);
                    norm = get_distance<double>(Am, RBn);
                    match = norm < tolerance;
                    all_equal = (i == j) && (j == k) && (k == l);
                    if (match && !all_equal)
                    {
                        std::vector<int> row = {i, j, k, l};
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
int2dvec_t find_unique_pairs(int2dvec_t &coincidences)
{
    int2dvec_t uniquePairs;

#pragma omp parallel for shared(uniquePairs) schedule(static) ordered collapse(2)
    for (int i = 0; i < coincidences.size(); i++)
    {
        for (int j = 0; j < coincidences.size(); j++)
        {
            int m1 = coincidences[i][0];
            int m2 = coincidences[i][1];
            int n1 = coincidences[i][2];
            int n2 = coincidences[i][3];

            int m3 = coincidences[j][0];
            int m4 = coincidences[j][1];
            int n3 = coincidences[j][2];
            int n4 = coincidences[j][3];

            int detM = m1 * m4 - m2 * m3;
            int detN = n1 * n4 - n2 * n3;

            if ((detM > 0) && (detN > 0) && (j != i))
            {
                std::vector<int> subvec{m1, m2, m3, m4, n1, n2, n3, n4};
                int gcd = find_gcd(subvec, 8);
                if (abs(gcd) == 1)
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
    double2dvec_t latticeA = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
    double2dvec_t latticeB = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};

#ifdef _OPENMP
    log_number_of_threads();
#endif

    double2dvec_t positions = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}};
    std::vector<int> atomic_numbers = {1, 1};
    int num_atom = 2, num_primitive_atom;
    double symprec = 1e-5;

    Atoms atoms;
    atoms.positions = positions;
    atoms.atomic_numbers = atomic_numbers;
    atoms.num_atom = num_atom;
    atoms.lattice = latticeA;

    std::vector<int> spins(atomic_numbers.size(), 1);
    std::vector<int> equivs(atomic_numbers.size(), 1);
    atoms.spins = spins;
    atoms.equivalent_atoms = equivs;

    int2dvec_t SuperCellMatrix = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};

    Atoms supercell = make_supercell(atoms, SuperCellMatrix);
    Atoms rotatoms = {};
    rotatoms = rotate_atoms_around_z(atoms, 90);

    // double2dvec_t basisA = {{latticeA[0][0], latticeA[0][1]}, {latticeA[1][0], latticeA[1][1]}};
    // double2dvec_t basisB = {{latticeB[0][0], latticeB[0][1]}, {latticeB[1][0], latticeB[1][1]}};

    // int Nmax = 3;
    // int Nmin = -Nmax;
    // double tolerance = 0.01;
    // double theta;
    // double thetaMin, thetaMax;
    // unsigned int numberOfAngles;
    // int2dvec_t coincidences;
    // int2dvec_t uniquePairs;

    // std::map<double, int2dvec_t> AnglesMN;

    // numberOfAngles = 2;
    // thetaMin = 0;
    // thetaMax = 90;
    // for (unsigned int i = 1; i <= numberOfAngles; i++)
    // {
    //     theta = (thetaMax - thetaMin) * (i / (double)numberOfAngles);
    //     printf("i %d  angle %.2f of %d \n", i, theta, numberOfAngles);
    //     coincidences = find_coincidences(basisA, basisB, theta, Nmin, Nmax, tolerance);
    //     std::cout << "coin size: " << coincidences.size() << std::endl;
    //     uniquePairs = find_unique_pairs(coincidences);
    //     std::cout << "unpair size: " << uniquePairs.size() << std::endl;
    //     if (uniquePairs.size() > 0)
    //     {
    //         AnglesMN[theta] = uniquePairs;
    //     };
    // };
    //print_map_key_2d_vector(AnglesMN);
}