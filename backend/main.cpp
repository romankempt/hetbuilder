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
#include "atom_class.h"
#include "atom_functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
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
                    int1dvec_t vecM = {i, j};
                    int1dvec_t vecN = {k, l};
                    double1dvec_t Am;
                    double1dvec_t Bn;
                    double1dvec_t RBn;
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
                        int1dvec_t row = {i, j, k, l};
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
 * First loop is OpenMP parallel. Second one cannot be collapsed because j > i to avoid repititions.
 * 
 * All pairs with an absolute greatest common divisor different from 1 are removed,
 * because they correspond to scalar multiples of other smaller super cells.
 */
int2dvec_t find_unique_pairs(int2dvec_t &coincidences)
{
    int2dvec_t uniquePairs;

#pragma omp parallel for shared(uniquePairs) schedule(static) ordered collapse(1)
    for (int i = 0; i < coincidences.size(); i++)
    {
        for (int j = i + 1; j < coincidences.size(); j++)
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

            if ((detM > 0) && (detN > 0))
            {
                int1dvec_t subvec{m1, m2, m3, m4, n1, n2, n3, n4};
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

struct Stack
{
    Atoms bottomLayer;
    Atoms topLayer;
    Atoms interface;
    double angle;
    int2dvec_t M;
    int2dvec_t N;
    //backtansform matrices somehow from spglib datasets?
};

std::vector<Stack> build_all_supercells(Atoms &bottom, Atoms &top, std::map<double, int2dvec_t> &AnglesMN, double &weight, double &distance, const int &no_idealize, const double &symprec, const double &angle_tolerance)
{
    std::vector<Stack> stacks;
    for (auto i = AnglesMN.begin(); i != AnglesMN.end(); ++i)
    {
        double theta = (*i).first;
        int2dvec_t pairs = (*i).second;
#pragma omp parallel for shared(stacks, AnglesMN, theta, pairs) schedule(static) ordered collapse(1)
        for (int j = 0; j < pairs.size(); j++)
        {
            int1dvec_t row = pairs[j];
            int2dvec_t M = {{row[0], row[1], 0}, {row[2], row[3], 0}, {0, 0, 1}};
            int2dvec_t N = {{row[4], row[5], 0}, {row[6], row[7], 0}, {0, 0, 1}};
            Atoms bottomLayer = make_supercell(bottom, M);
            Atoms topLayer = make_supercell(top, N);
            Atoms topLayerRot = rotate_atoms_around_z(topLayer, theta);
            Atoms interface = stack_atoms(bottomLayer, topLayerRot, weight, distance);
            int success = interface.standardize(1, no_idealize, symprec, angle_tolerance);
            if (success != 0)
            {
                Stack stack = {};
                stack.bottomLayer = bottomLayer;
                stack.topLayer = topLayerRot;
                stack.interface = interface;
                stack.angle = theta;
                stack.M = M;
                stack.N = N;
#pragma omp ordered
                stacks.push_back(stack);
            }
        };
    };
    return stacks;
};

int main()
{
    double2dvec_t latticeA = {{4, 1, 0}, {1, 4, 0}, {0, 0, 4}};
    double2dvec_t latticeB = {{4, 1, 0}, {1, 4, 0}, {0, 0, 4}};

#ifdef _OPENMP
    log_number_of_threads();
#endif

    double2dvec_t positions = {
        {0, 0, 0},
        {0.5, 0.5, 0.5}};
    int1dvec_t atomic_numbers = {1, 1};
    int num_atom = 2, num_primitive_atom;

    Atoms bottom(latticeA, positions, atomic_numbers);
    Atoms top(latticeB, positions, atomic_numbers);

    int2dvec_t SuperCellMatrix = {{2, -2, 0}, {-1, 2, 0}, {0, 0, 2}};

    double2dvec_t basisA = {{latticeA[0][0], latticeA[0][1]}, {latticeA[1][0], latticeA[1][1]}};
    double2dvec_t basisB = {{latticeB[0][0], latticeB[0][1]}, {latticeB[1][0], latticeB[1][1]}};

    int Nmax = 3;
    int Nmin = 1;
    double tolerance = 0.01;
    double weight = 0.5;
    double distance = 4.0;
    const int no_idealize = 0;
    const double symprec = 1e-5;
    const double angle_tolerance = 5.0;

    double1dvec_t thetas = {0.00, 45.0, 90.0};
    int2dvec_t coincidences;
    int2dvec_t uniquePairs;

    std::map<double, int2dvec_t> AnglesMN;

    for (int i = 0; i < thetas.size(); i++)
    {
        double theta = thetas[i];
        printf("i %d  angle %.2f of %d \n", i, theta, thetas.size());
        coincidences = find_coincidences(basisA, basisB, theta, Nmin, Nmax, tolerance);
        std::cout << "coin size: " << coincidences.size() << std::endl;
        uniquePairs = find_unique_pairs(coincidences);
        std::cout << "unpair size: " << uniquePairs.size() << std::endl;
        if (uniquePairs.size() > 0)
        {
            AnglesMN.insert(std::make_pair(theta, uniquePairs));
        };
    };
    //print_map_key_2d_vector(AnglesMN);

    std::vector<Stack> stacks;
    stacks = build_all_supercells(bottom, top, AnglesMN, weight, distance, no_idealize, symprec, angle_tolerance);
    std::cout << "standardized stacks size: " << stacks.size() << std::endl;
    for (int i = 0; i < stacks.size(); i++)
    {
        Atoms example = stacks[i].interface;
        example.print();
    }

    // next is filtering
}