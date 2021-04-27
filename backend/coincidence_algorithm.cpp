#include <set>

#include "logging_functions.h"
#include "math_functions.h"
#include "atom_class.h"
#include "atom_functions.h"
#include "helper_classes.h"
#include "interface_class.h"
#include "coincidence_algorithm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::map<double, std::vector<CoincidencePairs>> angle_dict_t;

/**
 * Solves the equation |Am - R(theta)Bn| < tolerance for a given angle theta.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t CoincidenceAlgorithm::find_coincidences(double2dvec_t &A, double2dvec_t &B, double &theta, int &Nmin, int &Nmax, double &tolerance)
{

    int2dvec_t coincidences;

#pragma omp parallel for default(none) shared(A, B, theta, Nmin, Nmax, tolerance, coincidences) schedule(static) ordered collapse(4)
    for (int i = -Nmax; i < (Nmax + 1); i++)
    {
        for (int j = -Nmax; j < (Nmax + 1); j++)
        {
            for (int k = -Nmax; k < (Nmax + 1); k++)
            {
                for (int l = -Nmax; l < (Nmax + 1); l++)
                {
                    if ((std::abs(i) > Nmin) && (std::abs(j) > Nmin) && (std::abs(k) > Nmin) && (std::abs(l) > Nmin))
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
                        }
                    }
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
int2dvec_t CoincidenceAlgorithm::find_unique_pairs(int2dvec_t &coincidences)
{
    int2dvec_t uniquePairs;

#pragma omp parallel for shared(uniquePairs, coincidences) schedule(static) ordered collapse(2)
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

/**
 * Reduces the number of unique pairs by filtering out the pairs with the same determinant.
 * 
 * Pairs with positive, symmetric entries are preferred.
 */
angle_dict_t CoincidenceAlgorithm::reduce_unique_pairs(std::map<double, int2dvec_t> &AnglesMN)
{
    angle_dict_t fAnglesMN;
    for (auto i = AnglesMN.begin(); i != AnglesMN.end(); ++i)
    {
        double theta = (*i).first;
        int2dvec_t pairs = (*i).second;
        std::vector<CoincidencePairs> CoinPairs;
        for (int j = 0; j < pairs.size(); j++)
        {
            int1dvec_t row = pairs[j];
            int2dvec_t M = {{row[0], row[1], 0},
                            {row[2], row[3], 0},
                            {0, 0, 1}};
            int2dvec_t N = {{row[4], row[5], 0},
                            {row[6], row[7], 0},
                            {0, 0, 1}};
            CoincidencePairs pair(M, N);
            CoinPairs.push_back(pair);
        }
        std::set<CoincidencePairs> s(CoinPairs.begin(), CoinPairs.end());
        std::vector<CoincidencePairs> v(s.begin(), s.end());
        fAnglesMN.insert(std::make_pair(theta, v));
    }
    return fAnglesMN;
};

/**
 * Builds all supercells, applying the supercell matrices M and N and the Rotation R(theta).
 * 
 * The unit cell of the stack (interface) is given bei C = A + weight * (B - A).
 * The interfaces are standardized via spglib for the given symprec and angle_tolerance.
 * The loop over the supecell generation and standardization is OpenMP parallel.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> CoincidenceAlgorithm::build_all_supercells(Atoms &bottom, Atoms &top, angle_dict_t &AnglesMN,
                                                                  double &weight, double &distance, int &no_idealize, double &symprec, double &angle_tolerance)
{
    std::vector<Interface> stacks;
    for (auto i = AnglesMN.begin(); i != AnglesMN.end(); ++i)
    {
        double theta = (*i).first;
        std::vector<CoincidencePairs> pairs = (*i).second;
#pragma omp parallel for shared(bottom, top, stacks, AnglesMN, theta, pairs) schedule(static) ordered collapse(1)
        for (int j = 0; j < pairs.size(); j++)
        {
            CoincidencePairs row = pairs[j];
            int2dvec_t M = row.M;
            int2dvec_t N = row.N;
            Atoms bottomLayer = make_supercell(bottom, M);
            Atoms topLayer = make_supercell(top, N);
            Atoms topLayerRot = rotate_atoms_around_z(topLayer, theta);
            Atoms interface = stack_atoms(bottomLayer, topLayerRot, weight, distance);
            int spacegroup = interface.standardize(1, no_idealize, symprec, angle_tolerance);
            if (spacegroup != 0)
            {
                Interface stack(bottomLayer, topLayerRot, interface, theta, M, N, spacegroup);
#pragma omp ordered
                stacks.push_back(stack);
            }
        };
    };
    return stacks;
};

/**
 * Filters the interfaces.
 * 
 * Interfaces are considered equal if their spacegroup, area and number of atoms matches.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> CoincidenceAlgorithm::filter_supercells(std::vector<Interface> &stacks)
{
    std::set<Interface> s(stacks.begin(), stacks.end());
    std::vector<Interface> v(s.begin(), s.end());

    return v;
};

/**
 * Executes the coincidence lattice search algorithm for given parameters.
 */
std::vector<Interface> CoincidenceAlgorithm::run(int Nmax,
                                                 int Nmin,
                                                 double1dvec_t angles,
                                                 double tolerance,
                                                 double weight,
                                                 double distance,
                                                 int no_idealize,
                                                 double symprec,
                                                 double angle_tolerance)
{
    int2dvec_t coincidences;
    std::map<double, int2dvec_t> AnglesMN;

    // basis is transposed
    double2dvec_t basisA = {{this->primitive_bottom.lattice[0][0], this->primitive_bottom.lattice[1][0]}, {this->primitive_bottom.lattice[0][1], this->primitive_bottom.lattice[1][1]}};
    double2dvec_t basisB = {{this->primitive_top.lattice[0][0], this->primitive_top.lattice[1][0]}, {this->primitive_top.lattice[0][1], this->primitive_top.lattice[1][1]}};

    for (int i = 0; i < angles.size(); i++)
    {
        double theta = angles[i];
        coincidences = find_coincidences(basisA, basisB, theta, Nmin, Nmax, tolerance);
        if (coincidences.size() > 0)
        {
            int2dvec_t uniquePairs;
            uniquePairs = find_unique_pairs(coincidences);
            if (uniquePairs.size() > 0)
            {
                AnglesMN.insert(std::make_pair(theta, uniquePairs));
            }
        };
    };

    std::vector<Interface> stacks;
    if (AnglesMN.size() > 0)
    {
        angle_dict_t fAnglesMN;
        fAnglesMN = reduce_unique_pairs(AnglesMN);
        stacks = build_all_supercells(this->primitive_bottom,
                                      this->primitive_top,
                                      fAnglesMN,
                                      weight,
                                      distance,
                                      no_idealize,
                                      symprec,
                                      angle_tolerance);
    }
    else
    {
        //std::cerr << "Could not find any coincidence pairs." << std::endl;
        return {};
    }

    std::vector<Interface> fstacks;
    if (stacks.size() > 0)
    {
        fstacks = filter_supercells(stacks);
    }

    return fstacks;
};