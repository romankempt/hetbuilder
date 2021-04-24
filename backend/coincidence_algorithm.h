#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <map>
#include <algorithm>
#include <set>

#include "spglib.h"

#include "logging_functions.h"
#include "math_functions.h"
#include "atom_class.h"
#include "atom_functions.h"

/**
 * Solves the equation |Am - R(theta)Bn| < tolerance for a given angle theta.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t find_coincidences(double2dvec_t &A, double2dvec_t &B, double &theta, int &Nmin, int &Nmax, double &tolerance);

/**
 * Constructs the independent pairs (m1,m2,m3,m4) and (n1,n2,n3,n4).
 * 
 * First loop is OpenMP parallel. Second one cannot be collapsed because j > i to avoid repititions.
 * 
 * All pairs with an absolute greatest common divisor different from 1 are removed,
 * because they correspond to scalar multiples of other smaller super cells.
 */
int2dvec_t find_unique_pairs(int2dvec_t &coincidences);

/**
 * Builds all supercells, applying the supercell matrices M and N and the Rotation R(theta).
 * 
 * The unit cell of the stack (interface) is given bei C = A + weight * (B - A).
 * The interfaces are standardized via spglib for the given symprec and angle_tolerance.
 * The loop over the supecell generation and standardization is OpenMP parallel.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> build_all_supercells(Atoms &bottom, Atoms &top, std::map<double, int2dvec_t> &AnglesMN, double &weight, double &distance, int &no_idealize, double &symprec, double &angle_tolerance);

/**
 * Filters the interfaces.
 * 
 * Interfaces are considered equal if their spacegroup, area and number of atoms matches.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> filter_supercells(std::vector<Interface> &stacks);

/**
 * Class definition of the lattice coincidence algorithm.
 * 
 * Executed by the run() method.
 */
class CoincidenceAlgorithm
{
public:
    Atoms primitive_bottom;
    Atoms primitive_top;
    int Nmax;
    int Nmin;
    double1dvec_t angles;
    double tolerance;
    double weight;
    double distance;
    int no_idealize;
    double symprec;
    double angle_tolerance;

    CoincidenceAlgorithm(Atoms cPrimitiveBottom,
                         Atoms cPrimitiveTop,
                         int cNmax = 10,
                         int cNmin = -10,
                         double1dvec_t cAngles = {0.0, 30.0, 60.0, 90.0},
                         double cTolerance = 0.01,
                         double cWeight = 0.5,
                         double cDistance = 4.0,
                         int cNoIdealize = 0,
                         double cSymPrec = 1e-5,
                         double cAngleTolerance = 5.0)
    {
        primitive_bottom = cPrimitiveBottom;
        primitive_top = cPrimitiveTop;
        Nmax = cNmax;
        Nmin = cNmin;
        angles = cAngles;
        tolerance = cTolerance;
        weight = cWeight;
        distance = cDistance;
        no_idealize = cNoIdealize;
        symprec = cSymPrec;
        angle_tolerance = cAngleTolerance;
    };

    std::vector<Interface> run()
    {
        int2dvec_t coincidences;
        std::map<double, int2dvec_t> AnglesMN;

        // basis is transposed
        double2dvec_t basisA = {{this->primitive_bottom.lattice[0][0], this->primitive_bottom.lattice[1][0]}, {this->primitive_bottom.lattice[0][1], this->primitive_bottom.lattice[1][1]}};
        double2dvec_t basisB = {{this->primitive_top.lattice[0][0], this->primitive_top.lattice[1][0]}, {this->primitive_top.lattice[0][1], this->primitive_top.lattice[1][1]}};

        for (int i = 0; i < this->angles.size(); i++)
        {
            double theta = this->angles[i];
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
            stacks = build_all_supercells(this->primitive_bottom,
                                          this->primitive_top,
                                          AnglesMN,
                                          this->weight,
                                          this->distance,
                                          this->no_idealize,
                                          this->symprec,
                                          this->angle_tolerance);
        }
        else
        {
            std::cerr << "Could not find any coincidence pairs." << std::endl;
            return {};
        }

        std::vector<Interface> fstacks;
        if (stacks.size() > 0)
        {
            fstacks = filter_supercells(stacks);
        }

        return fstacks;
    };
};