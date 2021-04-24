#pragma once
#include "math_functions.h"

class Atoms
{

public:
    double2dvec_t lattice;
    double2dvec_t positions;
    int1dvec_t atomic_numbers;
    int numAtom;

    Atoms(double2dvec_t cLattice = {},
          double2dvec_t cPositions = {},
          int1dvec_t cAtomicNumbers = {})
    {
        lattice = cLattice;
        positions = cPositions;
        atomic_numbers = cAtomicNumbers;
        numAtom = atomic_numbers.size();
    };

    void print(void);

    double2dvec_t get_scaled_positions(void);

    double2dvec_t scaled_positions_to_cartesian(double2dvec_t &scalpos);

    void scale_cell(double2dvec_t &newcell);

    Atoms operator+(const Atoms &b);

    void lattice_to_spglib_array(double arr[3][3]);

    void positions_to_spglib_array(double arr[][3]);

    void atomic_numbers_to_spglib_types(int arr[]);

    int standardize(int to_primitive = 1, int no_idealize = 0, double symprec = 1e-5, double angle_tolerance = 5.0);
};
