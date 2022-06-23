#pragma once
#include "math_functions.h"
#include "xtalcomp.h"

class Atoms
{

public:
    double2dvec_t lattice;
    double2dvec_t positions;
    int1dvec_t atomic_numbers;
    int numAtom;
    int1dvec_t indices;
    double1dvec_t magmoms;

    // constructor if neither indices nor magmoms were provided
    Atoms(double2dvec_t cLattice,
          double2dvec_t cPositions,
          int1dvec_t cAtomicNumbers)
    {
        lattice = cLattice;
        positions = cPositions;
        atomic_numbers = cAtomicNumbers;
        numAtom = atomic_numbers.size();
        for (int i = 0; i < numAtom; i++)
        {
            indices.push_back(i);
            magmoms.push_back(0.0);
        }
    };

    // constructor if indices were not provided
    Atoms(double2dvec_t cLattice,
          double2dvec_t cPositions,
          int1dvec_t cAtomicNumbers,
          double1dvec_t cMagmoms)
    {
        lattice = cLattice;
        positions = cPositions;
        atomic_numbers = cAtomicNumbers;
        magmoms = cMagmoms;
        numAtom = atomic_numbers.size();
        for (int i = 0; i < numAtom; i++)
            indices.push_back(i);
    };

    // constructor if magmoms were not provided
    Atoms(double2dvec_t cLattice,
          double2dvec_t cPositions,
          int1dvec_t cAtomicNumbers,
          int1dvec_t cIndices)
    {
        lattice = cLattice;
        positions = cPositions;
        atomic_numbers = cAtomicNumbers;
        numAtom = atomic_numbers.size();
        indices = cIndices;
        for (int i = 0; i < numAtom; i++)
        {
            magmoms.push_back(0.0);
        }
    };

    // constructor if magmoms and indices were provided
    Atoms(double2dvec_t cLattice,
          double2dvec_t cPositions,
          int1dvec_t cAtomicNumbers,
          int1dvec_t cIndices,
          double1dvec_t cMagmoms)
    {
        lattice = cLattice;
        positions = cPositions;
        atomic_numbers = cAtomicNumbers;
        numAtom = atomic_numbers.size();
        indices = cIndices;
        magmoms = cMagmoms;
    };

    // empty constructur
    Atoms(){};

    void print(void);

    int1dvec_t get_index_mapping(void);

    double2dvec_t get_scaled_positions(void);

    double2dvec_t scaled_positions_to_cartesian(double2dvec_t &scalpos);

    void scale_cell(double2dvec_t &newcell);

    Atoms operator+(const Atoms &b);

    void lattice_to_spglib_array(double arr[3][3]);

    void positions_to_spglib_array(double arr[][3]);

    void atomic_numbers_to_spglib_types(int arr[]);

    int standardize(int to_primitive = 1, int no_idealize = 0, double symprec = 1e-5, double angle_tolerance = 5.0);

    XcMatrix lattice_to_xtalcomp_cell();

    std::vector<unsigned int> atomic_numbers_to_xtalcomp_types();

    std::vector<XcVector> positions_to_xtalcomp_positions();

    bool xtalcomp_compare(Atoms &other);
};
