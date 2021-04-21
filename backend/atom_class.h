#pragma once
#include <vector>
#include <iostream>
#include <array>

#include "spglib.h"

#include "logging_functions.h"
#include "math_functions.h"

typedef std::vector<int> int1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

class Atoms
{

public:
    double2dvec_t lattice;
    double2dvec_t positions;
    int1dvec_t atomic_numbers;
    int num_atom;
    int1dvec_t spins;
    int1dvec_t equivalent_atoms;

    Atoms(double2dvec_t Lattice = {},
          double2dvec_t Positions = {},
          int1dvec_t AtomicNumbers = {},
          int1dvec_t Spins = {},
          int1dvec_t EquivalentAtoms = {})
    {
        lattice = Lattice;
        positions = Positions;
        atomic_numbers = AtomicNumbers;
        if (Spins.size() == 0)
        {
            int1dvec_t default_spins(AtomicNumbers.size(), 1);
            spins = default_spins;
        }
        else
        {
            spins = Spins;
        }
        if (EquivalentAtoms.size() == 0)
        {
            int1dvec_t default_equivs(AtomicNumbers.size(), 1);
            equivalent_atoms = default_equivs;
        }
        else
        {
            equivalent_atoms = EquivalentAtoms;
        }

        num_atom = atomic_numbers.size();
    };

    void print()
    {
        std::cout << "Lattice: " << std::endl;
        print_2d_vector(lattice);
        std::cout << "Positions: " << std::endl;
        print_2d_vector(positions);
        std::cout << "Scaled Positions: " << std::endl;
        double2dvec_t scalpos = get_scaled_positions();
        print_2d_vector(scalpos);
        std::cout << "Atomic numbers: " << std::endl;
        print_1d_vector(atomic_numbers);
        std::cout << "Atomic spins: " << std::endl;
        print_1d_vector(spins);
        std::cout << "Atomic equivalencies: " << std::endl;
        print_1d_vector(equivalent_atoms);
        std::cout << std::endl;
    };

    double2dvec_t get_scaled_positions()
    {
        double2dvec_t cell = this->lattice;
        double2dvec_t icell = invert_3x3_matrix(cell);
        double2dvec_t scaled_positions;
        for (int row = 0; row < positions.size(); row++)
        {
            double1dvec_t subvec = matrix3x3_dot_vec3x1(icell, positions[row]);
            scaled_positions.push_back(subvec);
        }
        return scaled_positions;
    };

    double2dvec_t scaled_positions_to_cartesian(double2dvec_t &scalpos)
    {
        double2dvec_t cell = this->lattice;
        double2dvec_t cart_pos;
        for (int row = 0; row < scalpos.size(); row++)
        {
            double1dvec_t subvec = matrix3x3_dot_vec3x1(cell, scalpos[row]);
            cart_pos.push_back(subvec);
        }
        return cart_pos;
    };

    void scale_cell(double2dvec_t &newcell)
    {
        double2dvec_t scal_pos = get_scaled_positions();
        this->lattice = newcell;
        double2dvec_t cart_pos = scaled_positions_to_cartesian(scal_pos);
        this->positions = cart_pos;
    };

    // Overload + operator to add two Atoms objects.
    Atoms operator+(const Atoms &b)
    {
        double2dvec_t pos1 = this->positions;
        double2dvec_t cell1 = this->lattice;
        int1dvec_t numbers1 = this->atomic_numbers;
        int1dvec_t spins1 = this->spins;
        int1dvec_t equiv1 = this->equivalent_atoms;

        for (int row = 0; row < b.num_atom; row++)
        {
            pos1.push_back(b.positions[row]);
            numbers1.push_back(b.atomic_numbers[row]);
            spins1.push_back(b.spins[row]);
            equiv1.push_back(b.equivalent_atoms[row]);
        }

        Atoms newAtoms(cell1, pos1, numbers1, spins1, equiv1);
        return newAtoms;
    };

    void lattice_to_spglib_array(double arr[3][3])
    {
        // note that the spglib lattice basis is transposed
        for (unsigned i = 0; (i < 3); i++)
        {
            for (unsigned j = 0; (j < 3); j++)
            {
                arr[j][i] = lattice[i][j];
            }
        }
    };

    void positions_to_spglib_array(double arr[][3])
    {
        double2dvec_t scalpos = get_scaled_positions();
        for (unsigned i = 0; i < scalpos.size(); i++)
        {
            for (unsigned j = 0; (j < 3); j++)
            {
                arr[i][j] = scalpos[i][j];
            }
        }
    };

    void atomic_numbers_to_spglib_types(int arr[])
    {
        for (unsigned i = 0; (i < num_atom); i++)
        {
            arr[i] = atomic_numbers[i];
        }
    };

    int standardize(const int to_primitive = 1, const int no_idealize = 0, const double symprec = 1e-5, const double angle_tolerance = 5.0)
    {
        double spglibBasis[3][3];
        lattice_to_spglib_array(spglibBasis);
        double spglibPos[num_atom][3];
        positions_to_spglib_array(spglibPos);
        int spglibTypes[num_atom];
        atomic_numbers_to_spglib_types(spglibTypes);

        int success = spgat_standardize_cell(spglibBasis,
                                             spglibPos,
                                             spglibTypes,
                                             num_atom,
                                             to_primitive,
                                             no_idealize,
                                             symprec,
                                             angle_tolerance);
        int spacegroup;
        char symbol[11];
        spacegroup = spgat_get_international(symbol,
                                             spglibBasis,
                                             spglibPos,
                                             spglibTypes,
                                             num_atom,
                                             symprec,
                                             angle_tolerance);
        if (success != 0)
        {
            for (unsigned i = 0; (i < 3); i++)
            {
                for (unsigned j = 0; (j < 3); j++)
                {
                    this->lattice[j][i] = spglibBasis[i][j];
                }
            }
            int arrSize = sizeof(spglibPos) / sizeof(spglibPos[0]);
            this->num_atom = arrSize;
            double2dvec_t spglibScalPos;
            int1dvec_t spglibNewTypes(num_atom, 0);
            for (unsigned i = 0; (i < num_atom); i++)
            {
                double1dvec_t subvec = {spglibPos[i][0], spglibPos[i][1], spglibPos[i][2]};
                spglibScalPos.push_back(subvec);
                spglibNewTypes[i] = spglibTypes[i];
            }

            double2dvec_t cart_pos = scaled_positions_to_cartesian(spglibScalPos);
            this->positions = cart_pos;
            this->atomic_numbers = spglibNewTypes;
        }

        return spacegroup;
    }
};

class Interface
{
public:
    Atoms BottomLayer;
    Atoms TopLayer;
    Atoms Stack;
    double Angle;
    int2dvec_t M;
    int2dvec_t N;
    int SpaceGroup;
    Interface(Atoms bottomLayer,
              Atoms topLayer,
              Atoms stack,
              double angle,
              int2dvec_t MatrixM,
              int2dvec_t MatrixN,
              int spaceGroup)
    {
        BottomLayer = bottomLayer;
        TopLayer = topLayer;
        Stack = stack;
        Angle = angle;
        M = MatrixM;
        N = MatrixN;
        SpaceGroup = spaceGroup;
    };

    bool operator==(Interface &other)
    {
        Atoms first = this->Stack;
        Atoms second = other.Stack;
        bool spgmatch = (this->SpaceGroup == other.SpaceGroup);
        bool nummatch = (first.num_atom == second.num_atom);
        double areaM1 = this->M[0][0] * this->M[1][1] - this->M[0][1] * this->M[1][0];
        double areaM2 = other.M[0][0] * other.M[1][1] - other.M[0][1] * other.M[1][0];
        bool areamatch = std::abs(areaM1 - areaM2) < 1e-6;
        bool is_equal = (spgmatch && nummatch && areamatch);
        return is_equal;
    };
};