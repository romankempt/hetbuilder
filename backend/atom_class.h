#pragma once
#include <vector>
#include <iostream>

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
    }
};