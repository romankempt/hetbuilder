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

    void print()
    {
        std::cout << "Lattice: " << std::endl;
        print_2d_vector(this->lattice);
        std::cout << "Positions: " << std::endl;
        print_2d_vector(this->positions);
        std::cout << "Scaled Positions: " << std::endl;
        double2dvec_t scalpos = get_scaled_positions();
        print_2d_vector(scalpos);
        std::cout << "Atomic Numbers: " << std::endl;
        for (auto j : this->atomic_numbers)
        {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    };

    double2dvec_t get_scaled_positions()
    {
        double2dvec_t cell = this->lattice;
        double2dvec_t icell = invert_3x3_matrix(cell);
        double2dvec_t icellT = transpose_matrix3x3(icell);
        double2dvec_t scaled_positions;
        for (int row = 0; row < this->positions.size(); row++)
        {
            double1dvec_t subvec = matrix3x3_dot_vec3x1(icellT, this->positions[row]);
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
            double1dvec_t subvec = vec1x3_dot_3x3_matrix(scalpos[row], cell);
            cart_pos.push_back(subvec);
        }
        return cart_pos;
    };

    void scale_cell(double2dvec_t &newcell)
    {
        double2dvec_t scal_pos = this->get_scaled_positions();
        this->lattice = newcell;
        double2dvec_t cart_pos = this->scaled_positions_to_cartesian(scal_pos);
        this->positions = cart_pos;
    };

    // Overload + operator to add two Atoms objects.
    Atoms operator+(const Atoms &b)
    {
        double2dvec_t pos1 = this->positions;
        double2dvec_t cell1 = this->lattice;
        int1dvec_t numbers1 = this->atomic_numbers;

        for (int row = 0; row < b.numAtom; row++)
        {
            pos1.push_back(b.positions[row]);
            numbers1.push_back(b.atomic_numbers[row]);
        }

        Atoms newAtoms(cell1, pos1, numbers1);
        return newAtoms;
    };

    void lattice_to_spglib_array(double arr[3][3])
    {
        // transposition
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
        double2dvec_t scalpos = this->get_scaled_positions();
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
        for (unsigned i = 0; (i < this->numAtom); i++)
        {
            arr[i] = this->atomic_numbers[i];
        }
    };

    int standardize(int to_primitive = 1, int no_idealize = 0, double symprec = 1e-5, double angle_tolerance = 5.0)
    {
        double spglibPos[this->numAtom][3];
        positions_to_spglib_array(spglibPos);

        double spglibBasis[3][3];
        lattice_to_spglib_array(spglibBasis);

        int spglibTypes[this->numAtom];
        atomic_numbers_to_spglib_types(spglibTypes);

        int newNumAtoms = spgat_standardize_cell(spglibBasis,
                                                 spglibPos,
                                                 spglibTypes,
                                                 this->numAtom,
                                                 to_primitive,
                                                 no_idealize,
                                                 symprec,
                                                 angle_tolerance);
        int spaceGroup;
        char symbol[11];
        spaceGroup = spgat_get_international(symbol,
                                             spglibBasis,
                                             spglibPos,
                                             spglibTypes,
                                             this->numAtom,
                                             symprec,
                                             angle_tolerance);
        if (newNumAtoms != 0)
        {
            // transposition
            for (unsigned i = 0; (i < 3); i++)
            {
                for (unsigned j = 0; (j < 3); j++)
                {
                    this->lattice[j][i] = spglibBasis[i][j];
                }
            }
            this->numAtom = newNumAtoms;
            double2dvec_t spglibScalPos;
            int1dvec_t spglibNewTypes(newNumAtoms, 0);
            for (unsigned i = 0; (i < newNumAtoms); i++)
            {
                double1dvec_t subvec = {spglibPos[i][0], spglibPos[i][1], spglibPos[i][2]};
                spglibScalPos.push_back(subvec);
                spglibNewTypes[i] = spglibTypes[i];
            }

            double2dvec_t cart_pos = scaled_positions_to_cartesian(spglibScalPos);
            this->positions = cart_pos;
            this->atomic_numbers = spglibNewTypes;
        }

        return spaceGroup;
    }
};

class Interface
{
public:
    Atoms bottomLayer;
    Atoms topLayer;
    Atoms stack;
    double angle;
    int2dvec_t M;
    int2dvec_t N;
    int spaceGroup;

    friend bool operator==(const Interface &c1, const Interface &c2);
    friend bool operator>(const Interface &c1, const Interface &c2);
    friend bool operator>=(const Interface &c1, const Interface &c2);
    friend bool operator<(const Interface &c1, const Interface &c2);
    friend bool operator<=(const Interface &c1, const Interface &c2);

    Interface(Atoms cBottomLayer,
              Atoms cTopLayer,
              Atoms cStack,
              double cAngle,
              int2dvec_t cMatrixM,
              int2dvec_t cMatrixN,
              int cSpaceGroup)
    {
        bottomLayer = cBottomLayer;
        topLayer = cTopLayer;
        stack = cStack;
        angle = cAngle;
        M = cMatrixM;
        N = cMatrixN;
        spaceGroup = cSpaceGroup;
    };
};

inline bool operator==(const Interface &c1, const Interface &c2)
{
    bool spgmatch = (c1.spaceGroup == c2.spaceGroup);
    bool nummatch = (c1.stack.numAtom == c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(std::abs(area1) - std::abs(area2)) < 1e-6;
    bool equals = (spgmatch && nummatch && areamatch);
    return equals;
}

inline bool operator>(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom > c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) > std::abs(area2);
    bool greater_than = (nummatch && areamatch);
    return greater_than;
}

inline bool operator>=(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom >= c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) >= std::abs(area2);
    bool greater_than = (nummatch && areamatch);
    return greater_than;
}

inline bool operator<(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom < c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) < std::abs(area2);
    bool less_than = (nummatch && areamatch);
    return less_than;
}

inline bool operator<=(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom <= c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) <= std::abs(area2);
    bool less_than = (nummatch && areamatch);
    return less_than;
}
