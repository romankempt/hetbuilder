#include "logging_functions.h"
#include "math_functions.h"
#include "atom_class.h"

#include "spglib.h"
#include "xtalcomp.h"

#include <algorithm>

// Prints lattice, positions, scaled positions and atomic numbers of atoms object.
void Atoms::print()
{
    std::cout << "Lattice: " << std::endl;
    print_2d_vector(this->lattice);
    std::cout << "Positions: " << std::endl;
    print_2d_vector(this->positions);
    std::cout << "Scaled Positions: " << std::endl;
    double2dvec_t scalpos = get_scaled_positions();
    print_2d_vector(scalpos);
    std::cout << "Magnetic Moments: '" << std::endl;
    for (auto j : this->magmoms)
    {
        std::cout << "\t" << j << " ";
    }
    std::cout << std::endl;
    std::cout << "Atomic Numbers: " << std::endl;
    for (auto j : this->atomic_numbers)
    {
        std::cout << "\t" << j << " ";
    }
    std::cout << std::endl;
};

// Returns the index mapping. Mostly useful for debugging the supercell to see which atoms were moved where.
int1dvec_t Atoms::get_index_mapping()
{
    return this->indices;
};

// Returns fractional coordinates.
double2dvec_t Atoms::get_scaled_positions()
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

// Converts fractional coordinates to cartesian coordinates.
double2dvec_t Atoms::scaled_positions_to_cartesian(double2dvec_t &scalpos)
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

// Scales cell of atoms object to size of newcell and adjusts cartesian atomic positions.
void Atoms::scale_cell(double2dvec_t &newcell)
{
    double2dvec_t scal_pos = this->get_scaled_positions();
    this->lattice = newcell;
    double2dvec_t cart_pos = this->scaled_positions_to_cartesian(scal_pos);
    this->positions = cart_pos;
};

// Overload + operator to add two Atoms objects.
Atoms Atoms::operator+(const Atoms &b)
{
    double2dvec_t pos1 = this->positions;
    double2dvec_t cell1 = this->lattice;
    int1dvec_t numbers1 = this->atomic_numbers;
    int1dvec_t indices1 = this->indices;
    double1dvec_t magmoms1 = this->magmoms;

    for (int row = 0; row < b.numAtom; row++)
    {
        pos1.push_back(b.positions[row]);
        numbers1.push_back(b.atomic_numbers[row]);
        indices1.push_back(b.indices[row]);
        magmoms1.push_back(b.magmoms[row]);
    }

    Atoms newAtoms(cell1, pos1, numbers1, indices1, magmoms1);
    return newAtoms;
};

// Helper function to convert double2dvec_t lattice to double array for spglib, which is transposed.
void Atoms::lattice_to_spglib_array(double arr[3][3])
{
    // transposition
    for (unsigned i = 0; (i < 3); i++)
    {
        for (unsigned j = 0; (j < 3); j++)
        {
            arr[j][i] = this->lattice[i][j];
        }
    }
};

// Helper function to convert double2dvec_t positions to double array for spglib.
void Atoms::positions_to_spglib_array(double arr[][3])
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

// Helper function to convert int1dvec_t atomic numbers to int array for spglib.
void Atoms::atomic_numbers_to_spglib_types(int arr[])
{
    for (unsigned i = 0; (i < this->numAtom); i++)
    {
        arr[i] = this->atomic_numbers[i];
    }
};

// Standardizes unit cell and positions via spglib. Returns spacegroup if successful, otherwise returns 0.
// Does not support magnetic moments. Magnetic moments are lost here.
int Atoms::standardize(int to_primitive, int no_idealize, double symprec, double angle_tolerance)
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
        double1dvec_t newMagneticMoments(newNumAtoms, 0.00);
        for (unsigned i = 0; (i < newNumAtoms); i++)
        {
            double1dvec_t subvec = {spglibPos[i][0], spglibPos[i][1], spglibPos[i][2]};
            spglibScalPos.push_back(subvec);
            spglibNewTypes[i] = spglibTypes[i];
            newMagneticMoments[i] = 0.00;
        }

        double2dvec_t cart_pos = scaled_positions_to_cartesian(spglibScalPos);
        this->positions = cart_pos;
        this->atomic_numbers = spglibNewTypes;
        this->magmoms = newMagneticMoments;
    }

    return spaceGroup;
}

// Helper function to convert double2dvec_t lattice to single vector for XtalComp.
XcMatrix Atoms::lattice_to_xtalcomp_cell()
{
    double arr[3][3];
    for (unsigned i = 0; (i < 3); i++)
    {
        for (unsigned j = 0; (j < 3); j++)
        {
            arr[i][j] = this->lattice[i][j];
        }
    }
    XcMatrix xtalcomp_cell(arr);
    return xtalcomp_cell;
};

// Helper function to convert int1dvec_t atomic numbers to int vector for XtalComp.
std::vector<unsigned int> Atoms::atomic_numbers_to_xtalcomp_types()
{
    std::vector<unsigned int> xtalcomp_types(std::begin(this->atomic_numbers), std::end(this->atomic_numbers));
    return xtalcomp_types;
};

// Helper function to convert double2dvec_t scaled positions to vector of XcVector positions for XtalComp.
std::vector<XcVector> Atoms::positions_to_xtalcomp_positions()
{
    std::vector<XcVector> xtalcomp_pos;
    XcVector subvector;
    double2dvec_t scalpos = this->get_scaled_positions();
    for (const auto &entries : scalpos)
    {
        subvector.set(entries[0], entries[1], entries[2]);
        xtalcomp_pos.push_back(subvector);
    }
    return xtalcomp_pos;
}

// Wrapper around the XtalComp Comparison Algorithm.
// Added additional check for magnetic moments.
bool Atoms::xtalcomp_compare(Atoms &other)
{
    XcMatrix cell1 = this->lattice_to_xtalcomp_cell();
    XcMatrix cell2 = other.lattice_to_xtalcomp_cell();
    std::vector<XcVector> pos1 = this->positions_to_xtalcomp_positions();
    std::vector<XcVector> pos2 = other.positions_to_xtalcomp_positions();
    std::vector<unsigned int> types1 = this->atomic_numbers_to_xtalcomp_types();
    std::vector<unsigned int> types2 = other.atomic_numbers_to_xtalcomp_types();
    // we compare if the sorted magmoms are identical
    double1dvec_t magmoms1 = this->magmoms;
    double1dvec_t magmoms2 = other.magmoms;

    std::sort(magmoms1.begin(), magmoms1.end());
    std::sort(magmoms2.begin(), magmoms2.end());
    bool match = false;
    if (magmoms1 == magmoms2)
    {
        match = XtalComp::compare(cell1, types1, pos1,
                                  cell2, types2, pos2,
                                  NULL, 0.05, 0.25, false);
        // std::cout << "XtalComp returns " << match << std::endl;
    }

    return match;
}