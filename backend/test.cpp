#include "coincidence_algorithm.h"

int test_coincidence_algorithm()
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

    //SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number)

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

    std::vector<Interface> stacks;
    std::vector<Interface> fstacks;
    stacks = build_all_supercells(bottom, top, AnglesMN, weight, distance, no_idealize, symprec, angle_tolerance);
    std::cout << "standardized stacks size: " << stacks.size() << std::endl;
    fstacks = filter_supercells(stacks);
    std::cout << "filtered stacks size: " << fstacks.size() << std::endl;
    for (int i = 0; i < fstacks.size(); i++)
    {
        Atoms example = fstacks[i].Stack;
        print_2d_vector(example.lattice);
        std::cout << std::endl;
    }

    return 0;
}