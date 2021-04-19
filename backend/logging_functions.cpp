#include <vector>
#include <iostream>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Prints 2d vector.
template <typename T>
void print_2d_vector(const std::vector<std::vector<T>> &vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        for (int j = 0; j < vec[i].size(); j++)
        {
            std::cout << vec[i][j] << ' ';
        }
        std::cout << std::endl;
    }
};

template void print_2d_vector<int>(const std::vector<std::vector<int>> &vec);
template void print_2d_vector<double>(const std::vector<std::vector<double>> &vec);

// Prints number of OpenMP threads.
void log_number_of_threads()
{
#ifdef _OPENMP
    int limit, maxthreads;
    limit = omp_get_thread_limit();
    maxthreads = omp_get_max_threads();
    //printf("Limit is %d.\n", limit);
    //printf("Max is %d.\n", maxthreads);
#pragma omp parallel
    printf("Using %d OpenMP threads.\n", omp_get_num_threads());
#endif
};

// Prints map of single-valued double key with a 2d vector of ints.
void print_map_key_2d_vector(std::map<double, int2dvec_t> const &m)
{
    for (auto i = m.begin(); i != m.end(); ++i)
    {
        std::cout << "Key :" << (*i).first << std::endl;
        int2dvec_t matrix = (*i).second;
        print_2d_vector<int>(matrix);
        std::cout << std::endl;
    };
};
