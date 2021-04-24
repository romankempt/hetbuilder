#include <vector>
#include <iostream>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "math_functions.h"
#include "logging_functions.h"

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Prints 1d vector.
template <typename T>
void print_1d_vector(std::vector<T> &vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << ' ';
    }
    std::cout << std::endl;
};

template void print_1d_vector<int>(int1dvec_t &vec);
template void print_1d_vector<double>(double1dvec_t &vec);

// Prints 2d vector.
template <typename T>
void print_2d_vector(const std::vector<std::vector<T>> &vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        for (int j = 0; j < vec[i].size(); j++)
        {
            std::cout << vec[i][j] << "  ";
        }
        std::cout << std::endl;
    }
};

template <typename T>
void print_2d_vector(std::vector<std::vector<T>> &vec)
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

template void print_2d_vector<int>(const int2dvec_t &vec);
template void print_2d_vector<double>(const double2dvec_t &vec);
template void print_2d_vector<int>(int2dvec_t &vec);
template void print_2d_vector<double>(double2dvec_t &vec);

// Prints number of OpenMP threads.
void log_number_of_threads()
{
#ifdef _OPENMP
    int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();

        /* Only master thread does this */
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Using %d OpenMP threads.\n", nthreads);
        }

    } /* All threads join master thread and disband */
#endif
};

// Get number of OpenMP threads.
int get_number_of_threads()
{
    int nthreads = 1;
#ifdef _OPENMP
    int tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel shared(nthreads, tid)
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();

        /* Only master thread does this */
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
        }

    } /* All threads join master thread and disband */
#endif
    return nthreads;
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
