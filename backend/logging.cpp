#include <vector>
#include <iostream>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;

template <typename T>
void print_2d_vector(const vector<vector<T>> &vec)
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

void log_number_of_threads()
{
    int limit, maxthreads;
    limit = omp_get_thread_limit();
    maxthreads = omp_get_max_threads();
    printf("Limit is %d.\n", limit);
    printf("Max is %d.\n", maxthreads);

#pragma omp parallel
    printf("Hello from thread %d of %d .\n", omp_get_thread_num(), omp_get_num_threads());
};

template <typename K, typename V>
void print_map_key_2d_vector(map<T, vector<vector<V>>> const &m)
{
    for (auto i = m.begin(); i != m.end(); ++i)
    {
        cout << "Key :" << (*i).first << endl;
        vector<vector<V>> matrix = (*i).second;
        print_2d_vector(matrix);
        cout << endl;
    };
};