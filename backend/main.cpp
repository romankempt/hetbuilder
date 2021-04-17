#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::cout, std::cin, std::endl;
using std::sin, std::cos, std::sqrt, std::pow, std::abs;
using std::vector;

typedef vector<vector<int>> int2dvec_t;

void BasisDotVector(const double (&basis)[2][2], const int (&vector)[2], double (&result)[2])
{
    result[0] = basis[0][0] * vector[0] + basis[0][1] * vector[1];
    result[1] = basis[1][0] * vector[0] + basis[1][1] * vector[1];
};

void RotateVector(const double (&vector)[2], const double &theta, double (&result)[2])
{
    double t = theta * 2 * M_PI / 180.0;
    double R[2][2] = {{cos(t), -sin(t)}, {sin(t), cos(t)}};
    result[0] = R[0][0] * vector[0] + R[0][1] * vector[1];
    result[1] = R[1][0] * vector[0] + R[1][1] * vector[1];
};

float getDistance(const double (&Am)[2], const double (&RBn)[2])
{
    double norm;
    norm = (Am[0] - RBn[0]) * (Am[0] - RBn[0]);
    norm += (Am[1] - RBn[1]) * (Am[1] - RBn[1]);
    norm = sqrt(norm);
    return norm;
};

// Function to return gcd of a and b
int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b % a, a);
}

// Function to find gcd of array of numbers
int findGCD(vector<int> &arr, int n)
{
    int result = arr[0];
    for (int i = 1; i < n; i++)
    {
        result = gcd(arr[i], result);

        if (result == 1)
        {
            return 1;
        }
    }
    return result;
}

/**
 * Solves the equation |Am - RBn| < tolerance.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t findCoincidences(const double (&A)[2][2], const double (&B)[2][2], const double &theta, const int &Nmin, const int &Nmax, const double &tolerance)
{
    double Am[2] = {}, Bn[2] = {};
    int vecM[2] = {}, vecN[2] = {};
    double RBn[2] = {};
    double norm;
    int match;
    bool all_equal;
    int2dvec_t coincidences;

    const int nCombinations = (int)pow((Nmax - Nmin + 1), 4);
    cout << "Doing " << nCombinations << " combinations." << endl;

    //omp_set_num_threads(2);
#pragma omp parallel for default(none) shared(A, B, theta, Nmin, Nmax, tolerance, Am, Bn, vecM, vecN, RBn, norm, match, all_equal, coincidences) schedule(dynamic) ordered collapse(4)
    for (int i = Nmin; i < (Nmax + 1); i++)
    {
        for (int j = Nmin; j < (Nmax + 1); j++)
        {
            for (int k = Nmin; k < (Nmax + 1); k++)
            {
                for (int l = Nmin; l < (Nmax + 1); l++)
                {

                    vecM[0] = i;
                    vecM[1] = j;
                    vecN[0] = k;
                    vecN[1] = l;
                    BasisDotVector(A, vecM, Am);
                    BasisDotVector(B, vecN, Bn);
                    RotateVector(Bn, theta, RBn);
                    norm = getDistance(Am, RBn);
                    match = norm < tolerance;
                    all_equal = (i == j) && (j == k) && (k == l);
                    if (match && !all_equal)
                    {
                        vector<int> row{i, j, k, l};
#pragma omp ordered
                        coincidences.push_back(row);
                    };
                }
            }
        }
    }
    if (coincidences.size() > 0)
    {
        return coincidences;
    }
    else
    {
        return {};
    };
};

int2dvec_t findUniquePairs(int2dvec_t &coincidences)
{
    int2dvec_t uniquePairs;
    int m1, m2, m3, m4, n1, n2, n3, n4;
    int detM, detN;
    int _gcd;

#pragma omp parallel for shared(m1, m2, m3, m4, n1, n2, n3, n4, detM, detN, _gcd) schedule(dynamic) ordered collapse(2)
    for (std::size_t i = 0; i < coincidences.size(); ++i)
    {
        for (std::size_t j = 0; j < coincidences.size(); ++j)
        {
            m1 = coincidences[i][0];
            m2 = coincidences[i][1];
            n1 = coincidences[i][3];
            n2 = coincidences[i][4];
            m3 = coincidences[j][0];
            m4 = coincidences[j][1];
            n3 = coincidences[j][3];
            n4 = coincidences[j][4];
            detM = m1 * m4 - m2 * m3;
            detN = n1 * n4 - n2 * n3;

            if ((detM > 0) && (detN > 0) && (j != i))
            {
                vector<int> subvec{m1, m2, m3, m4, n1, n2, n3, n4};
                _gcd = findGCD(subvec, 8);
                if (_gcd == 1)
                {
#pragma omp ordered
                    uniquePairs.push_back(subvec);
                };
            };
        };
    };

    // return {} if no unique pairs were found
    if (uniquePairs.size() > 0)
    {
        return uniquePairs;
    }
    else
    {
        return {};
    }
};

template <typename T>
void print_2d_vector(const vector<vector<T>> &vec)
{
    for (vector<vector<float>>::size_type i = 0; i < vec.size(); i++)
    {
        for (vector<float>::size_type j = 0; j < vec[i].size(); j++)
        {
            std::cout << vec[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

void log_numer_of_threads()
{
    int limit, maxthreads;
    limit = omp_get_thread_limit();
    maxthreads = omp_get_max_threads();
    printf("Limit is %d.\n", limit);
    printf("Max is %d.\n", maxthreads);

#pragma omp parallel
    printf("Hello from thread %d of %d .\n", omp_get_thread_num(), omp_get_num_threads());
};

int main()
{
    const double basisA[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
    const double basisB[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    const int Nmax = 10;
    const int Nmin = -Nmax;
    const double tolerance = 0.01;

    double theta;
    double thetaMin, thetaMax;
    unsigned int numberOfAngles;
    int2dvec_t coincidences;
    int2dvec_t uniquePairs;

    std::map<double, int2dvec_t> AnglesMN;

    log_numer_of_threads();

    numberOfAngles = 10;
    thetaMin = 0;
    thetaMax = 90;
    for (unsigned int i = 1; i <= numberOfAngles; i++)
    {
        theta = (thetaMax - thetaMin) * (i / (double)numberOfAngles);
        printf("i %d  angle %.2f of %d \n", i, theta, numberOfAngles);
        coincidences = findCoincidences(basisA, basisB, theta, Nmin, Nmax, tolerance);
        uniquePairs = findUniquePairs(coincidences);
        AnglesMN[theta] = uniquePairs;
    };

    //cout << "####################" << endl;
    //cout << "Coincidences: " << coincidences.size() << endl;
    //cout << "Unique Pairs " << uniquePairs.size() << endl;
    //print_2d_vector(results);
}