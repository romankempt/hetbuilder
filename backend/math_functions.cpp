#include <cmath>
#include <vector>
#include <array>

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