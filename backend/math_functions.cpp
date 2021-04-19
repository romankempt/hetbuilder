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

template <typename T>
std::vector<T> vec1x3_dot_3x3_matrix(std::vector<T> &a, std::vector<std::vector<T>> &matrix)
{
    std::vector<T> b(3, 0);
    for (int i = 0; i < a.size(); i++)
    {
        b[i] = a[0] * matrix[0][i] + a[1] * matrix[1][i] + a[2] * matrix[2][i];
    }
    return b;
};

template <typename T>
double get_3x3_matrix_determinant(vector<vector<T>> &mat)
{
    double determinant = 0;

    //finding determinant
    for (int i = 0; i < 3; i++)
        determinant = determinant + (mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] - mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]));

    return determinant;
};

template <typename T>
vector<vector<double>> invert_3x3_matrix(vector<vector<T>> &mat)
{
    double determinant = get_3x3_matrix_determinant(mat);
    vector<vector<double>> minv(3, vector<double>(3, 0)); // inverse of matrix m
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            minv[i][j] = ((mat[(j + 1) % 3][(i + 1) % 3] * mat[(j + 2) % 3][(i + 2) % 3]) - (mat[(j + 1) % 3][(i + 2) % 3] * mat[(j + 2) % 3][(i + 1) % 3])) / determinant;
    }
    return minv;
}