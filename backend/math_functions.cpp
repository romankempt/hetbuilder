#include <cmath>
#include <vector>
#include <array>
#include <set>
#include <map>

#include "math_functions.h"

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Function to return matrix vector product of basis A(2,2) with vector(2)
template <typename T1, typename T2>
std::vector<T1> basis_2x2_dot_2d_vector(const std::vector<std::vector<T1>> &basis, std::vector<T2> &vec)
{
    std::vector<T1> result;
    for (int i = 0; i < 2; i++)
    {
        result.push_back(basis[i][0] * vec[0] + basis[i][1] * vec[1]);
    }
    return result;
};

template int1dvec_t basis_2x2_dot_2d_vector<int, int>(const int2dvec_t &basis, int1dvec_t &vec);
template double1dvec_t basis_2x2_dot_2d_vector<double, int>(const double2dvec_t &basis, int1dvec_t &vec);
template double1dvec_t basis_2x2_dot_2d_vector<double, double>(const double2dvec_t &basis, double1dvec_t &vec);

// Function to rotate vector(2) by angle theta in degrees
template <typename T>
double1dvec_t rotate_2d_vector(std::vector<T> &vec, const double &theta)
{
    double1dvec_t result = {0.0, 0.0};
    double t = theta * M_PI / 180.0;
    double R[2][2] = {{cos(t), -sin(t)}, {sin(t), cos(t)}};
    result[0] = R[0][0] * vec[0] + R[0][1] * vec[1];
    result[1] = R[1][0] * vec[0] + R[1][1] * vec[1];
    return result;
};

template double1dvec_t rotate_2d_vector<int>(int1dvec_t &vec, const double &theta);
template double1dvec_t rotate_2d_vector<double>(double1dvec_t &vec, const double &theta);

// Returns distance |Am - RBn|
template <typename T>
double get_distance(std::vector<T> &Am, std::vector<T> &RBn)
{
    double norm;
    norm = (Am[0] - RBn[0]) * (Am[0] - RBn[0]);
    norm += (Am[1] - RBn[1]) * (Am[1] - RBn[1]);
    norm = sqrt(norm);
    return norm;
};

template double get_distance<int>(int1dvec_t &Am, int1dvec_t &RBn);
template double get_distance<double>(double1dvec_t &Am, double1dvec_t &RBn);

// Function to return gcd of a and b
int get_gcd(int a, int b)
{
    if (a == 0)
        return b;
    return get_gcd(b % a, a);
}

// Function to find gcd of array of numbers
int find_gcd(int1dvec_t &arr, int n)
{
    int result = arr[0];
    for (int i = 1; i < n; i++)
    {
        result = get_gcd(arr[i], result);

        if (result == 1)
        {
            return 1;
        }
    }
    return result;
}

// Function to perform dot product of row vector(3) times matrix(3,3)
template <typename T1, typename T2>
std::vector<T1> vec1x3_dot_3x3_matrix(std::vector<T1> &a, std::vector<std::vector<T2>> &matrix)
{
    std::vector<T1> b(3, 0);
    for (int i = 0; i < a.size(); i++)
    {
        b[i] = a[0] * matrix[0][i] + a[1] * matrix[1][i] + a[2] * matrix[2][i];
    }
    return b;
};

template int1dvec_t vec1x3_dot_3x3_matrix<int, int>(int1dvec_t &a, int2dvec_t &matrix);
template double1dvec_t vec1x3_dot_3x3_matrix<double, int>(double1dvec_t &a, int2dvec_t &matrix);
template double1dvec_t vec1x3_dot_3x3_matrix<double, double>(double1dvec_t &a, double2dvec_t &matrix);

// Function to perform dot product of matrix(3,3) times column vector
template <typename T1, typename T2>
std::vector<T1> matrix3x3_dot_vec3x1(std::vector<std::vector<T2>> &matrix, std::vector<T1> &a)
{
    std::vector<T1> b(3, 0);
    for (int i = 0; i < a.size(); i++)
    {
        b[i] = a[0] * matrix[i][0] + a[1] * matrix[i][1] + a[2] * matrix[i][2];
    }
    return b;
};

template int1dvec_t matrix3x3_dot_vec3x1<int, int>(int2dvec_t &matrix, int1dvec_t &a);
template double1dvec_t matrix3x3_dot_vec3x1<double, int>(int2dvec_t &matrix, double1dvec_t &a);
template double1dvec_t matrix3x3_dot_vec3x1<double, double>(double2dvec_t &matrix, double1dvec_t &a);

// Function to get determinant of 3x3 matrix
template <typename T>
double get_3x3_matrix_determinant(std::vector<std::vector<T>> &mat)
{
    double determinant = 0;

    //finding determinant
    for (int i = 0; i < 3; i++)
        determinant = determinant + (mat[0][i] * (mat[1][(i + 1) % 3] * mat[2][(i + 2) % 3] - mat[1][(i + 2) % 3] * mat[2][(i + 1) % 3]));

    return determinant;
};

template double get_3x3_matrix_determinant<int>(int2dvec_t &mat);
template double get_3x3_matrix_determinant<double>(double2dvec_t &mat);

// Function to get inverse of 3x3 matrix
template <typename T>
double2dvec_t invert_3x3_matrix(std::vector<std::vector<T>> &mat)
{
    double determinant = get_3x3_matrix_determinant(mat);
    double2dvec_t minv(3, std::vector<double>(3, 0)); // inverse of matrix m
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            minv[i][j] = ((mat[(j + 1) % 3][(i + 1) % 3] * mat[(j + 2) % 3][(i + 2) % 3]) - (mat[(j + 1) % 3][(i + 2) % 3] * mat[(j + 2) % 3][(i + 1) % 3])) / determinant;
    }
    return minv;
}

template double2dvec_t invert_3x3_matrix<int>(int2dvec_t &mat);
template double2dvec_t invert_3x3_matrix<double>(double2dvec_t &mat);

/**
* This function multiplies two 3x3 matrices and returns a 3x3 matrix.
*/
template <typename T1, typename T2>
std::vector<std::vector<T2>> matrix3x3_dot_matrix3x3(std::vector<std::vector<T1>> &mat1, std::vector<std::vector<T2>> &mat2)
{
    std::vector<std::vector<T2>> res(3, std::vector<T2>(3, 0));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
                res[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
    return res;
};

template int2dvec_t matrix3x3_dot_matrix3x3<int, int>(int2dvec_t &mat1, int2dvec_t &mat2);
template double2dvec_t matrix3x3_dot_matrix3x3<int, double>(int2dvec_t &mat1, double2dvec_t &mat2);
template double2dvec_t matrix3x3_dot_matrix3x3<double, double>(double2dvec_t &mat1, double2dvec_t &mat2);

// Returns transpose of 3x3 matrix mat.
template <typename T>
std::vector<std::vector<T>> transpose_matrix3x3(std::vector<std::vector<T>> &mat)
{
    std::vector<std::vector<T>> outtrans(mat[0].size(),
                                         std::vector<T>(mat.size()));
    for (int i = 0; i < mat.size(); i++)
    {
        for (int j = 0; j < mat[0].size(); j++)
        {
            outtrans[j][i] = mat[i][j];
        }
    }
    return outtrans;
};

template int2dvec_t transpose_matrix3x3(int2dvec_t &mat);
template double2dvec_t transpose_matrix3x3(double2dvec_t &mat);