#pragma once
#include <vector>
#include <array>

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Function to return matrix vector product of basis A(2,2) with vector(2)
template <typename T1, typename T2>
std::vector<T1> basis_2x2_dot_2d_vector(const std::vector<std::vector<T1>> &basis, std::vector<T2> &vec);

// Function to rotate vector(2) by angle theta in degrees
template <typename T>
double1dvec_t rotate_2d_vector(std::vector<T> &vec, const double &theta);

// Returns distance |Am - RBn|
template <typename T>
double get_distance(std::vector<T> &Am, std::vector<T> &RBn);

// Function to return gcd of a and b
int get_gcd(int a, int b);

// Function to find gcd of array of numbers
int find_gcd(std::vector<int> &arr, int n);

// Function to perform dot product of row vector(3) times matrix(3,3)
template <typename T1, typename T2>
std::vector<T1> vec1x3_dot_3x3_matrix(std::vector<T1> &a, std::vector<std::vector<T2>> &matrix);

// Function to perform dot product of matrix(3,3) times column vector(3)
template <typename T1, typename T2>
std::vector<T1> matrix3x3_dot_vec3x1(std::vector<std::vector<T2>> &matrix, std::vector<T1> &a);

// Function to get determinant of 3x3 matrix
template <typename T>
double get_3x3_matrix_determinant(std::vector<std::vector<T>> &mat);

// Function to get inverse of 3x3 matrix
template <typename T>
double2dvec_t invert_3x3_matrix(std::vector<std::vector<T>> &mat);

/**
* This function multiplies two 3x3 matrices and returns a 3x3 matrix.
*/
template <typename T1, typename T2>
std::vector<std::vector<T2>> matrix3x3_dot_matrix3x3(std::vector<std::vector<T1>> &mat1, std::vector<std::vector<T2>> &mat2);

template <typename T>
std::vector<std::vector<T>> transpose_matrix3x3(std::vector<std::vector<T>> &mat);