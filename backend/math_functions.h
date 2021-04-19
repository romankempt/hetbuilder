#pragma once
#include <vector>
#include <array>

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Function to return matrix vector product of basis A(2,2) with vector(2)
template <typename T1, typename T2>
std::vector<T1> basis_2x2_dot_2d_vector(const std::vector<std::vector<T1>> &basis, std::vector<T2> &vec);

// Function to rotate vector(2) by angle theta in degrees
template <typename T>
std::vector<double> rotate_2d_vector(std::vector<T> &vec, const double &theta);

// Returns distance |Am - RBn|
template <typename T>
double get_distance(std::vector<T> &Am, std::vector<T> &RBn);

// Function to return gcd of a and b
int get_gcd(int a, int b);

// Function to find gcd of array of numbers
int find_gcd(std::vector<int> &arr, int n);

// Function to perform dot product of row vector(3) times matrix(3,3)
template <typename T>
std::vector<T> vec1x3_dot_3x3_matrix(std::vector<T> &a, std::vector<std::vector<T>> &matrix);

// Function to get determinant of 3x3 matrix
template <typename T>
double get_3x3_matrix_determinant(std::vector<std::vector<T>> &mat);

// Function to get inverse of 3x3 matrix
template <typename T>
std::vector<std::vector<double>> invert_3x3_matrix(std::vector<std::vector<T>> &mat);