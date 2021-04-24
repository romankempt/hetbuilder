#pragma once
#include <vector>
#include <array>

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

template <typename T1, typename T2>
std::vector<T1> basis_2x2_dot_2d_vector(const std::vector<std::vector<T1>> &basis, std::vector<T2> &vec);

template <typename T>
double1dvec_t rotate_2d_vector(std::vector<T> &vec, const double &theta);

template <typename T>
double get_distance(std::vector<T> &Am, std::vector<T> &RBn);

int get_gcd(int a, int b);

int find_gcd(std::vector<int> &arr, int n);

template <typename T1, typename T2>
std::vector<T1> vec1x3_dot_3x3_matrix(std::vector<T1> &a, std::vector<std::vector<T2>> &matrix);

template <typename T1, typename T2>
std::vector<T1> matrix3x3_dot_vec3x1(std::vector<std::vector<T2>> &matrix, std::vector<T1> &a);

template <typename T>
double get_3x3_matrix_determinant(std::vector<std::vector<T>> &mat);

template <typename T>
double2dvec_t invert_3x3_matrix(std::vector<std::vector<T>> &mat);

template <typename T1, typename T2>
std::vector<std::vector<T2>> matrix3x3_dot_matrix3x3(std::vector<std::vector<T1>> &mat1, std::vector<std::vector<T2>> &mat2);

template <typename T>
std::vector<std::vector<T>> transpose_matrix3x3(std::vector<std::vector<T>> &mat);