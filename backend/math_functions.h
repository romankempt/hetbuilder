#pragma once
#include <cmath>
#include <vector>
#include <array>

using std::vector;

typedef vector<vector<int>> int2dvec_t;

void BasisDotVector(const double (&basis)[2][2], const int (&vector)[2], double (&result)[2]);

void RotateVector(const double (&vector)[2], const double &theta, double (&result)[2]);

float getDistance(const double (&Am)[2], const double (&RBn)[2]);

// Function to return gcd of a and b
int gcd(int a, int b);

// Function to find gcd of array of numbers
int findGCD(vector<int> &arr, int n);