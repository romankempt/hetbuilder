#pragma once

#include "math_functions.h"

class CoincidencePairs
{
public:
    int2dvec_t M;
    int2dvec_t N;

    CoincidencePairs(int2dvec_t &cM, int2dvec_t &cN)
    {
        M = cM;
        N = cN;
    };

    int score() const;
};

// Two coincidence matrices are qualitatively the same if they have yield the same area.
inline bool operator==(const CoincidencePairs &M1, const CoincidencePairs &M2)
{
    int2dvec_t m1 = M1.M;
    int2dvec_t m2 = M2.M;

    // I have to be careful to only compare M and M, not M and N
    int detm1 = m1[0][0] * m1[1][1] - m1[0][1] * m1[1][0];
    int detm2 = m2[0][0] * m2[1][1] - m2[0][1] * m2[1][0];
    bool same_area = std::abs(detm1 - detm2) < 1e-4;
    return same_area;
}

inline bool operator>(const CoincidencePairs &M1, const CoincidencePairs &M2)
{
    int score1 = M1.score();
    int score2 = M2.score();

    return (score1 < score2);
}

inline bool operator<(const CoincidencePairs &M1, const CoincidencePairs &M2)
{
    int score1 = M1.score();
    int score2 = M2.score();

    return (score1 > score2);
}