#pragma once

#include "math_functions.h"
#include <iostream>

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

class pBar
{
public:
    double neededProgress;
    std::string firstPartOfpBar;
    std::string lastPartOfpBar = "]",
                pBarFiller = "=",
                pBarUpdater = ">";
    bool show = false;

    pBar(double cneededProgress, std::string cfirstPartOfpBar, bool cshow)
    {
        neededProgress = cneededProgress;
        firstPartOfpBar = cfirstPartOfpBar;
        show = cshow;
    };

    pBar()
    {
        neededProgress = 100;
        firstPartOfpBar = "";
        show = false;
    };

    void update(double newProgress)
    {
        currentProgress += newProgress;
        amountOfFiller = (int)((currentProgress / neededProgress) * (double)pBarLength);
    }
    void print()
    {
        if (show)
        {
            currUpdateVal %= pBarUpdater.length();
            std::cout << "\r"             // Bring cursor to start of line
                      << firstPartOfpBar; // Print out first part of pBar
            for (int a = 0; a < amountOfFiller; a++)
            { // Print out current progress
                std::cout << pBarFiller;
            }
            std::cout << pBarUpdater[currUpdateVal];
            for (int b = 0; b < pBarLength - amountOfFiller; b++)
            { // Print out spaces
                std::cout << " ";
            }
            std::cout << lastPartOfpBar                                                  // Print out last part of progress bar
                      << " (" << (int)(100 * (currentProgress / neededProgress)) << "%)" // This just prints out the percent
                      << std::flush;
            currUpdateVal += 1;
        }
    }

private:
    int amountOfFiller,
        pBarLength = 50,        // I would recommend NOT changing this
        currUpdateVal = 0;      // Do not change
    double currentProgress = 0; // Do not change
                                // neededProgress = 100;   // I would recommend NOT changing this
};