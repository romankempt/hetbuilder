#include "helper_classes.h"

/**
 * Yields a score to allow for sorting coincidence pairs by how "arbitrarily nice" they are.
 * 
 * Prefers coincidence matrices that are symmetric and positive.
 */
int CoincidencePairs::score() const
{
    int sum = this->M[0][0] + this->M[0][1] + this->M[1][0] + this->M[1][1];
    int positivity = sum > 0;

    int offdiagsymm = (this->M[0][1] == this->M[1][0]);
    int diagsymm = (this->M[0][0] == this->M[1][1]);

    int score = positivity + offdiagsymm + diagsymm;
    return score;
};