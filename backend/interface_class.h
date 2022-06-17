#pragma once
#include "atom_class.h"

/**
 * Container class for an interface after the coincidence lattice search.
 *
 * Contains the bottom layer and top layer in supercell form, as well as the
 * stack thereof.
 *
 * Comparison operators are defined to allow for ordering and comparisons in a set.
 */
class Interface
{
public:
    Atoms bottomLayer;
    Atoms topLayer;
    Atoms stack;
    double angle;
    int2dvec_t M;
    int2dvec_t N;
    int spaceGroup;

    friend bool operator==(const Interface &c1, const Interface &c2);
    friend bool operator>(const Interface &c1, const Interface &c2);
    friend bool operator<(const Interface &c1, const Interface &c2);

    Interface(Atoms cBottomLayer,
              Atoms cTopLayer,
              Atoms cStack,
              double cAngle,
              int2dvec_t cMatrixM,
              int2dvec_t cMatrixN,
              int cspaceGroup)
    {
        bottomLayer = cBottomLayer;
        topLayer = cTopLayer;
        stack = cStack;
        angle = cAngle;
        M = cMatrixM;
        N = cMatrixN;
        spaceGroup = cspaceGroup;
    };

    void set_stack(Atoms &stack)
    {
        this->stack = stack;
    }

    void set_spacegroup(int &sg)
    {
        this->spaceGroup = sg;
    }
};

inline bool operator==(const Interface &c1, const Interface &c2)
{
    // bool spgmatch = (c1.spaceGroup == c2.spaceGroup);
    bool nummatch = (c1.stack.numAtom == c2.stack.numAtom);
    // bool anglematch = std::abs(c1.angle - c2.angle) < 1e-4;
    double area1 = std::abs(c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0]);
    double area2 = std::abs(c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0]);
    bool areamatch = std::abs(area1 - area2) < 1e-6;
    bool equals = (nummatch && areamatch);
    bool match = false;
    if (equals)
    {
        Atoms a1 = c1.stack;
        Atoms a2 = c2.stack;
        bool match = a1.xtalcomp_compare(a2);
    }
    return (equals && match);
}

inline bool operator==(Interface &c1, Interface &c2)
{
    bool match = c1.stack.xtalcomp_compare(c2.stack);
    return match;
}

inline bool operator>(const Interface &c1, const Interface &c2)
{
    return (c1.stack.numAtom > c2.stack.numAtom);
}

inline bool operator<(const Interface &c1, const Interface &c2)
{
    return (c1.stack.numAtom < c2.stack.numAtom);
}
