import ase .io
from ase.neighborlist import NeighborList
from ase.geometry.analysis import Analysis
import numpy as np
from itertools import combinations_with_replacement


def get_bond_lengths(atoms: "ase.atoms.AToms"):
    symbs = set(atoms.get_chemical_symbols())
    pairs = list(combinations_with_replacement(symbs, 2))
    ana = Analysis(atoms)
    unique_bonds = {}
    for p in pairs:
        bonds = ana.get_bonds(*p)
        if bonds == [[]] :
            continue
        bond_values = ana.get_values(bonds)
        avg, std = np.average(bond_values), np.std(bond_values)
        unique_bonds[p] = (avg, std)
    return unique_bonds


at1 = ase.io.read("C122Mo37S74_angle1.00_stress0.00.xyz")
at2 = ase.io.read("MoS2_2H_1L.xyz")
at3 = ase.io.read("graphene.xyz")

b1 = get_bond_lengths(at1)
print(b1)
b2 = get_bond_lengths(at2)
b3 = get_bond_lengths(at3)

print(b2, b3)

for key, item in b1.items():
    l, std = item
    for k2, i2 in b2.items():
        if (k2 == key) or (k2[::-1] == key):
            l2, std2 = i2
            print(std, std2)
            print((l2 - l)/l * 100)
            print( (std2 - std)/std * 100)
