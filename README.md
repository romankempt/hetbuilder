# Hetbuild

Builds 2D heterostructure interfaces via coincidence lattice theory.

## Installation

Requires a C++17 compiler and [cmake](https://cmake.org/) at least with version 3.18.4. A suitable compiler can be installed with Anaconda:

```bash
conda install -c conda-forge cxx-compiler git pip cmake
```

Then you can directly install from git:

```bash
pip install 2D-interface-builder.zip
```

## Documentation

Documentation is available at [ReadTheDocs](https://2d-interface-builder.readthedocs.io).

## Testing

Tests can be run in the project directory with

```bash
pytest -v tests
```

More tests will follow.

## Requirements

- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/) ase 3.19 or higher
- [Space Group Libary](https://atztogo.github.io/spglib/python-spglib.html) spglib
- [Scientific Python](https://www.scipy.org/) scipy (including numpy, matplotlib and pandas)
- [networkx](https://networkx.github.io/documentation/stable/install.html)
