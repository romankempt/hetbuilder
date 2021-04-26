# Hetbuilder - builds heterostructure interfaces

[![DOI](https://zenodo.org/badge/358881237.svg)](https://zenodo.org/badge/latestdoi/358881237)

Builds 2D heterostructure interfaces via coincidence lattice theory.

## Installation

Requires a C++17 compiler, [cmake](https://cmake.org/) at least with version 3.18.4, as well as [spglib](https://atztogo.github.io/spglib/python-spglib.html) and [pybind11](https://github.com/pybind/pybind11). 
All can be installed with Anaconda:

```bash
conda install -c conda-forge cxx-compiler git pip cmake spglib pybind11
```

Then you can directly install from git:

```bash
pip install git+https://github.com/romankempt/hetbuilder.git
```

## Documentation


## Testing

Tests can be run in the project directory with

```bash
pytest -v tests
```


## Requirements

- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
- [Space Group Libary](https://atztogo.github.io/spglib/python-spglib.html)
- [SciPy](https://www.scipy.org/)
- [pybind11](https://github.com/pybind/pybind11)
