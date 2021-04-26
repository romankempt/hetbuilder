# Hetbuilder - builds heterostructure interfaces

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
