# Hetbuilder - builds heterostructure interfaces

[![DOI](https://zenodo.org/badge/358881237.svg)](https://zenodo.org/badge/latestdoi/358881237)
[![Documentation Status](https://readthedocs.org/projects/hetbuilder/badge/?version=latest)](https://hetbuilder.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/hetbuilder.svg)](https://badge.fury.io/py/hetbuilder)

Builds 2D heterostructure interfaces via coincidence lattice theory.

## Installation

### Build-time dependencies

Requires a C++17 compiler and [cmake](https://cmake.org/).

It is also recommended to preinstall [spglib](https://atztogo.github.io/spglib/python-spglib.html) and [pybind11](https://github.com/pybind/pybind11).
Otherwise, these will be built during the installation from the submodules.

#### Installing with Anaconda

Create a clean conda environment:
```bash
conda env create -n hetbuilder python=3.9
```

Then install the build-time dependencies first:
```bash
conda install -c conda-forge cxx-compiler git pip cmake spglib pybind11
```

Then, you can install the project from pip:
```bash
pip install hetbuilder
```

If that does not work, try directly installing from git:
```bash
pip install git+https://github.com/romankempt/hetbuilder.git
```

#### Installing with pip

PyPI does not provide the library files of [spglib](https://atztogo.github.io/spglib/python-spglib.html). These will be built from the submodules at installation time, which might be time-consuming.

On Unix, you can install a `cxx-compiler` with:
```bash
sudo apt install build-essential
```


## First steps

The installation exposes a multi-level [typer](https://github.com/tiangolo/typer) CLI utility called `hetbuilder`:

```bash
hetbuilder --help
```

The `build` utility is fully implemented.
You can use any ASE-readable structure format to specify the lower and upper layer. They should be recognizable as two-dimensional, e.g., by having a zero vector in the *z*-direction.

```bash
hetbuilder build graphene.xyz MoS2.cif
```

This should open a [matplotlib](https://matplotlib.org/) interface looking like this:

![](pictures/interface.png)

## Documentation

Documentation is available at [Read the Docs](https://hetbuilder.readthedocs.io/en/latest/index.html).

## Testing

Tests can be run in the project directory with

```bash
pytest -v tests
```

## Citing

If you use this tool, please cite 10.5281/zenodo.4721346.

## Requirements

- [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
- [Space Group Libary](https://atztogo.github.io/spglib/python-spglib.html)
- [SciPy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/)
- [pybind11](https://github.com/pybind/pybind11)
- [typer](https://github.com/tiangolo/typer)
