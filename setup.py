#!/usr/bin/env python
import re
from setuptools import find_packages
from pathlib import Path
import os
import sys

VERSIONFILE = "hetbuilder/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


try:
    import skbuild
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

print("Scikit-build version:", skbuild.__version__)


# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="hetbuilder",
    version=version,
    author="Roman Kempt",
    author_email="roman.kempt@tu-dresden.de",
    description="A tool to build heterostructure interfaces based on coincidence lattice theory.",
    long_description=long_description,
    license="MIT",
    url="https://github.com/romankempt/hetbuilder.git",
    download_url="https://github.com/romankempt/hetbuilder.git",
    packages=find_packages(),
    package_data={"": ["*.xyz", "CMakeLists.txt"]},
    # package_dir={"": ""},
    cmake_install_dir="hetbuilder",
    include_package_data=True,
    # scripts=["./bin/hetbuilder"],
    install_requires=[
        "numpy",
        "scipy",
        "spglib",
        "matplotlib",
        "ase",
        "networkx",
        "pretty_errors",
        "rich",
        "typer",
        # "pybind11",
    ],
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    cmake_args=[
        "-DCMAKE_BUILD_TYPE={}".format("RELEASE"),  # not used on MSVC, but no harm,
    ],
    zip_safe=False,
)
