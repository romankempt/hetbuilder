#!/usr/bin/env python
import os
import re
import sys
import sysconfig
import platform
import subprocess
from pathlib import Path

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand

import re

VERSIONFILE = "hetbuilder/__init__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "TESTING" if self.debug else "RELEASE"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            "-DEXAMPLE_VERSION_INFO={}".format(self.distribution.get_version()),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
        ]
        build_args = []

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


setup(
    name="hetbuilder",
    version=version,
    author="Roman Kempt",
    author_email="roman.kempt@tu-dresden.de",
    description="A tool to build heterostructure interfaces based on coincidence lattice theory.",
    # long_description=open("README.md").read(),
    license="MIT",
    # url="https://github.com/AK-Heine/2D-Interface-Builder",
    # download_url="https://github.com/AK-Heine/2D-Interface-Builder",
    packages=find_packages(
        where="hetbuilder",
        exclude=[
            "*.tests",
            "*.tests.*",
            "tests.*",
            "tests",
            "WIP",
            "pictures",
            "examples",
            "docs",
            "*__pycache__*",
            "*vscode*",
        ],
    ),
    package_data={"": ["*.so"]},
    # scripts=["bin/build_interface"],
    install_requires=[
        "spglib",
        "numpy",
        "scipy",
        "matplotlib",
        "ase",
        "networkx",
        "cmake",
        "pybind11",
        "pretty_errors",
        "rich",
    ],
    classifiers=[
        "License :: OSI Approved :: MIT license",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    ext_modules=[CMakeExtension("hetbuilder_backend")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
