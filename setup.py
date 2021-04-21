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

import re

# VERSIONFILE = "interfacebuilder/__init__.py"
# verstrline = open(VERSIONFILE, "rt").read()
# VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
# mo = re.search(VSRE, verstrline, re.M)
# if mo:
#    version = mo.group(1)
# else:
#    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

version = "0.1"


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        build_directory = os.path.abspath(self.build_temp)

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + build_directory,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        # Assuming Makefiles
        build_args += ["--", "-j2"]

        self.build_args = build_args

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # CMakeLists.txt is in the same directory as this setup.py file
        cmake_list_dir = os.path.abspath(os.path.dirname(__file__))
        print("-" * 10, "Running CMake prepare", "-" * 40)
        subprocess.check_call(
            ["cmake", cmake_list_dir] + cmake_args, cwd=self.build_temp, env=env
        )

        print("-" * 10, "Building extensions", "-" * 40)
        cmake_cmd = ["cmake", "--build", "."] + self.build_args
        subprocess.check_call(cmake_cmd, cwd=self.build_temp)

        # Move from build temp to final position

        print("-" * 10, "Copying extensions", "-" * 40)
        for ext in self.extensions:
            self.move_output(ext)

    def move_output(self, ext):
        build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = build_temp / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file(source_path, dest_path)


ext_modules = [Extension("pybackend", sources=["backend/pybindings.cpp"])]

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
    packages=find_packages(),
    # scripts=["bin/build_interface"],
    install_requires=[
        "spglib",
        "numpy",
        "scipy",
        "matplotlib",
        "ase",
        "networkx",
        "cmake",
    ],
    classifiers=[
        "License :: OSI Approved :: MIT license",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    ext_modules=ext_modules,
    cmdclass=dict(build_ext=CMakeBuild),
    # zip_safe=False,
)
