#  ___________________________________________________________________________
#
#  PSELIB: A library of nonlinear optimization test problems from process
#  systems engineering
#  Copyright (c) 2024. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

from setuptools import setup

setup(
    name="pselib",
    version="0.1.0",
    description=(
        "A library of nonlinear optimization test problems from process systems"
        " engineering applications."
    ),
    author="Robert Parker",
    author_email="rbparker@lanl.gov",
    license="BSD-3",
    packages=["pselib"],
)
