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

"""Base class for test problems

Note that this file is named testproblem.py rather than test_problem.py because
otherwise it gets picked up by pytest, which tries to instantiate and run the
TestProblemBase abstract class.

"""
# ^ I think the solution here is to avoid naming classes TestXYZ. Even if their
# home file is not test_xyz, they will get treated as a test class when they
# are imported into a test module.
import abc


class PseTestProblemBase(abc.ABC):

    @abc.abstractmethod
    def uid(self):
        pass

    @abc.abstractmethod
    def create_instance(self):
        """Construct a Pyomo model of this test problem
        """
        pass
