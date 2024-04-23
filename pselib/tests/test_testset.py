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

import pselib
from pselib.testproblem import PseTestProblemBase
from pselib.testset import get_problem
import pytest

class TestTestSet:

    def test_problems(self):
        ts = pselib.TestSet()
        assert len(ts.problems) > 0
        for problem in ts.problems:
            print(type(problem))
            assert isinstance(problem, PseTestProblemBase)


if __name__ == "__main__":
    pytest.main([__file__])
