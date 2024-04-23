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

import pytest
import pyomo.environ as pyo
import pselib
from pselib.evaluation import assert_primal_feasibility


PRIMAL_FEASIBILITY_TOL = 1e-5


class TestSolveProblems:

    def test_solve_problems(self):
        ts = pselib.TestSet()
        solver = pyo.SolverFactory("ipopt")
        for problem in ts.problems:
            m = problem.create_instance()
            results = solver.solve(m, tee=True)
            pyo.assert_optimal_termination(results)
            assert_primal_feasibility(m, atol=PRIMAL_FEASIBILITY_TOL)


if __name__ == "__main__":
    pytest.main([__file__])
