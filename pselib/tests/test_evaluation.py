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
from pselib.evaluation import (
    assert_primal_feasibility,
    evaluate_objective_function,
)


def _make_simple_model():
    m = pyo.ConcreteModel()
    m.I = pyo.Set(initialize=range(1, 5))
    m.x = pyo.Var(m.I, initialize=1.0, bounds=(0, 5))
    m.eq = pyo.Constraint(pyo.Any)
    m.eq[1] = m.x[1] + m.x[2] == 3
    m.ineq = pyo.Constraint(pyo.Any)
    m.ineq[1] = m.x[3] + m.x[4] <= 7
    m.obj = pyo.Objective(expr=sum(m.x[i]**2 for i in m.I))
    return m


class TestPrimalFeasibility:

    def test_primal_feasibility(self):
        m = _make_simple_model()
        with pytest.raises(AssertionError):
            # eq[1] is not feasible
            assert_primal_feasibility(m)

        m.x[2] = 2
        # Now we are exactly feasible
        assert_primal_feasibility(m)

        m.x[2] = (2 + 1e-6)
        # Now we are only feasible within a modest tolerance
        with pytest.raises(AssertionError):
            assert_primal_feasibility(m)

        assert_primal_feasibility(m, atol=1e-5)

        m.x[3] = -1
        # Now we are infeasible due to a bound violation
        with pytest.raises(AssertionError):
            assert_primal_feasibility(m, atol=1e-5)

        m.x[3] = 6.5
        # Now we are infeasible due to an inequality violation
        with pytest.raises(AssertionError):
            assert_primal_feasibility(m, atol=1e-5)

        m.x[3] = 5
        m.x[4] = 2
        # An inequality and a bound are active, but we are still feasible
        assert_primal_feasibility(m, atol=1e-5)

    def test_evaluate_objective(self):
        m = _make_simple_model()
        obj_val = evaluate_objective_function(m)
        assert obj_val == 4
        m.x[1] = 0
        m.x[2] = 2
        m.x[3] = 1
        m.x[4] = 3
        obj_val = evaluate_objective_function(m)
        assert obj_val == 14
