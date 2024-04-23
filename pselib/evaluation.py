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

import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from idaes.core.util.model_statistics import large_residuals_set


def assert_primal_feasibility(model, atol=0.0):
    large_residuals = large_residuals_set(model, tol=atol)
    violated_bounds = ComponentSet()
    for var in model.component_data_objects(pyo.Var, active=True):
        if var.ub is None:
            ub_viol = 0.0
        else:
            ub_viol = var.value - var.ub
        if var.lb is None:
            lb_viol = 0.0
        else:
            lb_viol = var.lb - var.value
        if max(ub_viol, lb_viol) > atol:
            violated_bounds.add(var)
    assert not large_residuals and not violated_bounds


def evaluate_objective_function(model):
    objectives = list(model.component_data_objects(pyo.Objective, active=True))
    if len(objectives) > 1:
        obj_names = [obj.name for obj in objectives]
        raise RuntimeError(
            f"Only a single active objective function is supported."
            f" Found objectives: {obj_names}"
        )
    elif not objectives:
        raise RuntimeError("No active objective function found")
    return pyo.value(objectives[0].expr)
