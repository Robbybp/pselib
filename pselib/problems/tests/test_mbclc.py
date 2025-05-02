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
from pselib.problems.mbclc import (
    DynamicMbclcMethane,
    _get_setpoint_data,
)
import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface

from var_elim.elimination_callbacks import matching_elim_callback, d2_elim_callback
from idaes.core.util.model_statistics import large_residuals_set


class TestDynamicMbclcMethane:

    def test_create_instance(self):
        problem = DynamicMbclcMethane()
        m = problem.create_instance()
        igraph = IncidenceGraphInterface(m, include_inequality=False)
        ncon = len(igraph.constraints)
        nvar = len(igraph.variables)
        dof = nvar - ncon
        vdmp, cdmp = igraph.dulmage_mendelsohn()
        nv_uc = len(vdmp.unmatched)
        nc_oc = len(cdmp.unmatched)
        assert ncon > 0
        assert nvar > 0
        assert dof > 0
        assert nv_uc == dof
        assert nc_oc == 0
        objs = list(m.component_data_objects(pyo.Objective, active=True))
        assert len(objs) == 1

    def test_solve(self):
        problem = DynamicMbclcMethane()
        m = problem.create_instance(nxfe=20)
        parameters = {"fs.moving_bed.solid_inlet.temperature": 1900.0}
        problem.set_parameter_values(m, parameters)
        solver = pyo.SolverFactory("cyipopt")
        #solver.config.options["nlp_scaling_method"] = "user-scaling"
        #pyo.TransformationFactory("core.scale_model").apply_to(m)
        violated_cons = large_residuals_set(m, tol=1e-5)
        print("Violated constraints")
        print("--------------------")
        violated_cons = sorted(
            violated_cons,
            key=lambda con: abs(pyo.value(con.body - con.upper)),
            reverse=True,
        )
        for con in violated_cons:
            val = abs(pyo.value(con.body - con.upper))
            print(f"{val}  {con.name}")
        print("--------------------")
        #matching_elim_callback(m)
        #d2_elim_callback(m)
        res = solver.solve(m, tee=True)
        violated_cons = large_residuals_set(m, tol=1e-5)
        pyo.assert_optimal_termination(res)
        assert not violated_cons


if __name__ == "__main__":
    #pytest.main([__file__])
    TestDynamicMbclcMethane().test_solve()
