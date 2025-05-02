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

from pselib.testproblem import PseTestProblemBase
import pyomo.environ as pyo
import pyomo.contrib.mpc as mpc

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.models_extra.gas_solid_contactors.unit_models.moving_bed import (
    MBR as MovingBed,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction import (
    GasPhaseParameterBlock,
    SolidPhaseParameterBlock,
    HeteroReactionParameterBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import scaling as iscale


def _make_model(
    dynamic=False,
    ntfe=10,
    ntcp=1,
    tfe_width=15.0,
    nxfe=10,
):
    """Make a square reduction reactor model.

    Does not do any initialization, scaling, or set up anything necessary
    for the optimization problem.

    """
    if dynamic:
        # This is the *model* horizon
        horizon = ntfe * tfe_width

    # Create spatial grid
    # Why is this necessary? Shouldn't it be the default?
    xfe_list = [1.0*i/nxfe for i in range(nxfe+1)]

    # Create time grid
    if dynamic:
        time_list = [0.0, horizon]
    else:
        time_list = [0.0]

    # Create top-level and flowsheet models
    if dynamic:
        name = "Dynamic MB reduction CLC flowsheet model"
    else:
        name = "Steady-state MB reduction CLC flowsheet model"
    model = pyo.ConcreteModel(name=name)
    model.fs = FlowsheetBlock(
        dynamic=dynamic,
        time_set=time_list,
        time_units=pyo.units.s,
    )

    # Set up thermo props and reaction props
    model.fs.gas_properties = GasPhaseParameterBlock()
    model.fs.solid_properties = SolidPhaseParameterBlock()
    model.fs.hetero_reactions = HeteroReactionParameterBlock(
        solid_property_package=model.fs.solid_properties,
        gas_property_package=model.fs.gas_properties,
    )

    mb_config = {
        "finite_elements": nxfe,
        "has_holdup": True,
        "length_domain_set": xfe_list,
        "transformation_method": "dae.finite_difference",
        "gas_transformation_scheme": "BACKWARD",
        "solid_transformation_scheme": "FORWARD",
        "pressure_drop_type": "ergun_correlation",
        "gas_phase_config": {
            "property_package": model.fs.gas_properties,
        },
        "solid_phase_config": {
            "property_package": model.fs.solid_properties,
            "reaction_package": model.fs.hetero_reactions,
        },
    }

    # Create MovingBed unit model
    model.fs.moving_bed = MovingBed(**mb_config)

    if dynamic:
        # Apply time-discretization
        discretizer = pyo.TransformationFactory('dae.collocation')
        discretizer.apply_to(
            model, wrt=model.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU'
        )

    # Note that we need to fix inputs/disturbances *before* fixing initial,
    # conditions as we rely on having square subsystems to identify diff eqns.
    _fix_dof(model)

    if dynamic:
        raise RuntimeError()
    #    # Fix initial conditions
    #    dae = DAEInterface(model, model.fs.time)
    #    # Note that the "differential subsystem" is "not defined" at t0
    #    # as the discretization equations don't exist.
    #    t0 = model.fs.time.first()
    #    t1 = list(model.fs.time)[1]
    #    diff_deriv_disc_list = dae.get_valid_diff_deriv_disc_at_time(t1)
    #    # Note that these variables are already indexed.
    #    diff_vars = [var for var, _, _ in diff_deriv_disc_list]
    #    for var in diff_vars:
    #        var[t0].fix()

    #    # Initialize derivatives to zero
    #    model.fs.moving_bed.gas_phase.material_accumulation[...].set_value(0.0)
    #    model.fs.moving_bed.gas_phase.energy_accumulation[...].set_value(0.0)
    #    model.fs.moving_bed.solid_phase.material_accumulation[...].set_value(0.0)
    #    model.fs.moving_bed.solid_phase.energy_accumulation[...].set_value(0.0)

    return model


def _fix_dof(model):
    mb = model.fs.moving_bed
    # Fix geometry variables
    mb.bed_diameter.fix(6.5*pyo.units.m)
    mb.bed_height.fix(5*pyo.units.m)

    # Fix inlet variables for all t
    time = model.fs.time

    # The following code should be valid for steady state model as well
    # Fix inlets to nominal values
    for t in time:
        # Inlets fixed to ss values, no disturbances
        mb.gas_inlet.flow_mol[t].fix(128.20513*pyo.units.mol/pyo.units.s)
        #mb.gas_inlet.pressure[t].fix(2.00)
        mb.gas_inlet.pressure[t].fix(2.00*pyo.units.bar)

        mb.solid_inlet.flow_mass[t].fix(591.4*pyo.units.kg/pyo.units.s)
        mb.gas_inlet.temperature[t].fix(298.15*pyo.units.K)
        mb.gas_inlet.mole_frac_comp[t, "CO2"].fix(0.02499)
        mb.gas_inlet.mole_frac_comp[t, "H2O"].fix(0.00001)
        mb.gas_inlet.mole_frac_comp[t, "CH4"].fix(0.97500)

        mb.solid_inlet.temperature[t].fix(1183.15*pyo.units.K)
        mb.solid_inlet.particle_porosity[t].fix(0.27)
        mb.solid_inlet.mass_frac_comp[t, "Fe2O3"].fix(0.45)
        mb.solid_inlet.mass_frac_comp[t, "Fe3O4"].fix(1e-9)
        mb.solid_inlet.mass_frac_comp[t, "Al2O3"].fix(0.55)


def _initialize(model):
    mb = model.fs.moving_bed
    gas_phase_state_args = {
        "flow_mol": mb.gas_inlet.flow_mol[0].value,
        "temperature": mb.solid_inlet.temperature[0].value,
        "pressure": mb.gas_inlet.pressure[0].value,
        "mole_frac": {
            "CH4": mb.gas_inlet.mole_frac_comp[0, "CH4"].value,
            "CO2": mb.gas_inlet.mole_frac_comp[0, "CO2"].value,
            "H2O": mb.gas_inlet.mole_frac_comp[0, "H2O"].value,
        },
    }
    solid_phase_state_args = {
        "flow_mass": mb.solid_inlet.flow_mass[0].value,
        "particle_porosity": mb.solid_inlet.particle_porosity[0].value,
        "temperature": mb.solid_inlet.temperature[0].value,
        "mass_frac": {
            "Fe2O3": mb.solid_inlet.mass_frac_comp[0, "Fe2O3"].value,
            "Fe3O4": mb.solid_inlet.mass_frac_comp[0, "Fe3O4"].value,
            "Al2O3": mb.solid_inlet.mass_frac_comp[0, "Al2O3"].value,
        },
    }

    iscale.calculate_scaling_factors(model)

    model.fs.moving_bed.initialize(
        optarg={"tol": 1e-5},
        gas_phase_state_args=gas_phase_state_args,
        solid_phase_state_args=solid_phase_state_args,
    )


class SteadyMbclcMethane(PseTestProblemBase):

    uid = "MBCLC-METHANE-STEADY"

    def __init__(self):
        # TODO: allow specification of parameters as configuration arg?
        self._parameters = [
            pyo.ComponentUID("fs.moving_bed.solid_phase.properties[0,1].temperature"),
            pyo.ComponentUID("fs.moving_bed.solid_phase.properties[0,1].flow_mass"),
        ]
        _parameter_range_list = [
            (1000.0*pyo.units.K, 1400.0*pyo.units.K),
            (500.0*pyo.units.kg/pyo.units.s, 700.0*pyo.units.kg/pyo.units.s),
            #(1183.14, 1183.16),
            #(591.3, 592.5),
        ]
        self._parameter_ranges = dict(zip(self._parameters, _parameter_range_list))

    @property
    def parameters(self):
        return self._parameters

    @property
    def parameter_ranges(self):
        return self._parameter_ranges

    # Should this class be configured during instantiation?
    # Or at the call to construct?
    def create_instance(self, nfe=10):
        m = _make_model(nxfe=nfe)

        horizon = 60.0
        disturbance_dict = {"CO2": 0.5, "H2O": 0.0, "CH4": 0.5}
        disturbance = mpc.IntervalData(
            {
                "fs.moving_bed.gas_phase.properties[*,0.0].mole_frac_comp[%s]" % j: [val]
                for j, val in disturbance_dict.items()
            },
            [(0.0, horizon)],
        )

        setpoint_interface = mpc.DynamicModelInterface(m, m.fs.time)
        setpoint_inputs = mpc.ScalarData({
            "fs.moving_bed.gas_phase.properties[*,0.0].flow_mol": 128.2,
            "fs.moving_bed.solid_phase.properties[*,1.0].flow_mass": 591.4,
        })
        setpoint_interface.load_data(setpoint_inputs)
        iscale.calculate_scaling_factors(m)
        setpoint_interface.load_data(disturbance.get_data_at_time(horizon))

        _initialize(m)
        solver = pyo.SolverFactory("ipopt")
        solver.solve(m, tee=True)

        sp_target = {"fs.moving_bed.solid_phase.reactions[*,0.0].OC_conv": 0.95}
        m.fs.moving_bed.gas_inlet.flow_mol[:].unfix()
        (
            m.penalty_set, m.conv_penalty
        ) = setpoint_interface.get_penalty_from_target(sp_target)
        m.obj = pyo.Objective(expr=sum(m.conv_penalty.values()))

        for x in m.fs.moving_bed.length_domain:
            var = m.fs.moving_bed.gas_phase.properties[0, x].dens_mol
            m.fs.moving_bed.gas_phase.properties[0, x].scaling_factor[var] = 1.0

            var = m.fs.moving_bed.solid_phase.reactions[0, x].OC_conv
            m.fs.moving_bed.gas_phase.properties[0, x].scaling_factor[var] = 1.0

            var = m.fs.moving_bed.solid_phase.reactions[0, x].OC_conv_temp
            m.fs.moving_bed.gas_phase.properties[0, x].scaling_factor[var] = 1.0

            var = m.fs.moving_bed.Nu_particle[0, x]
            m.fs.moving_bed.scaling_factor[var] = 1.0
            var = m.fs.moving_bed.Pr_particle[0, x]
            m.fs.moving_bed.scaling_factor[var] = 1.0
            var = m.fs.moving_bed.Re_particle[0, x]
            m.fs.moving_bed.scaling_factor[var] = 1.0

        # The problem with scaling the model here is that all subsequent
        # parameters must be set with the scaling transformation in mind, which
        # is not practical.
        #pyo.TransformationFactory("core.scale_model").apply_to(m, rename=False)

        return m


if __name__ == "__main__":
    problem = SteadyMbclcMethane()
    m = problem.create_instance()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=True)
