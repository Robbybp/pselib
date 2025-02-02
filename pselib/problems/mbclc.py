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
from pselib.evaluation import assert_primal_feasibility
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
        # Fix initial conditions
        t0 = model.fs.time.first()
        for x in model.fs.moving_bed.length_domain:
            if x != model.fs.moving_bed.length_domain.first():
                model.fs.moving_bed.gas_phase.material_holdup[t0, x, ...].fix()
                model.fs.moving_bed.gas_phase.energy_holdup[t0, x, ...].fix()
            if x != model.fs.moving_bed.length_domain.last():
                model.fs.moving_bed.solid_phase.material_holdup[t0, x, ...].fix()
                model.fs.moving_bed.solid_phase.energy_holdup[t0, x, ...].fix()

        # Initialize derivatives to zero
        model.fs.moving_bed.gas_phase.material_accumulation[...].set_value(0.0)
        model.fs.moving_bed.gas_phase.energy_accumulation[...].set_value(0.0)
        model.fs.moving_bed.solid_phase.material_accumulation[...].set_value(0.0)
        model.fs.moving_bed.solid_phase.energy_accumulation[...].set_value(0.0)

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


def _get_disturbance_data(horizon=60.0):
    disturbance_dict = {"CO2": 0.5, "H2O": 0.0, "CH4": 0.5}
    disturbance = mpc.IntervalData(
        {
            "fs.moving_bed.gas_phase.properties[*,0.0].mole_frac_comp[%s]" % j: [val]
            for j, val in disturbance_dict.items()
        },
        [(0.0, horizon)],
    )
    return disturbance


def get_state_variable_names(space):
    setpoint_states = []
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].flow_mol" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].temperature" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].pressure" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].mole_frac_comp[CH4]" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].mole_frac_comp[H2O]" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.gas_phase.properties[*,%s].mole_frac_comp[CO2]" % x
        for x in space if x != space.first()
    )
    setpoint_states.extend(
        "fs.moving_bed.solid_phase.properties[*,%s].flow_mass" % x
        for x in space if x != space.last()
    )
    setpoint_states.extend(
        "fs.moving_bed.solid_phase.properties[*,%s].temperature" % x
        for x in space if x != space.last()
    )
    setpoint_states.extend(
        "fs.moving_bed.solid_phase.properties[*,%s].mass_frac_comp[Fe2O3]" % x
        for x in space if x != space.last()
    )
    setpoint_states.extend(
        "fs.moving_bed.solid_phase.properties[*,%s].mass_frac_comp[Fe3O4]" % x
        for x in space if x != space.last()
    )
    setpoint_states.extend(
        "fs.moving_bed.solid_phase.properties[*,%s].mass_frac_comp[Al2O3]" % x
        for x in space if x != space.last()
    )
    return setpoint_states


def _get_setpoint_model(
    nfe=10,
    target_conversion=0.95,
):
    m = _make_model(nxfe=nfe)

    horizon = 60.0
    disturbance = _get_disturbance_data(horizon=horizon)
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

    sp_target = {"fs.moving_bed.solid_phase.reactions[*,0.0].OC_conv": target_conversion}
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
    return m


def _get_setpoint_data(nfe=10, target_conversion=0.95):
    m = _get_setpoint_model(nfe=nfe, target_conversion=target_conversion)
    solver = pyo.SolverFactory("ipopt")
    solver.options["nlp_scaling_method"] = "user-scaling"
    solver.solve(m, tee=True)
    setpoint_interface = mpc.DynamicModelInterface(m, m.fs.time)
    setpoint_data = setpoint_interface.get_data_at_time()
    return setpoint_data


class SteadyMbclcMethane(PseTestProblemBase):

    @classmethod
    @property
    def uid(cls):
        return "MBCLC-METHANE-STEADY"

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
        m = _get_setpoint_model(nfe=nfe)
        return m


class DynamicMbclcMethane(PseTestProblemBase):

    @classmethod
    @property
    def uid(cls):
        return "MBCLC-METHANE-DYNAMIC"

    def __init__(self):
        # TODO: allow specification of parameters as configuration arg?
        self._parameters = [
            # We will set solid inlet temperature for the first minute of operation
            "fs.moving_bed.solid_inlet.temperature",
            # We multiply initial material holdup by a factor for the top 20% of the
            # reactor.
            "initial_material_holdup_factor",
        ]
        _parameter_range_list = [
            (800.0*pyo.units.K, 2000.0*pyo.units.K),
            # These are just factors that multiply the initial steady state values
            # This is to avoid having to specify a different holdup for each component.
            (0.8, 2.0),
        ]
        self._parameter_ranges = dict(zip(self._parameters, _parameter_range_list))

    @property
    def parameters(self):
        return self._parameters

    @property
    def parameter_ranges(self):
        return self._parameter_ranges

    def set_parameter_values(self, model, parameters):
        for key, val in parameters.items():
            # Accept string or ComponentUID as a key
            if key == "fs.moving_bed.solid_inlet.temperature":
                for t in model.fs.time:
                    # TODO: set this time threshold as an instance argument?
                    # TODO: Use find_nearest_index instead of float tolerances
                    if t <= 60.0+1e-6:
                        model.fs.moving_bed.solid_inlet.temperature[t].fix(val)
            elif key == "initial_material_holdup_factor":
                t0 = model.fs.time.first()
                for x in model.fs.moving_bed.length_domain:
                    # TODO: Set this threshold as well
                    if x >= 0.8-1e-6:
                        for var in model.fs.moving_bed.solid_phase.material_holdup[t0, x, ...]:
                            var.set_value(var.value * val)

    def create_instance(self, nxfe=20, ntfe=20, tfe_width=15.0):
        m = _make_model(dynamic=True, nxfe=nxfe, ntfe=ntfe, tfe_width=tfe_width)
        horizon = tfe_width * ntfe
        # Add objective function
        dyn = mpc.DynamicModelInterface(m, m.fs.time)
        setpoint = _get_setpoint_data(nfe=nxfe, target_conversion=0.97)

        # Solve a model where the setpoint is the initial state
        #m_sp = _make_model(dynamic=False, nxfe=nxfe)
        #_initialize(m_sp)
        #solver = pyo.SolverFactory("ipopt")
        #solver.solve(m_sp, tee=True)
        #sp_dyn = mpc.DynamicModelInterface(m_sp, m_sp.fs.time)
        #setpoint = sp_dyn.get_data_at_time()

        # Extract variable that I actually want to use in the setpoint
        setpoint_varnames = get_state_variable_names(m.fs.moving_bed.length_domain)
        setpoint_vars = [m.find_component(name) for name in setpoint_varnames]
        weight_data = {}
        for name in setpoint_varnames:
            if "temperature" in name:
                weight_data[name] = (1 / 100.0)**2
            elif "flow_mass" in name:
                weight_data[name] = (1 / 1000.0)**2
            elif "flow_mol" in name:
                weight_data[name] = (1 / 100.0)**2
            else:
                weight_data[name] = 1.0
        m.setpoint_set, m.setpoint_expr = dyn.get_penalty_from_target(
            setpoint, variables=setpoint_vars, weight_data=weight_data
        )
        m.setpoint_objective = pyo.Objective(
            expr=sum(m.setpoint_expr[i, t] for i in m.setpoint_set for t in m.fs.time)
        )

        # Unfix control variables
        for t in m.fs.time:
            if t != m.fs.time.first():
                m.fs.moving_bed.gas_inlet.flow_mol[t].unfix()
                m.fs.moving_bed.solid_inlet.flow_mass[t].unfix()
        # Add piecewise-constant constraints
        sample_points = [t for i, t in enumerate(m.fs.time) if (i % 2 == 0)]
        dof_vars = [
            m.fs.moving_bed.gas_inlet.flow_mol,
            m.fs.moving_bed.solid_inlet.flow_mass,
        ]
        m.pwc_set, m.pwc_con = dyn.get_piecewise_constant_constraints(dof_vars, sample_points)

        # Set bounds
        # Some "operational" bounds on control inputs
        m.fs.moving_bed.gas_inlet.flow_mol.setlb(80.0)
        m.fs.moving_bed.gas_inlet.flow_mol.setub(200.0)
        m.fs.moving_bed.solid_inlet.flow_mass.setlb(400.0)
        #m.fs.moving_bed.solid_inlet.flow_mass.setub(1200.0)
        m.fs.moving_bed.solid_inlet.flow_mass.setub(1000.0)
        # Some safety/physical limit bounds
        m.fs.moving_bed.gas_phase.properties[:, :].temperature.setlb(200.0)
        m.fs.moving_bed.gas_phase.properties[:, :].temperature.setub(2200.0)
        m.fs.moving_bed.solid_phase.properties[:, :].temperature.setlb(200.0)
        m.fs.moving_bed.solid_phase.properties[:, :].temperature.setub(2200.0)
        m.fs.moving_bed.gas_phase.properties[:, :].flow_mol.setlb(100.0)
        m.fs.moving_bed.gas_phase.properties[:, :].flow_mol.setub(1000.0)
        m.fs.moving_bed.solid_phase.properties[:, :].flow_mass.setlb(100.0)
        m.fs.moving_bed.solid_phase.properties[:, :].flow_mass.setub(2000.0)

        # Initialize
        # Initialization should come after the problem is all set up, because it should
        # be equivalent to call initialize after this method instead.
        m_init = _make_model(dynamic=False, nxfe=nxfe)
        _initialize(m_init)
        solver = pyo.SolverFactory("ipopt")
        solver.solve(m_init, tee=True)
        # TODO: Raise error if time is not a set in the model?
        init_dyn = mpc.DynamicModelInterface(m_init, m_init.fs.time)
        init_data = init_dyn.get_data_at_time()
        # Initialize dynamic model to an initial steady state solution
        dyn.load_data(init_data)
        # At this point, non-time-indexed variables in the dynamic model have not been
        # initialized.
        scalar_data = init_dyn.get_scalar_variable_data()
        dyn.load_data(scalar_data)

        # Apply disturbance
        # Disturbance needs to be applied after initialization. Otherwise initialization
        # with the initial steady state will override the disturbance.
        disturbance = _get_disturbance_data(horizon=horizon)
        dyn.load_data(disturbance)

        ## Setting an additional disturbance for the first minute of operation makes
        ## the problem more difficult to solve. We stop converging at about 1800 K.
        #for t in m.fs.time:
        #    if t <= 60.0+1e-6: # Account for float roundoff in time points
        #        m.fs.moving_bed.solid_inlet.temperature[t].fix(1800.0*pyo.units.K)

        ## Apply a perturbation to the initial condition of the dynamic model.
        ## This doesn't even converge with a 1% perturbation for all x. Maybe I need
        ## to only perturb a limited x?
        #t0 = m.fs.time.first()
        #for x in m.fs.moving_bed.length_domain:
        #    # I don't think we need to skip the inlet here.
        #    if x >= 0.8-1e-6:
        #        for var in m.fs.moving_bed.solid_phase.material_holdup[t0, x, ...]:
        #            var.set_value(var.value * 0.9)

        assert_primal_feasibility(m_init, atol=1e-5)
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        # Set some dummy scaling factors to avoid warnings
        # ... these aren't getting rid of the warnings
        var = m.fs.moving_bed.solid_phase.heat
        m.scaling_factor[var] = 1.0
        var = m.fs.moving_bed.solid_phase.rate_reaction_extent
        m.scaling_factor[var] = 1.0
        var = m.fs.moving_bed.solid_phase.area
        m.scaling_factor[var] = 1.0
        var = m.fs.moving_bed.gas_phase.heat
        m.scaling_factor[var] = 1.0
        var = m.fs.moving_bed.gas_phase.area
        m.scaling_factor[var] = 1.0

        # Scale model
        iscale.calculate_scaling_factors(m)

        # Set any model parameters we know about at this point
        return m


if __name__ == "__main__":
    # TODO: CLI argument for steady vs dynamic
    steady = False
    if steady:
        problem = SteadyMbclcMethane()
        m = problem.create_instance()
        solver = pyo.SolverFactory("ipopt")
        solver.solve(m, tee=True)
    else:
        problem = DynamicMbclcMethane()
        m = problem.create_instance()
        solver = pyo.SolverFactory("ipopt")
        solver.solve(m, tee=True)
        dyn = mpc.DynamicModelInterface(m, m.fs.time)
