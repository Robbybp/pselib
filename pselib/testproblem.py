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


"""
TODO:
__init__ method: What should be required to define an "instance of a problem"?
    - parameters
    - parameter ranges
    How are parameters (the keys themselves, not the values) validated?
    We may want to vary things that are not explicit mutable parameters in the
    model. Maybe they are arguments to some model construction function...

How to handle parameters that are not "model parameters" (i.e., they are not
accessible via find_component)?
    - E.g. a parameter of the setpoint in a control problem
    - We definitely want to support these, right?
    - Some API consistency should be supported:
      - Print name of parameter
      - Get range of values
      - Create a model with particular values of these parameters. "Instance parameters"
        (like nfe) must be provided at model construction time. "Model parameters"
        can be provided here. An advantage of setting model parameters here is that
        we can scale them if necessary.
        By convention, setting model parameters here should be equivalent to setting
        them directly on the model immediately following construction (as long as
        scaling is done properly).

set_model_parameter function that automatically applies scaling to parameters if
necessary?

initialize_model method. I'd like this library to be useful for testing custom
initialization strategies, so I need to separate model initialization.

scale_model method. I'd like this library to be useful for testing custom scaling
methods, so I need a separate scaling method.

initialize=True and scale=True should probably be arguments to create_instance,
so the uninterested user doesn't have to think about these.

create_instance: This name is a bit confusing, since "instance" could mean an
instance of this class. create_model is probably better.

"""
class PseTestProblemBase(abc.ABC):

    @abc.abstractmethod
    def uid(self):
        pass

    @abc.abstractmethod
    def create_instance(self):
        """Construct a Pyomo model of this test problem
        """
        pass
