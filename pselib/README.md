## IDAES Test Set

A test set of process systems engineering-related models and optimization problems.
Intended to benchmark nonlinear equation and optimization solvers.

## Core API

The most common way to interact with IDAES Test Set is via the `TestSet` class.
A `TestSet` defines a set of test problems. By default, a `TestSet` object
contains every problem in the IDAES Test Set.

Test problems are represented as `TestProblem` objects. `TestProblem`s provide
a consistent interface for interacting with the underlying Pyomo models that
make up the test problems. Every `TestProblem` has a unique string identifier
that can be accessed by the `uid` field, and can construct a Pyomo model
representing the test problem with the `create_instance` method.

The following example iterates over all test problems, creates a Pyomo model
of the problem, and inspects a few properties.
```python
import pyomo.environ as pyo
import idaes_testset as its
ts = its.TestSet
for problem in ts.problems:
    print(type(problem))
    print(problem.uid)
    model = problem.create_instance()
    print(type(problem))

    nvar = len(list(model.component_objects(pyo.Var)))
    print(f"N. variables: {nvar}")
```

## Adding a model to IDAES Test Set
