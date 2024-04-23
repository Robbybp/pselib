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

from pselib.problems.mbclc import SteadyMbclcMethane


ALL_TEST_PROBLEM_CLASSES = [
    SteadyMbclcMethane,
]

TEST_PROBLEM_LOOKUP = {
    problem.uid: problem for problem in ALL_TEST_PROBLEM_CLASSES
}


def get_problem(uid):
    return TEST_PROBLEM_LOOKUP[uid]()


class TestSet:
    """Set of test problems
    """

    def __init__(self, test_problem_uids=None):
        # Cache instantiated TestProblems on the TestSet instance.
        # Note that there may be some reason to defer instantiation, for example
        # if we want to instantiate multiple times with different parameters
        # at some point.
        if test_problem_uids is None:
            self._problems = list(
                problemclass() for problemclass in TEST_PROBLEM_LOOKUP.values()
            )
        else:
            self._problems = [
                TEST_PROBLEM_LOOKUP[uid]() for uid in test_problem_uids
            ]

    @property
    def problems(self):
        return self._problems
