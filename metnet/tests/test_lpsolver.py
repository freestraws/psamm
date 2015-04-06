#!/usr/bin/env python

import unittest

from metnet.lpsolver import lp

try:
    from metnet.lpsolver import cplex
except ImportError:
    cplex = None


class TestCplexProblem(unittest.TestCase):
    def setUp(self):
        self.solver = cplex.Solver()

    @unittest.skipIf(cplex is None, 'solver not available')
    def test_objective_reset_on_set_linear_objective(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)

        # Solve first time, maximize x
        prob.set_linear_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('x'), 10)

        # Solve second time, maximize y
        # If the objective is not properly cleared,
        # the second solve will still maximize x.
        prob.set_linear_objective(prob.var('y'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('y'), 10)

    @unittest.skipIf(cplex is None, 'solver not available')
    def test_result_to_bool_conversion_on_optimal(self):
        '''Run a feasible LP problem and check that the result evaluates to True'''
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)

        prob.set_linear_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertTrue(result)

    @unittest.skipIf(cplex is None, 'solver not available')
    def test_result_to_bool_conversion_on_infeasible(self):
        '''Run an infeasible LP problem and check that the result evaluates to False'''
        prob = self.solver.create_problem()
        prob.define('x', 'y', 'z', lower=0, upper=10)
        prob.add_linear_constraints(2*prob.var('x') == -prob.var('y'),
                                    prob.var('x') + prob.var('z') >= 6,
                                    prob.var('z') <= 3)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)
        prob.set_linear_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertFalse(result)


if __name__ == '__main__':
    unittest.main()
