# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
'''
   .. module:: Mosaic
   :synopsis:  This class implements the epitope selection functionality
    to construct so called mosaic vaccines

    This module builds upon Coopr's Pyomo, an embedded algebraic modeling
    languages [2].

    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice


.. moduleauthor:: dorigatti

'''

from __future__ import division, print_function

import copy
import itertools as itr
import logging
import math
import multiprocessing as mp
import StringIO
import sys
import time

import numpy as np
import pyomo.environ as aml
import pyomo.kernel as pmo
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Utility import generate_overlap_graph
from Fred2.Utility import solve_TSP_LKH as _solve_tsp
from pyomo.core.expr.numeric_expr import SumExpression
from pyomo.opt import SolverFactory, TerminationCondition


class TeamOrienteeringIlp:

    # TODO eventually should allow minimum coverage of certain subtypes
    def __init__(self, num_teams, vertex_reward, edge_cost, max_edge_cost, max_vertices,
                 model_type='dfj', solver='gurobi_persistent'):

        self.logger = logging.getLogger(self.__class__.__name__)

        if num_teams < 1:
            raise ValueError('at least one tema needed')

        self._model = None
        self._result = None
        self._solver = SolverFactory(solver)
        self._vertex_reward = vertex_reward
        self._num_teams = num_teams
        self._edge_cost = edge_cost
        self._max_edge_cost = max_edge_cost
        self._max_vertices = max_vertices

    def build_model(self):
        self.logger.info('Building model...')

        self._model = aml.ConcreteModel()

        # model parameters
        self._model.TeamCount = aml.Param(initialize=self._num_teams)
        self._model.MaxEdgeCost = aml.Param(initialize=self._max_edge_cost)
        self._model.MaxVertexCount = aml.Param(initialize=self._max_vertices)

        # graph objects
        self._model.Teams = aml.RangeSet(0, self._num_teams - 1)
        self._model.Nodes = aml.RangeSet(0, len(self._vertex_reward) - 1)
        self._model.Edges = aml.Set(initialize=[
            (i, j)
            for i in xrange(len(self._vertex_reward))
            for j in xrange(len(self._vertex_reward))
            if i != j
        ])

        # reward of nodes and cost of edges
        self._model.r = aml.Param(self._model.Nodes, initialize=lambda model, n: self._vertex_reward[n])
        self._model.d = aml.Param(self._model.Edges, initialize=lambda model, u, v: self._edge_cost[u][v])

        # indicator variables for nodes and arcs
        self._model.x = aml.Var(self._model.Nodes * self._model.Nodes * self._model.Teams, domain=aml.Binary, initialize=0)
        self._model.y = aml.Var(self._model.Nodes * self._model.Teams, domain=aml.Binary, initialize=0)

        # objective of the model: maximize reward collected from visited nodes
        self._model.Obj = aml.Objective(
            rule=lambda model: sum(model.y[n, t] * model.r[n] for n in model.Nodes for t in model.Teams),
            sense=aml.maximize
        )

        # every team must leave and come back
        self._model.TeamsLeave = aml.Constraint(
            rule=lambda model: sum(model.x[0, n, t] for n in model.Nodes for t in model.Teams if n > 0) == model.TeamCount
        )
        self._model.TeamsReturn = aml.Constraint(
            rule=lambda model: sum(model.x[n, 0, t] for n in model.Nodes for t in model.Teams if n > 0) == model.TeamCount
        )

        # every vertex must be visited at most once
        self._model.VertexVisit = aml.Constraint(
            (n for n in self._model.Nodes if n != 0),
            rule=lambda model, node: sum(model.y[node, t] for t in model.Teams) <= 1
        )
        
        # incoming connnections = outgoing connections = node selected
        # i.e. no sources or sinks (implies path is connected)
        #      enforces consistency between x and y (i.e. node touched by arcs if and only if it is selected)
        #      and at most one path passes from the node
        self._model.Incoming = aml.Constraint(
            ((n, t) for n in self._model.Nodes for t in self._model.Teams),
            rule=lambda model, node, team: sum(
                model.x[node, v, team] for v in model.Nodes if v != node
            ) == model.y[node, team]
        )
        self._model.Outgoing = aml.Constraint(
            ((n, t) for n in self._model.Nodes for t in self._model.Teams),
            rule=lambda model, node, team: sum(
                model.x[v, node, team] for v in model.Nodes if v != node
            ) == model.y[node, team]
        )

        # each path must not exceed the specified edge cost
        if self._max_edge_cost > 0:
            self.logger.info('Maximum edge cost for each tour is %f', self._max_edge_cost)
            self._model.MaxPathEdgeCost = aml.Constraint(
                self._model.Teams,
                rule=lambda model, team: sum(model.x[u, v, team] * model.d[u, v] for (u, v) in model.Edges) <= model.MaxEdgeCost
            )
        else:
            self.logger.info('No maximum edge cost enforced.')

        # each path must not pass through more vertices than specified
        if self._max_vertices > 0:
            self.logger.info('Maximum vertex count for each tour is %d', self._max_vertices)
            self._model.MaxPathVertices = aml.Constraint(
                self._model.Teams,
                rule=lambda model, team: sum(model.y[n, team] for n in model.Nodes if n != 0) <= model.MaxVertexCount
            )
        else:
            self.logger.info('No maximum vertex count enforced.')

        self.logger.info('Model build!')
        return self

    def solve(self, options=None):
        # if logging is configured, gurobipy will print messages there *and* on stdout
        # so we silence its logger and redirect all stdout to our own logger
        logging.getLogger('gurobipy.gurobipy').disabled = True

        class LoggingStdOut:
            def __init__(self):
                self.logger = logging.getLogger('stdout')

            def write(self, message):
                self.logger.debug(message.rstrip())

            def flush(self, *args, **kwargs):
                pass

        sys.stdout = LoggingStdOut()

        try:
            return self._solve(options)
        except Exception:
            # restore stdout so that handlers can print normally
            # https://docs.python.org/3/library/sys.html#sys.__stdout__
            sys.stdout = sys.__stdout__
            raise
        finally:
            sys.stdout = sys.__stdout__

    def _solve(self, options=None):
        ''' solves the model optimally
        '''
        options = options if options is not None else {
            'MIPFocus': 3,  # focus on improving the bound
            'Method': 3,    # use the concurrent solver for root relaxation
        }

        self.logger.info('Solving started')
        if self._model is None:
            raise RuntimeError('must call build_model before solve')
        self._solver.set_instance(self._model)
        self._subtour_constraints = 0

        while True:
            res = self._solver.solve(options=options, tee=1, save_results=False, report_timing=True)
            if res.solver.termination_condition != TerminationCondition.optimal:
                raise RuntimeError('Could not solve problem - %s . Please check your settings' % res.Solution.status)

            self._solver.load_vars()
            team_tours_dict = self._extract_solution_from_model()
            self.logger.debug('Solution contains the following tours: %s', team_tours_dict)

            valid_solution = True
            all_tours_edges = []
            for i, tour in enumerate(team_tours_dict):
                this_tour_edges = self._extract_tours_from_arcs(tour)
                self.logger.debug('Solution for team %d contains the following tour(s):', i)
                for tt in this_tour_edges:
                    self.logger.debug('   %s', tt)

                if len(this_tour_edges) > 1 or not any(u == 0 for u, v in this_tour_edges[0]):
                    valid_solution = False
                    self._eliminate_subtours_dfj(this_tour_edges)
                    self.logger.debug('Subtour elimination constraints updated (%d inserted so far)' % (
                        self._subtour_constraints
                    ))
                else:
                    all_tours_edges.append(this_tour_edges[0])
            
            if valid_solution:
                break

        res.write(num=1)

        self.logger.info('Solved successfully')
        self._result = all_tours_edges
        return self._result

    def _eliminate_subtours_dfj(self, tours):
        ''' adds DFJ subtour elimination constraints
        '''
        for tour in tours:
            tour_nodes = set(i for i, _ in tour)
            for team in xrange(self._num_teams):
                self._subtour_constraints += 1
                name = 'Subtour_%d' % self._subtour_constraints
                constraint = pmo.constraint(
                    body=sum(
                        self._model.x[i, j, team]
                        for i in tour_nodes
                        for j in tour_nodes
                        if i != j
                    ),
                    ub=len(tour_nodes) - 1
                )
                setattr(self._model, name, constraint)
                self._solver.add_constraint(getattr(self._model, name))

    def _extract_solution_from_model(self):
        ''' returns a list of dictionaries i -> list of j containing the tours found by the model
        '''
        tours = []
        for t in range(self._num_teams):
            vertices = [
                n for n in self._model.Nodes
                if 0.98 <= self._model.y[n, t].value
            ]
            self.logger.debug('Team %d selected nodes %s', t, vertices)

            edges = {}
            for i, j in self._model.Edges:
                if 0.98 <= self._model.x[i, j, t].value <= 1.5:
                    edges[i] = j
            tours.append(edges)

            self.logger.debug('Team %d selected edges %s', t, edges)
            assert set(vertices) == set(edges)
        return tours

    @staticmethod
    def _extract_tours_from_arcs(arcs):
        ''' given a dictionary of arcs, returns a list of tours, where every tour is a list of arcs
        '''
        assert set(arcs.keys()) == set(arcs.values()), 'arcs do not form a set of tours: %r' % arcs
        not_assigned, assigned = set(arcs.keys()), set()

        tours = []
        while not_assigned:
            tour, cursor = [], not_assigned.pop()
            while not tour or cursor not in assigned:
                tour.append((cursor, arcs[cursor]))
                assigned.add(cursor)
                not_assigned.discard(cursor)
                cursor = arcs[cursor]
            tours.append(tour)
        return tours


def test():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.DEBUG)

    rewards = [0, 100, 100, 0, 100, 100, 0]
    costs = [[1] * 7 for _ in range(7)]

    top = TeamOrienteeringIlp(3, rewards, costs, 3, 3)
    top.build_model()
    print(top.solve())


if __name__ == '__main__':
    test()