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


import copy
import itertools as itr
import logging
import math
import multiprocessing as mp
import sys
import time
from collections import defaultdict
from io import StringIO

import numpy as np
import pyomo.environ as aml
import pyomo.kernel as pmo
from pyomo.core.expr.numeric_expr import SumExpression
from pyomo.opt import SolverFactory, TerminationCondition

from epytope.Core.Result import EpitopePredictionResult


class TeamOrienteeringIlp:

    def __init__(self, num_teams, vertex_reward, edge_cost, max_edge_cost, max_vertices,
                 type_coverage=None, min_type_coverage=None, min_avg_type_conservation=None,
                 lazy_subtour_elimination=False, solver='gurobi_persistent'):

        self.logger = logging.getLogger(self.__class__.__name__)

        if num_teams < 1:
            raise ValueError('at least one team needed')

        self._model = None
        self._result = None
        self._solver = None
        self._team_max_vertices_constraints = []
        self._team_max_edge_cost_constraints = []
        self._lazy_subtour_elimination = lazy_subtour_elimination
        self._solver_type = solver
        self._vertex_reward = vertex_reward
        self._num_teams = num_teams

        self._max_edge_cost = max_edge_cost
        self._max_vertices = max_vertices
        self._type_coverage = type_coverage
        self._min_type_coverage = min_type_coverage
        self._min_avg_type_conservation = min_avg_type_conservation

        if isinstance(edge_cost, dict):
            # in this case, the edge costs is a dictionary (u, v) -> cost
            self.logger.debug('Using sparse mode to build model')
            self._is_graph_sparse = True
        else:
            # in this case, we have a matrix (array of arrays)
            self.logger.debug('Using dense mode to build model')
            self._is_graph_sparse = False
        self._edge_cost = edge_cost

    def build_model(self):
        self.logger.info('Building model...')

        self._model = aml.ConcreteModel()

        # model parameters
        self.logger.debug('Adding graph objects...')
        self._model.TeamCount = aml.Param(initialize=self._num_teams, mutable=True)
        self._model.MaxEdgeCost = aml.Param(initialize=self._max_edge_cost, mutable=True)
        self._model.MaxVertexCount = aml.Param(initialize=self._max_vertices, mutable=True)

        self._model.Teams = aml.RangeSet(0, self._num_teams - 1)
        self._model.Nodes = aml.RangeSet(0, len(self._vertex_reward) - 1)
        self._model.r = aml.Param(self._model.Nodes, initialize=lambda model, n: self._vertex_reward[n])

        if self._is_graph_sparse:
            nodes_in, nodes_out = defaultdict(list), defaultdict(list)
            for u, v in self._edge_cost:
                nodes_in[v].append(u)
                nodes_out[u].append(v)

            self._model.Edges = aml.Set(initialize=self._edge_cost.keys())
            self._model.NodesIn = aml.Set(self._model.Nodes, initialize=lambda model, node: nodes_in[node])
            self._model.NodesOut = aml.Set(self._model.Nodes, initialize=lambda model, node: nodes_out[node])
            self._model.d = aml.Param(self._model.Edges, initialize=lambda model, u, v: self._edge_cost[(u, v)])
        else:
            self._model.Edges = aml.Set(initialize=self._model.Nodes * self._model.Nodes,
                                        filter=lambda model, u, v: u != v)
            self._model.d = aml.Param(self._model.Edges, initialize=lambda model, u, v: self._edge_cost[u][v])

        # indicator variables for nodes and arcs
        self._model.x = aml.Var(self._model.Edges * self._model.Teams, domain=aml.Binary, initialize=0)
        self._model.y = aml.Var(self._model.Nodes * self._model.Teams, domain=aml.Binary, initialize=0)

        # objective of the model: maximize reward collected from visited nodes
        self._model.Objective = aml.Objective(
            rule=lambda model: sum(model.y[n, t] * model.r[n] for n in model.Nodes for t in model.Teams),
            sense=aml.maximize
        )

        # every team must leave and come back
        self.logger.debug('Adding leave and return constraints...')
        self._model.TeamsLeave = aml.Constraint(
            rule=lambda model: sum(
                model.x[(u, v), t] for u, v in model.Edges for t in model.Teams if u == 0
            ) == model.TeamCount
        )
        self._model.TeamsReturn = aml.Constraint(
            rule=lambda model: sum(
                model.x[(u, v), t] for u, v in model.Edges for t in model.Teams if v == 0
            ) == model.TeamCount
        )

        # every vertex must be visited at most once
        self.logger.debug('Adding visit count constraint...')
        self._model.VertexVisit = aml.Constraint(
            (n for n in self._model.Nodes if n != 0),
            rule=lambda model, node: sum(model.y[node, t] for t in model.Teams) <= 1
        )

        # incoming connnections = outgoing connections = node selected
        # i.e. no sources or sinks (implies path is connected)
        #      enforces consistency between x and y (i.e. node touched by arcs if and only if it is selected)
        #      and at most one path passes from the node
        self.logger.debug('Adding consistency and connectedness constraints...')
        if self._is_graph_sparse:
            def in_rule(model, node, team): return sum(
                model.x[(v, node), team] for v in model.NodesIn[node]
            ) == model.y[node, team]

            def out_rule(model, node, team): return sum(
                model.x[(node, v), team] for v in model.NodesOut[node]
            ) == model.y[node, team]
        else:
            def in_rule(model, node, team): return sum(
                model.x[(node, v), team] for v in model.Nodes if v != node
            ) == model.y[node, team]

            def out_rule(model, node, team): return sum(
                model.x[(v, node), team] for v in model.Nodes if v != node
            ) == model.y[node, team]

        self._model.Incoming = aml.Constraint(
            ((n, t) for n in self._model.Nodes for t in self._model.Teams),
            rule=in_rule
        )
        self._model.Outgoing = aml.Constraint(
            ((n, t) for n in self._model.Nodes for t in self._model.Teams),
            rule=out_rule
        )

        # subtour elimination, if required
        if not self._lazy_subtour_elimination:
            self.logger.debug('Adding subtour elimination constraints...')
            self._model.u = aml.Var(self._model.Nodes * self._model.Teams, bounds=(1.0, len(self._model.Nodes) - 1))
            self._model.SubTour = aml.Constraint((
                (u, v, t)
                for u, v in self._model.Edges
                for t in self._model.Teams
                if u != 0 and v != 0
            ), rule=lambda model, u, v, t: (
                model.u[u, t] - model.u[v, t] + 1 <= (len(model.Nodes) - 1) * (1 - model.x[(u, v), t])
            ))

        if self._type_coverage and (self._min_type_coverage or self._min_avg_type_conservation):
            assert (not self._min_type_coverage
                    or not self._min_avg_type_conservation
                    or len(self._min_avg_type_conservation) == len(self._min_type_coverage))

            self.logger.debug('Adding coverage information...')
            # type_coverage is a binary tensor s.t. C_ijk = 1 iff vertex j covers option k of type i
            # min_type_coverage is a vector s.t. c_i is the minimum options of type i that the vaccine must cover
            self._model.Types = aml.RangeSet(0, len(self._type_coverage) - 1)
            self._model.Options = aml.RangeSet(0, len(self._type_coverage[0][0]) - 1)
            self._model.TypeCoverage = aml.Param(self._model.Types * self._model.Nodes * self._model.Options,
                                                 initialize=lambda model, t, n, o: self._type_coverage[t][n][o])

            # indicator variable 1 iff at least one team visits at least one vertex of option i of type j
            self._model.OptionCovered = aml.Var(
                self._model.Types * self._model.Options, domain=aml.Binary, initialize=0)
            self._model.OptionCoveredConstraint = aml.Constraint(
                self._model.Types * self._model.Options,
                rule=lambda model, typ, option: sum(
                    model.y[n, t] * model.TypeCoverage[typ, n, option]
                    for n in model.Nodes for t in model.Teams
                ) >= model.OptionCovered[typ, option]
            )

            self.logger.debug('Enforcing minimum coverage and/or conservation with %d types and %d options per type',
                              len(self._type_coverage), len(self._type_coverage[0][0]))

            if self._min_type_coverage:
                self.logger.debug('Adding minimum coverage for each type...')
                # sum of above indicator variables must be at least the minimum option coverage for each type
                self._model.MinOptionCoverage = aml.Param(
                    self._model.Types, initialize=lambda model, typ: self._min_type_coverage[typ]
                )
                self._model.MinOptionCoverageConstraint = aml.Constraint(
                    self._model.Types,
                    rule=lambda model, typ: (model.MinOptionCoverage[typ], sum(
                        model.OptionCovered[typ, o] for o in model.Options
                    ), None)
                )

            if self._min_avg_type_conservation:
                # every vertex must cover a minimum number of different options
                self.logger.debug('Adding average vertex conservation...')
                self._model.MinOptionConservation = aml.Param(
                    self._model.Types, initialize=lambda model, typ: self._min_avg_type_conservation[typ]
                )
                self._model.MinOptionConservationConstraint = aml.Constraint(
                    self._model.Types, rule=lambda model, typ: sum(
                        model.y[n, t] * (
                            sum(model.TypeCoverage[typ, n, o] for o in model.Options)
                            - model.MinOptionConservation[typ]
                        )
                        for t in model.Teams
                        for n in model.Nodes
                    ) >= 0)
        else:
            self.logger.info('No coverage enforced')

        self._solver = SolverFactory(self._solver_type)
        self._solver.set_instance(self._model)

        self.update_max_vertices(self._max_vertices)
        self.update_max_edge_cost(self._max_edge_cost)

        self.logger.info('Model build!')
        return self

    def update_max_vertices(self, max_vertices):
        def get_constraint(team):
            if max_vertices > 0:
                return pmo.constraint(expr=sum(
                    self._model.y[n, team] for n in self._model.Nodes if n != 0
                ) <= self._model.MaxVertexCount)
            else:
                return None

        self._max_vertices = max_vertices
        self._model.MaxVertexCount.set_value(max_vertices)
        self._team_max_vertices_constraints = self._update_constraint_for_all_teams(
            self._team_max_vertices_constraints, get_constraint, 'MaxVerticesForTeam%d'
        )

        if max_vertices < 0:
            self.logger.info('No maximum vertex count enforced.')
        else:
            self.logger.info('Maximum vertex count for each tour is %d', self._max_vertices)

    def update_max_edge_cost(self, max_edge_cost):
        def get_constraint(team):
            if max_edge_cost > 0:
                return pmo.constraint(expr=sum(
                    self._model.x[(u, v), team] * self._model.d[u, v] for u, v in self._model.Edges
                ) <= self._model.MaxEdgeCost)
            else:
                return None

        self._max_edge_cost = max_edge_cost
        self._model.MaxEdgeCost.set_value(max_edge_cost)
        self._team_max_edge_cost_constraints = self._update_constraint_for_all_teams(
            self._team_max_edge_cost_constraints, get_constraint, 'MaxEdgeCostForTeam%d'
        )

        if max_edge_cost > 0:
            self.logger.info('Maximum edge cost for each tour is %f', self._max_edge_cost)
        else:
            self.logger.info('No maximum edge cost enforced.')

    def _update_constraint_for_all_teams(self, current_constraints, constraint_fn, name_fmt):
        for name in current_constraints:
            constr = getattr(self._model, name)
            try:
                self._solver.remove_constraint(constr)
            except KeyError:
                # happens after model is built, but not solved. not sure why
                pass

            setattr(self._model, name, None)
            del constr

        new_constraints = []
        for team in range(self._num_teams):
            name = name_fmt % team
            constr = constraint_fn(team)
            if constr is None:
                continue
            setattr(self._model, name, constr)
            new_constraints.append(name)
            self._solver.add_constraint(constr)
        return new_constraints

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

        self.logger.info('Solving started')
        if self._model is None:
            raise RuntimeError('must call build_model before solve')
        self._subtour_constraints = 0

        while True:
            res = self._solver.solve(options=options or {}, tee=1, save_results=False, report_timing=True)
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
                    assert self._lazy_subtour_elimination, 'subtour elimination failed'
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
            for team in range(self._num_teams):
                self._subtour_constraints += 1
                name = 'Subtour_%d' % self._subtour_constraints
                constraint = pmo.constraint(
                    body=sum(
                        self._model.x[i, j, team]
                        for i in tour_nodes
                        for j in tour_nodes
                        if i != j and (i, j) in self._model.Edges
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

    def explore_edge_cost_vertex_reward_tradeoff(self, steps):
        # introduce variables for vertex reward and edge cost
        self._model.VertexReward = aml.Var()
        self._model.AssignVertexReward = pmo.constraint(expr=sum(
            self._model.y[n, t] * self._model.r[n] for n in self._model.Nodes for t in self._model.Teams
        ) == self._model.VertexReward)

        self._model.EdgeCost = aml.Var()
        self._model.AssignEdgeCost = aml.Constraint(expr=sum(
            self._model.x[u, v, t] * self._model.d[u, v] for (u, v) in self._model.Edges for t in self._model.Teams
        ) == self._model.EdgeCost)

        self._solver.add_var(self._model.VertexReward)
        self._solver.add_constraint(self._model.AssignVertexReward)
        self._solver.add_var(self._model.EdgeCost)
        self._solver.add_constraint(self._model.AssignEdgeCost)

        # step 1: obtain maximum vertex reward
        self.logger.info('Obtaining maximum vertex reward...')
        del self._model.Objective
        self._model.Objective = aml.Objective(expr=self._model.VertexReward, sense=aml.maximize)
        self._solver.set_objective(self._model.Objective)
        self.solve()
        max_reward = aml.value(self._model.VertexReward)
        self.logger.info('Maximum reward is %f with cost %f', max_reward, aml.value(self._model.EdgeCost))

        # step 2: obtain minimum edge cost
        self.logger.info('Obtaining minumum edge cost...')
        del self._model.Objective
        self._model.Objective = aml.Objective(expr=self._model.EdgeCost, sense=aml.minimize)
        self._solver.set_objective(self._model.Objective)
        self.solve()
        max_cost = aml.value(self._model.EdgeCost)
        self.logger.info('Minimum cost is %f with reward %f', max_cost, aml.value(self._model.VertexReward))

        # step 3: obtain minumum edge cost, conditioned on maximum vertex reward
        self.logger.info('Obtaining minimum cost conditioned on maximum reward...')
        self._model.ForcedReward = aml.Param(initialize=max_reward)
        self._model.MaxVertexReward = pmo.constraint(expr=self._model.VertexReward == self._model.ForcedReward)
        self._solver.add_constraint(self._model.MaxVertexReward)
        self.solve()
        max_cost_max_reward = aml.value(self._model.EdgeCost)
        self.logger.info('Cost is %f with reward %f', max_cost_max_reward, aml.value(self._model.VertexReward))

        # step 4: iterate between these two values
        self._solver.remove_constraint(self._model.MaxVertexReward)
        del self._model.MaxVertexReward
        del self._model.ForcedReward
        self._model.EdgeCostSlackValue = aml.Param(initialize=0.0, mutable=True)
        self._model.EdgeCostSlack = aml.Var(within=aml.NonNegativeReals)
        self._model.Epsilon = aml.Param(initialize=1e-4)
        del self._model.Objective
        self._model.Objective = aml.Objective(
            rule=lambda model: model.VertexReward + model.Epsilon * model.EdgeCostSlack,
            sense=aml.maximize
        )
        self._solver.add_var(self._model.EdgeCostSlack)
        self._solver.set_objective(self._model.Objective)

        for i in range(steps):
            value = max_cost_max_reward + i * (max_cost - max_cost_max_reward) / float(steps - 1)
            self.logger.info('======')
            self.logger.info('Iteration %d - Cost bound is %f', i + 1, value)
            self._model.EdgeCostSlackValue.set_value(value)
            self._model.EdgeCostConstr = pmo.constraint(
                expr=self._model.EdgeCost + self._model.EdgeCostSlack == self._model.EdgeCostSlackValue
            )
            self._solver.add_constraint(self._model.EdgeCostConstr)

            vaccine = self.solve()
            reward, cost = aml.value(self._model.VertexReward), aml.value(self._model.EdgeCost)
            self.logger.info('Obtained reward %f with cost %f', reward, cost)
            yield vaccine, reward, cost

            self._solver.remove_constraint(self._model.EdgeCostConstr)
            del self._model.EdgeCostConstr


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
