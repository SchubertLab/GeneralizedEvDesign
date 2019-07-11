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

import itertools as itr
import copy
import math
import time
import numpy as np

import multiprocessing as mp
import logging

import pyomo.environ as aml
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.core.expr.numeric_expr import SumExpression
import pyomo.kernel as pmo

from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Utility import generate_overlap_graph
from Fred2.Utility import solve_TSP_LKH as _solve_tsp


class MosaicVaccineLazyILP(object):

    def __init__(self, predicted_affinities, threshold=None, max_vaccine_aminoacids=100,
                 max_vaccine_epitopes=999999999999, verbosity=0,
                 subtour_elimination='dfj'):
        if not isinstance(predicted_affinities, EpitopePredictionResult):
            raise ValueError('first input parameter is not of type EpitopePredictionResult')
        elif subtour_elimination not in ['dfj', 'mtz']:
            raise ValueError('subtour elimination method must be either "dfj" or "mtz"')
        
        self._subtour_elimination_method = subtour_elimination
        self.__raw_affinities = predicted_affinities
        self.__alleles = copy.deepcopy(self.__raw_affinities.columns.values.tolist())
        self.__allele_probs = self.__fill_allele_probs(self.__alleles)
        self.__solver = SolverFactory('gurobi_persistent')
        self.__verbosity = verbosity
        self.__changed = True
        self.__max_vaccine_aminoacids = max_vaccine_aminoacids
        self.__result = None
        self.__thresh = threshold or {}
        self.__pep_to_index = {}
        self.__max_vaccine_epitopes = max_vaccine_epitopes
        self.__subtour_constraints = 0

        self.__process_parameters()

        self.__arc_cost = generate_overlap_graph(self.__peptides[1:])
        self._initial_peptides, self._initial_arcs = self._find_initial_solution()

        self.build_model()

    def __process_parameters(self):
        self.__peptides = ['start']
        self.__immunogenicities = [0]
        self.__variations = set()
        self.__conservations = {0: 0}
        self.__epitope_variations = {}
        self.__allele_bindings = {}

        method = self.__raw_affinities.index.values[0][1]  # TODO find better way of filtering by method
        res_df = self.__raw_affinities.xs(self.__raw_affinities.index.values[0][1], level='Method')
        threshold_mask = res_df.apply(lambda x: any(  # FIXME shouldn't we log-transform first if method is one of .... ?
            x[allele] > self.__thresh.get(allele.name, -float('inf'))
            for allele in res_df.columns
        ), axis=1)

        if threshold_mask.sum() == 0:
            raise ValueError('binding affinity threshold too high, no peptides selected')
        else:
            print(threshold_mask.sum(), 'peptides above threshold, breakdown:')
            for col in res_df.columns:
                thr = self.__thresh[col.name]
                print('   ', col, np.sum(res_df[col] > thr))

        res_df = res_df[threshold_mask]
        for i, tup in enumerate(res_df.itertuples()):
            peptide = tup[0]
            i += 1  # skip first dummy peptide
            self.__peptides.append(peptide)
            self.__pep_to_index[peptide] = i
            immunogen = 0
            for allele, bind_affinity in itr.izip(self.__alleles, tup[1:]):
                if method in ['smm', 'smmpmbec', 'arb', 'comblibsidney']:
                    # log-transform binding strenghts and thresholds
                    bind_affinity = min(1., max(0.0, 1.0 - math.log(bind_affinity, 50000)))
                    if allele.name in self.__thresh:
                        bind_thr = min(1., max(0.0, 1.0 - math.log(self.__thresh.get(allele.name), 50000)))
                    else:
                        bind_thr = float('-inf')
                else:
                    bind_thr = self.__thresh.get(allele.name, float('-inf'))

                if bind_affinity >= bind_thr:
                    self.__allele_bindings.setdefault(allele.name, set()).add(i)
                    immunogen += self.__allele_probs[allele.name] * bind_affinity
            self.__immunogenicities.append(immunogen)

            proteins = set(peptide.get_all_proteins())
            self.__conservations[i] = len(proteins)
            for prot in proteins:
                self.__variations.add(prot.gene_id)
                self.__epitope_variations.setdefault(prot.gene_id, set()).add(i)

        #calculate conservation
        total = len(self.__variations)
        self.__conservations = {
            e: (v / total) if total > 0 else 1
            for e, v in self.__conservations.iteritems()
        }

    def __fill_allele_probs(self, alleles):
        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        if len(no_prob) > 0:
            #group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            for g, v in no_prob_grouped.iteritems():
                total_loc_a = len(v)
                if g in prob_grouped:
                    remaining_mass = 1.0 - sum(a.prob for a in prob_grouped[g])
                    for a in v:
                        a.prob = remaining_mass / total_loc_a
                else:
                    for a in v:
                        a.prob = 1.0 / total_loc_a

        return {a.name: a.prob for a in alleles}
    
    def _find_initial_solution(self):
        ''' find an initial solution to initialize the variables
            greedily select the peptides with highest immunogenicity
            until the maximum length of the vaccine is reached

            returns the set of epitopes and the set of arcs
        '''
        sorted_by_immunogenicity = sorted(
            self.__peptides, key=lambda pep: self.__immunogenicities[self.__pep_to_index.get(pep, 0)], reverse=True
        )

        solution_peptides, solution_arcs = [0], []
        aminoacids_length = 0
        for i, pep in enumerate(sorted_by_immunogenicity):
            if i >= self.__max_vaccine_epitopes:
                break
            
            idx = self.__pep_to_index[pep]
            aminoacids_length += self.__arc_cost[solution_peptides[-1], idx]
            if aminoacids_length > self.__max_vaccine_aminoacids:
                break
            
            solution_arcs.append((solution_peptides[-1], idx))
            solution_peptides.append(idx)
    
        solution_arcs.append((solution_peptides[-1], 0))       
        if self.__verbosity:
            print('The model will be initialized with the following feasible solution:')
            for u, v in solution_arcs:
                print('    %d (%s) -> %d (%s) (cost: %d, immunogenicity: %.2f)' % (
                    u, self.__peptides[u], v, self.__peptides[v], 
                    self.__arc_cost[u, v], self.__immunogenicities[v]
                ))
    
        return set(solution_peptides), set(solution_arcs)
    
    def build_model(self):
        if self.__verbosity:
            print('Building model...')
        self.model = aml.ConcreteModel()
        self.__build_model_sets()
        self.__build_model_params()
        self.__build_model_variables()

        self.__build_model_constraint_connectivity()
        self.__build_model_constraint_length()
        #self.__build_model_optionals()
        
        self.__build_model_objective()

        if self.__verbosity:
            print('Model built!')

    def __build_model_params(self):
        self.model.i         = aml.Param(self.model.Nodes, initialize=lambda model, i: self.__immunogenicities[i])
        self.model.d         = aml.Param(self.model.Arcs, initialize=lambda model, i, j: self.__arc_cost[i][j])
        self.model.k         = aml.Param(initialize=self.__max_vaccine_epitopes, within=aml.PositiveIntegers, mutable=True)
        self.model.TMAX      = aml.Param(initialize=self.__max_vaccine_aminoacids, within=aml.PositiveIntegers, mutable=True)
    
    def __build_model_sets(self):
        self.model.Nodes     = aml.RangeSet(0, len(self.__peptides) - 1)
        self.model.Arcs      = self.model.Nodes * self.model.Nodes # FIXME dont know if this includes self references as well....

    def __build_model_variables(self):
        self.model.x         = aml.Var(self.model.Arcs, domain=aml.Binary, bounds=(0, 1),
                                       initialize=lambda model, u, v: (u, v) in self._initial_arcs)
        self.model.y         = aml.Var(self.model.Nodes, domain=aml.Binary,
                                       initialize=lambda model, u: u in self._initial_peptides)

        if self._subtour_elimination_method == 'mtz':
            self.model.u     = aml.Var(self.model.Nodes - set([0]), bounds=(1.0, len(self.__peptides) - 1))

    def __build_model_objective(self):
        self.model.Obj       = aml.Objective(rule=lambda model: sum(model.y[i] * model.i[i] for i in model.Nodes),
                                             sense=aml.maximize)

    def __build_model_constraint_connectivity(self):
        ''' number of selected incoming connections must equal number of selected outgoing connections
            any selected node must have exactly one outgoing connection (so, implicitly, exactly one incoming as well)
        '''
        def InsideOutRule(model, node):
            outgoing = sum(model.x[node, j] for j in model.Nodes if j != node)
            incoming = sum(model.x[i, node] for i in model.Nodes if i != node)
            return outgoing == incoming

        self.model.IncomingOutgoing = aml.Constraint(self.model.Nodes, rule=InsideOutRule)
        self.model.Selection = aml.Constraint(self.model.Nodes, rule=lambda model, node: (
            sum(model.x[node, i] + model.x[i, node] for i in model.Nodes if i != node) == 2 * model.y[node]
        ))
    
    def __build_model_constraint_length(self):  # aka the Knapsack constraint
        ''' do not exceed desired number of aminoacids and peptides in the vaccine 
        '''
        self.model.maxNofAminoacids = aml.Constraint(rule=lambda model: (
            None, sum(model.d[i, j] * model.x[i, j] for i, j in model.Arcs), model.TMAX))
        self.model.maxNofEpitopes = aml.Constraint(rule=lambda model: (
            None, sum(model.y[n] for n in model.Nodes), model.k))

    def __build_model_optionals(self):
        self.model.c         = aml.Param(self.model.Nodes, initialize=lambda model, e: self.__conservations[e], mutable=True)
        self.model.t_allele  = aml.Param(initialize=0, within=aml.NonNegativeIntegers, mutable=True)
        self.model.t_var     = aml.Param(initialize=0, within=aml.NonNegativeIntegers, mutable=True)
        self.model.t_c       = aml.Param(initialize=0.0, within=aml.NonNegativeReals, mutable=True)

        self.model.Q         = aml.Set(initialize=self.__variations)
        self.model.A         = aml.Set(initialize=self.__allele_bindings.keys())
        self.model.E_var     = aml.Set(self.model.Q, initialize=lambda mode, v: self.__epitope_variations[v])
        self.model.A_I       = aml.Set(self.model.A, initialize=lambda model, allele: self.__allele_bindings[allele])

        self.model.z         = aml.Var(self.model.A, domain=aml.Binary, initialize=0)
        self.model.w         = aml.Var(self.model.Q, domain=aml.Binary, initialize=0)

        self.model.IsAlleleCovConst = aml.Constraint(self.model.A, rule=lambda model, allele:
            sum(model.y[e] for e in model.A_I[allele]) >= model.z[allele])
        self.model.MinAlleleCovConst = aml.Constraint(rule=lambda model:
            sum(model.z[allele] for allele in model.A) >= model.t_allele)
        self.model.IsAntigenCovConst = aml.Constraint(self.model.Q, rule=lambda model, q:
            sum(model.y[e] for e in model.E_var[q]) >= model.w[q])
        self.model.MinAntigenCovConst = aml.Constraint(rule=lambda model:
            sum(model.w[q] for q in model.Q) >= model.t_var)
        self.model.EpitopeConsConst = aml.Constraint(self.model.Nodes, rule=lambda model, e:
            (1 - model.c[e]) * model.y[e] <= 1 - model.t_c)

        # TODO need a way to re-activate them
        self.model.IsAlleleCovConst.deactivate()
        self.model.MinAlleleCovConst.deactivate()
        self.model.IsAntigenCovConst.deactivate()
        self.model.MinAntigenCovConst.deactivate()
        self.model.EpitopeConsConst.deactivate()

    def solve(self, options=None):
        ''' solves the model optimally
        '''
        if not self.__changed:
            return self.__result

        options = options if options is not None else {
            'MIPFocus': 3,  # focus on improving the bound
            'Method': 3,    # use the concurrent solver for root relaxation
        }
        self.__solver.set_instance(self.model)
        while True:
            res = self.__solver.solve(options=options, tee=bool(self.__verbosity), save_results=False,
                                      report_timing=True)
            if res.solver.termination_condition != TerminationCondition.optimal:
                raise RuntimeError('Could not solve problem - %s . Please check your settings' % res.Solution.status)

            self.__solver.load_vars()
            arcs = self.__extract_solution_from_model()
            if self.__verbosity:
                print('Solution contains the following arcs:', arcs)
            
            tours = self._extract_tours_from_arcs(arcs)
            if self.__verbosity:
                print('Solution contains the following tour(s):')
                for tt in tours:
                    print('   ', tt)

            if len(tours) > 1:
                if self._subtour_elimination_method == 'mtz':
                    self._eliminate_subtours_mtz(arcs, tours)
                else:
                    self._eliminate_subtours_dfj(arcs, tours)

                if self.__verbosity:
                    print('========================')
                    print('Invalid solution returned!')
                    print('%s subtour elimination constraints updated (%d inserted so far)' % (
                        self._subtour_elimination_method.upper(), self.__subtour_constraints
                    ))
                
                self._reinitialize_model(tours)
            else:
                break

        if self.__verbosity > 0:
            res.write(num=1)

        self.__result = tours[0]
        vaccine_peptides = [self.__peptides[j] for i, j in self.__result[:-1]]
        return self.__result, vaccine_peptides

    def _reinitialize_model(self, tours):    
        ''' re-initializes the model with a new feasible solution
            first, the current (invalid) solution is erased
            then, a new solution is generated and inserted
        '''

        for tour in tours:
            for u, v in tour:
                self.model.x[u, v] = self.model.y[u] = 0
        
        new_peptides, new_arcs = self._find_initial_solution()
        for u, v in new_arcs:
            self.model.x[u, v] = self.model.y[u] = 1
        
    def _eliminate_subtours_dfj(self, arcs, tours):
        ''' adds DFJ subtour elimination constraints
        '''
        for tour in tours:
            tour_nodes = set(i for i, _ in tour)
            self.__subtour_constraints += 1
            name = 'Subtour_%d' % self.__subtour_constraints
            constraint = pmo.constraint(
                body=sum(self.model.x[i, j] for i in tour_nodes for j in tour_nodes if i != j),
                ub=len(tour_nodes) - 1
            )
            setattr(self.model, name, constraint)
            self.__solver.add_constraint(getattr(self.model, name))

    def _eliminate_subtours_mtz(self, arcs, tours):
        ''' adds MTZ subtour elimination constraints
        '''
        def add_mtz_constraint(i, j):
            self.__subtour_constraints += 1
            name = 'Subtour_%d' % self.__subtour_constraints
            constraint = pmo.constraint(
                body=self.model.u[i] - self.model.u[j] + (len(self.__peptides) - 1) * self.model.x[i, j],
                ub=len(self.__peptides) - 2
            )
            setattr(self.model, name, constraint)
            self.__solver.add_constraint(getattr(self.model, name))

        for i, j in arcs.iteritems():
            if i != 0 and j != 0:
                add_mtz_constraint(j, i)
                add_mtz_constraint(i, j)

    def __extract_solution_from_model(self):
        ''' returns a dictionary i -> j containing the tour found by the model
        '''
        tour_arcs = {i: j for i, j in self.model.Arcs
                     if 0.98 <= self.model.x[i, j].value <= 1.5}
        return tour_arcs
    
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