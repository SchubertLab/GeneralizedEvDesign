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
import StringIO
import sys

import multiprocessing as mp
import logging

import pyomo.environ as aml
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.core.expr.numeric_expr import SumExpression
import pyomo.kernel as pmo

from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Utility import generate_overlap_graph
from Fred2.Utility import solve_TSP_LKH as _solve_tsp


class MosaicVaccineLazyILP:

    def __init__(self, predicted_affinities, affinity_threshold=None, 
                 max_vaccine_aminoacids=100, max_vaccine_epitopes=None,
                 min_allele_coverage=None, min_antigen_coverage=None,
                 min_epitope_conservation=None, verbosity=0,
                 subtour_elimination='dfj'):

        self.logger = logging.getLogger(self.__class__.__name__)

        if not isinstance(predicted_affinities, EpitopePredictionResult):
            raise ValueError('first input parameter is not of type EpitopePredictionResult')
        elif subtour_elimination not in ['dfj', 'mtz']:
            raise ValueError('subtour elimination method must be either "dfj" or "mtz"')
        
        self._model = None
        self._raw_affinities = predicted_affinities
        self._alleles = copy.deepcopy(self._raw_affinities.columns.values.tolist())
        self._allele_probs = self._fill_allele_probs(self._alleles)
        self._solver = SolverFactory('gurobi_persistent')
        self._verbosity = verbosity
        self._result = None
        self._thresh = affinity_threshold or {}
        self._subtour_elimination_method = subtour_elimination

        self._process_parameters()

        # desired vaccine properties
        self._max_vaccine_aminoacids = max_vaccine_aminoacids
        self._max_vaccine_epitopes = max_vaccine_epitopes
        self._min_allele_coverage = min_allele_coverage
        self._min_antigen_coverage = min_antigen_coverage
        self._min_epitope_conservation = min_epitope_conservation

    def _process_parameters(self):
        self._peptides = ['start']
        self._pep_to_index = {}
        self._immunogenicities = [0]
        self._all_genes = set()
        self._conservations = {0: 0}
        self._gene_to_epitope = {}
        self._allele_bindings = {}
        self._peptide_bindings = {}

        method = self._raw_affinities.index.values[0][1]  # TODO find better way of filtering by method
        self.logger.info('Using binding affinities from method "%s"', method)
        res_df = self._raw_affinities.xs(self._raw_affinities.index.values[0][1], level='Method')
        threshold_mask = res_df.apply(lambda x: any(  # FIXME shouldn't we log-transform first if method is one of .... ?
            x[allele] > self._thresh.get(allele, -float('inf'))
            for allele in res_df.columns
        ), axis=1)

        if threshold_mask.sum() == 0:
            raise ValueError('Binding affinity threshold too high, no peptides selected')

        self.logger.info('%d peptides above threshold for at least one allele', threshold_mask.sum())
        self.logger.debug('Breakdown:')
        for col in res_df.columns:
            self.logger.debug('    %s %d', col, np.sum(res_df[col] > self._thresh[col]))

        res_df = res_df[threshold_mask]
        for i, tup in enumerate(res_df.itertuples()):
            peptide = tup[0]
            i += 1  # skip first dummy peptide
            self._peptides.append(peptide)
            self._pep_to_index[peptide] = i
            immunogen = 0
            for allele, bind_affinity in itr.izip(self._alleles, tup[1:]):
                if method in ['smm', 'smmpmbec', 'arb', 'comblibsidney']:
                    # log-transform binding strenghts and thresholds and clip in [0-1]
                    bind_affinity = min(1., max(0.0, 1.0 - math.log(bind_affinity, 50000)))
                    if allele in self._thresh:
                        bind_thr = min(1., max(0.0, 1.0 - math.log(self._thresh.get(allele), 50000)))
                    else:
                        bind_thr = float('-inf')
                else:
                    bind_thr = self._thresh.get(allele)

                if bind_affinity >= bind_thr:
                    self._allele_bindings.setdefault(allele.name, set()).add(i)
                    self._peptide_bindings.setdefault(i, set()).add(allele)
                    immunogen += self._allele_probs[allele.name] * bind_affinity

            self._immunogenicities.append(immunogen)

            proteins = set(peptide.get_all_proteins())
            self._conservations[i] = len(proteins)
            for prot in proteins:
                self._all_genes.add(prot.gene_id)
                self._gene_to_epitope.setdefault(prot.gene_id, set()).add(i)

        # calculate epitope conservation, i.e. number of proteins that contain it over number of genes
        # TODO this does not make any sense ?!
        total = len(self._all_genes)
        self._conservations = {
            e: (v / total) if total > 0 else 1
            for e, v in self._conservations.iteritems()
        }

    def _fill_allele_probs(self, alleles):
        # test if allele prob is set, if not set allele prob uniform
        # if only partly set infer missing values (assuming uniformity of missing value)
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
            self._peptides, key=lambda pep: self._immunogenicities[self._pep_to_index.get(pep, 0)], reverse=True
        )

        solution_peptides, solution_arcs = [0], []
        aminoacids_length = total_immunogenicity = 0
        for i, pep in enumerate(sorted_by_immunogenicity):
            if self._max_vaccine_epitopes > 0 and i >= self._max_vaccine_epitopes:
                break
            
            idx = self._pep_to_index[pep]
            arc_cost = self._arc_cost[solution_peptides[-1], idx]
            if self._max_vaccine_aminoacids > 0 and aminoacids_length + arc_cost > self._max_vaccine_aminoacids:
                break
            
            aminoacids_length += arc_cost
            total_immunogenicity += self._immunogenicities[idx]
            solution_arcs.append((solution_peptides[-1], idx))
            solution_peptides.append(idx)
    
        solution_arcs.append((solution_peptides[-1], 0))       
        self.logger.debug('The model will be initialized with the following feasible '
                          'solution with cost %d and immunogenicity %.2f:', aminoacids_length, total_immunogenicity)
        for u, v in solution_arcs:
            self.logger.debug('    %5d (%9s) -> %5d (%9s) (cost: %2d, immunogenicity: %4.2f)',
                              u, self._peptides[u], v, self._peptides[v],
                              self._arc_cost[u, v], self._immunogenicities[v])

        return set(solution_peptides), set(solution_arcs)

    def build_model(self):
        self._compute_arcs_cost()
        self._initial_peptides, self._initial_arcs = self._find_initial_solution()

        self.logger.info('Building model...')

        self.model = aml.ConcreteModel()
        self._build_base_tsp_model()
        self._build_model_vaccine_constraints()

        self.logger.info('Model built!')
        
    def _compute_arcs_cost(self):
        self.logger.debug('Computing arc costs...')
        # FIXME following solution based on suffix trees gives the wrong answer
        #self._arc_cost = generate_overlap_graph(self._peptides[1:])

        self._arc_cost = np.zeros((len(self._peptides), len(self._peptides)))

        for i, peptide_from in enumerate(self._peptides):
            for j, peptide_to in enumerate(self._peptides):
                cost = None
                if j == 0 or i == j:
                    cost = 0
                elif i == 0:
                    cost = len(peptide_to)
                else:  # compute longest suffix-prefix
                    k = 1
                    peptide_to, peptide_from = str(peptide_to), str(peptide_from)
                    while k < len(peptide_from) and k < len(peptide_to) and peptide_from[-k:] == peptide_to[:k]:
                        k += 1
                    cost = len(peptide_to) - k + 1
                
                self._arc_cost[i, j] = cost

        # following code checks that the overlaps are computed correctly
        for i, pfrom in enumerate(self._peptides):
            for j, pto in enumerate(self._peptides):
                overlap = int(9 - self._arc_cost[i, j])
                if overlap == 0 or i == 0 or j == 0:
                    continue

                pfrom, pto = str(pfrom), str(pto)
                assert pfrom[-overlap:] == pto[:overlap], (i, pfrom, j, pto, overlap)
        
    def _build_base_tsp_model(self):
        ''' builds a MIP model suitable for solving the orienteering problem
        '''
        # graph objects
        self.model.Nodes = aml.RangeSet(0, len(self._peptides) - 1)
        self.model.Arcs = aml.Set(initialize=[(i, j)
                                  for i in xrange(len(self._peptides))
                                  for j in xrange(len(self._peptides))
                                  if i != j])

        # reward of nodes and cost of arcs
        self.model.i = aml.Param(self.model.Nodes, initialize=lambda model, i: self._immunogenicities[i])
        self.model.d = aml.Param(self.model.Arcs, initialize=lambda model, i, j: self._arc_cost[i][j])

        # indicator variables for nodes and arcs
        self.model.x = aml.Var(self.model.Arcs, domain=aml.Binary, bounds=(0, 1),
                                       initialize=lambda model, u, v: (u, v) in self._initial_arcs)
        self.model.y = aml.Var(self.model.Nodes, domain=aml.Binary,
                                       initialize=lambda model, u: u in self._initial_peptides)

        # node potentials necessary for MTZ subtour elimination constraints
        if self._subtour_elimination_method == 'mtz':
            self.model.u = aml.Var(self.model.Nodes - set([0]), bounds=(1.0, len(self._peptides) - 1))

        # objective of the model
        self.model.Obj = aml.Objective(rule=lambda model: sum(model.y[i] * model.i[i] for i in model.Nodes),
                                             sense=aml.maximize)

        # ensures that every node has the same number of incoming and outgoing connections
        def TransitRule(model, node):
            outgoing = sum(model.x[node, j] for j in model.Nodes if j != node)
            incoming = sum(model.x[i, node] for i in model.Nodes if i != node)
            return outgoing == incoming
        self.model.IncomingOutgoing = aml.Constraint(self.model.Nodes, rule=TransitRule)

        # ensures consistency between x and y, and that only two connections (one in, one out) are allowed
        self.model.Selection = aml.Constraint(self.model.Nodes, rule=lambda model, node: (
            sum(model.x[node, i] + model.x[i, node] for i in model.Nodes if i != node) == 2 * model.y[node]
        ))

    def _build_model_vaccine_constraints(self):
        ''' inserts the constraints necessary to enforce the desired properties of the vaccine
        '''
        # include at most tmax aminoacids in the vaccine
        if self._max_vaccine_aminoacids > 0:
            self.logger.info('Vaccine will have at most %d aminoacids', self._max_vaccine_aminoacids)
            self.model.TMAX = aml.Param(initialize=self._max_vaccine_aminoacids, within=aml.PositiveIntegers)
            self.model.maxNofAminoacids = aml.Constraint(rule=lambda model: (
                None, sum(model.d[i, j] * model.x[i, j] for i, j in model.Arcs), model.TMAX))
        else:
            self.logger.info('No limit on number of aminoacids in the vaccine')
        
        # include at most k epitopes in the vaccine
        if self._max_vaccine_epitopes > 0:
            if self._max_vaccine_epitopes <= 1:
                max_epitopes = int(len(self._peptides) * self._max_vaccine_epitopes)
            else:
                max_epitopes = self._max_vaccine_epitopes

            self.logger.info('Vaccine will have at most %d epitopes', max_epitopes)
            self.model.k = aml.Param(initialize=max_epitopes, within=aml.PositiveIntegers)
            self.model.maxNofEpitopes = aml.Constraint(rule=lambda model: (
                None, sum(model.y[n] for n in model.Nodes), model.k))
        else:
            self.logger.info('No limit on number of epitopes in the vaccine')

        # cover at least t_allele different alleles
        if self._min_allele_coverage > 0:
            if self._min_allele_coverage <= 1:
                min_alleles = int(len(self._alleles) * self._min_allele_coverage)
            else:
                min_alleles = self._min_allele_coverage

            self.logger.info('Vaccine will cover at least %d alleles', min_alleles)
            self.model.t_allele = aml.Param(initialize=min_alleles, within=aml.NonNegativeIntegers)
            self.model.A = aml.Set(initialize=self._allele_bindings.keys())
            self.model.A_I = aml.Set(self.model.A, initialize=lambda model, allele: self._allele_bindings[allele])
            self.model.z = aml.Var(self.model.A, domain=aml.Binary, initialize=0)

            self.model.MinAlleleCovConst = aml.Constraint(rule=lambda model:
                sum(model.z[allele] for allele in model.A) >= model.t_allele)
            self.model.IsAlleleCovConst = aml.Constraint(self.model.A, rule=lambda model, allele:
                sum(model.y[e] for e in model.A_I[allele]) >= model.z[allele])
        else:
            self.logger.info('No minimum allele coverage enforced')

        # cover at least t_var different antigens
        if self._min_antigen_coverage > 0:
            if self._min_antigen_coverage <= 1:
                min_antigens = int(len(self._all_genes) * self._min_antigen_coverage)
            else:
                min_antigens = self._min_antigen_coverage

            self.logger.info('Vaccine will cover at least %d antigens', min_antigens)
            self.model.t_var = aml.Param(initialize=min_antigens, within=aml.NonNegativeIntegers)
            self.model.Q = aml.Set(initialize=self._all_genes)
            self.model.E_var = aml.Set(self.model.Q, initialize=lambda mode, v: self._gene_to_epitope[v])
            self.model.w = aml.Var(self.model.Q, domain=aml.Binary, initialize=0)

            self.model.MinAntigenCovConst = aml.Constraint(rule=lambda model:
                sum(model.w[q] for q in model.Q) >= model.t_var)
            self.model.IsAntigenCovConst = aml.Constraint(self.model.Q, rule=lambda model, q:
                sum(model.y[e] for e in model.E_var[q]) >= model.w[q])
        else:
            self.logger.info('No minimum antigen coverage enforced')

        # select only epitopes with conservation larger than t_c
        # FIXME currently unused since I don't understand how we compute the conservation
        #if self._min_epitope_conservation > 0:
        #    print(self._conservations)
        #    self.model.c = aml.Param(self.model.Nodes, initialize=lambda model, e: self._conservations[e])
        #    self.model.t_c = aml.Param(initialize=self._min_epitope_conservation, within=aml.NonNegativeReals)
        #    self.model.EpitopeConsConst = aml.Constraint(self.model.Nodes, rule=lambda model, e: (
        #        None, (1 - model.c[e]) * model.y[e] - (1 - model.t_c), 0.0))

    def solve(self, options=None):
        # if logging is configured, gurobipy will print messages there *and* on stdout
        # so we replace stdout with a fake file that ignores everything
        class DevNull:
            def write(self, *args, **kwargs):
                pass
            
            def flush(self, *args, **kwargs):
                pass
        
        sys.stdout = DevNull()

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
        self.logger.info('Using %s subtour elimination constraints', self._subtour_elimination_method.upper())
        if self.model is None:
            raise RuntimeError('must call build_model before solve')
        self._solver.set_instance(self.model)
        self._subtour_constraints = 0

        while True:
            res = self._solver.solve(options=options, tee=1, save_results=False, report_timing=True)
            if res.solver.termination_condition != TerminationCondition.optimal:
                raise RuntimeError('Could not solve problem - %s . Please check your settings' % res.Solution.status)

            self._solver.load_vars()
            arcs = self._extract_solution_from_model()
            self.logger.debug('Solution contains the following arcs:', arcs)
            
            tours = self._extract_tours_from_arcs(arcs)
            self.logger.debug('Solution contains the following tour(s):')
            for tt in tours:
                self.logger.debug('   %s', tt)

            if len(tours) > 1 or 0 not in (u for u, v in tours[0]):
                if self._subtour_elimination_method == 'mtz':
                    self._eliminate_subtours_mtz(arcs, tours)
                else:
                    self._eliminate_subtours_dfj(arcs, tours)

                self.logger.debug('========================')
                self.logger.debug('Invalid solution returned!')
                self.logger.debug('%s subtour elimination constraints updated (%d inserted so far)' % (
                    self._subtour_elimination_method.upper(), self._subtour_constraints
                ))
                
                self._reinitialize_model(tours)
            else:
                break
        
        if self._verbosity > 0:
            res.write(num=1)

        self.logger.info('Solved successfully')
        self._result = self._solution_summary(tours[0])
        return self._result

    def _solution_summary(self, tour):
        res = EvaluationResult()
        res.vaccine = []
        res.epitopes = []
        res.peptide_to_gene = {}
        res.peptide_to_hla = {}
        res.immunogenicities = {}

        for i, j in tour[:-1]:
            pep = self._peptides[j]
            cost = self._arc_cost[i, j]
            idx = self._pep_to_index[pep]

            # check overlaps are correct
            assert i == 0 or cost == 9 or str(self._peptides[i])[int(cost)-9:] == str(pep)[:-int(cost)], (
                'wrong overlap', str(self._peptides[i]), str(self._peptides[j]), int(9 - cost)
            )

            res.epitopes.append(pep)
            res.peptide_to_gene[pep] = set([p.gene_id for p in pep.get_all_proteins()])
            res.peptide_to_hla[pep] = set(self._peptide_bindings[idx])
            res.vaccine.append(str(pep)[-int(cost):])
            res.immunogenicities[pep] = self._immunogenicities[idx]
        
        res.vaccine = ''.join(res.vaccine)
        
        return res

    def _reinitialize_model(self, tours):    
        ''' removes the specified tours and re-initializes
            the model with the initial feasible solution
        '''

        for tour in tours:
            for u, v in tour:
                self.model.x[u, v] = self.model.y[u] = 0
        
        for u, v in self._initial_arcs:
            self.model.x[u, v] = self.model.y[u] = 1
        
    def _eliminate_subtours_dfj(self, arcs, tours):
        ''' adds DFJ subtour elimination constraints
        '''
        for tour in tours:
            tour_nodes = set(i for i, _ in tour)
            self._subtour_constraints += 1
            name = 'Subtour_%d' % self._subtour_constraints
            constraint = pmo.constraint(
                body=sum(self.model.x[i, j] for i in tour_nodes for j in tour_nodes if i != j),
                ub=len(tour_nodes) - 1
            )
            setattr(self.model, name, constraint)
            self._solver.add_constraint(getattr(self.model, name))

    def _eliminate_subtours_mtz(self, arcs, tours):
        ''' adds MTZ subtour elimination constraints
        '''
        def add_mtz_constraint(i, j):
            self._subtour_constraints += 1
            name = 'Subtour_%d' % self._subtour_constraints
            constraint = pmo.constraint(
                body=self.model.u[i] - self.model.u[j] + (len(self._peptides) - 1) * self.model.x[i, j],
                ub=len(self._peptides) - 2
            )
            setattr(self.model, name, constraint)
            self._solver.add_constraint(getattr(self.model, name))

        for i, j in arcs.iteritems():
            if i != 0 and j != 0:
                add_mtz_constraint(j, i)
                add_mtz_constraint(i, j)

    def _extract_solution_from_model(self):
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


class EvaluationResult:
    vaccine = None
    epitopes = None
    genes = None
    alleles = None
    peptide_to_hla = None
    peptide_to_gene = None
    immunogenicities = None

    def pretty_print(self, print_fn=None):
        print_fn = print_fn or print

        print_fn('')
        print_fn('==============================')
        print_fn('=       Vaccine Summary      =')
        print_fn('==============================')
        for i in range(0, len(self.vaccine), 50):
            print_fn('%s%s' % (
                'Sequence : ' if i == 0 else ' ' * 11,
                self.vaccine[i:i + 50]
            ))
        
        print_fn('')
        genes_covered = sorted(reduce(lambda r, s: r | s, self.peptide_to_gene.values()))
        alleles_covered = sorted(map(str, reduce(lambda r, s: r | s, self.peptide_to_hla.values())))
        compression = 1 - float(len(self.vaccine)) / sum(len(str(e)) for e in self.epitopes)
        print_fn('  Number of Aminoacids : %d' % len(self.vaccine))
        print_fn('    Number of Epitopes : %d'% len(self.epitopes))
        print_fn('           Compression : %.2f %%' % (100 * compression))
        print_fn('  Total Immunogenicity : %.2f' % sum(self.immunogenicities.values()))
        print_fn('      Covers %3d Genes : ' % len(genes_covered) + ', '.join(genes_covered))
        print_fn('    Covers %3d Alleles : ' % len(alleles_covered) + ', '.join(alleles_covered))

        print_fn('')
        print_fn('')
        print_fn('Epitope Breakdown')
        print_fn('=================')
        print_fn('')
        for i, ep in enumerate(self.epitopes):
            print_fn('       Epitope %03d : %s' % (i + 1, ep))
            print_fn('    Immunogenicity : %.2f' % self.immunogenicities[ep])
            print_fn('    Covers Gene(s) : ' + ', '.join(sorted(self.peptide_to_gene[ep])))
            print_fn('  Covers Allele(s) : ' + ', '.join(sorted(map(str, self.peptide_to_hla[ep]))))
            print_fn('')