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
from pyomo.core.expr.numeric_expr import SumExpression
from pyomo.opt import SolverFactory, TerminationCondition

from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Utility import generate_overlap_graph
from Fred2.Utility import solve_TSP_LKH as _solve_tsp


class DataContainer:
    ''' contains the data structures and common pre-processing logic
        needed by the other classes
    '''

    peptides = None
    pep_to_index = None
    immunogenicities = None
    all_genes = None
    conservations = None
    gene_to_epitope = None
    allele_bindings = None
    peptide_bindings = None
    raw_affinities = None
    alleles = None
    allele_probs = None
    thresh = None
    arc_cost = None

    def __init__(self, predicted_affinities, affinity_threshold=None):
        self.logger = logging.getLogger(self.__class__.__name__)

        if not isinstance(predicted_affinities, EpitopePredictionResult):
            raise ValueError('first input parameter is not of type EpitopePredictionResult')

        self.raw_affinities = predicted_affinities
        self.alleles = copy.deepcopy(self.raw_affinities.columns.values.tolist())
        self.allele_probs = self._fill_allele_probs(self.alleles)
        self.thresh = affinity_threshold or {}
        self._process_parameters()
        self._compute_arcs_cost()

    def _process_parameters(self):
        self.peptides = ['start']
        self.pep_to_index = {}
        self.immunogenicities = [0]
        self.all_genes = set()
        self.conservations = {0: 0}
        self.gene_to_epitope = {}
        self.allele_bindings = {}
        self.peptide_bindings = {}

        method = self.raw_affinities.index.values[0][1]  # TODO find better way of filtering by method
        self.logger.info('Using binding affinities from method "%s"', method)
        res_df = self.raw_affinities.xs(self.raw_affinities.index.values[0][1], level='Method')
        threshold_mask = res_df.apply(lambda x: any(  # FIXME shouldn't we log-transform first if method is one of .... ?
            x[allele] > self.thresh.get(allele, -float('inf'))
            for allele in res_df.columns
        ), axis=1)

        if threshold_mask.sum() == 0:
            raise ValueError('Binding affinity threshold too high, no peptides selected')

        self.logger.info('%d peptides above threshold for at least one allele', threshold_mask.sum())
        self.logger.debug('Breakdown:')
        for col in res_df.columns:
            self.logger.debug('    %s %d', col, np.sum(res_df[col] > self.thresh[col]))

        res_df = res_df[threshold_mask]
        for i, tup in enumerate(res_df.itertuples()):
            peptide = tup[0]
            i += 1  # skip first dummy peptide
            self.peptides.append(peptide)
            self.pep_to_index[peptide] = i
            immunogen = 0
            for allele, bind_affinity in itr.izip(self.alleles, tup[1:]):
                if method in ['smm', 'smmpmbec', 'arb', 'comblibsidney']:
                    # log-transform binding strenghts and thresholds and clip in [0-1]
                    bind_affinity = min(1., max(0.0, 1.0 - math.log(bind_affinity, 50000)))
                    if allele in self.thresh:
                        bind_thr = min(1., max(0.0, 1.0 - math.log(self.thresh.get(allele), 50000)))
                    else:
                        bind_thr = float('-inf')
                else:
                    bind_thr = self.thresh.get(allele)

                if bind_affinity >= bind_thr:
                    self.allele_bindings.setdefault(allele.name, set()).add(i)
                    self.peptide_bindings.setdefault(i, set()).add(allele)
                    immunogen += self.allele_probs[allele.name] * bind_affinity

            self.immunogenicities.append(immunogen)

            proteins = set(peptide.get_all_proteins())
            self.conservations[i] = len(proteins)
            for prot in proteins:
                self.all_genes.add(prot.gene_id)
                self.gene_to_epitope.setdefault(prot.gene_id, set()).add(i)

        # calculate epitope conservation, i.e. number of proteins that contain it over number of genes
        # TODO this does not make any sense ?!
        total = len(self.all_genes)
        self.conservations = {
            e: (v / total) if total > 0 else 1
            for e, v in self.conservations.iteritems()
        }

    @staticmethod
    def _fill_allele_probs(alleles):
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

    def _compute_arcs_cost(self):
        self.logger.debug('Computing arc costs...')
        # FIXME following solution based on suffix trees gives the wrong answer
        #self._arc_cost = generate_overlap_graph(self._peptides[1:])

        self.arc_cost = np.zeros((len(self.peptides), len(self.peptides)))

        for i, peptide_from in enumerate(self.peptides):
            for j, peptide_to in enumerate(self.peptides):
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

                self.arc_cost[i, j] = cost

        # following code checks that the overlaps are computed correctly
        for i, pfrom in enumerate(self.peptides):
            for j, pto in enumerate(self.peptides):
                overlap = int(9 - self.arc_cost[i, j])
                if overlap == 0 or i == 0 or j == 0:
                    continue

                pfrom, pto = str(pfrom), str(pto)
                assert pfrom[-overlap:] == pto[:overlap], (i, pfrom, j, pto, overlap)


class MosaicVaccineILP:

    def __init__(self, predicted_affinities, affinity_threshold=None,
                 max_vaccine_aminoacids=100, max_vaccine_epitopes=None,
                 min_allele_coverage=None, min_antigen_coverage=None,
                 min_epitope_conservation=None, verbosity=0,
                 model_type='dfj', solver='gurobi_persistent'):

        self.logger = logging.getLogger(self.__class__.__name__)

        if not isinstance(predicted_affinities, EpitopePredictionResult):
            raise ValueError('first input parameter is not of type EpitopePredictionResult')
        elif model_type not in ['dfj', 'mtz-g', 'mtz-l']:
            raise ValueError('model type must be either "dfj", "mtz-g" or "mtz-l"')
        elif (model_type == 'mtz-l' or model_type == 'dfj') and not solver.endswith('_persistent'):
            raise ValueError('persistent solver required for lazy constraints')

        self._data = DataContainer(predicted_affinities, affinity_threshold)

        self._model_type = model_type
        self._model = None
        self._solver = SolverFactory(solver)
        self._verbosity = verbosity
        self._result = None

        # desired vaccine properties
        self._max_vaccine_aminoacids = max_vaccine_aminoacids
        self._max_vaccine_epitopes = max_vaccine_epitopes
        self._min_allele_coverage = min_allele_coverage
        self._min_antigen_coverage = min_antigen_coverage
        self._min_epitope_conservation = min_epitope_conservation

    def _find_initial_solution(self):
        ''' find an initial solution to initialize the variables
            greedily select the peptides with highest immunogenicity
            until the maximum length of the vaccine is reached

            returns the set of epitopes and the set of arcs
        '''
        sorted_by_immunogenicity = sorted(
            self._data.peptides,
            key=lambda pep: self._data.immunogenicities[self._data.pep_to_index.get(pep, 0)],
            reverse=True
        )

        solution_peptides, solution_arcs = [0], []
        aminoacids_length = total_immunogenicity = 0
        for i, pep in enumerate(sorted_by_immunogenicity):
            if self._max_vaccine_epitopes > 0 and i >= self._max_vaccine_epitopes:
                break

            idx = self._data.pep_to_index[pep]
            arc_cost = self._data.arc_cost[solution_peptides[-1], idx]
            if self._max_vaccine_aminoacids > 0 and aminoacids_length + arc_cost > self._max_vaccine_aminoacids:
                break

            aminoacids_length += arc_cost
            total_immunogenicity += self._data.immunogenicities[idx]
            solution_arcs.append((solution_peptides[-1], idx))
            solution_peptides.append(idx)

        solution_arcs.append((solution_peptides[-1], 0))
        self.logger.debug('The model will be initialized with the following feasible '
                          'solution with cost %d and immunogenicity %.2f:', aminoacids_length, total_immunogenicity)
        for u, v in solution_arcs:
            self.logger.debug('    %5d (%9s) -> %5d (%9s) (cost: %2d, immunogenicity: %4.2f)',
                              u, self._data.peptides[u], v, self._data.peptides[v],
                              self._data.arc_cost[u, v], self._data.immunogenicities[v])

        return set(solution_peptides), set(solution_arcs)

    def build_model(self):
        self._initial_peptides, self._initial_arcs = self._find_initial_solution()

        self.logger.info('Building model...')

        self._model = aml.ConcreteModel()
        self._build_base_tsp_model()
        self._build_model_vaccine_constraints()

        self.logger.info('Model built!')

    def _build_base_tsp_model(self):
        ''' builds a MIP model suitable for solving the orienteering problem
        '''
        # graph objects
        self._model.Nodes = aml.RangeSet(0, len(self._data.peptides) - 1)
        self._model.Arcs = aml.Set(initialize=[
            (i, j)
            for i in xrange(len(self._data.peptides))
            for j in xrange(len(self._data.peptides))
            if i != j
        ])

        # reward of nodes and cost of arcs
        self._model.i = aml.Param(self._model.Nodes, initialize=lambda model, i: self._data.immunogenicities[i])
        self._model.d = aml.Param(self._model.Arcs, initialize=lambda model, i, j: self._data.arc_cost[i][j])

        # indicator variables for nodes and arcs
        self._model.x = aml.Var(self._model.Arcs, domain=aml.Binary, bounds=(0, 1),
                                initialize=lambda model, u, v: (u, v) in self._initial_arcs)
        self._model.y = aml.Var(self._model.Nodes, domain=aml.Binary,
                                initialize=lambda model, u: u in self._initial_peptides)

        # node potentials necessary for MTZ subtour elimination constraints
        if self._model_type.startswith('mtz'):
            self._model.u = aml.Var(self._model.Nodes - set([0]), bounds=(1.0, len(self._data.peptides) - 1))

        # if required, add subtour elimination constraints now
        if self._model_type == 'mtz-g':
            self._model.SubTour = aml.Constraint((
                (i, j)
                for i in xrange(1, len(self._data.peptides))
                for j in xrange(1, len(self._data.peptides))
                if i != j
            ), rule=lambda model, i, j: (
                None,
                model.u[i] - model.u[j] + (len(self._data.peptides) - 1) * model.x[i, j],
                len(self._data.peptides) - 2
            ))

        # objective of the model
        self._model.Obj = aml.Objective(rule=lambda model: sum(model.y[i] * model.i[i] for i in model.Nodes),
                                        sense=aml.maximize)

        # ensures that every node has the same number of incoming and outgoing connections
        def TransitRule(model, node):
            outgoing = sum(model.x[node, j] for j in model.Nodes if j != node)
            incoming = sum(model.x[i, node] for i in model.Nodes if i != node)
            return outgoing == incoming
        self._model.IncomingOutgoing = aml.Constraint(self._model.Nodes, rule=TransitRule)

        # ensures consistency between x and y, and that only two connections (one in, one out) are allowed
        self._model.Selection = aml.Constraint(self._model.Nodes, rule=lambda model, node: (
            sum(model.x[node, i] + model.x[i, node] for i in model.Nodes if i != node) == 2 * model.y[node]
        ))

    def _build_model_vaccine_constraints(self):
        ''' inserts the constraints necessary to enforce the desired properties of the vaccine
        '''
        # include at most tmax aminoacids in the vaccine
        if self._max_vaccine_aminoacids > 0:
            self.logger.info('Vaccine will have at most %d aminoacids', self._max_vaccine_aminoacids)
            self._model.TMAX = aml.Param(initialize=self._max_vaccine_aminoacids, within=aml.PositiveIntegers)
            self._model.maxNofAminoacids = aml.Constraint(rule=lambda model: (
                None, sum(model.d[i, j] * model.x[i, j] for i, j in model.Arcs), model.TMAX))
        else:
            self.logger.info('No limit on number of aminoacids in the vaccine')

        # include at most k epitopes in the vaccine
        if self._max_vaccine_epitopes > 0:
            if self._max_vaccine_epitopes <= 1:
                max_epitopes = int(len(self._data.peptides) * self._max_vaccine_epitopes)
            else:
                max_epitopes = self._max_vaccine_epitopes

            self.logger.info('Vaccine will have at most %d epitopes', max_epitopes)
            self._model.k = aml.Param(initialize=max_epitopes, within=aml.PositiveIntegers)
            self._model.maxNofEpitopes = aml.Constraint(rule=lambda model: (
                None, sum(model.y[n] for n in model.Nodes), model.k))
        else:
            self.logger.info('No limit on number of epitopes in the vaccine')

        # cover at least t_allele different alleles
        if self._min_allele_coverage > 0:
            if self._min_allele_coverage <= 1:
                min_alleles = int(len(self._data.alleles) * self._min_allele_coverage)
            else:
                min_alleles = self._min_allele_coverage

            self.logger.info('Vaccine will cover at least %d alleles', min_alleles)
            self._model.t_allele = aml.Param(initialize=min_alleles, within=aml.NonNegativeIntegers)
            self._model.A = aml.Set(initialize=self._data.allele_bindings.keys())
            self._model.A_I = aml.Set(self._model.A, initialize=lambda model, allele: self._data.allele_bindings[allele])
            self._model.z = aml.Var(self._model.A, domain=aml.Binary, initialize=0)

            self._model.MinAlleleCovConst = aml.Constraint(
                rule=lambda model: sum(model.z[allele] for allele in model.A) >= model.t_allele)
            self._model.IsAlleleCovConst = aml.Constraint(
                self._model.A, rule=lambda model, allele: sum(model.y[e] for e in model.A_I[allele]) >= model.z[allele])
        else:
            self.logger.info('No minimum allele coverage enforced')

        # cover at least t_var different antigens
        if self._min_antigen_coverage > 0:
            if self._min_antigen_coverage <= 1:
                min_antigens = int(len(self._data.all_genes) * self._min_antigen_coverage)
            else:
                min_antigens = self._min_antigen_coverage

            self.logger.info('Vaccine will cover at least %d antigens', min_antigens)
            self._model.t_var = aml.Param(initialize=min_antigens, within=aml.NonNegativeIntegers)
            self._model.Q = aml.Set(initialize=self._data.all_genes)
            self._model.E_var = aml.Set(self._model.Q, initialize=lambda mode, v: self._data.gene_to_epitope[v])
            self._model.w = aml.Var(self._model.Q, domain=aml.Binary, initialize=0)

            self._model.MinAntigenCovConst = aml.Constraint(
                rule=lambda model: sum(model.w[q] for q in model.Q) >= model.t_var)
            self._model.IsAntigenCovConst = aml.Constraint(
                self._model.Q, rule=lambda model, q: sum(model.y[e] for e in model.E_var[q]) >= model.w[q])
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
            arcs = self._extract_solution_from_model()
            self.logger.debug('Solution contains the following arcs:', arcs)

            tours = self._extract_tours_from_arcs(arcs)
            self.logger.debug('Solution contains the following tour(s):')
            for tt in tours:
                self.logger.debug('   %s', tt)

            if len(tours) > 1 or 0 not in (u for u, v in tours[0]):
                if self._model_type == 'mtz-l':
                    self._eliminate_subtours_mtz(arcs, tours)
                elif self._model_type == 'dfj':
                    self._eliminate_subtours_dfj(arcs, tours)
                else:
                    raise RuntimeError('This is a bug: I should not have reached this point. '
                                       'This means that the MTZ model is mis-specified.')

                self.logger.debug('========================')
                self.logger.debug('Invalid solution returned!')
                self.logger.debug('%s subtour elimination constraints updated (%d inserted so far)' % (
                    self._model_type.upper(), self._subtour_constraints
                ))
                self.logger.debug('========================')

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
            pep = self._data.peptides[j]
            cost = self._data.arc_cost[i, j]
            idx = self._data.pep_to_index[pep]

            # check overlaps are correct
            assert i == 0 or cost == 9 or str(self._data.peptides[i])[int(cost)-9:] == str(pep)[:-int(cost)], (
                'wrong overlap', str(self._data.peptides[i]), str(self._data.peptides[j]), int(9 - cost)
            )

            res.epitopes.append(pep)
            res.peptide_to_gene[pep] = set([p.gene_id for p in pep.get_all_proteins()])
            res.peptide_to_hla[pep] = set(self._data.peptide_bindings[idx])
            res.vaccine.append(str(pep)[-int(cost):])
            res.immunogenicities[pep] = self._data.immunogenicities[idx]

        res.vaccine = ''.join(res.vaccine)

        return res

    def _reinitialize_model(self, tours):
        ''' removes the specified tours and re-initializes
            the model with the initial feasible solution
        '''

        for tour in tours:
            for u, v in tour:
                self._model.x[u, v] = self._model.y[u] = 0

        for u, v in self._initial_arcs:
            self._model.x[u, v] = self._model.y[u] = 1

    def _eliminate_subtours_dfj(self, arcs, tours):
        ''' adds DFJ subtour elimination constraints
        '''
        for tour in tours:
            tour_nodes = set(i for i, _ in tour)
            self._subtour_constraints += 1
            name = 'Subtour_%d' % self._subtour_constraints
            constraint = pmo.constraint(
                body=sum(self._model.x[i, j] for i in tour_nodes for j in tour_nodes if i != j),
                ub=len(tour_nodes) - 1
            )
            setattr(self._model, name, constraint)
            self._solver.add_constraint(getattr(self._model, name))

    def _eliminate_subtours_mtz(self, arcs, tours):
        ''' adds MTZ subtour elimination constraints
        '''
        def add_mtz_constraint(i, j):
            self._subtour_constraints += 1
            name = 'Subtour_%d' % self._subtour_constraints
            constraint = pmo.constraint(
                body=self._model.u[i] - self._model.u[j] + (len(self._data.peptides) - 1) * self._model.x[i, j],
                ub=len(self._data.peptides) - 2
            )
            setattr(self._model, name, constraint)
            self._solver.add_constraint(getattr(self._model, name))

        for i, j in arcs.iteritems():
            if i != 0 and j != 0:
                add_mtz_constraint(j, i)
                add_mtz_constraint(i, j)

    def _extract_solution_from_model(self):
        ''' returns a dictionary i -> j containing the tour found by the model
        '''
        tour_arcs = {i: j for i, j in self._model.Arcs
                     if 0.98 <= self._model.x[i, j].value <= 1.5}
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
