# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
   .. module:: Mosaic
   :synopsis:  This class implements the epitope selection functionality
    to construct so called mosaic vaccines

    This module builds upon Coopr's Pyomo, an embedded algebraic modeling
    languages [2].

    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice


.. moduleauthor:: schubert

"""

from __future__ import division, print_function

import itertools as itr
import copy
import math
import time
import numpy as np

import multiprocessing as mp

import pyomo.environ as aml
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.core.expr.numeric_expr import SumExpression

from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Utility import generate_overlap_graph
from Fred2.Utility import solve_TSP_LKH as _solve_tsp


class MosaicVaccineILP(object):

    def __init__(self, aff, threshold=None, t_max=100, k=999999999999, solver="cplex", verbosity=0):
        if not isinstance(aff, EpitopePredictionResult):
            raise ValueError("first input parameter is not of type EpitopePredictionResult")

        self.__raw_affinities = aff
        self.__alleles = copy.deepcopy(self.__raw_affinities.columns.values.tolist())
        self.__allele_probs = self.__fill_allele_probs(self.__alleles)
        self.__solver = SolverFactory(solver)
        self.__verbosity = verbosity
        self.__changed = True
        self.__t_max = t_max
        self.__result = None
        self.__thresh = {} if threshold is None else threshold
        self.__pep_to_index = {}
        self.__k = k

        self.__init_model_parameters()

        arc_cost = generate_overlap_graph(self.__peptides[1:])

        ################################################################################################################
        # initializes MIP model for OR with MTZ sub tour elimination
        ################################################################################################################
        model = aml.ConcreteModel()

        ################################################################################################################
        print("Sets")
        ################################################################################################################
        model.Nodes     = aml.RangeSet(0, len(self.__peptides) - 1)
        model.Q         = aml.Set(initialize=self.__variations)
        model.A         = aml.Set(initialize=self.__allele_bindings.keys())
        model.E_var     = aml.Set(model.Q, initialize=lambda mode, v: self.__epitope_variations[v])
        model.A_I       = aml.Set(model.A, initialize=lambda model, allele: self.__allele_bindings[allele])
        model.Arcs      = model.Nodes * model.Nodes # FIXME dont know if this includes self references as well....

        ################################################################################################################
        print("Params")
        ################################################################################################################
        model.i         = aml.Param(model.Nodes, initialize=lambda model, i: self.__immunogenicities[i])
        model.d         = aml.Param(model.Arcs, initialize=lambda model, i, j: arc_cost[i][j])
        model.k         = aml.Param(initialize=self.__k, within=aml.PositiveIntegers, mutable=True)
        model.c         = aml.Param(model.Nodes, initialize=lambda model, e: self.__conservations[e], mutable=True)
        model.TMAX      = aml.Param(initialize=self.__t_max, within=aml.PositiveIntegers, mutable=True)
        model.t_allele  = aml.Param(initialize=0, within=aml.NonNegativeIntegers, mutable=True)
        model.t_var     = aml.Param(initialize=0, within=aml.NonNegativeIntegers, mutable=True)
        model.t_c       = aml.Param(initialize=0.0, within=aml.NonNegativeReals, mutable=True)

        ################################################################################################################
        print("Variables")
        ################################################################################################################
        model.x         = aml.Var(model.Arcs, domain=aml.Binary, bounds=(0, 1), initialize=0)
        model.z         = aml.Var(model.A, domain=aml.Binary, initialize=0)
        model.y         = aml.Var(model.Nodes, domain=aml.Binary, initialize=0)
        model.w         = aml.Var(model.Q, domain=aml.Binary, initialize=0)
        model.u         = aml.Var(model.Nodes - set([0]), bounds=(1.0, len(self.__peptides) - 1))

        ################################################################################################################
        model.Obj = aml.Objective(rule=lambda model: sum(model.y[i] * model.i[i] for i in model.Nodes), sense=aml.maximize)

        ################################################################################################################
        print("Constraints")
        ################################################################################################################
        # conecitvity constraint
        model.Conn1 = aml.Constraint(model.Nodes, rule=lambda model, node: (
            None, sum(model.x[node, j] for j in model.Nodes if j != node), 1.0
        ))
        model.Conn2 = aml.Constraint(model.Nodes, rule=lambda model, node: (
            None, sum(model.x[i, node] for i in model.Nodes if i != node), 1.0
        ))
        # Equality constraint
        def Equal_rule(model, node):
            aa = sum(model.x[node, j] for j in model.Nodes if j != node)
            bb = sum(model.x[i, node] for i in model.Nodes if i != node)
            if isinstance(aa, SumExpression) or isinstance(bb, SumExpression):
                return aa == bb
            else:
                return Constraint.Feasible

        model.Equal = aml.Constraint(model.Nodes, rule=Equal_rule)
        # Knapsack Constriant
        model.Knapsack = aml.Constraint(rule=lambda model: (
            None, sum(model.d[i, j] * model.x[i, j] for i, j in model.Arcs), self.__t_max
        ))

        print("selection specific start")
        model.IsNodeSelected        = aml.Constraint(model.Nodes, rule=lambda model, n: sum(model.x[n, m] for m in model.Nodes if n != m) >= model.y[n])
        model.nofEpitopes           = aml.Constraint(rule=lambda model: sum(model.y[n] for n in model.Nodes) <= model.k)
        model.IsAlleleCovConst      = aml.Constraint(model.A, rule=lambda model, allele: sum(model.y[e] for e in model.A_I[allele]) >= model.z[allele])
        model.MinAlleleCovConst     = aml.Constraint(rule=lambda model: sum(model.z[allele] for allele in model.A) >= model.t_allele)
        model.IsAntigenCovConst     = aml.Constraint(model.Q, rule=lambda model, q: sum(model.y[e] for e in model.E_var[q]) >= model.w[q])
        model.MinAntigenCovConst    = aml.Constraint(rule=lambda model: sum(model.w[q] for q in model.Q) >= model.t_var)
        model.EpitopeConsConst      = aml.Constraint(model.Nodes, rule=lambda model, e: (1 - model.c[e]) * model.y[e] <= 1 - model.t_c)

        #constraints
        model.IsAlleleCovConst.deactivate()
        model.MinAlleleCovConst.deactivate()
        model.IsAntigenCovConst.deactivate()
        model.MinAntigenCovConst.deactivate()
        model.EpitopeConsConst.deactivate()
        print("selection specific end")

        print("MTZ start")
        # Subout Elimination MTZ
        model.SubTour = aml.Constraint((
            (i, j)
            for i in xrange(1, len(self.__peptides))
            for j in xrange(1, len(self.__peptides))
            if i != j
        ), rule=lambda model, i, j: (
            None, model.u[i] - model.u[j] + (len(self.__peptides) - 1) * model.x[i, j], len(self.__peptides) - 2
        ))
        print("MTZ end")

        self.instance = model
        if verbosity:
            model.pprint()

    def __init_model_parameters(self):
        self.__peptides = ["start"]
        self.__immunogenicities = [0]
        self.__variations = set()
        self.__conservations = {0: 0}
        self.__epitope_variations = {}
        self.__allele_bindings = {}

        method = self.__raw_affinities.index.values[0][1]  # TODO find better way of filtering by method
        res_df = self.__raw_affinities.xs(self.__raw_affinities.index.values[0][1], level="Method")
        res_df = res_df[res_df.apply(lambda x: any(  # FIXME shouldn't we log-transform first if method is one of .... ?
            x[allele] > self.__thresh.get(allele.name, -float("inf"))
            for allele in res_df.columns
        ), axis=1)]

        for i, tup in enumerate(res_df.itertuples()):
            peptide = tup[0]
            i += 1  # skip first dummy peptide
            self.__peptides.append(peptide)
            self.__pep_to_index[peptide] = i
            immunogen = 0
            for allele, bind_affinity in itr.izip(self.__alleles, tup[1:]):
                if method in ["smm", "smmpmbec", "arb", "comblibsidney"]:
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

    def set_k(self, k):
        """
            Sets the number of epitopes to select

            :param int k: The number of epitopes
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        tmp = self.instance.k.value
        try:
            getattr(self.instance, str(self.instance.k)).set_value(int(k))
            self.__changed = True
        except ValueError:
            self.__changed = False
            getattr(self.instance, str(self.instance.k)).set_value(int(tmp))
            raise ValueError('set_k', 'An error has occurred during setting parameter k. Please check if k is integer.')


    def set_Tmax(self, t_max):
        """
            Sets the number of epitopes to select

            :param int t_max: The maximum length of the mosaic vaccine
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        old_tmax = self.instance.TMAX.value
        try:
            getattr(self.instance, str(self.instance.TMAX)).set_value(int(t_max))
            self.__changed = True
        except ValueError:
            self.__changed = False
            getattr(self.instance, str(self.instance.TMAX)).set_value(int(old_tmax))
            raise ValueError('set_t_max',
                             'An error has occurred during setting parameter Tmax. Please check if t_max is integer.')

    def activate_allele_coverage_const(self, minCoverage):
        """
            Enables the allele coverage constraint

            :param float minCoverage: Percentage of alleles which have to be covered [0,1]
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        # parameter
        mc = self.instance.t_allele.value

        try:
            getattr(self.instance, str(self.instance.t_allele)).set_value(int(len(self.__alleles) * minCoverage))
            #variables

            #constraints
            self.instance.IsAlleleCovConst.activate()
            self.instance.MinAlleleCovConst.activate()
            self.__changed = True
        except ValueError:
            getattr(self.instance, str(self.instance.t_allele)).set_value(mc)
            self.__changed = False
            raise ValueError(
                'activate_allele_coverage_const","An error occurred during activation of of the allele coverage constraint. ' +
                'Please check your specified minimum coverage parameter to be in the range of 0.0 and 1.0.')

    def deactivate_allele_coverage_const(self):
        """
            Deactivates the allele coverage constraint
        """

        # parameter
        self.__changed = True

        #constraints
        self.instance.IsAlleleCovConst.deactivate()
        self.instance.MinAlleleCovConst.deactivate()

    def activate_antigen_coverage_const(self, t_var):
        """
            Activates the variation coverage constraint

            :param int t_var: The number of epitopes which have to come from each variation
            :raises ValueError: If the input variable is not in the same domain as the parameter

        """
        tmp = self.instance.t_var.value
        try:
            getattr(self.instance, str(self.instance.t_var)).set_value(int(len(self.instance.Q)*t_var))
            self.instance.IsAntigenCovConst.activate()
            self.instance.MinAntigenCovConst.activate()
            self.__changed = True
        except ValueError:
            getattr(self.instance, str(self.instance.t_var)).set_value(int(tmp))
            self.instance.IsAntigenCovConst.deactivate()
            self.instance.MinAntigenCovConst.deactivate()
            self.__changed = False
            raise ValueError("activate_antigen_coverage_const",
                            "An error has occurred during activation of the coverage constraint. Please make sure your input is an integer.")

    def deactivate_antigen_coverage_const(self):
        """
            Deactivates the variation coverage constraint
        """
        self.__changed = True
        self.instance.IsAntigenCovConst.deactivate()
        self.instance.MinAntigenCovConst.deactivate()

    def activate_epitope_conservation_const(self, t_c, conservation=None):
        """
            Activates the epitope conservation constraint

            :param float t_c: The percentage of conservation an epitope has to have [0.0,1.0].
            :param: conservation: A dict with key=:class:`~Fred2.Core.Peptide.Peptide` specifying a different
                                  conservation score for each :class:`~Fred2.Core.Peptide.Peptide`
            :type conservation: dict(:class:`~Fred2.Core.Peptide.Peptide`,float)
            :raises ValueError: If the input variable is not in the same domain as the parameter
        """
        if t_c < 0 or t_c > 1:
            raise ValueError("activate_epitope_conservation_const",
                            "The conservation threshold is out of its numerical bound. It has to be between 0.0 and 1.0.")

        self.__changed = True
        getattr(self.instance, str(self.instance.t_c)).set_value(float(t_c))
        if conservation is not None:
            for i,e in enumerate(self.__peptides):
                if e in conservation:
                    getattr(self.instance, str(self.instance.c))[i] = conservation[e]
                else:
                    getattr(self.instance, str(self.instance.c))[i] = 0.0

        self.instance.EpitopeConsConst.activate()

    def deactivate_epitope_conservation_const(self):
        """
            Deactivates epitope conservation constraint
        """
        self.__changed = True
        self.instance.EpitopeConsConst.deactivate()

    def solve(self, options=None):
        """
            solves the model optimally
        """
        if not self.__changed:
            return self.__result

        options = dict() if options is None else options
        instance = self.instance

        res = self.__solver.solve(instance, options=options, tee=1)
        instance.solutions.load_from(res)
        if self.__verbosity > 0:
            res.write(num=1)

        if res.solver.termination_condition != TerminationCondition.optimal:
            raise RuntimeError("Could not solve problem - " + str(res.Solution.status) + ". Please check your settings")

        sol = []
        unsorted = dict([(i, j) for i, j in instance.Arcs if 0.98 <= instance.x[i, j].value <= 1.5])
        i = 0

        while unsorted:
            j = unsorted[i]
            sol.append((i, j))
            del unsorted[i]
            i = j

        self.__result = sol
        seqs = [self.__peptides[i] for i,j in sol if i != 0]
        return sol, seqs