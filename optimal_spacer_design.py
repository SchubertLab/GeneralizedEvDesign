# as part of this package.
"""
.. module:: EpitopeAssembly.EpitopeAssembly
   :synopsis: This module contains all classes for EpitopeAssembly.
.. moduleauthor:: schubert & dorigatti
"""

# this module was modified as follows:
#  - we only kept the optimal spacer design and removed everything else
#  - updated documentation and renamed EpitopeAssemblyWithSpacer
#  - tidied-up code


import copy
import itertools as itr
import logging
import math
import multiprocessing as mp
import os
import subprocess
import warnings
from collections import defaultdict
from tempfile import NamedTemporaryFile

from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
from tqdm import tqdm

import utilities
from epytope.CleavagePrediction.PSSM import APSSMCleavageSitePredictor
from epytope.Core import Allele, Peptide, Protein
from epytope.Core.Base import ACleavageSitePrediction, AEpitopePrediction
from epytope.Core.Generator import generate_peptides_from_proteins
from epytope.EpitopePrediction.PSSM import APSSMEpitopePrediction


class OptimalSpacerDesign(object):

    def __init__(self, peptides, cleav_pred, epi_pred, alleles, k=5, en=9,
                 threshold=None, solver="glpk", alpha=0.99, beta=0,
                 verbosity=0):
        """

        :param peptides: A list of :class:`~epytope.Core.Peptide.Peptide` which shell be arranged
        :type peptides: list(:class:`~epytope.Core.Peptide.Peptide`)
        :param cleav_pred: A :class:`~epytope.CleavagePrediction.PSSM.APSSMCleavageSitePredictor` (PSSM only)
        :type cleav_pred: :class:`~epytope.Core.Base.ACleavageSitePredictor`
        :param epi_pred: A :class:`~epytope.EpitopePrediction.PSSM.APSSMEpitopePrediction` (PSSM only)
        :type epi_pred: :class:`~epytope.Core.Base.AEpitopePredictor`
        :param alleles: A list of :class:`~epytope.Core.Allele.Allele` for which predictions should be made
        :type alleles: list(:class:`~epytope.Core.Allele.Allele`)
        :param int k: List of spacer lengths to optimize on
        :param int en: Length of epitopes
        :param dict(str,float) threshold: A dictionary specifying the epitope prediction threshold for each
                                          :class:`~epytope.Core.Allele.Allele`
        :param str solver: Specifies the solver to use (must be callable by pyomo)
        :param float alpha: Specifies how how much junction-cleavage score can be sacrificed  to gain lower
                            neo-immunogenicity
        :param float beta: Specifies how how much noe-immunogenicity score can be sacrificed to gain lower non-junction
                           cleavage score
        :param int verbosity: Specifies how verbos the class will be, 0 means normal, >0 debug mode
        """

        self.logger = logging.getLogger(self.__class__.__name__)

        # test input
        if not isinstance(cleav_pred, APSSMCleavageSitePredictor):
            raise ValueError("Second input must be a PSSM-based cleavage site predictor.")

        if not isinstance(epi_pred, APSSMEpitopePrediction):
            raise ValueError("Third input must be a PSSM-based epitope predictor.")

        if en not in epi_pred.supportedLength:
            raise ValueError("Specified epitope length of en=%i is not supported by %s" % (en, epi_pred.name))

        _alleles = [copy.deepcopy(a) for a in alleles if a in epi_pred.supportedAlleles]

        if not _alleles:
            raise ValueError("Specified alleles are not supported by %s" % epi_pred.name)

        # infere probability if not already set

        # test if allele prob is set, if not set allele prob uniform
        # if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
            if a.prob is None:
                no_prob.append(a)
            else:
                prob.append(a)

        if len(no_prob) > 0:
            # group by locus
            no_prob_grouped = {}
            prob_grouped = {}
            for a in no_prob:
                no_prob_grouped.setdefault(a.locus, []).append(a)
            for a in prob:
                prob_grouped.setdefault(a.locus, []).append(a)

            for g, v in no_prob_grouped.items():
                total_loc_a = len(v)
                if g in prob_grouped:
                    remaining_mass = 1.0 - sum(a.prob for a in prob_grouped[g])
                    for a in v:
                        a.prob = remaining_mass/total_loc_a
                else:
                    for a in v:
                        a.prob = 1.0/total_loc_a
        probs = {a.name: a.prob for a in _alleles}
        if verbosity:
            for a in _alleles:
                print(a.name, a.prob)

        self.spacer = {}

        # start constructing model
        self.__solver = solver
        self.__verbosity = verbosity
        self.__changed = True
        self.__k = k
        self.__result = None
        self.__thresh = {a.name: 0 for a in alleles} if threshold is None else threshold
        self.__alleles = _alleles
        self.__epi_pred = epi_pred
        self.__clev_pred = cleav_pred
        self.__en = en
        self.__alpha = alpha
        self.__beta = beta
        self.__peptides = list(peptides)

    def solve(self, threads=None, options=None):
        """
        Solve the epitope assembly problem with spacers optimally using integer linear programming.

        .. note::

            This can take quite long and should not be done for more and 30 epitopes max!
            Also, one has to disable pre-solving steps in order to use this model.

        :param int threads: Number of threads used for spacer design.
                            Be careful, if options contain solver threads it will allocate threads*solver_threads cores!
        :param dict(str,str) options: Solver specific options as keys and parameters as values
        :return: A list of ordered :class:`~epytope.Core.Peptide.Peptide`
        :rtype: list(:class:`~epytope.Core.Peptide.Peptide`)
        """
        def __load_model(name, model):
            return getattr(__import__("epytope.Data.pssms."+name+".mat."+model, fromlist=[model]), model)

        options = dict() if options is None else options

        self.logger.debug('Preparing parameters...')
        cn = min(self.__clev_pred.supportedLength)
        cl_pssm = __load_model(self.__clev_pred.name, self.__clev_pred.name+"_"+str(cn))
        cleav_pos = self.__clev_pred.cleavagePos
        en = self.__en
        epi_pssms = {}
        allele_prob = {}
        delete_alleles = []
        if self.__epi_pred.name in ["smm", "smmpmbec", "comblibsidney"]:
            self.__thresh = {k: (1-math.log(v, 50000) if v != 0 else 0) for k, v in self.__thresh.items()}
        for a in self.__alleles:
            allele_prob[a.name] = a.prob
            try:
                pssm = __load_model(self.__epi_pred.name, "%s_%i" % (self.__epi_pred.convert_alleles([a])[0], en))
                if self.__epi_pred.name in ["smm", "smmpmbec", "comblibsidney"]:
                    for j, v in pssm.items():
                        for aa, score in v.items():
                            epi_pssms[j, aa, a.name] = 1/10. - math.log(math.pow(10, score), 50000)
                else:
                    for j, v in pssm.items():
                        for aa, score in v.items():
                            epi_pssms[j, aa, a.name] = score
            except ImportError:
                delete_alleles.append(a)

        # delete alleles from model that generated an error while loading matrices
        self.logger.debug('Deleted %d alleles', len(delete_alleles))
        for a in delete_alleles:
            del allele_prob[a.name]
            del self.__thresh[a.name]

        if not epi_pssms:
            raise ValueError("Selected alleles with epitope length are not supported by the prediction method.")

        # run spacer designs in parallel using multiprocessing
        task_count = len(self.__peptides) * (len(self.__peptides) - 1) * len(self.__k)
        task_params = (
            (
                str(ei), str(ej), i, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob,
                self.__alpha, self.__thresh, self.__solver, self.__beta, options
            )
            for ei, ej in itr.product(self.__peptides, repeat=2) if ei != ej
            for i in self.__k
        )
        res = utilities.parallel_apply(_runs_lexmin, task_params, threads)

        self.logger.debug('Designing spacers between all pairs of epitopes...')
        opt_spacer = {}
        adj_matrix = {}
        inf = float("inf")
        for i, (ei, ej, score, epi, spacer, c1, c2, non_c) in tqdm(enumerate(res), total=task_count):
            if adj_matrix.get((ei, ej), inf) > -min(c1, c2):
                adj_matrix[(ei, ej)] = -min(c1, c2)
                opt_spacer[(ei, ej)] = spacer

        self.logger.info('Optimal spacers found!')
        self.spacer = opt_spacer
        self.adj_matrix = adj_matrix
        return self


def _runs_lexmin(*a):
    """
    private used to unpack arguments send to processes
    :param a:
    :return: ei,ej,cleavage_score,imm_score,c1_score,c2_score,non-junction_score
    """
    spacer, cleav, epi, good_cleav, bad_cleav, non_c = _spacer_design(*a)
    return a[0], a[1], cleav, epi, spacer, good_cleav, bad_cleav, non_c


def _spacer_design(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, alpha,
                   thresh, solver, beta=0, options=None):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: start epitope
        :param str ej: end epitope
        :param int k: length of spacer
        :param int en: epitope length
        :param int cn: cleavage-site string length
        :param dict(int,dict(string,float)) cl_pssm: a cleavage site prediction PSSM as dict-of-dicts
        :param dict(int,dict(string,float)) epi_pssm: a epitope prediction PSSM as dict-of-dicts
        :param int cleav_pos: integer specifying at which AA within the epitope of length cn the cleave is predicted
        :param dict(strin,float) allele_prob: a dict of HLA alleles as string (i.e. A*02:01) and probabilities [0,1]
        :param float alpha: specifies the first-order influence on the objectives [0,1]
        :param float thresh: specifies at which score a peptide is considered as epitope
        :param string solver: string specifying which ILP solver should be used
        :param dict(str,str) options: solver specific options as keys and parameters as values
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """
    options = dict() if options is None else options

    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g = len(ei)+k-cleav_pos
        c1 = sum(cl_pssm[j][seq[i + j]] for j in range(cn))+cl_pssm.get(-1, {}).get("con", 0)
        c2 = sum(cl_pssm[j][seq[g + j]] for j in range(cn))+cl_pssm.get(-1, {}).get("con", 0)
        non_c = sum(
            sum(
                cl_pssm[j][seq[k + j]] for j in range(cn) if k != i and k != g
            ) + cl_pssm.get(-1, {}).get("con", 0)
            for k in range(len(seq) - (cn - 1))
        )

        imm = sum(prob * sum(
            max(
                sum(
                    epi_pssms[j, seq[i + j], a] for j in range(en)
                ) + epi_pssms.get((-1, "con", a), 0) - thresh[a],
                0
            ) for i in range(len(seq)-en)
        ) for a, prob in allele_prob.items())

        return "", (c1+c2)/2, imm, c1, c2, non_c

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i, v in p.items():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i, v in p.items():
            for a, score in v.items():
                norm.setdefault(i, {}).update({a: (score-min_p)/(max_p - min_p)})
        return norm

    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = list('ACDEFGHIKLMNPQRSTVWY')
    model = ConcreteModel()
    le = len(ei) + len(ej) + k
    neg_inf = -float("inf")
    ep = {}

    def sequence_set(model, i):
        if i < len(ei):
            return [ei[i]]
        elif i < len(ei)+k:
            return alphabet
        else:
            return [ej[i-len(ei)-k]]

    model.A = Set(initialize=list(allele_prob.keys()))
    model.C = Set(initialize=list(range(cn)))
    model.EN = Set(initialize=list(range(en)))
    model.L = Set(initialize=list(range(le)))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i, a) for i in range(le) for a in model.S[i]])
    model.R = Set(initialize=list(range(le-(en-1))))

    # param
    model.f = Param(model.C, model.Sigma, initialize=lambda model, i, a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1, "con", m), 0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1, {}).get("con", 0))
    # epitope part
    model.i = Param(model.EN, model.Sigma, model.A, initialize=lambda model, i, a, m: epi_pssms[i, a, m])
    model.tau_epi = Param(initialize=10**6, mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, 0))

    # Variables
    model.x = Var(model.AUX, domain=Binary)
    model.y = Var(model.R, model.A, domain=NonNegativeReals)

    # objective linear
    model.obj_cleav = Objective(
        rule=lambda model: 0.5 * (
            sum(model.f[i, a]*model.x[model.ci+i, a] for i in model.C for a in model.S[model.ci+i])
            + sum(model.f[j, a]*model.x[model.cj+j, a] for j in model.C for a in model.S[model.cj+j])
            + 2 * model.bc
        ),
        sense=maximize)

    model.obj_epi = Objective(rule=lambda model: sum(
        model.y[i, a]*model.p[a] for a in model.A for i in model.R
    ), sense=minimize)

    model.obj_non_cleav = Objective(
        rule=lambda model: sum(
            model.f[j, a]*model.x[j+i, a] for i in range(le-(cn-1))
            for j in model.C for a in model.S[i+j]
            if i != model.ci and i != model.cj
        ) + (le - (cn - 1) - 2) * model.bc,
        sense=minimize
    )

    # constraints
    model.cons = Constraint(model.L, rule=lambda model, i: 1 == sum(
        model.x[i, a] for a in model.S[i]
    ))

    model.max_imm_c = Constraint(
        model.R, model.A, rule=lambda model, i, m: model.y[i, m] >= sum(
            model.x[i+j, a]*model.i[j, a, m]
            for j in model.EN for a in model.S[i+j]
        ) + model.bi[m] - model.t_a[m]
    )

    # neo-epitope constraint
    model.c_epi = Constraint(rule=lambda model: model.tau_epi >= sum(
        model.y[i, a]*model.p[a] for a in model.A for i in model.R
    ))

    # cleavage constraint
    model.c_cleavage = Constraint(rule=lambda model: 0.5 * model.tau_cleav <= (
        sum(model.f[i, a]*model.x[model.ci+i, a] for i in model.C for a in model.S[model.ci+i])
        + sum(model.f[j, a]*model.x[model.cj+j, a] for j in model.C for a in model.S[model.cj+j])
        + 2 * model.bc
    ))

    instance = model
    solver = SolverFactory(solver)

    instance.obj_epi.deactivate()
    instance.obj_non_cleav.deactivate()
    instance.c_epi.deactivate()
    instance.c_cleavage.deactivate()

    res = solver.solve(instance, options=options)  # , tee=True)

    if (
            res.solver.status == SolverStatus.ok
            and res.solver.termination_condition == TerminationCondition.optimal
    ):
        instance.solutions.load_from(res)
        obj_cleav = instance.obj_cleav()

        instance.obj_cleav.deactivate()
        instance.obj_epi.activate()
        instance.c_cleavage.activate()

        # set bound of now inactive objective
        getattr(instance, "tau_cleav").set_value(alpha*obj_cleav)
        res2 = solver.solve(instance, options=options)  # , tee=True)
        if (
                res2.solver.status == SolverStatus.ok
                and res2.solver.termination_condition == TerminationCondition.optimal
        ):
            instance.solutions.load_from(res2)

            if beta:

                # print "In thrid objective"
                obj_imm = instance.obj_epi()

                instance.obj_epi.deactivate()
                instance.obj_non_cleav.activate()
                instance.c_epi.activate()

                getattr(instance, "tau_epi").set_value((2-beta)*obj_imm)

                res3 = solver.solve(instance, options=options)  # , tee=True)
                if (
                        res3.solver.status == SolverStatus.ok
                        and res3.solver.termination_condition == TerminationCondition.optimal
                ):
                    instance.solutions.load_from(res3)
                    ci = float(sum(
                        cl_pssm[i][a] * instance.x[model.ci + i, a].value
                        for i in instance.C for a in instance.S[instance.ci+i])
                    ) + cl_pssm.get(-1, {}).get("con", 0)
                    cj = float(sum(
                        cl_pssm[j][a] * instance.x[model.cj + j, a].value
                        for j in instance.C for a in instance.S[instance.cj + j]
                    )) + cl_pssm.get(-1, {}).get("con", 0)
                    imm = float(sum(
                        instance.y[i, a].value * instance.p[a]
                        for a in instance.A for i in instance.R
                    ))
                    non_c = float(sum(
                        cl_pssm[j][a]*instance.x[j+i, a].value
                        for i in range(le-(cn-1)) for j in instance.C for a in instance.S[i+j]
                        if i != instance.ci and i != instance.cj
                    ))

                    return "".join(
                        a for i in range(len(ei), len(ei) + k) for a in instance.S[i]
                        if instance.x[i, a].value
                    ), float(ci + cj)/2, imm, float(ci), float(cj), non_c
                else:
                    raise RuntimeError("Problem could not be solved. Please check your input.")
            else:
                ci = float(sum(
                    cl_pssm[i][a] * instance.x[model.ci + i, a].value
                    for i in instance.C for a in instance.S[instance.ci + i]
                )) + cl_pssm.get(-1, {}).get("con", 0)
                cj = float(sum(
                    cl_pssm[j][a] * instance.x[model.cj + j, a].value
                    for j in instance.C for a in instance.S[instance.cj + j]
                )) + cl_pssm.get(-1, {}).get("con", 0)
                non_c = float(sum(
                    cl_pssm[j][a] * instance.x[j + i, a].value
                    for i in range(le - (cn - 1)) for j in instance.C for a in instance.S[i+j]
                    if i != instance.ci and i != instance.cj
                ))

                return "".join(
                    a for i in range(len(ei), len(ei) + k) for a in instance.S[i]
                    if instance.x[i, a].value
                ), float(ci + cj)/2, instance.obj_epi(), float(ci), float(cj), non_c
        else:
            raise RuntimeError("Problem could not be solved. Please check your input.")
    else:
        raise RuntimeError("Problem could not be solved. Please check your input.")
