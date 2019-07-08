
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


def _calc_set_popcover(args):
    '''
    calculates the PopCover objective for a given set of peptides
    :param args: (selection, beta, p_h, p_g, R)
    :return: the PopCover objective
    '''
    print(args)
    selection, beta, p_h, p_g, R, n_allele, n_vars = args
    score = 0
    for h in xrange(n_allele):
        for g in xrange(n_vars):
            E = sum(R.get((h, g, e), 0) for e in selection)
            r = E > 0
            score += (r*p_h[h]*p_g[g])/(beta+E)
    return score


def _calc_gcb_popcover_update(args):
    '''
    calculates the gcb update scores by solving a TSP for the current selection
    and calculating the difference in score and path length with the previous solution

    :param args: tuple of
    :return:
    '''
    selection, beta, p_h, p_g, R, O, prev_obj, prev_length, n_alleles, n_vars = args
    curr_popcover = _calc_set_popcover((selection, beta, p_h, p_g, R, n_alleles, n_vars))
    d_popcover = curr_popcover - prev_obj
    order = _solve_tsp(O, warmstart=None)

    curr_length = _calculate_length(order, O)
    d_length = curr_length - prev_length
    print(selection,'previous pop', prev_obj, 'diff', d_popcover, 'prev length',prev_length, 'diff', d_length)
    return d_popcover/d_length, curr_popcover, curr_length, order, selection[-1]


def _calculate_length(seqs, overlaps):
    '''
    calculates the length of the current mosaic vaccine

    :param list(int) seqs: a list of peptides indices
    :param overlaps: an overlap graph
    :type: numpy.array
    :return: length of the mosaic vaccine
    '''
    return sum(overlaps[seqs[i], seqs[i+1]] for i in xrange(len(seqs)-1))


def _popcover_score(args):
    '''
    calculates the greedy PopCover score
    :param args: a tuple of arguments ()
    :return: the popcover score for peptide i
    '''
    i, beta, p_h, p_g, E, R, n_allele, n_var = args

    return sum((R.get((h, g, i), 0)*p_g[g]*p_h[h])/(beta+E.get((h, g), 0))
               for h in xrange(n_allele)
               for g in xrange(n_var))


class MosaicVaccineGreedy(object):
    '''
    Greedy implementation of the mosaic problem basted on the generalized cost-benefit algorithm
    PopCover's objective function is used instead of OptiType's as it models implicitly the anitgen
    and HLA allele coverage constraints in it objective function.

    PopCover's immunogenicity function is defined as:


    S_j^{H,G} = \sum_{h \in H} \sum_{g \in G} \frac{R_{h,g}^j \codt p_h \codt p_g}{\beta + E_{h,g}}

    The sum is over all genomes g and HLA alleles h. Rjki is 1 if epitope j is present in genome g and is presented
    by allele h, and Ehg is the number of times allele h has been targeted by epitopes in genome g by the already
    selected set of epitopes, p_h is the frequency of allele h in a given population and p_g is the genomes frequency

    ::Note::
    p_g is assumed to be 1 if not specified otherwise.

    '''

    def __init__(self, _results, threshold=None, beta=0.1, t_max=100, p_g=None, verbosity=0, processes=1):
        # check input data
        if not isinstance(_results, EpitopePredictionResult):
            raise ValueError('first input parameter is not of type EpitopePredictionResult')

        # start constructing model
        _alleles = copy.deepcopy(_results.columns.values.tolist())
        self.__m = mp.Manager()
        #test if allele prob is set, if not set allele prob uniform
        #if only partly set infer missing values (assuming uniformity of missing value)
        prob = []
        no_prob = []
        for a in _alleles:
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
                        a.prob = remaining_mass/total_loc_a
                else:
                    for a in v:
                        a.prob = 1.0/total_loc_a

        method = _results.index.values[0][1]
        res_df = _results.xs(_results.index.values[0][1], level='Method')
        res_df = res_df[res_df.apply(lambda x: any(x[a] > threshold.get(a.name, -float('inf'))
                                                   for a in res_df.columns), axis=1)]


        #construct structures needed for fast calculation
        variations = list(set(prot for p in res_df.index for prot in p.get_all_proteins()))
        v_to_idx = {p:i for i,p in enumerate(variations)}
        a_to_idx = {a:i for i,a in enumerate(_alleles)}
        if not variations:
            raise ValueError('Epitopes do not have protein origin annotated which is needed for algorithm.')

        peps = ['start']
        alleles_epi = {}
        epi_alleles = {}
        epi_var = {}
        var_epi = {}
        R = self.__m.dict()
        imm = self.__m.dict()
        for i,tup in enumerate(res_df.itertuples()):
            k = i+1
            p = tup[0]
            peps.append(p)
            for a, s in itr.izip(_alleles, tup[1:]):
                if method in ['smm', 'smmpmbec', 'arb', 'comblibsidney']:
                    try:
                        thr = min(1., max(0.0, 1.0 - math.log(threshold.get(a.name),
                                                              50000))) if a.name in threshold else -float('inf')
                    except:
                        thr = 0

                    a_idx = a_to_idx[a]
                    if s >= thr:

                        alleles_epi.setdefault(a_idx, set()).add(i)
                        epi_alleles.setdefault(i, set()).add(a_idx)
                        for pr in p.get_all_proteins():
                            pr_idx = v_to_idx[pr]
                            R[a_idx, pr_idx, k] = 1
                    imm[a_idx,k] = s if s >= thr else 0

                else:
                    a_idx = a_to_idx[a]
                    if s >= threshold.get(a.name, -float('inf')):
                        alleles_epi.setdefault(a_idx, set()).add(i)
                        epi_alleles.setdefault(i, set()).add(a_idx)
                        for pr in p.get_all_proteins():
                            pr_idx = v_to_idx[pr]
                            R[a_idx, pr_idx, k] = 1
                    imm[a_idx, k] = s if s >= threshold.get(a.name, -float('inf')) else 0

            for pr in p.get_all_proteins():
                pr_idx = v_to_idx[pr]
                epi_var.setdefault(i, set()).add(pr_idx)
                var_epi.setdefault(pr_idx, set()).add(i)

        del v_to_idx
        del a_to_idx

        self.__verbosity = verbosity
        self.__processes = processes
        self.__beta = beta
        self.__tMax = t_max
        self.__p_h = np.array([a.prob for a in _alleles])
        if p_g is None:
            self.__p_g = np.ones(shape=len(variations))
        else:
            raise NotImplementedError
        self.__peps = peps
        self.__alleles = _alleles
        self.__variaitons = variations
        self.__R = R
        self.__alleles_epi = alleles_epi
        self.__epi_alleles = epi_alleles
        self.__var_epi = var_epi
        self.__epi_var = epi_var
        self.__overlap = generate_overlap_graph(peps[1:])
        #self.__imm = imm

    def solve(self):
        '''
        Uses a greedy approximation algorithm to construct a mosaic vaccine


        :return:
        '''

        def __find_update_idx(i, selection):
            '''
            calculates the set of peptides who's score has to be updated

            :param int i: the index of the current selection
            :return: list of peptide indices for update
            '''
            effected_alleles = self.__epi_alleles[i]
            effected_vars = self.__epi_var[i]

            epi_all = []
            for a in effected_alleles:
                epi_all.extend(self.__alleles_epi[a])

            epi_va = []
            for v in effected_vars:
                epi_va.extend(self.__var_epi[v])

            return list(set(epi_va).intersection(set(epi_all)) - set(selection))

        pool = mp.Pool(processes=self.__processes)

        n_peps = len(self.__peps)
        n_alleles = len(self.__alleles)
        n_vars = len(self.__variaitons)
        R = self.__R
        selection_gcb = [0]
        selection_naive = [0]
        beta = self.__beta
        scoring = np.zeros(shape=n_peps)
        p_h = self.__p_h
        p_g = self.__p_g
        E = mp.Manager().dict()

        # a list of peptides that need to be updated each round
        to_update = range(1, n_peps-1, 1)
        scoring[0] = -float('inf')

        #initialize the score
        print('Initialization')
        scoring[to_update] = pool.map(_popcover_score,
                                      [(i, beta, p_h, p_g, E, R, n_alleles, n_vars) for i in
                                       to_update])

        print('First Selection')
        curr_selection = np.argmax(scoring)
        selection_naive.append(curr_selection)

        # first gcb selection is the peptide that has the largest score normalized by the average suffix prefix matches
        print('GCB first selection')
        curr_selection_gcb = np.argmax(scoring / self.__overlap.mean(axis=1))
        curr_length_gbc = len(self.__peps[curr_selection_gcb-1])
        selection_gcb.append(curr_selection_gcb)
        curr_length_naive = len(self.__peps[curr_selection-1])

        scoring[curr_selection] = -float('inf')


        # simple selection
        print('Naive Optimization')
        while curr_length_naive < self.__tMax:
            print(curr_selection, curr_length_naive, selection_naive)
            # find indices that have to be updated
            for i in self.__epi_alleles[curr_selection]:
                for j in self.__epi_var[curr_selection]:
                    E[i, j] = E.get((i, j), 0) + 1

            to_update = __find_update_idx(curr_selection, selection_naive)
            print(to_update)
            scoring[to_update] = pool.map(_popcover_score,
                                          [(i, beta, p_h, p_g, E, R, n_alleles, n_vars) for i in
                                           to_update])

            # search for the element with max score that still fits:
            best_score = -float('inf')
            best_element = None
            for i in xrange(n_peps):
                if scoring[i] >= best_score:
                    best_element = i
                    best_score = scoring[i]
            if best_element is None:
                break

            curr_selection = np.argmax(scoring)
            scoring[curr_selection] = -float('inf')
            selection_naive.append(curr_selection)
            curr_length_naive = _calculate_length(selection_naive, self.__overlap)


        # GCB part
        print('GCB optimization')
        pep_idx = set(range(n_peps))
        pep_idx.remove(0)
        pep_idx.remove(selection_naive[1])
        O = self.__overlap
        prev_obj = _calc_set_popcover((selection_gcb, beta, p_h, p_g, R, n_alleles, n_vars))
        while curr_length_gbc < self.__tMax and pep_idx:

            tasks = []
            for i in pep_idx:
                tmp_idx = selection_gcb+[i]
                tmp_idx.sort()
                task = (selection_gcb + [i], beta, p_h, p_g, R,
                        O[np.ix_(tmp_idx, tmp_idx)], prev_obj,
                        curr_length_gbc, n_alleles, n_vars)
                tasks.append(task)

            diff = pool.map(_calc_gcb_popcover_update, tasks)
            print(diff)
            (diff, imm, length, selection, i) = max(diff)
            print('Difference:', diff, 'Imm:', imm, 'length', length, 'Selected:',i)
            selection_gcb = selection
            pep_idx.remove(i)
            curr_length_gbc = length
            prev_obj = imm

        if _calc_set_popcover((selection_gcb, beta, p_h, p_g, R, n_alleles, n_vars)) >= _calc_set_popcover(
                (selection_naive, beta, p_h, p_g, R, n_alleles, n_vars)):
            return [self.__peps[i] for i in selection_gcb]
        else:
            return [self.__peps[i] for i in selection_naive]