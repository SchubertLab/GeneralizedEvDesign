from __future__ import division, print_function
import utilities

import logging
import os
import time
from builtins import map
from collections import defaultdict
from random import sample as random_sample

import click
import Fred2
import numpy as np
import pandas as pd
from Fred2.Core import (Allele, Peptide, Protein,
                        generate_peptides_from_proteins)
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import (EpitopePredictionResult,
                                     EpitopePredictorFactory)
from Fred2.EpitopeSelection.PopCover import PopCover
from Fred2.IO import FileReader
from Fred2.Utility import generate_overlap_graph

from mosaic_vaccine_ilp import (DataContainer, EvaluationResult,
                                MosaicVaccineILP)
from MosaicVaccineGreedy import MosaicVaccineGreedy
from team_orienteering_ilp import TeamOrienteeringIlp


LOGGER = None


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
def main(verbose):
    global LOGGER
    level = (logging.DEBUG) if verbose else logging.INFO
    LOGGER = logging.getLogger()
    LOGGER.setLevel(level)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    LOGGER.addHandler(sh)

    fh = logging.FileHandler('dev/last-run.log', 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    LOGGER.addHandler(fh)


def compute_immunogenicities(bindings):
    pass


@main.command()
@click.argument('alleles-file', type=click.Path())
@click.argument('bindings-file', type=click.Path())
@click.option('--mosaics', '-m', default=1, help='How many different mosaics to produce')
@click.option('--max-aminoacids', '-a', default=0, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=3, help='Maximum length of the vaccine in epitopes')
def mosaic(alleles_file, bindings_file, mosaics, max_aminoacids, max_epitopes):
    ''' Produces a mosaic vaccine starting from the desired epitopes and their binding strength
    '''

    program_start_time = time.time()

    # load alleles
    thresh = utilities.get_alleles_and_thresholds(alleles_file)
    LOGGER.info('Loaded %d alleles', len(thresh))

    # load bindings
    bindings = pd.read_csv(bindings_file, index_col=['Seq'])
    if list(bindings.Method.unique()) != ['netmhcpan']:
        raise ValueError('Wrong binding prediction method used. Please use netmhcpan!')
    bindings = bindings.drop('Method', axis=1)
    LOGGER.info('Loaded %d peptides', len(bindings))
    
    # compute immunogenicities
    immunogens = bindings.apply(lambda row: sum(
        (bind * thresh.T[allele].frequency / 100 if bind > thresh.T[allele].threshold else 0)
        for allele, bind in zip(row.index, row)
    ), axis=1)
    immunogens = immunogens[immunogens > 0]
    if len(immunogens) == 0:
        raise ValueError('Thresholds too high, no peptides bind!')
    LOGGER.info('%d peptides bind to at least one allele', len(immunogens))

    # compute edge cost and create solver
    epitopes = [''] + immunogens.index.to_list()
    vertex_rewards = [0] + immunogens.to_list()
    edge_cost = utilities.compute_suffix_prefix_cost(epitopes)
    solver = TeamOrienteeringIlp(
        num_teams=mosaics, vertex_reward=vertex_rewards, edge_cost=edge_cost,
        max_edge_cost=max_aminoacids, max_vertices=max_epitopes,
    )

    # find optimal design
    solver_build_time = time.time()
    solver.build_model()
    solver_start_time = time.time()
    result = solver.solve()
    solver_end_time = time.time()

    # print info
    LOGGER.info('Vaccine info:')
    for i, mosaic in enumerate(result):
        LOGGER.info('Mosaic #%d', i + 1)
        for _, vertex in mosaic[:-1]:
            LOGGER.info('    %s - IG: %.2f' % (epitopes[vertex], immunogens[vertex]))

    LOGGER.info('==== Stopwatch')
    LOGGER.info('          Total time : %.2f s', solver_end_time - program_start_time)
    LOGGER.info('      Pre-processing : %.2f s', solver_build_time - program_start_time)
    LOGGER.info(' Model creation time : %.2f s', solver_start_time - solver_build_time)
    LOGGER.info('        Solving time : %.2f s', solver_end_time - solver_start_time)


if __name__ == '__main__':
    main()