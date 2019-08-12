
from __future__ import division, print_function
import multiprocessing as mp
import csv
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
from Fred2.EpitopeSelection import OptiTope, PopCover
from Fred2.IO import FileReader
from Fred2.Utility import generate_overlap_graph

from mosaic_vaccine_ilp import (DataContainer, EvaluationResult,
                                MosaicVaccineILP)
from team_orienteering_ilp import TeamOrienteeringIlp

from pyomo.core.expr.numeric_expr import SumExpression
from pyomo.opt import SolverFactory, TerminationCondition
import pyomo.environ as aml
import pyomo.kernel as pmo


LOGGER = None


@click.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-cleavages', type=click.Path())
@click.argument('output-frontier', type=click.Path())
@click.option('--pareto-steps', '-s', default=10, help='How many points to discover on the pareto frontier')
@click.option('--cocktail', '-c', default=1, help='How many strains to include in the vaccine cocktail')
@click.option('--max-aminoacids', '-a', default=0, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=10, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', default=0.0, help='Vaccine must cover at least this many alleles')
@click.option('--min-proteins', default=0.0, help='Vaccine must cover at least this many proteins')
def main(input_epitopes, input_cleavages, output_frontier, cocktail, pareto_steps,
         max_aminoacids, max_epitopes, min_alleles, min_proteins):
    ''' Explore the trade-off between the immunogenicity and the cleavage likelihood for string-of-beads vaccines.
        Outputs the Pareto frontier of the two objectives.
    '''

    global LOGGER
    LOGGER = utilities.init_logging(verbose=False)

    # load epitopes
    epitopes = utilities.load_epitopes(input_epitopes)
    LOGGER.info('Loaded %d epitopes', len(epitopes))

    # read cleavage scores
    cleavage_epitopes = set()
    with open(input_cleavages) as f:
        cleavages = {}
        for row in csv.DictReader(f):
            cleavages[(row['from'], row['to'])] = float(row['score'])
            cleavage_epitopes.add(row['from'])
            cleavage_epitopes.add(row['to'])
    LOGGER.info('Loaded %d cleavage scores', len(cleavages))

    # compute edge cost
    edge_cost, vertices, vertices_rewards = [], [], []
    vertex_to_epitope = [''] + list(cleavage_epitopes)
    for ep_from in vertex_to_epitope:
        vertices.append(ep_from)
        vertices_rewards.append(0 if ep_from == '' else epitopes[ep_from]['immunogen'])
        edge_cost.append([
            cleavages[(ep_from, ep_to)] if ep_from != '' and ep_to != '' else 0.0
            for ep_to in vertex_to_epitope
        ])
    
    type_coverage, min_type_coverage = utilities.compute_coverage_matrix(epitopes, min_alleles, min_proteins)

    # find optimal design
    solver_build_time = time.time()
    solver = TeamOrienteeringIlp(
        num_teams=cocktail, vertex_reward=vertices_rewards, edge_cost=edge_cost,
        type_coverage=type_coverage, min_type_coverage=min_type_coverage,
        max_edge_cost=0, max_vertices=max_epitopes, 
    )
    solver.build_model()
    reward_cost = solver.explore_edge_cost_vertex_reward_tradeoff(pareto_steps)
    with open(output_frontier, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('immunogenicity', 'cleavage'))
        for ig, cl in reward_cost:
            writer.writerow((ig, cl))
            f.flush()  # write progress immediately
            LOGGER.info('Immunogenicity: %.3f - Cleavage: %.3f', ig, cl)


if __name__ == '__main__':
    main()