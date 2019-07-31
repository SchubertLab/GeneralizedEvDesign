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


LOGGER = None


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
def main(verbose):
    global LOGGER
    LOGGER = utilities.init_logging(verbose)


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--mosaics', '-m', default=1, help='How many different mosaics to produce')
@click.option('--max-aminoacids', '-a', default=0, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=10, help='Maximum length of the vaccine in epitopes')
def mosaic(input_epitopes, output_vaccine, mosaics, max_aminoacids, max_epitopes):
    program_start_time = time.time()

    # load bindings
    with open(input_epitopes) as f:
        bindings = []
        for row in csv.DictReader(f):
            row['immunogen'] = float(row['immunogen'])
            if row['immunogen'] > 0:
                bindings.append(row)
    LOGGER.info('Loaded %d epitopes', len(bindings))

    # compute edge cost and create solver
    epitopes = [''] + [b['epitope'] for b in bindings]
    vertex_rewards = [0] + [b['immunogen'] for b in bindings]
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

    # print info and save
    with open(output_vaccine, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('cocktail', 'index', 'epitope'))
        for i, mosaic in enumerate(result):
            LOGGER.info('Mosaic #%d', i + 1)
            for j, (_, vertex) in enumerate(mosaic[:-1]):
                writer.writerow((i, j, bindings[vertex - 1]['epitope']))
                LOGGER.info('    %s - IG: %.2f', bindings[vertex - 1]['epitope'], bindings[vertex - 1]['immunogen'])

    LOGGER.info('==== Stopwatch')
    LOGGER.info('          Total time : %.2f s', solver_end_time - program_start_time)
    LOGGER.info('      Pre-processing : %.2f s', solver_build_time - program_start_time)
    LOGGER.info(' Model creation time : %.2f s', solver_start_time - solver_build_time)
    LOGGER.info('        Solving time : %.2f s', solver_end_time - solver_start_time)


@main.command()
@click.argument('input-affinities', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--epitopes', '-e', default=10, help='Number of epitopes to include in the vaccine')
def optitope(input_affinities, input_alleles, output_vaccine, epitopes):
    allele_data = utilities.get_alleles_and_thresholds(input_alleles).to_dict('index')
    thresholds = {allele.replace('HLA-', ''): data['threshold'] for allele, data in allele_data.iteritems()}
    LOGGER.info('Loaded %d alleles', len(thresholds))

    affinities = utilities.affinities_from_csv(input_affinities, allele_data)
    LOGGER.info('Loaded %d affinities', len(affinities))

    LOGGER.info("Creating vaccine...")
    model = OptiTope(affinities, thresholds, k=epitopes, solver='gurobi')
    vaccine = model.solve()

    LOGGER.info('Vaccine summary:')
    with open(output_vaccine, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('cocktail', 'index', 'epitope'))
        total_ig = 0.0
        for i, epitope in enumerate(vaccine):
            writer.writerow((0, i, epitope))
            epitope_immunog = sum(model.instance.p[a] * model.instance.i[epitope, a]
                                  for a in model.instance.A)
            total_ig += epitope_immunog
            LOGGER.info('    %s - IG: %.2f', epitope, epitope_immunog)
        LOGGER.info('Total immunogenicity: %.2f', total_ig)


@main.command()
@click.argument('input-peptides', type=click.Path())
@click.argument('input-affinities', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--epitopes', '-e', default=10, help='Number of epitopes to include in the vaccine')
@click.option('--processes', '-p', default=-1, help='Number of processes to use for parallel computation')
def popcover(input_peptides, input_affinities, input_alleles, output_vaccine, processes, epitopes):
    with open(input_peptides) as f:
        reader = csv.DictReader(f)
        peptides = {Peptide(r['peptide']): set(r['proteins'].split(';')) for r in reader}
    LOGGER.info('Loaded %d peptides', len(peptides))

    allele_data = utilities.get_alleles_and_thresholds(input_alleles).to_dict('index')
    thresholds = {allele.replace('HLA-', ''): data['threshold'] for allele, data in allele_data.iteritems()}
    LOGGER.info('Loaded %d alleles', len(thresholds))

    affinities = utilities.affinities_from_csv(input_affinities, allele_data, peptides)
    LOGGER.info('Loaded %d affinities', len(affinities))

    LOGGER.info("Creating vaccine...")
    model = PopCover(affinities, thresholds, k=epitopes,
                     processes=processes if processes > 0 else (mp.cpu_count() + processes))
    vaccine = model.solve()

    with open(output_vaccine, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('cocktail', 'index', 'epitope'))
        for i, epitope in enumerate(vaccine):
            writer.writerow((0, i, epitope))
            LOGGER.info('    %s', epitope)


if __name__ == '__main__':
    main()