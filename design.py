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
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages on the console')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(verbose, log_file):
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-overlaps', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--top-proteins', help='Only consider the top epitopes by protein coverage', type=float)
@click.option('--top-immunogen', help='Only consider the top epitopes by immunogenicity', type=float)
@click.option('--top-alleles', help='Only consider the top epitopes by allele coverage', type=float)
@click.option('--cocktail', '-c', default=1, help='How many strains to include in the vaccine cocktail')
@click.option('--max-aminoacids', '-a', default=[0], multiple=True, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=[10], multiple=True, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', default=0.0, help='Vaccine must cover at least this many alleles')
@click.option('--min-proteins', default=0.0, help='Vaccine must cover at least this many proteins')
@click.option('--greedy-subtour', '-g', is_flag=True, help='Insert MTZ subtour elimination at the beginning')
@click.option('--min-overlap', '-o', default=0, help='Minimum epitope overlap')
def mosaic(input_epitopes, input_overlaps, output_vaccine, cocktail, max_aminoacids, max_epitopes, greedy_subtour,
           top_proteins, top_immunogen, top_alleles, min_alleles, min_proteins, min_overlap):

    # load epitopes
    epitope_data = utilities.load_epitopes(input_epitopes, top_immunogen, top_alleles, top_proteins).values()
    LOGGER.info('Loaded %d epitopes', len(epitope_data))

    # load edge cost
    vertex_rewards = [0] + [b['immunogen'] for b in epitope_data]
    epitopes = [''] + [b['epitope'] for b in epitope_data]
    epitope_index = {e: i for i, e in enumerate(epitopes)}

    edges = {}
    for i, e in enumerate(epitope_data):
        edges[(0, i + 1)] = len(e['epitope'])
        edges[(i + 1, 0)] = 0

    with open(input_overlaps) as f:
        for row in csv.DictReader(f):
            i, j = epitope_index[row['from']], epitope_index[row['to']]
            cost = float(row['cost'])
            if i != j and cost <= 9 - min_overlap:
                edges[(i, j)] = cost

    # the overlap file does not contain pairs that do not overlap, so we have to add them manually if needed
    if min_overlap <= 0:
        for i in range(1, len(epitopes)):
            for j in range(1, len(epitopes)):
                if i != j and (i, j) not in edges:
                    edges[(i, j)] = 9

    LOGGER.info('Kept %d edges (from %d)', len(edges), len(epitopes) * (len(epitopes) - 1))

    # compute hla and protein coverage
    type_coverage, min_type_coverage = utilities.compute_coverage_matrix(epitope_data, min_alleles, min_proteins)

    # find optimal design
    solver = TeamOrienteeringIlp(
        num_teams=cocktail, vertex_reward=vertex_rewards, edge_cost=edges,
        max_edge_cost=0, max_vertices=0, lazy_subtour_elimination=not greedy_subtour,
        type_coverage=type_coverage, min_type_coverage=min_type_coverage
    )
    solver.build_model()

    # sort ascending so that the previous solution is still feasible
    for ep in sorted(map(int, max_epitopes)):
        for am in sorted(map(int, max_aminoacids)):
            LOGGER.info('Solving with at most %d epitopes and most %d aminoacids', ep, am)
            solver.update_max_vertices(ep)
            solver.update_max_edge_cost(am)

            result = solver.solve()

            # print info and save
            fname = output_vaccine.format(a=am, e=ep)
            LOGGER.info('Saving result to %s', fname)
            total_ig = 0.0
            with open(fname, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(('cocktail', 'index', 'epitope'))
                for i, mosaic in enumerate(result):
                    LOGGER.info('Mosaic #%d', i + 1)
                    for j, (_, vertex) in enumerate(mosaic[:-1]):
                        writer.writerow((i, j, epitope_data[vertex - 1]['epitope']))
                        total_ig += epitope_data[vertex - 1]['immunogen']
                        LOGGER.info(
                            '    %s - IG: %.2f',
                            epitope_data[vertex - 1]['epitope'],
                            epitope_data[vertex - 1]['immunogen']
                        )
                LOGGER.info('Total immunogenicity: %.3f', total_ig)


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-cleavages', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--cocktail', '-c', default=1, help='How many strains to include in the vaccine cocktail')
@click.option('--max-aminoacids', '-a', default=0, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=10, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', default=0.0, help='Vaccine must cover at least this many alleles')
@click.option('--min-proteins', default=0.0, help='Vaccine must cover at least this many proteins')
def string_of_beads(input_epitopes, input_cleavages, output_vaccine, cocktail,
                    max_aminoacids, max_epitopes, min_alleles, min_proteins):
    program_start_time = time.time()

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
        max_edge_cost=max_aminoacids, max_vertices=max_epitopes,
    )
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
                epitope = epitopes[vertex_to_epitope[vertex]]
                writer.writerow((i, j, epitope['epitope']))
                LOGGER.info('    %s - IG: %.2f', epitope['epitope'], epitope['immunogen'])

    LOGGER.info('==== Stopwatch')
    LOGGER.info('          Total time : %.2f s', solver_end_time - program_start_time)
    LOGGER.info('      Pre-processing : %.2f s', solver_build_time - program_start_time)
    LOGGER.info(' Model creation time : %.2f s', solver_start_time - solver_build_time)
    LOGGER.info('        Solving time : %.2f s', solver_end_time - solver_start_time)


@main.command()
@click.argument('input-peptides', type=click.Path())
@click.argument('input-affinities', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--epitopes', '-e', default=10, help='Number of epitopes to include in the vaccine')
@click.option('--min-alleles', default=0.0, help='Vaccine must cover at least this many alleles')
@click.option('--min-proteins', default=0.0, help='Vaccine must cover at least this many proteins')
def optitope(input_affinities, input_peptides, input_alleles, output_vaccine, epitopes, min_alleles, min_proteins):
    with open(input_peptides) as f:
        reader = csv.DictReader(f)
        peptides = {
            # we don't really need the actual protein sequence, just fill it with the id to make it unique
            Peptide(r['peptide']): set(Protein(gid, gene_id=gid) for gid in r['proteins'].split(';'))
            for r in reader
        }
    LOGGER.info('Loaded %d peptides', len(peptides))

    allele_data = utilities.get_alleles_and_thresholds(input_alleles).to_dict('index')
    thresholds = {allele.replace('HLA-', ''): data['threshold'] for allele, data in allele_data.iteritems()}
    LOGGER.info('Loaded %d alleles', len(thresholds))

    affinities = utilities.affinities_from_csv(input_affinities, allele_data, peptide_coverage=peptides)
    LOGGER.info('Loaded %d affinities', len(affinities))

    LOGGER.info("Creating vaccine...")
    model = OptiTope(affinities, thresholds, k=epitopes, solver='gurobi')
    if min_alleles > 0:
        model.activate_allele_coverage_const(min_alleles)
        LOGGER.info('Vaccine will cover at least %f alleles', min_alleles)
    if min_proteins > 0:
        model.activate_antigen_coverage_const(min_proteins)
        LOGGER.info('Vaccine will cover at least %f proteins', min_proteins)
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