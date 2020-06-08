from __future__ import division, print_function
import multiprocessing as mp
import csv
import utilities

import logging
import os
import time
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
from Fred2.EpitopeSelection import OptiTope #, PopCover
from Fred2.IO import FileReader

from team_orienteering_ilp import TeamOrienteeringIlp


LOGGER = None


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages on the console')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(verbose, log_file):
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)


@main.command()
@click.argument('input-proteins', type=click.Path())
@click.argument('input-alleles', type=click.Path())
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
@click.option('--min-avg-prot-conservation', default=0.0, help='On average, epitopes in the vaccine must cover at least this many proteins')
@click.option('--min-avg-alle-conservation', default=0.0, help='On average, epitopes in the vaccine must cover at least this many alleles')
@click.option('--greedy-subtour', '-g', is_flag=True, help='Insert MTZ subtour elimination at the beginning')
@click.option('--min-overlap', '-o', default=0, help='Minimum epitope overlap')
def mosaic(max_epitopes, max_aminoacids, output_vaccine, **kwargs):
    # get model instance
    solver, data = utilities.get_mosaic_solver_instance(LOGGER, **kwargs)
    solver.build_model()

    # preserve user sorting
    # initializing a larger problem with the optimal (and still feasible) solution of a smaller one
    # usually makes it much faster to solve. however, too small problems could be infeasible, in which
    # case it usually takes a very long time to prove them so
    for ep in map(int, max_epitopes):
        for am in map(int, max_aminoacids):
            LOGGER.info('Solving with at most %d epitopes and most %d aminoacids', ep, am)
            solver.update_max_vertices(ep)
            solver.update_max_edge_cost(am)

            try:
                result = solver.solve()
            except RuntimeError:
                LOGGER.info('Problem was found infeasible with %d epitopes and %d aminoacids', ep, am)
                continue

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
                        writer.writerow((i, j, data['epitope_data'][vertex - 1]['epitope']))
                        total_ig += data['epitope_data'][vertex - 1]['immunogen']
                        LOGGER.info(
                            '    %s - IG: %.2f',
                            data['epitope_data'][vertex - 1]['epitope'],
                            data['epitope_data'][vertex - 1]['immunogen']
                        )
                LOGGER.info('Total immunogenicity: %.3f', total_ig)


@main.command()
@click.argument('input-proteins', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-cleavages', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--cocktail', '-c', default=1, help='How many strains to include in the vaccine cocktail')
@click.option('--max-aminoacids', '-a', default=0, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=10, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', default=0.0, help='Vaccine must cover at least this many alleles')
@click.option('--min-proteins', default=0.0, help='Vaccine must cover at least this many proteins')
@click.option('--min-avg-prot-conservation', default=0.0, help='On average, epitopes in the vaccine must cover at least this many proteins')
@click.option('--min-avg-alle-conservation', default=0.0, help='On average, epitopes in the vaccine must cover at least this many alleles')
@click.option('--greedy-subtour', '-g', is_flag=True, help='Insert MTZ subtour elimination at the beginning')
def string_of_beads(input_proteins, input_alleles, input_epitopes, input_cleavages, output_vaccine,
                    cocktail, greedy_subtour, max_aminoacids, max_epitopes, min_alleles, min_proteins,
                    min_avg_prot_conservation, min_avg_alle_conservation):
    program_start_time = time.time()

    # load proteins
    LOGGER.info('Reading sequences...')
    proteins = FileReader.read_fasta(input_proteins, in_type=Protein)
    LOGGER.info('%d proteins read', len(proteins))

    # load alleles
    alleles = [Allele(a) for a in utilities.get_alleles_and_thresholds(input_alleles).index]
    LOGGER.info('Loaded %d alleles', len(alleles))

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
    LOGGER.info('Kept %d epitopes with available clevages', len(vertices) - 1)

    type_coverage, min_type_coverage, min_avg_type_conservation = utilities.compute_coverage_matrix([
        epitopes[e] for e in vertex_to_epitope[1:]
    ], min_alleles, min_proteins, min_avg_prot_conservation, min_avg_alle_conservation, len(proteins), len(alleles))

    # find optimal design
    solver_build_time = time.time()
    solver = TeamOrienteeringIlp(
        num_teams=cocktail, vertex_reward=vertices_rewards, edge_cost=edge_cost,
        type_coverage=type_coverage, min_type_coverage=min_type_coverage,
        min_avg_type_conservation=min_avg_type_conservation, max_edge_cost=max_aminoacids,
        max_vertices=max_epitopes, lazy_subtour_elimination=not greedy_subtour
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
            LOGGER.info('    %s - %.2f', epitope, epitope_immunog)
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
