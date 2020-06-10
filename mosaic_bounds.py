

import csv
import logging
import multiprocessing as mp
import os
import time
from collections import defaultdict
from random import sample as random_sample

import click
import numpy as np
import pandas as pd
import pyomo.environ as aml
import pyomo.kernel as pmo

import Fred2
import utilities
from Fred2.Core import (Allele, Peptide, Protein,
                        generate_peptides_from_proteins)
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import (EpitopePredictionResult,
                                     EpitopePredictorFactory)
from Fred2.EpitopeSelection import OptiTope  # , PopCover
from Fred2.IO import FileReader
from team_orienteering_ilp import TeamOrienteeringIlp

LOGGER = None


def optimize_coverage(solver):
    # update objective to maximize coverage
    # WARNING: only for the first type

    del solver._model.Objective
    solver._model.Objective = aml.Objective(
        rule=lambda model: sum(
            model.OptionCovered[0, o] for o in model.Options
        ) + sum(
            model.y[n, t] * model.r[n] for n in model.Nodes for t in model.Teams
        ),
        sense=aml.maximize
    )
    solver._solver.set_objective(solver._model.Objective)


def optimize_conservation(solver):
    # update objective to maximize coverage
    # WARNING: only for the first type
    del solver._model.Objective
    solver._model.Objective = aml.Objective(
        rule=lambda model: sum(
            model.y[n, t] * (sum(model.TypeCoverage[0, n, o] for o in model.Options))
            for t in model.Teams
            for n in model.Nodes
        ) + sum(
            model.y[n, t] * model.r[n] for n in model.Nodes for t in model.Teams
        ),
        sense=aml.maximize
    )
    solver._solver.set_objective(solver._model.Objective)


@click.command()
@click.argument('bound-what', type=click.Choice(['immunogen', 'coverage', 'conservation']))
@click.argument('input-proteins', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-overlaps', type=click.Path())
@click.argument('output-vaccine', type=click.Path())
@click.option('--top-proteins', default=0.0, help='Only consider the top epitopes by protein coverage')
@click.option('--top-immunogen', default=0.0, help='Only consider the top epitopes by immunogenicity')
@click.option('--top-alleles', default=0.0, help='Only consider the top epitopes by allele coverage')
@click.option('--cocktail', '-c', default=1, help='How many strains to include in the vaccine cocktail')
@click.option('--max-aminoacids', '-a', default=[0], multiple=True, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=[10], multiple=True, help='Maximum length of the vaccine in epitopes')
@click.option('--greedy-subtour', '-g', is_flag=True, help='Insert MTZ subtour elimination at the beginning')
@click.option('--min-overlap', '-o', default=0, help='Minimum epitope overlap')
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages on the console')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(verbose, log_file, bound_what, max_epitopes, max_aminoacids, output_vaccine, **kwargs):
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)

    # add coverage information to the model by enforcing at least 2 proteins covered
    if bound_what != 'immunogen':
        kwargs['min_proteins'] = 2
    solver, data = utilities.get_mosaic_solver_instance(LOGGER, **kwargs)
    solver.build_model()

    # update objective
    if bound_what == 'coverage':
        optimize_coverage(solver)
        LOGGER.info('Optimizing pathogen coverage')
    elif bound_what == 'conservation':
        optimize_conservation(solver)
        LOGGER.info('Optimizing average conservation')
    else:
        LOGGER.info('Optimizing immunogenicity')

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


if __name__ == '__main__':
    main()
