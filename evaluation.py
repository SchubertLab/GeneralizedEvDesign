from __future__ import division, print_function
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
from Fred2.EpitopeSelection.PopCover import PopCover
from Fred2.IO import FileReader
from Fred2.Utility import generate_overlap_graph

from mosaic_vaccine_ilp import (DataContainer, EvaluationResult,
                                MosaicVaccineILP)
from team_orienteering_ilp import TeamOrienteeringIlp


LOGGER = None


@click.command()
@click.argument('input-sequences', type=click.Path())
@click.argument('input-peptides', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-vaccine', type=click.Path())
@click.argument('output-summary', type=click.Path())
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
def main(input_sequences, input_peptides, input_alleles, input_epitopes, input_vaccine, output_summary, verbose):
    global LOGGER
    LOGGER = utilities.init_logging(verbose)

    # load vaccine
    with open(input_vaccine) as f:
        vaccine = {}
        for row in csv.DictReader(f):
            if row['cocktail'] not in vaccine:
                vaccine[row['cocktail']] = {}
            vaccine[row['cocktail']][int(row['index'])] = row['epitope']

        cocktail = []
        for mosaic in vaccine.values():
            ordered = sorted(mosaic.items(), key=lambda x: x[0])
            cocktail.append([e for _, e in ordered])
    LOGGER.info('Vaccine loaded')

    # load alleles
    allele_data = utilities.get_alleles_and_thresholds(
        input_alleles).to_dict('index')
    LOGGER.info('Loaded %d alleles', len(allele_data))

    # load peptides coverage
    max_coverage = {}
    with open(input_peptides) as f:
        for row in csv.DictReader(f):
            max_coverage[row['peptide']] = row
    LOGGER.info('Loaded %d peptides with coverage', len(max_coverage))

    # load bindings
    epitopes = {}
    with open(input_epitopes) as f:
        for row in csv.DictReader(f):
            row['immunogen'] = float(row['immunogen'])
            row['alleles'] = row['alleles'].split(';')
            row['proteins'] = row['proteins'].split(';')
            if row['immunogen'] > 0:
                epitopes[row['epitope']] = row
    LOGGER.info('Loaded %d epitopes', len(epitopes))

    # load sequences
    proteins = FileReader.read_fasta(input_sequences, in_type=Protein)
    LOGGER.info('Loaded %d proteins', len(proteins))

    # print vaccine
    for i, mosaic in enumerate(cocktail):
        LOGGER.info('Mosaic #%d: %s', i + 1, ', '.join(mosaic))\

    # compute immunogenicity
    immunogen = sum(
        epitopes[epitope]['immunogen']
        for mosaic in cocktail
        for epitope in mosaic
    )
    LOGGER.info('Vaccine has immunogenicity %.3f', immunogen)

    # compute population coverage
    def compute_allele_coverage(alleles):
        prob_locus_covered = {'HLA-A': 0.0, 'HLA-B': 0.0, 'HLA-C': 0.0}
        for allele in set(alleles):
            prob_locus_covered[allele[:5]
                               ] += allele_data[allele]['frequency'] / 100.0
        coverage = 1 - reduce(lambda p, q: p * q,
                              ((1 - p)**2 for p in prob_locus_covered.values()))
        return coverage

    max_coverage = compute_allele_coverage(allele_data.keys())
    LOGGER.info(
        'The given set of alleles covers %.2f%% of the population', 100 * max_coverage)

    vaccine_coverage = compute_allele_coverage(
        allele
        for mosaic in cocktail
        for epitope in mosaic
        for allele in epitopes[epitope]['alleles']
    )
    LOGGER.info('The vaccine covers %.2f%% of the population (%.2f%% of the maximum theoretical coverage)',
                100 * vaccine_coverage, 100 * vaccine_coverage / max_coverage)

    # compute pathogen coverage
    proteins_covered = reduce(lambda s, t: s | t, (
        set(epitopes[epitope]['proteins'])
        for mosaic in cocktail
        for epitope in mosaic
    ))
    LOGGER.info('Vaccine covers %d proteins (%.2f%% of the total)',
                len(proteins_covered), 100 * len(proteins_covered) / len(proteins))
    
    # compute epitope conservation
    conservations = [
        len(epitopes[epitope]['proteins'])
        for mosaic in cocktail
        for epitope in mosaic
    ]
    epitope_conservation = sum(conservations) / len(conservations) / len(proteins)
    LOGGER.info('Average epitope conservation: %.2f%%', 100 * epitope_conservation)

    vaccine_stats = {
        'immunogen': immunogen,
        'pop_coverage': vaccine_coverage,
        'max_pop_coverage': max_coverage,
        'rel_pop_coverage': vaccine_coverage / max_coverage,
        'prot_coverage': len(proteins_covered),
        'norm_prot_coverage': len(proteins_covered) / len(proteins),
        'conservation': epitope_conservation,
    }
    with open(output_summary, 'w') as f:
        writer = csv.DictWriter(f, vaccine_stats.keys())
        writer.writeheader()
        writer.writerow(vaccine_stats)


if __name__ == '__main__':
    main()
