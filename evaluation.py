from __future__ import division, print_function
import re
import glob
from collections import defaultdict
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


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages on the console')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(verbose, log_file):
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)


def compute_allele_coverage(alleles, allele_data):
    prob_locus_covered = {'HLA-A': 0.0, 'HLA-B': 0.0, 'HLA-C': 0.0}
    for allele in set(alleles):
        prob_locus_covered[allele[:5]
                            ] += allele_data[allele]['frequency'] / 100.0
    coverage = 1 - reduce(lambda p, q: p * q,
                            ((1 - p)**2 for p in prob_locus_covered.values()))
    return coverage


def evaluate_epitopes(epitopes, epitope_data, allele_data, protein_count):
    # compute immunogenicity
    immunogen = sum(epitope_data[epi]['immunogen'] for epi in epitopes)

    # compute population coverage
    max_coverage = compute_allele_coverage(allele_data.keys(), allele_data)

    vaccine_alleles = set([
        allele
        for epi in epitopes
        for allele in epitope_data[epi]['alleles']
    ])
    vaccine_coverage = compute_allele_coverage(vaccine_alleles, allele_data)

    # compute pathogen coverage
    proteins_covered = reduce(lambda s, t: s | t, (
        set(epitope_data[epi]['proteins'])
        for epi in epitopes
    ))
    
    # compute epitope conservation
    conservations = [
        len(epitope_data[epi]['proteins'])
        for epi in epitopes
    ]
    epitope_conservation = sum(conservations) / len(conservations) / protein_count

    stats = {
        'num_epitopes': len(epitopes),
        'immunogen': immunogen,
        'alleles': len(vaccine_alleles),
        'pop_coverage': vaccine_coverage,
        'max_pop_coverage': max_coverage,
        'rel_pop_coverage': vaccine_coverage / max_coverage,
        'prot_coverage': len(proteins_covered),
        'norm_prot_coverage': len(proteins_covered) / protein_count,
        'conservation': epitope_conservation,
    }
    LOGGER.info('The polypeptide has %d epitopes', stats['num_epitopes'])
    LOGGER.info('The epitopes have immunogenicity %.3f', stats['immunogen'])
    LOGGER.info('The epitopes cover %d alleles', stats['alleles'])
    LOGGER.info('The maximum population coverage is %.2f%%', 100 * stats['max_pop_coverage'])
    LOGGER.info('The epitopes cover %.2f%% of the population (%.2f%% of the maximum)',
                100 * stats['pop_coverage'], 100 * stats['rel_pop_coverage'])
    LOGGER.info('The epitopes cover %d proteins (%.2f%% of the total)', stats['prot_coverage'],
                100 * stats['norm_prot_coverage'])
    LOGGER.info('The average epitope conservation is %.2f%%', 100 * stats['conservation'])

    return stats


@main.command()
@click.argument('input-sequences', type=click.Path())
@click.argument('input-peptides', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-vaccine', type=click.Path())
@click.argument('output-summary', type=click.Path())
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
def vaccine(input_sequences, input_peptides, input_alleles, input_epitopes, input_vaccine, output_summary, verbose):
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
    peptides = {}
    with open(input_peptides) as f:
        for row in csv.DictReader(f):
            peptides[row['peptide']] = row['proteins'].split(';')
    LOGGER.info('Loaded %d peptides with coverage', len(peptides))

    # load epitopes (also fill peptides since some design methods do not use epitopes)
    epitope_data = {
        pep: {'immunogen': 0.0, 'alleles': [], 'proteins': prots}
        for pep, prots in peptides.iteritems()
    }
    with open(input_epitopes) as f:
        for row in csv.DictReader(f):
            row['immunogen'] = float(row['immunogen'])
            row['alleles'] = row['alleles'].split(';')
            row['proteins'] = row['proteins'].split(';')
            if row['immunogen'] > 0:
                epitope_data[row['epitope']] = row
    LOGGER.info('Loaded %d epitopes', len(epitope_data))

    # load sequences
    proteins = FileReader.read_fasta(input_sequences, in_type=Protein)
    LOGGER.info('Loaded %d proteins', len(proteins))

    # print stats for each mosaic
    for i, mosaic in enumerate(cocktail):
        LOGGER.info('---')
        LOGGER.info('Mosaic #%d - %d epitopes', i + 1, len(mosaic))
        for epi in mosaic:
            LOGGER.info('    %s', epi)
        evaluate_epitopes(mosaic, epitope_data, allele_data, len(proteins))

    # write csv
    LOGGER.info('---')
    vaccine_stats = evaluate_epitopes([
        epi for mosaic in cocktail for epi in mosaic
    ], epitope_data, allele_data, len(proteins))
    with open(output_summary, 'w') as f:
        writer = csv.DictWriter(f, vaccine_stats.keys())
        writer.writeheader()
        writer.writerow(vaccine_stats)


@main.command()
@click.argument('path-spec')
@click.option('--output-aggregate', '-a', help='Where to save the aggregate results')
@click.option('--path-format', help='Regex that specifies which parts of the path go to which column of the result')
@click.option('--summary-by', '-s', multiple=True, help='Log summarized evaluation metrics after grouping by these columns')
@click.option('--output-summary', '-S', help='In addition to logging, save the summary to this file')
def aggregate(path_spec, output_aggregate, path_format, summary_by, output_summary):
    ''' Aggregates several evaluations, optionally summarizing them.
    '''
    fnames = glob.glob(path_spec)
    if not fnames:
        LOGGER.error('Path specification did not match any files!')
        return -1
    
    dataframes = []
    for f in fnames:
        LOGGER.debug('Parsing %s...', f)
        if path_format:
            match = re.match(path_format, f)
            if match is None:
                LOGGER.error('File "%s" did not match the given pattern, quitting', f)
                return -2

            groups = match.groupdict()
            if not groups:
                groups = dict(zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ', match.groups()))
                if not groups:
                    LOGGER.error('No capturing groups specified in the regex')
                    return -3
        else:
            groups = {}
        
        df = pd.read_csv(f)
        df['source'] = f
        for col, val in groups.iteritems():
            df[col] = val
        
        dataframes.append(df)
    LOGGER.info('Parsed %d result files', len(dataframes))

    res_df = pd.concat(dataframes, ignore_index=True)

    for row in res_df.to_string().split('\n'):
        LOGGER.info(row)

    if output_aggregate:
        res_df.to_csv(output_aggregate, index=False)
        LOGGER.info('Saved raw results to %s', output_aggregate)

    if path_format:
        LOGGER.info('Summary of the results, grouped by %s', ', '.join(summary_by))
        if summary_by:
            summary = res_df.groupby(list(summary_by)).apply(lambda g: g.describe().T)
            # bring the columns at the outermost level to facilitate comparing the same metric among different evaluations
            summary.index = summary.index.reorder_levels([len(summary_by)] + range(len(summary_by)))
        else:
            summary = res_df.describe().T

        summary.sort_index(inplace=True)
        for row in summary.to_string().split('\n'):
            LOGGER.info(row)

        if output_summary:
            summary.to_csv(output_summary)
            LOGGER.info('Saved summary to %s', output_summary)


if __name__ == '__main__':
    main()
