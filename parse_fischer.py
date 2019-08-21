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
import pyomo.environ as aml
import pyomo.kernel as pmo


LOGGER = None


@click.command()
@click.argument('input-peptides', type=click.Path())
@click.argument('input-vaccine-sequences', type=click.Path())
@click.argument('output-vaccine-epitopes', type=click.Path())
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(input_peptides, input_vaccine_sequences, output_vaccine_epitopes, verbose, log_file):
    ''' Reads the vaccine produced by Fischer's online tool and converts it into epitopes
    '''
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)

    LOGGER.info('Reading peptides...')
    with open(input_peptides) as f:
        peptides = set(r['peptide'] for r in csv.DictReader(f))
    LOGGER.info('Read %d peptides', len(peptides))

    LOGGER.info('Reading vaccine...')
    mosaics = FileReader.read_fasta(input_vaccine_sequences, in_type=Protein)
    LOGGER.info('Vaccine has %d mosaic(s)', len(mosaics))

    with open(output_vaccine_epitopes, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('cocktail', 'index', 'epitope'))

        for c, mos in enumerate(mosaics):
            pep_count = unk_count = 0
            for i in range(0, len(mos) - 8):
                pep = mos[i:i + 9]
                assert len(pep) == 9

                if pep in peptides:
                    writer.writerow((c, pep_count, pep))
                    pep_count += 1
                else:
                    unk_count += 1

            LOGGER.info('Mosaic %d - Recognized: %d Unknown %d', c + 1, pep_count, unk_count)


if __name__ == '__main__':
    main()