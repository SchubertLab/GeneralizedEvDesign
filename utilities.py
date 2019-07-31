import numpy as np
import logging
import pandas as pd
from Fred2.Core import Allele, Peptide
from Fred2.EpitopePrediction import EpitopePredictionResult
import csv


def get_alleles_and_thresholds(allele_file):
    df = pd.read_csv(allele_file, index_col=['allele'])
    return df


def bindings_from_csv(bindings_file, allele_data):
    df = pd.read_csv(bindings_file)
    df['Seq'] = df.Seq.apply(Peptide)
    df = df.set_index(['Seq', 'Method'])
    df.columns = [Allele(c, allele_data[c]['frequency'] / 100) for c in df.columns]
    return EpitopePredictionResult(df)


def init_logging(verbose):
    level = (logging.DEBUG) if verbose else logging.INFO
    logger = logging.getLogger()
    logger.setLevel(level)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    fh = logging.FileHandler('dev/last-run.log', 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


def compute_suffix_prefix_cost(strings):
    all_costs = np.zeros((len(strings), len(strings)))
    for i, string_from in enumerate(strings):
        for j, string_to in enumerate(strings):
            cost = None
            if j == 0 or i == j:
                cost = 0
            elif i == 0:
                cost = len(string_to)
            else:  # compute longest suffix-prefix
                k = 1
                string_to, string_from = str(string_to), str(string_from)
                while k < len(string_from) and k < len(string_to) and string_from[-k:] == string_to[:k]:
                    k += 1
                cost = len(string_to) - k + 1

            all_costs[i, j] = cost
    return all_costs