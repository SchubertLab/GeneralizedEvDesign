import numpy as np
import pandas as pd
from Fred2.Core import Allele, Peptide
from Fred2.EpitopePrediction import EpitopePredictionResult
import csv


def get_alleles_and_thresholds(allele_file):
    df = pd.read_csv(allele_file, index_col=['allele'])
    return df


def bindings_from_csv(bindings_file):
    df = pd.read_csv(bindings_file, index_col=['Seq', 'Method'])
    df['Seq'] = df.Seq.apply(Peptide)
    df = df.set_index(['Seq', 'Method'])
    df.columns = [Allele(c) for c in df.column]
    bindings = EpitopePredictionResult(df)


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