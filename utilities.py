import numpy as np
import logging
import pandas as pd
from Fred2.Core import Allele, Peptide, Protein
from Fred2.IO import FileReader
from Fred2.EpitopePrediction import EpitopePredictionResult
import csv


def compute_coverage_matrix(epitope_data, min_alleles, min_proteins):
    # compute allele coverage matrix
    type_coverage, min_type_coverage = [], []
    if min_alleles > 0:
        alleles = [''] + list(set(a for e in epitope_data for a in e['alleles']))
        type_coverage.append(np.array([[0] * len(alleles)] + [
            [int(a in e['alleles']) for a in alleles]
            for e in epitope_data
        ]))
        count = int(min_alleles) if min_alleles > 1 else int(min_alleles * len(alleles))
        min_type_coverage.append(count)
        LOGGER.info('Vaccine will cover at least %d alleles', count)

    # compute protein coverage matrix
    if min_proteins:
        proteins = [''] + list(set(p for e in epitope_data for p in e['proteins']))
        type_coverage.append(np.array([[0] * len(proteins)] + [
            [int(p in e['proteins']) for p in proteins]
            for e in epitope_data
        ]))
        count = int(min_proteins) if min_proteins > 1 else int(min_proteins * len(proteins))
        min_type_coverage.append(count)
        LOGGER.info('Vaccine will cover at least %d proteins', count)

    # must pad all matrices to the same size
    if len(type_coverage) > 1:
        max_rows = max(a.shape[0] for a in type_coverage)
        max_cols = max(a.shape[1] for a in type_coverage)

        type_coverage = [
            np.pad(arr, ((0, max_rows - arr.shape[0]), (0, max_cols - arr.shape[1])), 'constant', constant_values=0)
            for arr in type_coverage
        ]

    return type_coverage, min_type_coverage


def load_epitopes(epitopes_file, top_immunogen=None, top_alleles=None, top_proteins=None):
    ''' loads the epitopes from the given file, returning a dictionary mapping the epitope string to its data
        optionally filters the epitopes by only taking the top N with the highest immunogenicity,
        or with the largest allele/protein coverage. if multiple options are given, the union of the
        matching epitopes is returned.
    '''
    with open(epitopes_file) as f:
        epitope_data = {}
        for row in csv.DictReader(f):
            row['immunogen'] = float(row['immunogen'])
            row['proteins'] = set(row['proteins'].split(';'))
            row['alleles'] = set(row['alleles'].split(';'))
            epitope_data[row['epitope']] = row

    if top_immunogen is None and top_alleles is None and top_proteins is None:
        return epitope_data

    def filter_epitopes(epitopes, top_count, top_key):
        assert top_count > 0
        count = int(top_count) if top_count > 1 else int(top_count * len(epitopes))
        best = sorted(epitopes, key=lambda e: top_key(epitopes[e]), reverse=True)
        return set(best[:count])

    top_epitopes = set()
    if top_immunogen:
        top_epitopes.update(filter_epitopes(epitope_data, top_immunogen, lambda e: e['immunogen']))
    if top_alleles:
        top_epitopes.update(filter_epitopes(epitope_data, top_alleles, lambda e: len(e['alleles'])))
    if top_proteins:
        top_epitopes.update(filter_epitopes(epitope_data, top_proteins, lambda e: len(e['proteins'])))
    
    return {e: epitope_data[e] for e in top_epitopes}


def get_alleles_and_thresholds(allele_file):
    df = pd.read_csv(allele_file, index_col=['allele'])
    return df


def read_annotated_proteins(proteins_file):
    ''' Reads proteins from a fasta file and extracts their metadata from the header.
        Currently follows the format of the HIV database
    '''
    proteins = FileReader.read_fasta(proteins_file, in_type=Protein)
    for prot in proteins:
        parts = prot.transcript_id.split('.')
        prot.transcript_id = parts[-1]
    return proteins


def affinities_from_csv(bindings_file, allele_data=None, peptide_coverage=None, proteins=None):
    ''' Loads binding affinities from a csv file. Optionally, augments alleles with probability
        and peptides with protein coverage.
    '''
    df = pd.read_csv(bindings_file)

    df['Seq'] = df.Seq.apply(Peptide)
    if peptide_coverage is not None:
        for pep in df.Seq:
            for prot in peptide_coverage[str(pep)]:
                pep.proteins[prot] = prot

    df = df.set_index(['Seq', 'Method'])

    if allele_data is not None:
        df.columns = [Allele(c, allele_data[c]['frequency'] / 100) for c in df.columns]
    else:
        df.columns = [Allele(c) for c in df.columns]

    return EpitopePredictionResult(df)


def init_logging(verbose):
    level = (logging.DEBUG) if verbose else logging.INFO

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    sh = logging.StreamHandler()
    sh.setLevel(level)
    sh.setFormatter(fmt)

    fh = logging.FileHandler('dev/last-run.log', 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    logger.addHandler(sh)
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