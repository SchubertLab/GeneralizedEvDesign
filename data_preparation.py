
import csv
import heapq
import itertools
import json
import logging
import multiprocessing as mp
import os
import random
import subprocess as sp
import tempfile
import time
import traceback
from collections import defaultdict
from queue import Queue
from random import sample as random_sample

import click
import numpy as np
import pandas as pd

import Fred2
import utilities
from Fred2.CleavagePrediction import (PCM, CleavageFragmentPredictorFactory,
                                      CleavageSitePredictorFactory)
from Fred2.Core import (Allele, Peptide, Protein,
                        generate_peptides_from_proteins)
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import (BIMAS, EpitopePredictionResult,
                                     EpitopePredictorFactory)
from Fred2.IO import FileReader
from optimal_spacer_design import OptimalSpacerDesign
from team_orienteering_ilp import TeamOrienteeringIlp

LOGGER = None


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages on the console')
@click.option('--log-file', '-l', type=click.Path(), help='Where to store the logs')
def main(verbose, log_file):
    global LOGGER
    LOGGER = utilities.init_logging(verbose, log_file, log_append=False)


@main.command()
@click.argument('input-sequences', type=click.Path())
@click.argument('count')
@click.argument('output-sequences', type=click.Path())
@click.option('--seed', '-s', help='Seed to use for the random selection')
def random_sequences(input_sequences, count, output_sequences, seed):
    with open(input_sequences) as f:
        prots = []
        for row in f:
            if row.startswith('>'):
                prots.append([row])
            else:
                prots[-1].append(row)
    LOGGER.info('Read %d sequences', len(prots))

    count = float(count)
    if count < 0:
        count = len(prots) + count
    elif count < 1:
        count = int(count * len(prots))
    else:
        count = int(count)

    random.seed(seed)
    sample = random.sample(prots, count)
    with open(output_sequences, 'w') as f:
        for prot in sample:
            f.writelines(prot)

    LOGGER.info('Randomly selected %d sequences', count)


@main.command()
@click.argument('input-sequences', type=click.Path(exists=True))
@click.argument('output-peptides', type=click.Path())
@click.option('--max-edits', '-e', default=0, help='Maximum edits allowed')
@click.option('--top-n', '-n', default=-1, help='Only keep the top N peptides by coverage')
def extract_peptides(input_sequences, max_edits, output_peptides, top_n):
    ''' Extract peptides from the given sequences and computes protein coverage for each peptide.
        Coverage can be computed allowing for inexact matching.

        In other words, it first generates all peptides that appear in the input proteins,
        and stores which proteins each peptide appears in. Then, for every peptide, it
        finds all peptides that can be obtained by changing at most max-edits aminoacids,
        and counts the proteins that contain the edited peptides.
    '''
    LOGGER.info('Reading sequences...')
    proteins = FileReader.read_fasta(input_sequences, in_type=Protein)
    LOGGER.info('%d proteins read', len(proteins))

    LOGGER.info('Extracting protein coverage for each peptide...')
    all_peptides = utilities.Trie()
    proteins_by_peptide = {}

    for i, prot in enumerate(proteins):
        aminoacids = ''.join(c for c in prot._data if c.isalpha())  # remove non-aminoacids from alignment
        peptides_in_this_protein = set()  # make sure we only count peptides once per protein
        for j in range(len(aminoacids) - 8):
            seq = str(aminoacids[j:j+9])
            if seq not in peptides_in_this_protein:
                peptides_in_this_protein.add(seq)
                all_peptides.insert(seq)
                if seq not in proteins_by_peptide:
                    proteins_by_peptide[seq] = set()
                proteins_by_peptide[seq].add(i)

        if utilities.is_percent_barrier(i, len(proteins), 5):
            LOGGER.debug('%d proteins analyzed (%.2f%%) and %d peptides extracted...',
                         i + 1, 100 * (i + 1) / len(proteins), len(proteins_by_peptide))

    LOGGER.info('Computing reachability...')
    top_peptides = []
    with open(output_peptides, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('peptide', 'proteins'))

        for i, peptide in enumerate(proteins_by_peptide):
            # find reachable peptides and which proteins they belong to
            reachable_proteins = set()
            for reachable, edits in all_peptides.reachable_strings(peptide, max_edits):
                reachable_proteins.update(proteins_by_peptide[reachable])

            # now either update the top N or save reachability to file
            if top_n > 0:
                heapq.heappush(top_peptides, (len(reachable_proteins), peptide, reachable_proteins))
                if len(top_peptides) > top_n:
                    heapq.heappop(top_peptides)
            else:
                writer.writerow((peptide, ';'.join(list(map(str, reachable_proteins)))))

            if utilities.is_percent_barrier(i, len(proteins_by_peptide), 2.5):
                LOGGER.debug('%d peptides analyzed (%.2f%%)...', i + 1, 100 * (i + 1) / len(proteins_by_peptide))

        # save the top N to file
        if top_n > 0:
            LOGGER.info('Saving top peptides to file')
            for _, peptide, proteins in top_peptides:
                writer.writerow((peptide, ','.join(list(map(str, reachable_proteins)))))


def predict_rank_with_netmhcpan(batch, alleles):
    tmp = tempfile.mktemp(prefix='netmhcpan-rank-')
    with open(tmp, 'w') as f:
        for pep in batch:
            f.write(f'{pep}\n')

    # call netmhcpan to predict
    output = sp.check_output([
        'netMHCpan', '-f', tmp, '-inptype', '1', '-a', ','.join(
            str(a).replace('*', '') for a in alleles
        )
    ])

    # write output to temp file
    with open(tmp + '-out', 'wb') as f:
        f.write(output)

    # parse netmhcpan output
    num_hr = 0
    binding = []
    header = [
        'Pos', 'HLA', 'Peptide', 'Core', 'Of', 'Gp', 'Gl', 'Ip',
        'Il', 'Icore', 'Identity', 'Score', '%Rank', 'BindLevel'
    ]
    for i, row in enumerate(output.decode().split('\n')):
        if row.startswith('---------'):
            num_hr += 1
            continue

        parts = row.split()
        if num_hr == 1:
            assert parts == header  # sanity check
        elif num_hr == 2:
            data = dict(zip(header, parts))
            binding.append((data['Peptide'], data['HLA'], 100 - float(data['%Rank'])))
        elif num_hr == 4:
            num_hr = 0

    # convert to dataframe following fred's format
    res = pd.DataFrame(
        binding, columns=['Seq', 'HLA', 'rank']
    ).pivot(
        index='Seq', columns='HLA', values='rank'
    )
    res['Method'] = 'netmhcpan-rank'
    cols = ['Method'] + res.columns[:-1].tolist()
    return res[cols]


def get_binding_affinity_process(predictor, batch, alleles):
    if predictor == 'netmhcpan-rank':
        return predict_rank_with_netmhcpan(batch, alleles)
    elif predictor == 'netmhcpan':
        from Fred2.EpitopePrediction.External import NetMHCpan_4_0
        predictor = NetMHCpan_4_0()
        return predictor.predict(batch, alleles)
    elif predictor == 'pickpocket':
        from Fred2.EpitopePrediction.External import PickPocket_1_1
        predictor = PickPocket_1_1()
        return predictor.predict(batch, alleles)
    elif predictor == 'mhcflurry':
        from Fred2.EpitopePrediction.ANN import MHCFlurryPredictor_1_4_3
        predictor = MHCFlurryPredictor_1_4_3()
        ic50 = predictor.predict(batch, alleles)
        return 1 - np.log(ic50) / np.log(50000)
    else:
        raise ValueError(f'unknown predictor: "{predictor}"')


@main.command()
@click.argument('input-alleles', type=click.Path())
@click.argument('input-peptides', type=click.Path())
@click.argument('output-affinities', type=click.Path())
@click.option('--processes', '-p', default=-1)
@click.option('--predictor', '-P', default='netmhcpan')
def compute_affinities(input_alleles, input_peptides, output_affinities, processes, predictor):
    ''' Computes the binding affinities between the given peptides and HLA alleles
    '''
    alleles = [
        Allele(a.replace('HLA-', ''))
        for a in utilities.get_alleles_and_thresholds(input_alleles).index
    ]
    LOGGER.info('Loaded %d alleles', len(alleles))

    with open(input_peptides) as f:
        reader = csv.DictReader(f)
        peptides = [(
            Peptide(r['peptide']), len(r['proteins'].split(';'))
        ) for r in reader]

    peptides.sort(key=lambda p: p[1], reverse=True)
    LOGGER.info('Loaded %d peptides', len(peptides))

    results = utilities.parallel_apply(get_binding_affinity_process, (
        (predictor.lower(), batch, alleles)
        for batch in utilities.batches((p for p, _ in peptides), bsize=256)
    ), processes)

    count = 0
    for bindings in results:
        bindings.to_csv(output_affinities, header=(count == 0), mode=('w' if count == 0 else 'a'))
        count += len(bindings)
        LOGGER.debug('Processed %d peptides (%.2f%%)...', count, 100 * count / len(peptides))


@main.command()
@click.argument('input-alleles', type=click.Path())
@click.argument('input-peptides', type=click.Path())
@click.argument('input-affinities', type=click.Path())
@click.argument('output-bindings', type=click.Path())
def extract_epitopes(input_alleles, input_peptides, input_affinities, output_bindings):
    ''' Extract epitopes, their immunogenicity and their coverage.
    '''

    # load alleles
    alleles = utilities.get_alleles_and_thresholds(input_alleles).to_dict('index')
    LOGGER.info('Loaded %d alleles', len(alleles))

    # load affinities, compute bindings and immunogenicities
    epitopes = {}
    with open(input_affinities) as f:
        for row in csv.DictReader(f):
            row.pop('Method')

            bindings, immunogen = [], 0.0
            for col, val in row.items():
                if not col.startswith('HLA'):
                    continue

                val = float(val)
                immunogen += val * alleles[col]['frequency'] / 100
                if val >= alleles[col]['threshold']:
                    bindings.append(col)

            epitopes[row['Seq']] = {
                'alleles': ';'.join(bindings),
                'immunogen': immunogen,
            }

    if not epitopes:
        LOGGER.error('No epitopes found!')
        return
    else:
        LOGGER.info('Found %d epitopes', len(epitopes))

    # load protein coverage
    coverage = {}
    all_proteins = set()
    with open(input_peptides) as f:
        for row in csv.DictReader(f):
            epitope = row.pop('peptide')
            coverage[epitope] = row
            all_proteins.update(set(row['proteins'].split(';')))
    LOGGER.info('Loaded %d proteins and %d peptides', len(all_proteins), len(coverage))

    # merge epitopes and coverage
    merged = []
    for epitope, data in epitopes.items():
        data.update(coverage[epitope])
        data['epitope'] = epitope
        merged.append(data)
    LOGGER.info('Merged coverage and affinities')

    with open(output_bindings, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=merged[0].keys())
        writer.writeheader()
        writer.writerows(merged)
    LOGGER.info('Saved %d epitopes', len(merged))


def get_cleavage_score_process(penalty, cleavage_model, window_size, epitopes):
    #predictor = CleavageSitePredictorFactory(cleavage_model)
    assert cleavage_model.lower() == 'pcm'
    from Fred2.CleavagePrediction import PCM
    predictor = PCM()

    results = []
    for ep_from, ep_to in epitopes:
        preds = predictor.predict(Peptide(ep_from + ep_to))
        score = 0.0
        join_pos = len(ep_from) - 1
        half_size = int((window_size - 1) / 2)
        for i, (_, lik) in enumerate(preds.values):
            if i - half_size <= join_pos <= i + half_size:
                weight = -1 if i == join_pos else penalty
                score += weight * lik
        results.append((ep_from, ep_to, score))
    return results


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('output-cleavages', type=click.Path())
@click.option('--top-proteins', default=0.0, help='Only consider the top epitopes by protein coverage')
@click.option('--top-immunogen', default=0.0, help='Only consider the top epitopes by immunogenicity')
@click.option('--top-alleles', default=0.0, help='Only consider the top epitopes by allele coverage')
@click.option('--penalty', '-P', default=0.1, help='How much to penalize wrong cleavages around the desired cleavage site')
@click.option('--cleavage-window', '-w', default=5, help='Size of the window to consider for wrong cleavages')
@click.option('--cleavage-model', '-c', default='PCM', help='Which model to use to predict cleavage sites')
@click.option('--processes', '-p', default=-1, help='Number of processes to use for parallel computation')
def compute_cleavages(input_epitopes, output_cleavages, cleavage_model, penalty, processes,
                      cleavage_window, top_proteins, top_immunogen, top_alleles):

    epitopes = list(utilities.load_epitopes(input_epitopes, top_immunogen, top_alleles, top_proteins).keys())
    LOGGER.info('Loaded %d epitopes', len(epitopes))

    LOGGER.info('Predicting cleavage sites of all pairs...')
    results = utilities.parallel_apply(get_cleavage_score_process, (
        (penalty, cleavage_model, cleavage_window, [(e, f) for f in epitopes])
        for e in epitopes
    ), processes)

    with open(output_cleavages, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('from', 'to', 'score'))
        for i, res in enumerate(results):
            for e, f, score in res:
                writer.writerow((e, f, score))

            if utilities.is_percent_barrier(i, len(epitopes), 1):
                LOGGER.debug('Processed %d cleavage pairs (%.2f%%)...',
                             len(epitopes) * (i + 1), 100 * (i + 1) / len(epitopes))


def compute_overlaps_process(epitope, other_epitopes):
    all_costs = []
    for other in other_epitopes:
        cost = utilities.compute_suffix_prefix_cost(epitope, other)
        if 0 < cost < 9:
            all_costs.append((epitope, other, cost))
    return all_costs


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('output-overlaps', type=click.Path())
@click.option('--processes', '-p', default=-1, help='Number of processes to use for parallel computation')
def compute_overlaps(input_epitopes, output_overlaps, processes):
    ''' Compute the all-pairs overlap cost for the epitopes
    '''
    epitopes = list(utilities.load_epitopes(input_epitopes).keys())
    LOGGER.info('Loaded %d epitopes', len(epitopes))

    LOGGER.info('Computing overlaps of all pairs...')
    results = utilities.parallel_apply(compute_overlaps_process, (
        (e, epitopes) for e in epitopes
    ), processes)

    by_cost = defaultdict(list)
    for i, res in enumerate(results):
        for epi_from, epi_to, cost in res:
            by_cost[cost].append((epi_from, epi_to))

        if utilities.is_percent_barrier(i, len(epitopes), 1):
            LOGGER.debug('Processed %d overlap pairs (%.2f%%)...',
                         len(epitopes) * (i + 1), 100 * (i + 1) / len(epitopes))

    LOGGER.info('Writing to file')
    with open(output_overlaps, 'w') as f:
        f.write('from,to,cost\n')
        for ov in sorted(by_cost.keys()):
            for epi_from, epi_to in by_cost[ov]:
                f.write('%s,%s,%d\n' % (epi_from, epi_to, ov))


@main.command()
@click.argument('input-epitopes', type=click.Path())
@click.argument('input-alleles', type=click.Path())
@click.argument('output-spacers', type=click.Path())
@click.option('--top-proteins', default=0.0, help='Only consider the top epitopes by protein coverage')
@click.option('--top-immunogen', default=0.0, help='Only consider the top epitopes by immunogenicity')
@click.option('--top-alleles', default=0.0, help='Only consider the top epitopes by allele coverage')
@click.option('--solver', default='gurobi', help='ILP solver to use')
@click.option('--pssm-cleavage', default='PCM', help='PSSM-based cleavage site predictor')
@click.option('--pssm-epitope', default='BIMAS', help='PSSM-based epitope predictor')
@click.option('--spacer-length', '-l', type=int, default=[0, 4], multiple=True, help='Spacer length(s) to consider')
@click.option('--alpha', '-a', default=0.99, help='Specifies how how much junction-cleavage score can be '
              'sacrificed  to gain lower neo-immunogenicity')
@click.option('--beta', '-b', default=0.0, help='Specifies how how much noe-immunogenicity score can be '
              'sacrificed to gain lower non-junction cleavage score')
@click.option('--processes', '-p', default=-1, help='Number of processes to use for parallel computation')
def design_spacers(input_epitopes, input_alleles, top_proteins, top_immunogen,
                   top_alleles, solver, pssm_cleavage, alpha, beta,
                   spacer_length, pssm_epitope, processes, output_spacers):

    all_epitopes = list(utilities.load_epitopes(input_epitopes, top_immunogen, top_alleles, top_proteins).keys())
    epitopes = [e for e in all_epitopes if 'X' not in e]
    LOGGER.debug('Removed %d epitopes with unknown amino acids', len(all_epitopes) - len(epitopes))
    LOGGER.info('Loaded %d epitopes', len(epitopes))

    alleles_df = utilities.get_alleles_and_thresholds(input_alleles)
    allele_list = [
        Allele(a.replace('HLA-', ''), prob=row.frequency / 100)
        for a, row in alleles_df.iterrows()
    ]
    threshold = {a.replace('HLA-', ''): row.threshold for a, row in alleles_df.iterrows()}
    LOGGER.info('Loaded %d alleles', len(allele_list))

    if pssm_cleavage != 'PCM':
        raise ValueError('Only PCM supported as cleavage predictor')
    cleavage_predictor = PCM()  # TODO use factory when it works
    if pssm_epitope != 'BIMAS':
        raise ValueError('Only BIMAS supported as epitope predictor')
    epitope_predictor = BIMAS()  # TODO use factory when it works

    designer = OptimalSpacerDesign(
        epitopes, cleavage_predictor, epitope_predictor,
        allele_list, threshold=threshold, solver=solver,
        k=spacer_length, alpha=alpha, beta=beta,
    ).solve(threads=processes)

    LOGGER.info('Writing results...')
    with open(output_spacers, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('from', 'to', 'score', 'spacer'))
        writer.writerows(
            (ei, ej, designer.adj_matrix[ei, ej], designer.spacer[ei, ej])
            for ei in epitopes for ej in epitopes if ei != ej
        )


if __name__ == '__main__':
    main()
