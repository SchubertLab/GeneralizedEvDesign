from __future__ import division, print_function
import itertools
import random
import numpy as np
from Fred2.CleavagePrediction import CleavageSitePredictorFactory, CleavageFragmentPredictorFactory
import utilities

import heapq
import json
import logging
import multiprocessing as mp
import os
import Queue
import time
import traceback
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

import csv


LOGGER = None


class Trie:
    def __init__(self):
        self.children = {}
    
    def _get_child(self, letter, create=True):
        if create and self.children.get(letter) is None:
            self.children[letter] = Trie()
        return self.children.get(letter)
    
    def insert(self, string, pos_in_string=0):
        if pos_in_string >= len(string):
            return
        
        child = self._get_child(string[pos_in_string], create=True)
        child.insert(string, pos_in_string + 1)
    
    def reachable_strings(self, string, mistakes_allowed, pos_in_string=0, mistakes_done=0):
        ''' yields all strings in the trie that can be reached from the given strings
            by changing at most `mistakes_allowed` characters, and the number of characters changed
        '''
        if not isinstance(string, list):
            string = list(string)

        if pos_in_string >= len(string):
            yield ''.join(string), mistakes_done
            return
        
        if mistakes_allowed - mistakes_done <= 0:
            child = self._get_child(string[pos_in_string], create=False)
            if child is not None:
                reachable = child.reachable_strings(string, mistakes_allowed,
                                                    pos_in_string + 1, mistakes_done)
                for s in reachable:
                    yield s
        else:
            for letter, child in self.children.iteritems():
                if letter == string[pos_in_string]:
                    reachable = child.reachable_strings(string, mistakes_allowed,
                                                        pos_in_string + 1, mistakes_done)
                    for s in reachable:
                        yield s
                else:
                    correct = string[pos_in_string]
                    string[pos_in_string] = letter
                    reachable = child.reachable_strings(string, mistakes_allowed,
                                                        pos_in_string + 1, mistakes_done + 1)
                    for s in reachable:
                        yield s
                    string[pos_in_string] = correct


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
def main(verbose):
    global LOGGER
    LOGGER = utilities.init_logging(verbose)


def is_percent_barrier(i, n, p):
    ''' returns true if i is on the boundary between two p% blocks of n '''
    return int(100.0 / p * (i + 1) / n) > int(100 / p * i / n)


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
    all_peptides = Trie()
    proteins_by_peptide = {}
    for i, prot in enumerate(proteins):
        aminoacids = ''.join(c for c in prot._data if c.isalpha())  # remove non-aminoacids from alignment
        peptides_in_this_protein = set()  # make sure we only count peptides once per protein
        for j in xrange(len(aminoacids) - 8):
            seq = str(aminoacids[j:j+9])
            if seq not in peptides_in_this_protein:
                peptides_in_this_protein.add(seq)
                all_peptides.insert(seq)
                if seq not in proteins_by_peptide:
                    proteins_by_peptide[seq] = set()
                proteins_by_peptide[seq].add(i)

        if is_percent_barrier(i, len(proteins), 5):
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

            if is_percent_barrier(i, len(proteins_by_peptide), 2.5):
                LOGGER.debug('%d peptides analyzed (%.2f%%)...', i + 1, 100 * (i + 1) / len(proteins_by_peptide))

        # save the top N to file
        if top_n > 0:
            LOGGER.info('Saving top peptides to file')
            for _, peptide, proteins in top_peptides:
                writer.writerow((peptide, ','.join(list(map(str, reachable_proteins)))))


def get_binding_affinity_process(batch, alleles):
    return EpitopePredictorFactory('netmhcpan').predict(batch, alleles)


@main.command()
@click.argument('input-alleles', type=click.Path())
@click.argument('input-peptides', type=click.Path())
@click.argument('output-affinities', type=click.Path())
@click.option('--processes', '-p', default=-1)
def compute_affinities(input_alleles, input_peptides, output_affinities, processes):
    ''' Computes the binding affinities between the given peptides and HLA alleles
    '''
    def iterbatches(it, bsize):
        batch = []
        while True:
            try:
                batch.append(next(it))
            except StopIteration:
                break

            if len(batch) >= bsize:
                yield batch
                batch = []

        if batch:
            yield batch

    alleles = [Allele(a) for a in utilities.get_alleles_and_thresholds(input_alleles).index]
    LOGGER.info('Loaded %d alleles', len(alleles))

    with open(input_peptides) as f:
        reader = csv.DictReader(f)
        peptides = [(
            Peptide(r['peptide']), len(r['proteins'].split(';'))
        ) for r in reader]
    
    peptides.sort(key=lambda p: p[1], reverse=True)
    LOGGER.info('Loaded %d peptides', len(peptides))

    pool = mp.Pool(processes=processes if processes > 0 else (mp.cpu_count() + processes))

    try:
        tasks = []
        for batch in iterbatches((p for p, _ in peptides), bsize=2048):
            tasks.append(pool.apply_async(get_binding_affinity_process, (batch, alleles)))
        
        count = 0
        for result in tasks:
            bindings = result.get(999999)
            bindings.to_csv(output_affinities, header=(count == 0), mode=('w' if count == 0 else 'a'))
            count += len(bindings)
            LOGGER.debug('Processed %d peptides (%.2f%%)...', count, 100 * count / len(peptides))
    except:
        pool.terminate()
        pool.join()
        raise
    else:
        pool.close()


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
        warned = False
        for row in csv.DictReader(f):
            if row['Method'] != 'netmhcpan':
                if not warned:
                    LOGGER.warn('Wrong affinity prediction method used; please use netmhcpan! Some rows will be skipped')
                    warned = True
                continue
            row.pop('Method')

            bindings, immunogen = [], 0.0
            for col, val in row.iteritems():
                if not col.startswith('HLA'):
                    continue
                
                val = float(val)
                immunogen += val * alleles[col]['frequency'] / 100
                if val > alleles[col]['threshold']:
                    bindings.append(col)

            if bindings:
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
    for epitope, data in epitopes.iteritems():
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
    predictor = CleavageSitePredictorFactory(cleavage_model)

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
@click.argument('input-epitopes')
@click.argument('output-cleavages')
@click.option('--top-conservation', help='Only consider the top epitopes by protein coverage', type=float)
@click.option('--penalty', '-P', default=0.1, help='How much to penalize wrong cleavages around the desired cleavage site')
@click.option('--cleavage-window', '-w', default=5, help='Size of the window to consider for wrong cleavages')
@click.option('--cleavage-model', '-c', default='PCM', help='Which model to use to predict cleavage sites')
@click.option('--processes', '-p', default=-1, help='Number of processes to use for parallel computation')
def compute_cleavages(input_epitopes, output_cleavages, cleavage_model, penalty, processes, cleavage_window, top_conservation):
    with open(input_epitopes) as f:
        epitope_data = list(csv.DictReader(f))
    LOGGER.info('Loaded %d epitopes', len(epitope_data))

    if top_conservation > 0:
        count = int(top_conservation) if top_conservation > 1 else int(top_conservation * len(epitope_data))
        epitope_data.sort(key=lambda e: e['proteins'].count(';'), reverse=True)
        epitope_data = epitope_data[:count]
        LOGGER.info('Kept %d epitopes after filtering by conservation', count)
    epitopes = [e['epitope'] for e in epitope_data]
    
    LOGGER.info('Predicting cleavage sites of all pairs...')
    pool = mp.Pool(processes=processes if processes > 0 else (mp.cpu_count() + processes))

    try:
        tasks = []
        for e in epitopes:
            tasks.append(pool.apply_async(
                get_cleavage_score_process,
                (penalty, cleavage_model, cleavage_window, [(e, f) for f in epitopes])
            ))
        
        with open(output_cleavages, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(('from', 'to', 'score'))
            for i, result in enumerate(tasks):
                for e, f, score in result.get(99999):
                    writer.writerow((e, f, score))
                if is_percent_barrier(i, len(tasks), 1):
                    LOGGER.debug('Processed %d cleavage pairs (%.2f%%)...',
                                 len(epitopes) * (i + 1), 100 * (i + 1) / len(tasks))
    except:
        pool.terminate()
        pool.join()
        raise
    else:
        pool.close()


if __name__ == '__main__':
    main()
