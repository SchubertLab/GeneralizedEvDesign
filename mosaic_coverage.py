from __future__ import print_function, division
import traceback
import Queue
import json
import multiprocessing as mp

import pandas as pd
import click

from mosaic_tools import group_peptides_by_gene

try:
    import cPickle as pickle
except ImportError:
    import pickle

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
from MosaicVaccineGreedy import MosaicVaccineGreedy


LOGGER = None



class Trie:
    def __init__(self, alphabet='ACEDGFIHKMLNQPSRTWVYX'):
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
            for letter in self.children:
                child = self._get_child(letter, create=False)
                if child is None:
                    continue

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
    level = (logging.DEBUG) if verbose else logging.INFO
    LOGGER = logging.getLogger()
    LOGGER.setLevel(level)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    LOGGER.addHandler(sh)

    fh = logging.FileHandler('dev/last-run.log', 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    LOGGER.addHandler(fh)


def is_percent_barrier(i, n, p):
    ''' returns true if i is on the boundary between two 100p% blocks of n '''
    return int(100.0 / p * (i + 1) / n) > int(100 / p * i / n)


class Process(mp.Process):
    def __init__(self, in_queue, out_queue, all_peptides, max_edits, protein_count_by_peptide):
        super(Process, self).__init__()

        self.in_queue = in_queue
        self.out_queue = out_queue

        self.all_peptides = all_peptides
        self.max_edits = max_edits
        self.protein_count_by_peptide = protein_count_by_peptide
    
    def run(self):
        timeouts = 0
        while True:
            try:
                obj = self.in_queue.get(timeout=1)
                timeouts = 0
                if obj is None:
                    return
                
                try:
                    result = self.compute_reachable_count(obj)
                except:
                    traceback.print_exc()
                    continue

                while True:
                    try:
                        self.out_queue.put(result, timeout=1)
                        break
                    except Queue.Full:
                        timeouts += 1
                        if timeouts >= 30:
                            return
            except Queue.Empty:
                timeouts += 1
                if timeouts >= 30:
                    return
    
    def compute_reachable_count(self, peptide):
        protein_count_by_edits = {}
        for reachable, edits in self.all_peptides.reachable_strings(peptide, self.max_edits):
            if edits not in protein_count_by_edits:
                protein_count_by_edits[edits] = 0
            protein_count_by_edits[edits] += self.protein_count_by_peptide[reachable]

        res = {'match-%d' % (len(peptide) - e): c for e, c in protein_count_by_edits.iteritems()}
        res['peptide'] = peptide
        return '%s\n' % json.dumps(res)


@main.command()
@click.argument('input-sequences', type=click.Path())
@click.argument('output-file', type=click.Path())
@click.option('--max-edits', '-e', default=2, help='Maximum edits allowed')
@click.option('--processes', '-p', default=-2, help='How many processes to use')
def compute_coverage(input_sequences, max_edits, output_file, processes):
    ''' Computes protein coverage for each peptide, allowing for inexact matching

        In other words, it first generates all peptides that appear in the input proteins,
        and counts how many proteins each peptide appears in. Then, for every peptide, it
        finds all peptides that can be obtained by changing at most max-edits aminoacids,
        and counts the proteins that contain the edited peptides.
    '''
    LOGGER.info('Reading sequences...')
    proteins = FileReader.read_fasta(input_sequences, in_type=Protein)
    LOGGER.info('%d proteins read', len(proteins))

    LOGGER.info('Extracting peptides and counting how many proteins each peptide covers...')
    all_peptides = Trie()
    protein_count_by_peptide = {}
    for i, prot in enumerate(proteins):
        aminoacids = ''.join(c for c in prot._data if c.isalpha())  # remove non-aminoacids from alignment
        peptides_in_this_protein = set()
        for j in xrange(len(aminoacids) - 8):
            seq = str(aminoacids[j:j+9])
            if seq not in peptides_in_this_protein:
                peptides_in_this_protein.add(seq)
                all_peptides.insert(seq)
                if seq not in protein_count_by_peptide:
                    protein_count_by_peptide[seq] = 0
                protein_count_by_peptide[seq] = 1

        if is_percent_barrier(i, len(proteins), 5):
            LOGGER.debug('%d proteins analyzed (%.2f%%) and %d peptides extracted...',
                         i + 1, 100 * (i + 1) / len(proteins), len(protein_count_by_peptide))

    LOGGER.info('Computing reachability')
    if processes < 0:
        processes = mp.cpu_count() + processes - 1

    in_queue, out_queue = mp.Queue(), mp.Queue()
    processes = [
        Process(in_queue, out_queue, all_peptides, max_edits, protein_count_by_peptide)
        for _ in range(processes)
    ]
    with open(output_file, 'w') as f:
        peptides = list(protein_count_by_peptide.keys())

        # warm up queue
        cursor = 0
        for i in range(min(len(peptides), 100)):
            in_queue.put(peptides[cursor])
            cursor += 1
        
        for p in processes:
            p.start()
        
        # process results, insert a new task every time we get a new result
        # this works as long as we consume quicker than the producers
        done = 0
        while done < len(peptides):
            res = out_queue.get()
            f.write(res)
            done += 1

            if cursor < len(peptides):
                in_queue.put(peptides[cursor])
                cursor += 1
            else:
                in_queue.put(None)

            if is_percent_barrier(done, len(peptides), 5):
                LOGGER.debug('%d peptides analyzed (%.2f%%)', done + 1, 100 * (done + 1) / len(protein_count_by_peptide))
    
    for p in processes:
        p.join()


if __name__ == '__main__':
    main()