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


@main.command()
@click.argument('input-sequences', type=click.Path())
@click.argument('output-file', type=click.Path())
@click.option('--max-edits', '-e', default=2, help='Maximum edits allowed')
def compute_coverage(input_sequences, max_edits, output_file):
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
    proteins_by_peptide = {}
    for i, prot in enumerate(proteins):
        aminoacids = ''.join(c for c in prot._data if c.isalpha())  # remove non-aminoacids from alignment
        peptides_in_this_protein = set()
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

    LOGGER.info('Computing reachability')
    with open(output_file, 'w') as f:
        for i, peptide in enumerate(proteins_by_peptide):
            proteins_by_edits = {}
            for reachable, edits in all_peptides.reachable_strings(peptide, max_edits):
                if edits not in proteins_by_edits:
                    proteins_by_edits[edits] = set()
                proteins_by_edits[edits].update(proteins_by_peptide[reachable])

            res = {'match-%d' % (len(peptide) - e): len(c) for e, c in proteins_by_edits.iteritems()}
            res['peptide'] = peptide
            f.write('%s\n' % json.dumps(res))

            if is_percent_barrier(i, len(proteins_by_peptide), 2.5):
                LOGGER.debug('%d peptides analyzed (%.2f%%)...', i + 1, 100 * (i + 1) / len(proteins_by_peptide))


if __name__ == '__main__':
    main()