from __future__ import division, print_function

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


def get_alleles_and_thresholds():
    # threshold in ic50 (lower is better)
    return {
        Allele('A*01:01',  4.498520 / 100.0): 625,
        Allele('A*02:01', 10.693200 / 100.0): 625,
        Allele('A*02:05',  0.884956 / 100.0): 625,
        Allele('A*03:01',  3.687320 / 100.0): 625,
        Allele('A*11:01',  7.522120 / 100.0): 625,
        Allele('A*24:02', 12.905600 / 100.0): 625,
        Allele('A*31:01',  2.433630 / 100.0): 625,
        Allele('A*68:01',  1.769910 / 100.0): 625,
        Allele('B*07:02',  3.613570 / 100.0): 625,
        Allele('B*08:01',  2.949850 / 100.0): 625,
        Allele('B*15:01',  2.064900 / 100.0): 625,
        Allele('B*27:02',  0.147493 / 100.0): 625,
        Allele('B*27:05',  1.106190 / 100.0): 625,
        Allele('B*35:01',  3.244840 / 100.0): 625,
        Allele('B*37:01',  0.442478 / 100.0): 625,
        Allele('B*38:01',  0.663717 / 100.0): 625,
        Allele('B*39:01',  1.769910 / 100.0): 625,
        Allele('B*40:01',  5.309730 / 100.0): 625,
        Allele('B*40:06',  0.516224 / 100.0): 625,
        Allele('B*44:03',  2.212390 / 100.0): 625,
        Allele('B*51:01',  3.244840 / 100.0): 625,
        Allele('B*51:02',  0.221239 / 100.0): 625,
        Allele('B*52:01',  0.884956 / 100.0): 625,
        Allele('B*58:01',  2.654870 / 100.0): 625,
        Allele('C*04:01',  8.259590 / 100.0): 625,
        Allele('C*06:02',  5.088500 / 100.0): 625,
        Allele('C*07:02',  9.660770 / 100.0): 625,
    }


def group_peptides_by_gene(input_file, peptides_path):
    ''' reads the input proteins and generates 9mers
        then counts pairs of (gene, peptide) and how many sequences for each gene
    '''
    proteins = FileReader.read_fasta(input_file, in_type=Protein)
    if peptides_path and os.path.exists(peptides_path):  # load from cache (assumes cache is up to date)
        LOGGER.info('Loading peptides from %s...', peptides_path)
        with open(peptides_path) as f:
            sequence_count_by_gene, sequence_count_by_gene_and_peptide = pickle.load(f)
    else:
        LOGGER.info('Generating peptides...')
        peptides = {}
        sequence_count_by_gene_and_peptide = {}
        sequence_count_by_gene = {}
        for prot in proteins:
            prot._data = ''.join(c for c in prot._data if c.isalpha())  # remove non-aminoacids from alignment

            # parse header
            parts = prot.transcript_id.split('.')
            if len(parts) > 4 and parts[0] == 'X':  # specific to hiv1_2017_aligned_sequences.fasta
                subtype, country, accession, gene = parts[1], parts[2], parts[-2], parts[-1]

                # update protein
                prot.gene_id = gene
                prot.transcript_id = accession + '.' + gene
            else:
                gene = ''

            if gene not in sequence_count_by_gene_and_peptide:
                sequence_count_by_gene_and_peptide[gene] = {}
                sequence_count_by_gene[gene] = 0
            sequence_count_by_gene[gene] += 1

            # update counts
            counted = set()
            for i in xrange(len(prot) - 8):
                seq = str(prot[i:i+9])

                if seq not in peptides:
                    peptides[seq] = Peptide(seq)

                pep = peptides[seq]
                if pep not in sequence_count_by_gene_and_peptide[gene]:
                    sequence_count_by_gene_and_peptide[gene][pep] = 0

                pep.proteins[prot.transcript_id] = prot
                pep.proteinPos[prot.transcript_id].append(i)
                if seq not in counted:  # only count peptides once per sequence
                    sequence_count_by_gene_and_peptide[gene][pep] += 1
                    counted.add(pep)

        if peptides_path:
            with open(peptides_path, 'w') as f:
                pickle.dump((sequence_count_by_gene, sequence_count_by_gene_and_peptide), f)
            LOGGER.info('All peptides saved to %s', peptides_path)
    
    return sequence_count_by_gene, sequence_count_by_gene_and_peptide


def get_peptides(input_file, min_conservation, peptides_path):
    ''' peptide conservation is computed separately by gene
    '''
    sequence_count_by_gene, count_by_gene_and_peptide = group_peptides_by_gene(input_file, peptides_path)
    LOGGER.info('Loaded %d proteins and %d peptides', sum(sequence_count_by_gene.values()),
                sum(map(len, count_by_gene_and_peptide.values())))
    
    LOGGER.info('Computing conserved peptides...')
    conserved_peptides = set()
    if 0 < min_conservation <= 1:   # proportion if between 0 and 1
        for gene, peptide_counts in count_by_gene_and_peptide.iteritems():
            for peptide, count in peptide_counts.iteritems():
                assert count <= sequence_count_by_gene[gene]
                cons = count / sequence_count_by_gene[gene]
                if cons >= min_conservation:
                    conserved_peptides.add(peptide)
    elif min_conservation > 1:      # top-n if larger than 1
        for gene, peptide_counts in count_by_gene_and_peptide.iteritems():
            conserved_peptides.update(sorted(
                peptide_counts.keys(), key=peptide_counts.get, reverse=True
            )[:int(min_conservation)])
    
    LOGGER.info('Found %d peptides above conservation threshold', len(conserved_peptides))
    if not conserved_peptides:
        raise RuntimeError('Improper conservation threshold, no peptides selected')
    return conserved_peptides


def get_binding_affinities_and_thresholds(peptides, bindings_path):
    ''' can either provide realistic binding affinities and thresholds, or
        use randomized values (for benchmarking purposes)
    '''
    LOGGER.info('Generating binding affinities...')
    allele_thresholds = get_alleles_and_thresholds()

    if not bindings_path or not os.path.exists(bindings_path):
        bindings = EpitopePredictorFactory('netmhcpan').predict(peptides, allele_thresholds.keys())
        if bindings_path:
            bindings.to_csv(bindings_path)
            LOGGER.info('Binding affinities saved to %s', bindings_path)
    else:
        seq_to_peptide = {str(pep): pep for pep in peptides}
        bindings = pd.read_csv(bindings_path)
        bindings['Seq'] = [seq_to_peptide[seq] for seq in bindings['Seq']]
        bindings = bindings.set_index(['Seq', 'Method'])
        bindings.columns = list(map(Allele, bindings.columns))
        bindings = EpitopePredictionResult(bindings)
        LOGGER.info('Binding affinities loaded from %s', bindings_path)

    return bindings, allele_thresholds


def get_solver_class(solver):
    solver_kwargs = {}

    if solver == 'gcb':
        solver_cls = MosaicVaccineGreedy
        solver_kwargs['processes'] = 1
    else:
        solver_cls = MosaicVaccineILP
        solver_kwargs['model_type'] = solver

    return solver_cls, solver_kwargs


@main.command()
@click.argument('input-sequences', type=click.Path('r'))
@click.argument('solver', type=click.Choice([
    'gcb',      # Generalized cost benefit heuristic
    'dfj',      # ILP solver with Dantzig-Fulkerson-Johnson subtour elimination constraints added lazily
    'mtz-l',    # ILP solver with Miller-Tucker-Zemlin subtour elimination constraints added lazily
    'mtz-g',    # ILP solver with all Miller-Tucker-Zemlin suboutr elimination constraints added at the beginning
]))
@click.option('--binding-affinities', '-b', type=click.Path(), help='Load (if exists) or save (if not) HLA-peptide binding affinities')
@click.option('--peptides', '-p', type=click.Path(), help='Load (if exists) or save (if not) computed peptides')
@click.option('--max-aminoacids', '-a', default=15, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=0.0, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', '-A', default=0.0, help='Minimum number of alleles to cover with the vaccine')
@click.option('--min-antigens', '-g', default=0.0, help='Minimum antigens to cover with the vaccine')
@click.option('--min-conservation', '-c', default=0.0, help='Minimum conservation of selected epitopes')
def design_vaccine(input_sequences, solver, binding_affinities, max_aminoacids, max_epitopes,
                   min_alleles, min_antigens, min_conservation, peptides):

    program_start_time = time.time()

    peptides = get_peptides(input_sequences, min_conservation, peptides)
    bindings, thresh = get_binding_affinities_and_thresholds(peptides, binding_affinities)

    solver_cls, solver_kwargs = get_solver_class(solver)
    LOGGER.debug('Using solver %s', solver_cls)

    solver_creation_time = time.time()
    solver = solver_cls(
        bindings, thresh, max_vaccine_aminoacids=max_aminoacids, max_vaccine_epitopes=max_epitopes,
        min_allele_coverage=min_alleles, min_antigen_coverage=min_antigens,
        min_epitope_conservation=None, **solver_kwargs
    )

    solver_build_time = time.time()
    solver.build_model()
    
    solver_start_time = time.time()
    result = solver.solve()
    solver_end_time = time.time()

    result.pretty_print(LOGGER.info)

    LOGGER.info('==== Stopwatch')
    LOGGER.info('          Total time : %.2f s', solver_end_time - program_start_time)
    LOGGER.info('  Inputs preparation : %.2f s', solver_creation_time - program_start_time)
    LOGGER.info('      Pre-processing : %.2f s', solver_build_time - solver_creation_time)
    LOGGER.info(' Model creation time : %.2f s', solver_start_time - solver_build_time)
    LOGGER.info('        Solving time : %.2f s', solver_end_time - solver_start_time)


@main.command()
@click.argument('input-sequences', type=click.Path('r'))
@click.argument('epitopes', nargs=-1)
@click.option('--min-conservation', '-c', default=0.0, help='Minimum conservation of selected epitopes')
@click.option('--random', '-r', default=0, help='Create a vaccine with a random number of epitopes')
@click.option('--greedy', '-g', default=0, help='Create a vaccine with the most immunogenic epitopes')
def inspect_vaccine(input_sequences, epitopes, min_conservation, random, greedy):
    if random > 0 and (epitopes or greedy):
        LOGGER.warning('Both epitopes and a vaccine generation strategy specified! I will use the given epitopes')
    elif random <= 0 and not epitopes and not greedy:
        LOGGER.error('Neither a vaccine strategy or the epitopes were specified.')
        return

    peptides = get_peptides(input_sequences, min_conservation)
    LOGGER.info('%d peptides generated', len(peptides))
    bindings, thresh = get_binding_affinities_and_thresholds(peptides, randomize=None, bindings_path=None)

    data = DataContainer(bindings, thresh)

    if not epitopes:
        if random > 0:
            LOGGER.info('Selecting %d random peptides', random)
            epitopes = list(map(str, random_sample(data.peptides, random)))
        elif greedy > 0:
            LOGGER.info('Selecting the %d most immunogenic peptides', greedy)
            epitopes = list(sorted(data.peptides[1:], key=lambda i: data.immunogenicities[data.pep_to_index[i]],
                                   reverse=True))[:greedy]

    tour = []
    last = 0
    for ep in epitopes:
        if ep not in data.pep_to_index:
            LOGGER.error('Peptide not found: %s', ep)
            return
        idx = data.pep_to_index[ep]
        tour.append((last, idx))
        last = idx
    tour.append((last, 0))

    result = EvaluationResult.build(data, tour)
    result.pretty_print(LOGGER.info)


if __name__ == '__main__':
    main()