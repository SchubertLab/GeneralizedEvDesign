from __future__ import division, print_function
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
    LOGGER = logging.getLogger('mosaic-tools')
    LOGGER.setLevel((logging.DEBUG) if verbose else logging.INFO)

    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    LOGGER.addHandler(sh)

    fh = logging.FileHandler('dev/last-run.log', 'w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    LOGGER.addHandler(fh)


def get_alleles_and_thresholds():
    return {     # According to Table 1 in Vider-Shalit et al. 2007
        #                  population cov.  :   threshold
        Allele('A*01:01',  4.498520 / 100.0):   0.08710, 
        Allele('A*02:01', 10.693200 / 100.0):   1.25720, 
        Allele('A*02:05',  0.884956 / 100.0):   0.41790, 
        Allele('A*03:01',  3.687320 / 100.0):   0.05527, 
        Allele('A*11:01',  7.522120 / 100.0):   0.04356, 
        Allele('A*24:02', 12.905600 / 100.0):   1.15830, 
        Allele('A*31:01',  2.433630 / 100.0):   0.09604, 
        Allele('A*68:01',  1.769910 / 100.0):   2.39910, 
        #Allele('B*04:01',  1.548670 / 100.0):   5.22380,   # not supported by netMHCpan
        Allele('B*07:02',  3.613570 / 100.0):   1.62320, 
        Allele('B*08:01',  2.949850 / 100.0):   0.07150, 
        Allele('B*15:01',  2.064900 / 100.0):   0.21666, 
        Allele('B*27:02',  0.147493 / 100.0):  50.04300, 
        Allele('B*27:05',  1.106190 / 100.0): 163.86800, 
        Allele('B*35:01',  3.244840 / 100.0):   3.20170, 
        Allele('B*37:01',  0.442478 / 100.0):   1.42390, 
        Allele('B*38:01',  0.663717 / 100.0):   7.92910, 
        Allele('B*39:01',  1.769910 / 100.0):   4.73410, 
        Allele('B*40:01',  5.309730 / 100.0):  21.17820, 
        Allele('B*40:06',  0.516224 / 100.0):   1.04370, 
        Allele('B*44:03',  2.212390 / 100.0):   6.20900, 
        Allele('B*51:01',  3.244840 / 100.0):  39.20670, 
        Allele('B*51:02',  0.221239 / 100.0): 155.05000, 
        Allele('B*52:01',  0.884956 / 100.0):   8.58850, 
        Allele('B*58:01',  2.654870 / 100.0):  13.35620, 
        Allele('C*04:01',  8.259590 / 100.0):   4.06100, 
        Allele('C*06:02',  5.088500 / 100.0):   2.16860, 
        Allele('C*07:02',  9.660770 / 100.0):   2.18310, 
    }


def get_peptides(input_file, min_conservation, peptides_path):
    if os.path.exists(peptides_path):
        with open(peptides_path) as f:
            conserved_peptides = pickle.load(f)
        LOGGER.info('Loaded %d conserved peptides from %s', len(conserved_peptides), peptides_path)
        return conserved_peptides
    
    proteins = FileReader.read_fasta(input_file, in_type=Protein)
    if len(proteins[0].transcript_id.split('.')) != 6:
        return list(generate_peptides_from_proteins(proteins, 9))

    # the following code is specific to HIV1_2017_aligned_sequences.fasta 's header format
    LOGGER.info('Computing peptide conservation...')
    peptide_counts = {}  # TODO should probably cache this stuff
    protein_count = 0
    for prot in proteins:
        # parse header
        parts = prot.transcript_id.split('.')
        subtype, country, accession, gene = parts[0], parts[1], parts[-2], parts[-1]

        # update protein
        prot._data = prot._data.replace('-', '').replace('*', '')  # remove non-aminoacids from alignment
        prot.gene_id = gene
        prot.transcript_id = accession + '.' + gene
        protein_count += 1

        # update conservations
        # TODO divide by antigen?
        for i in xrange(len(prot) - 8):
            seq = str(prot[i:i+9])
            if seq not in peptide_counts:
                pep = Peptide(seq)
                peptide_counts[pep] = 0

            pep.proteins[prot.transcript_id] = prot
            pep.proteinPos[prot.transcript_id].append(i)
            peptide_counts[pep] += 1

    LOGGER.info('Loaded %d proteins and %d peptides', protein_count, len(peptide_counts))
    
    # find sufficiently conserved peptides
    conserved_peptides = []
    if 0 < min_conservation <= 1:   # proportion if between 0 and 1
        for peptide, count in peptide_counts.iteritems():
            cons = count / protein_count
            conserved_peptides.append(peptide)
        LOGGER.info('%d peptides above conservation threshold', len(conserved_peptides))
    elif min_conservation > 1:      # top-n if larger than 1
        conserved_peptides = sorted(
            peptide_counts.keys(), key=peptide_counts.get, reverse=True
        )[:int(min_conservation)]
        LOGGER.info('%d peptides have a conservation above %.2f %%',
                    min_conservation, 100 * peptide_counts[conserved_peptides[-1]] / protein_count)

    if not conserved_peptides:
        raise RuntimeError('Improper conservation threshold, no peptides selected')
    
    if peptides_path:
        with open(peptides_path, 'w') as f:
            pickle.dump(conserved_peptides, f)
        LOGGER.info('Conserved peptides saved to %s', peptides_path)

    return conserved_peptides


def get_binding_affinities_and_thresholds(peptides, randomize, bindings_path):
    ''' can either provide realistic binding affinities and thresholds, or
        use randomized values (for benchmarking purposes)
    '''
    LOGGER.info('Generating binding affinities...')
    allele_thresholds = get_alleles_and_thresholds()
    if not bindings_path:
        bindings = EpitopePredictorFactory('netmhcpan').predict(peptides, allele_thresholds.keys())
    else:
        bindings = pd.read_csv(bindings_path)
        bindings['Seq'] = list(map(Peptide, bindings['Seq']))  # we don't have source protein
        bindings = bindings.set_index(['Seq', 'Method'])
        bindings.columns = list(map(Allele, bindings.columns))
        bindings = EpitopePredictionResult(bindings)
        LOGGER.info('Binding affinities loaded from %s', bindings_path)

    if randomize > 0:
        if randomize > 1:
            randomize = randomize / 100
        if bindings_path:
            LOGGER.warn('randomize and binding affinities both specified! The specified affinities will not')
            LOGGER.warn('be overwritten, but randominze will not work as intended if the affinities are not')
            LOGGER.warn('suitably randomized (i.e. uniformly between 0 and 1).')
            LOGGER.warn('')
            LOGGER.warn('To generate randomized affinities, run once with randomize and *without* binding affinities.')
            LOGGER.warn('Randomized affinities will be generated and stored in "./resources/bindings.csv"')
            LOGGER.warn('and can be loaded in subsequent runs.')
            LOGGER.warn('')
            LOGGER.warn('If you already did this, this warning is safe to ignore.')
            LOGGER.warn('')

        # chosen so that 100*randomize % of peptides have at least one allele
        # with larger binding strength, assuming that these are uniform(0, 1)
        tt = (1 - randomize)**(1.0 / len(allele_thresholds))
        thresh = {}
        for col in bindings.columns:
            if not bindings_path:
                bindings[col] = np.random.random(size=len(bindings))
            thresh[col] = tt
    else:
        thresh = allele_thresholds

    if randomize <= 0 and bindings_path:
        bindings.to_csv(bindings_path)
        LOGGER.info('Binding affinities saved to %s', bindings_path)

    return bindings, thresh


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
@click.option('--randomize', '-r', default=0.0, help='Randomly assign affinities and select a given portion of epitopes')
def design_vaccine(input_sequences, solver, randomize, binding_affinities, max_aminoacids, max_epitopes,
                   min_alleles, min_antigens, min_conservation, peptides):

    program_start_time = time.time()

    peptides = get_peptides(input_sequences, min_conservation, peptides)
    bindings, thresh = get_binding_affinities_and_thresholds(peptides, randomize, binding_affinities)

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
