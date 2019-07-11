from __future__ import print_function
import numpy as np
from builtins import map

import time
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from MosaicVaccineILP import MosaicVaccineILP
from MosaicVaccineLazyILP import MosaicVaccineLazyILP
from MosaicVaccineGreedy import MosaicVaccineGreedy
from Fred2.Core import Protein, Allele
from Fred2.Core import generate_peptides_from_proteins
from Fred2.IO import FileReader
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopeSelection.PopCover import PopCover
from Fred2.Utility import generate_overlap_graph
import click
import Fred2


THRESHOLD_PRESETS = [  # counts based on hivgen.fasta
    {'A*01:01': 0.1,  'B*07:02': 1.0, 'C*03:01': 0.3},   # 11 epitopes
    {'A*01:01': 100,  'B*07:02': 100, 'C*03:01': 0.23},  # 528 epitopes
    {'A*01:01': 0.04, 'B*07:02': 100, 'C*03:01': 100},   # 666 epitopes
    {'A*01:01': 100,  'B*07:02': 0.8, 'C*03:01': 100},   # 817 epitopes
    {'A*01:01': 0.04, 'B*07:02': 100, 'C*03:01': 0.23},  # 1194 epitopes
    {'A*01:01': 0.04, 'B*07:02': 0.8, 'C*03:01': 100},   # 1483 epitopes
    {'A*01:01': 0.04, 'B*07:02': 0.8, 'C*03:01': 0.23},  # 1841 epitopes
    {'A*01:01': 0.01, 'B*07:02': 0.01, 'C*03:01': 0.01}, # 8463 epitopes
]


def get_binding_affinities_and_thresholds(peptides, alleles, randomize, threshold):
    bindings = EpitopePredictorFactory('BIMAS').predict(peptides, alleles)
    if randomize > 0:
        if randomize > 1:
            randomize = randomize / 100
        thresh = {}
        affinities = np.random.random(size=len(bindings))
        for col in bindings.columns:
            bindings[col] = affinities
            thresh[col.name] = 1 - randomize
    else:
        thresh = THRESHOLD_PRESETS[threshold] if threshold < len(THRESHOLD_PRESETS) else THRESHOLD_PRESETS[-1]

    bindings.to_csv('resources/bindings.csv')
    return bindings, thresh


def get_solver(solver, verbose, max_epitopes, max_aminoacids):
    solver_kwargs = {
        'verbosity': int(verbose),
        'max_vaccine_epitopes': max_epitopes,
        'max_vaccine_aminoacids': max_aminoacids,
    }

    if solver == 'gcb':
        solver_cls = MosaicVaccineGreedy
    elif solver == 'dfj':
        solver_cls = MosaicVaccineLazyILP
        solver_kwargs['subtour_elimination'] = 'dfj'
    elif solver == 'mtz-l':
        solver_cls = MosaicVaccineLazyILP
        solver_kwargs['subtour_elimination'] = 'mtz'
    elif solver == 'mtz-g':
        solver_cls = MosaicVaccineILP
    else:
        raise ValueError('unknown solver')
        
    return solver_cls, solver_kwargs


@click.command()
@click.argument('input-file', type=click.Path('r'))
@click.argument('solver', type=click.Choice([
    'gcb',      # Generalized cost benefit heuristic
    'dfj',      # ILP solver with Dantzig-Fulkerson-Johnson subtour elimination constraints added lazily
    'mtz-l',    # ILP solver with Miller-Tucker-Zemlin subtour elimination constraints added lazily
    'mtz-g',    # ILP solver with all Miller-Tucker-Zemlin suboutr elimination constraints added at the beginning
]))
@click.option('--max-aminoacids', '-a', default=15, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=999, help='Maximum length of the vaccine in epitopes')
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
@click.option('--threshold', '-t', default=0, help='Select one of the threshdld presets (0-7)')
@click.option('--measure-time', '-T', is_flag=True, help='Print the total time at the end') 
@click.option('--randomize', '-r', default=-1.0, help='Randomly assign affinities and select a given portion of epitopes')
def main(input_file, solver, max_aminoacids, max_epitopes, verbose, threshold, measure_time, randomize):
    program_start_time = time.time()

    alleles = [Allele('HLA-A*01:01'), Allele('HLA-B*07:02'), Allele('HLA-C*03:01')]
    prot_seqs = FileReader.read_fasta(input_file, in_type=Protein)
    peptides = list(generate_peptides_from_proteins(prot_seqs, 9))
    print(len(peptides), 'peptides generated')
    bindings, thresh = get_binding_affinities_and_thresholds(peptides, alleles, randomize, threshold)

    if verbose:
        print('Threshold and binding affinity distribution per allele')
        for al in alleles:
            print('   ', al, thresh[al.name])
            print('   ', np.percentile(bindings[al], [0, 25, 50, 75, 90, 95, 98, 99, 100]))

    solver_cls, solver_kwargs = get_solver(solver, verbose, max_epitopes, max_aminoacids)
    if verbose:
        print('Using solver', solver_cls)

    solver_creation_time = time.time()
    solver = solver_cls(bindings, thresh, **solver_kwargs)
    
    solver_start_time = time.time()
    tour, peptides = solver.solve()
    solver_end_time = time.time()

    print(tour)
    print(peptides)

    if measure_time:
        print('==== Stopwatch')
        print('Total time           : %.2f s' % (solver_end_time - program_start_time))
        print('Pre-processing time  : %.2f s' % (solver_creation_time - program_start_time))
        print('Model creation time  : %.2f s' % (solver_start_time - solver_creation_time))
        print('Solving time         : %.2f s' % (solver_end_time - solver_start_time))


if __name__ == '__main__':
    main()