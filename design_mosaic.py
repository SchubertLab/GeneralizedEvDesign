from __future__ import print_function
import pandas as pd
import numpy as np
from builtins import map

import time
from Fred2.Core.Peptide import Peptide
from MosaicVaccineILP import MosaicVaccineILP
from MosaicVaccineLazyILP import MosaicVaccineLazyILP
from MosaicVaccineGreedy import MosaicVaccineGreedy
from Fred2.Core import Protein, Allele, Peptide
from Fred2.Core import generate_peptides_from_proteins
from Fred2.IO import FileReader
from Fred2.EpitopePrediction import EpitopePredictorFactory, EpitopePredictionResult
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


def get_binding_affinities_and_thresholds(peptides, alleles, randomize, threshold, bindings_path):
    if not bindings_path:
        bindings = EpitopePredictorFactory('BIMAS').predict(peptides, alleles)
    else:
        bindings = pd.read_csv(bindings_path)
        bindings['Seq'] = list(map(Peptide, bindings['Seq']))  # we don't have source protein
        bindings = bindings.set_index(['Seq', 'Method'])
        bindings.columns = list(map(Allele, bindings.columns))
        bindings = EpitopePredictionResult(bindings)
        print('Binding affinities loaded from', bindings_path)

    if randomize > 0:
        if randomize > 1:
            randomize = randomize / 100
        if bindings_path:
            print('WARNING: randomize and binding affinities both specified! The specified affinities will not')
            print('WARNING: be overwritten, but randominze will not work as intended if the affinities are not')
            print('WARNING: suitably randomized (i.e. uniformly between 0 and 1).')
            print('WARNING: ')
            print('WARNING: To generate randomized affinities, run once with randomize and *without* binding affinities.')
            print('WARNING: Randomized affinities will be generated and stored in "./resources/bindings.csv"')
            print('WARNING: and can be loaded in subsequent runs.')
            print('WARNING: ')
            print('WARNING: If you already did this, this warning is safe to ignore.')
            print()

        # chosen so that 100*randomize % of peptides have at least one allele
        # with larger binding strength, assuming that these are uniform(0, 1)
        tt = (1 - randomize)**(1.0 / len(alleles))
        thresh = {}
        for col in bindings.columns:
            if not bindings_path:
                bindings[col] = np.random.random(size=len(bindings))
            thresh[col.name] = tt
    else:
        thresh = THRESHOLD_PRESETS[threshold] if threshold < len(THRESHOLD_PRESETS) else THRESHOLD_PRESETS[-1]

    if not bindings_path:
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
@click.option('--binding-affinities', '-b', type=click.Path('r'), help='Pre-computed binding HLA-peptide binding affinities')
@click.option('--max-aminoacids', '-a', default=15, help='Maximum length of the vaccine in aminoacids')
@click.option('--max-epitopes', '-e', default=999, help='Maximum length of the vaccine in epitopes')
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
@click.option('--threshold', '-t', default=0, help='Select one of the threshdld presets (0-7)')
@click.option('--measure-time', '-T', is_flag=True, help='Print the total time at the end') 
@click.option('--randomize', '-r', default=-1.0, help='Randomly assign affinities and select a given portion of epitopes')
def main(input_file, solver, max_aminoacids, max_epitopes, verbose, threshold, measure_time, randomize, binding_affinities):
    program_start_time = time.time()

    alleles = [Allele('HLA-A*01:01'), Allele('HLA-B*07:02'), Allele('HLA-C*03:01')]
    prot_seqs = FileReader.read_fasta(input_file, in_type=Protein)
    peptides = list(generate_peptides_from_proteins(prot_seqs, 9))
    print(len(peptides), 'peptides generated')

    bindings, thresh = get_binding_affinities_and_thresholds(peptides, alleles, randomize,
                                                             threshold, binding_affinities)

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