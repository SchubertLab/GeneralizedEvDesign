from __future__ import print_function
import numpy as np
from builtins import map

import time
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from MosaicVaccineILP import MosaicVaccineILP
from MosaicVaccineLazyILP import MosaicVaccineLazyILP
from Fred2.Core import Protein, Allele
from Fred2.Core import generate_peptides_from_proteins
from Fred2.IO import FileReader
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopeSelection.PopCover import PopCover
from Fred2.Utility import generate_overlap_graph
import click
import Fred2


@click.command()
@click.argument('input-file', type=click.Path('r'))
@click.option('--max-aminoacids', '-a', default=15)
@click.option('--max-epitopes', '-e', default=999)
@click.option('--lazy', '-l', is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
@click.option('--threshold', '-t', default=0)
@click.option('--measure-time', '-T', is_flag=True)
def main(input_file, max_aminoacids, max_epitopes, lazy, verbose, threshold, measure_time):
    program_start_time = time.time()
    proteins = []
    alleles = [Allele('HLA-A*01:01'), Allele('HLA-B*07:02'), Allele('HLA-C*03:01')]
    prot_seqs = FileReader.read_fasta(input_file, in_type=Protein)
    peptides = list(generate_peptides_from_proteins(prot_seqs, 9))
    print(len(peptides), 'peptides generated')
    bindings = EpitopePredictorFactory('BIMAS').predict(peptides, alleles)
    bindings.to_csv('resources/bindings.csv')

    threshold_presets = [  # counts based on hivgen.fasta
        {'A*01:01': 0.1,  'B*07:02': 1.0, 'C*03:01': 0.3},   # 11 epitopes
        {'A*01:01': 100,  'B*07:02': 100, 'C*03:01': 0.23},  # 528 epitopes
        {'A*01:01': 0.04, 'B*07:02': 100, 'C*03:01': 100},   # 666 epitopes
        {'A*01:01': 100,  'B*07:02': 0.8, 'C*03:01': 100},   # 817 epitopes
        {'A*01:01': 0.04, 'B*07:02': 100, 'C*03:01': 0.23},  # 1194 epitopes
        {'A*01:01': 0.04, 'B*07:02': 0.8, 'C*03:01': 100},   # 1483 epitopes
        {'A*01:01': 0.04, 'B*07:02': 0.8, 'C*03:01': 0.23},  # 1841 epitopes
    ]
    thresh = threshold_presets[threshold] if threshold < len(threshold_presets) else threshold_presets[-1]

    if verbose:
        print('Binding affinity distribution per allele')
        for al in alleles:
            print('   ', al, np.percentile(bindings[al], [0, 25, 50, 75, 90, 95, 98, 99, 100]))

    solver_cls = MosaicVaccineLazyILP if lazy else MosaicVaccineILP
    solver_creation_time = time.time()
    n = solver_cls(bindings, thresh, max_vaccine_aminoacids=max_aminoacids,
                   max_vaccine_epitopes=max_epitopes, verbosity=int(verbose))
    solver_start_time = time.time()
    tour, peptides = n.solve()
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