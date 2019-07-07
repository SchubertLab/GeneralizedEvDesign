from __future__ import print_function
from builtins import map
__author__ = 'Schubert'

import unittest
from pprint import pprint


from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
#from Fred2.EpitopeSelection.Mosaic import MosaicVaccineILP, MosaicVaccineGreedy, _calculate_length
from Mosaic import MosaicVaccineILP, MosaicVaccineGreedy, _calculate_length

import numpy as np
import unittest
import os
import pandas
from Fred2.Core import Protein, Allele
from Fred2.Core import generate_peptides_from_proteins
from Fred2.IO import FileReader
from Fred2.EpitopePrediction import EpitopePredictorFactory, EpitopePredictionResult
from Fred2.EpitopeSelection.PopCover import PopCover
from Fred2.Utility import generate_overlap_graph
import tempfile
import inspect
import Fred2
import click


@click.command()
@click.argument('sequences', type=click.Path(exists=True), default='testSequences.fasta')
@click.option('--trim', '-t', default=-1)
@click.option('--verbose', '-v', default=0)
@click.option('--threshold', '-T', default=0.01)
def main(sequences, trim, verbose, threshold):
    # some common alleles
    alleles = [Allele("HLA-A*01:01"), Allele("HLA-B*07:02"), Allele("HLA-C*03:01")] 
    if trim > 0:
        tempf = tempfile.mktemp()
        with open(tempf, 'w') as f:
            with open(sequences, 'r') as g:
                for i, row in enumerate(g):
                    if i <= trim:
                        f.write(row)
                    else:
                        break
        sequences = tempf
    prot_seqs = FileReader.read_fasta(sequences, in_type=Protein)
    peptides = list(generate_peptides_from_proteins(prot_seqs, 9))
    print(len(peptides), 'peptides generated')
    bindings = EpitopePredictorFactory("BIMAS").predict(peptides, alleles)
    bindings.to_csv('bindings.csv')
    top_k = int(threshold * len(bindings))
    assert top_k > 0
    top_peptides = set([
        pep for a in alleles
        for pep, _ in bindings.sort_values('HLA-' + a.name, ascending=False).head(top_k).index.values
    ])
    vaccine_peptides = EpitopePredictionResult(bindings[bindings.index.isin(top_peptides, level=0)])
    print(len(vaccine_peptides), 'peptides selected')

    #quants = bindings.quantile(threshold)
    #thresh = {a.name: q for a, q in quants.to_dict().items()}
    #thresh = {'A*10:01': 10.0, 'B*07:02': 10.0, 'C*03:01': 10.0}
    #print('using the following thresholds', thresh)

    print('Building model')
    n = MosaicVaccineILP(vaccine_peptides, threshold=None, t_max=15, solver="cbc", verbosity=verbose)
    print('Solving model')
    print(n.solve())


if __name__ == '__main__':
    main()
