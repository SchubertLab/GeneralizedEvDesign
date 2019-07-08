from __future__ import print_function
import numpy as np
from builtins import map

from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.EpitopePrediction import EpitopePredictorFactory
from MosaicVaccineILP import MosaicVaccineILP
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
def main(input_file, max_aminoacids, max_epitopes):
    proteins = []
    alleles = [Allele("HLA-A*01:01"), Allele("HLA-B*07:02"), Allele("HLA-C*03:01")]
    prot_seqs = FileReader.read_fasta(input_file, in_type=Protein)
    peptides = list(generate_peptides_from_proteins(prot_seqs, 9))
    print(len(peptides), 'peptides generated')
    result = EpitopePredictorFactory("BIMAS").predict(peptides, alleles)
    result.to_csv('resources/bindings.csv')
    thresh = {
        "A*01:01": 0.03,
        "B*07:02": 1,
        "C*03:01": 1,
    }

    for al in alleles:
        print(al, np.percentile(result[al], [0, 1, 2, 5, 10, 25, 50, 75, 90, 95, 98, 99, 100]))

    n = MosaicVaccineILP(result, thresh, max_vaccine_aminoacids=max_aminoacids, max_vaccine_epitopes=max_epitopes,
                         solver="cbc")
    print(n.solve())


if __name__ == '__main__':
    main()