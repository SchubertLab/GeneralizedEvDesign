from __future__ import print_function

import pandas as pd
import numpy as np
from builtins import map

import time
import logging

from collections import defaultdict
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


LOGGER = None


def get_alleles_and_thresholds():
    return {     # According to Table 1 in Toussaint et al. (2011)
        #                  population cov.  :   threshold
        Allele('A*01:01',  4.498520 / 100.0):   0.08710, 
        Allele('A*02:01', 10.693200 / 100.0):   1.25720, 
        Allele('A*02:05',  0.884956 / 100.0):   0.41790, 
        Allele('A*03:01',  3.687320 / 100.0):   0.05527, 
        Allele('A*11:01',  7.522120 / 100.0):   0.04356, 
        Allele('A*24:02', 12.905600 / 100.0):   1.15830, 
        Allele('A*31:01',  2.433630 / 100.0):   0.09604, 
        Allele('A*68:01',  1.769910 / 100.0):   2.39910, 
        Allele('B*04:01',  1.548670 / 100.0):   5.22380, 
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


def compute_consensus_and_conservation(sequences):
    # warning: removes non-aminoacids (e.g. gaps) from the consensus sequence
    length = max(map(len, sequences))

    consensus, conservation = [], []
    for i in range(length):
        freqs = defaultdict(int)
        for p in sequences:
            if i < len(p):
                freqs[p[i]] += 1
        
        cons = None
        for k, v in freqs.items():
            if cons is None or v > cons[1]:
                cons = k, v
        
        if cons[0].isalpha():
            consensus.append(cons[0])
            conservation.append(float(cons[1]) / len(sequences))
    
    return ''.join(consensus), conservation


def get_peptides(input_file, min_conservation):
    proteins = FileReader.read_fasta(input_file, in_type=Protein)
    if len(proteins[0].transcript_id.split('.')) != 6:
        return list(generate_peptides_from_proteins(proteins, 9))

    # the following code is specific to HIV1_2017_aligned_sequences.fasta 's format
    # we first group proteins by gene and subtype
    proteins_by_subtype_and_gene = defaultdict(list)
    for prot in proteins:
        try:
            subtype, country, strain1, strain2, accession, gene = prot.transcript_id.split('.')
        except ValueError:
            LOGGER.debug('Cannot parse %s', prot.transcript_id)
            continue

        prot.gene_id = gene
        prot.transcript_id = accession
        prot.vars = {
            'subtype': subtype,
            'country': country,
            'strain1': strain1,
            'strain2': strain2,
            'accession': accession,
            'gene': gene,
        }

        proteins_by_subtype_and_gene[(subtype, gene)].append(prot)

    # then for each gene and subtype we compute the consensus sequence and aminoacid conservation
    # and extract conserved peptides from the consensus sequence 
    peptides = {}
    for (subtype, gene), pros in proteins_by_subtype_and_gene.iteritems():
        if len(pros) > 10:
            LOGGER.debug('Generating conserved consensus peptides for gene %s of subtype %s (%d proteins)',
                         gene, subtype, len(pros))
        else:
            LOGGER.debug('Not enough proteins (%d) to generate a reliable consensus for gene %s of subtype %s',
                         len(pros), gene, subtype)
            continue

        consensus, conservation = compute_consensus_and_conservation(pros)
        protein = Protein(consensus, gene_id=gene, transcript_id='%s-%s' % (subtype, gene))
        for i in xrange(len(consensus) - 8):
            if all(c >= min_conservation for c in conservation[i:i + 9]):
                seq = consensus[i:i + 9]
                if seq not in peptides:
                    peptides[seq] = Peptide(seq)
                peptides[seq].proteins[protein.transcript_id] = protein
                peptides[seq].proteinPos[protein.transcript_id].append(i)

    return peptides.values()


def get_binding_affinities_and_thresholds(peptides, randomize, bindings_path):
    ''' can either provide realistic binding affinities and thresholds, or
        use randomized values (for benchmarking purposes)
    '''
    allele_thresholds = get_alleles_and_thresholds()
    if not bindings_path:
        bindings = EpitopePredictorFactory('BIMAS').predict(peptides, allele_thresholds.keys())
    else:
        bindings = pd.read_csv(bindings_path)
        bindings['Seq'] = list(map(Peptide, bindings['Seq']))  # we don't have source protein
        bindings = bindings.set_index(['Seq', 'Method'])
        bindings.columns = list(map(Allele, bindings.columns))
        bindings = EpitopePredictionResult(bindings)
        LOGGER.info('Binding affinities loaded from', bindings_path)

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

    if not bindings_path:
        bindings.to_csv('resources/bindings.csv')

    return bindings, thresh


def get_solver_class(solver):
    solver_kwargs = {}

    if solver == 'gcb':
        solver_cls = MosaicVaccineGreedy
        raise NotImplementedError()
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
@click.option('--max-epitopes', '-e', default=0.0, help='Maximum length of the vaccine in epitopes')
@click.option('--min-alleles', '-A', default=0.0, help='Minimum number of alleles to cover with the vaccine')
@click.option('--min-antigens', '-g', default=0.0, help='Minimum antigens to cover with the vaccine')
@click.option('--min-conservation', '-c', default=0.0, help='Minimum conservation of selected epitopes')
@click.option('--verbose', '-v', is_flag=True, help='Print debug messages')
@click.option('--randomize', '-r', default=0.0, help='Randomly assign affinities and select a given portion of epitopes')
def main(input_file, solver, verbose, randomize, binding_affinities,
         max_aminoacids, max_epitopes, min_alleles, min_antigens, min_conservation):

    logging.basicConfig(level=(logging.DEBUG) if verbose else logging.INFO,
                        format='%(asctime)s %(name)s %(levelname)s: %(message)s')

    global LOGGER
    LOGGER = logging.getLogger('design-mosaic')

    program_start_time = time.time()

    peptides = get_peptides(input_file, min_conservation)
    LOGGER.info('%d peptides generated', len(peptides))
    bindings, thresh = get_binding_affinities_and_thresholds(peptides, randomize, binding_affinities)

    solver_cls, solver_kwargs = get_solver_class(solver)
    LOGGER.debug('Using solver %s', solver_cls)

    solver_creation_time = time.time()
    solver = solver_cls(
        bindings, thresh, max_vaccine_aminoacids=max_aminoacids, max_vaccine_epitopes=max_epitopes,
        min_allele_coverage=min_alleles, min_antigen_coverage=min_antigens,
        min_epitope_conservation=None, verbosity=verbose, **solver_kwargs
    )

    solver_build_time = time.time()
    solver.build_model()
    
    solver_start_time = time.time()
    result = solver.solve()
    solver_end_time = time.time()

    result.pretty_print(LOGGER.info)

    LOGGER.info('==== Stopwatch')
    LOGGER.info('Total time           : %.2f s', solver_end_time - program_start_time)
    LOGGER.info('Inputs preparation   : %.2f s', solver_creation_time - program_start_time)
    LOGGER.info('Pre-processing       : %.2f s', solver_build_time - solver_creation_time)
    LOGGER.info('Model creation time  : %.2f s', solver_start_time - solver_build_time)
    LOGGER.info('Solving time         : %.2f s', solver_end_time - solver_start_time)


if __name__ == '__main__':
    main()