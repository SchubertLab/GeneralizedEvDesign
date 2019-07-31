Generalized Epitope-based Vaccine Design
====

## Data Preparation
Because of the long time required to prepare the inputs necessary to design a vaccine, this process is split in three stages, with intermediate results saved in files.

1. Extract the peptides and their coverage from the proteins you want to target:

```
python data_preparation.py -v extract-peptides resources/hiv1-bc-env-small.fasta dev/hiv1-bc-env-coverage-small.csv
```

2. Compute the binding affinities between these peptides and the HLA alleles of interest:

```
python data_preparation.py -v compute-affinities resources/alleles-small.csv dev/hiv1-bc-env-small-coverage.csv dev/hiv1-bc-env-small-affinities.csv
```

3. Extract the epitopes, their immunogenicity and their protein and HLA coverage.

```
python data_preparation.py -v extract-epitopes resources/alleles-small.csv dev/hiv1-bc-env-small-coverage.csv dev/hiv1-bc-env-small-affinities.csv dev/hiv1-bc-env-small-epitopes.csv
```

## Vaccine Design
 - Mosaic

```
python design.py -v mosaic dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-vaccine.csv
```

## Vaccine Evaluation

Evaluation computes the following metrics: total immunogenicity, allele coverage, pathogen coverage, average epitope conservation and population coverage.

```
python evaluation.py dev/hiv1-bc-env-small-vaccine.csv dev/hiv1-bc-env-small-bindings.csv
```