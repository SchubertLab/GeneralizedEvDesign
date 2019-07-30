Generalized Epitope-based Vaccine Design
====

## Data Preparation
1. Extract the peptides from the proteins you want to target, as well as the proteins covered by each peptide:

```
python data_preparation.py -v extract-peptides resources/hiv1-bc-env.fasta dev/hiv1-bc-env-coverage.csv
```

2. Compute the binding affinities between these peptides and the HLA alleles of interest:

```
python data_preparation.py -v compute-bindings dev/hiv1-bc-env-coverage.csv dev/hiv1-bc-env-bindings.csv
```