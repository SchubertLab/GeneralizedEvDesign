Generalized Epitope-based Vaccine Design
====

## Data Preparation
Because of the long time required to prepare the inputs necessary to design a vaccine, this process is split in three stages, with intermediate results saved in CSV files.

The inputs required are the proteins to target, in fasta format, e.g.

```
>B.JP.2004.04JPDR6075B.AB221125
MRVTGIRKNCQLLWKWGTMLLGMLMICSALEQLWVTVYYGVPVWKDATTT
LFCASDAKAYDTEMHNVWATHACVPTDPNPQEVVLVNVTEEFNMWKNNMV
EQMHEDIINLWDQSLKPCVKLTPLCVTLHCTDLKNVTANSTASNSSTNWK
HGIRPVVSTQLLLNGSLAEEEVVLRSENFTNNAKTIIVQLNKTVVINCTR
PNNNTRKGIHIGAGRAIYATGAIIGDIRQAHCNLSETDWKNTLNKTVRKL
...
```

and the HLA alleles of interest, along with their percent frequency in the population and the binding threshold:

| allele      | frequency | threshold |
| ----------- | --------- | --------- |
| HLA-A*02:01 |   10.693… |     0.405 |
| HLA-A*24:02 |   12.906… |     0.405 |
| HLA-B*40:01 |    5.310… |     0.405 |
| HLA-B*51:01 |    3.245… |     0.405 |

The resources folder contains small samples that you can play with.

Following are the commands required to prepare the data; the last file is the output file. Use `--help` to see more information about each command.

1. Extract the peptides and their coverage from the proteins you want to target:

   ```
   python data_preparation.py -v extract-peptides resources/hiv1-bc-env-small.fasta dev/hiv1-bc-env-small-coverage.csv
   ```

   Sample output:

   | peptide   | proteins    |
   | --------- | ----------- |
   | IKQACPKVT | 8;13        |
   | RNLCLFGYH | 15          |
   | KEYALFYTL | 19          |
   | KTLEQIAEK | 9           |
   | RSSLRGLQR | 0;8;3;13;16 |
   | ...       | ...         |

2. Compute the binding affinities between these peptides and the HLA alleles of interest. These affinities are ic50 values scaled as follows: `affinity = 1 - log(ic50) / log(50000)`, so that 50, 500, 5000 and 50000 nM are scaled to 0.638, 0.426, 0.213 and 0.000 respectively.

   ```
   python data_preparation.py -v compute-affinities resources/alleles-small.csv dev/hiv1-bc-env-small-coverage.csv dev/hiv1-bc-env-small-affinities.csv
   ```
   
   Sample output:
   
   | Seq       | Method    | HLA-A*02:01 | HLA-A*24:02 | HLA-B*40:01 | ... |
   | --------- | --------- | ----------- | ----------- | ----------- | --- |
   | AADKLWVTV | netmhcpan |      0.256… |      0.070… |      0.138… | ... |
   | AAELLGRSS | netmhcpan |      0.016… |      0.004… |      0.028… | ... |
   | AAEQLWVTV | netmhcpan |      0.152… |      0.064… |      0.127… | ... |
   | AAGSTMGAA | netmhcpan |      0.093… |      0.011… |      0.037… | ... |
   | AAHCNISEG | netmhcpan |      0.024… |      0.016… |      0.043… | ... |
   | ...       | ...       |         ... |         ... |         ... | ... |

3. Extract the epitopes, their immunogenicity and their protein and HLA coverage:

   ```
   python data_preparation.py -v extract-epitopes resources/alleles-small.csv dev/hiv1-bc-env-small-coverage.csv dev/hiv1-bc-env-small-affinities.csv dev/hiv1-bc-env-small-epitopes.csv
   ```
   
   Sample output:
   
   | immunogen | alleles                 | proteins | epitope   |
   | --------- | ----------------------- | -------- | --------- |
   |    0.124… | HLA-A*02:01             | 9;18     | QMQEDIISL |
   |    0.142… | HLA-A*24:02             | 8;13     | SWFSITNWL |
   |    0.200… | HLA-A\*24:02;HLA-A\*02:01 | 8;13     | YQRWWIWSI |
   |    0.154… | HLA-C\*07:02;HLA-A\*24:02 | 19;14    | MYAPPIEGL |
   |    0.057… | HLA-A*02:01             | 6;15     | LLALDSWAS |

## Vaccine Design
 - Mosaic

   ```
   python design.py -v mosaic dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-vaccine-mosaic.csv
   ```

 - OptiTope

   ```
   python design.py -v optitope dev/hiv1-bc-env-small-affinities.csv resources/alleles-small.csv dev/hiv1-bc-env-small-vaccine-optitope.csv
   ```

Sample output:

| cocktail | index | epitope   |
| -------- | ----- | --------- |
|  0       |     0 | YQRWWIWSI |
|  0       |     1 | YTDTIYWLL |
|  0       |     2 | LLQYWSQEL |
|  0       |     3 | YFPNKTMNF |
|      ... |   ... | ...       |  

## Vaccine Evaluation
Evaluation computes the following metrics: total immunogenicity, allele coverage, pathogen coverage, average epitope conservation and population coverage. The population coverage is also computed relative to the maximum theoretical coverage that can be achieved with the given alleles.

```
python evaluation.py -v resources/hiv1-bc-env-small.fasta dev/hiv1-bc-env-small-coverage.csv resources/alleles-small.csv dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-vaccine-mosaic.csv dev/hiv1-bc-env-small-evaluation-mosaic.csv
```

Sample output:

| norm_prot_coverage | prot_coverage | pop_coverage | conservation | rel_pop_coverage | immunogen | max_pop_coverage |
| ------------------ | ------------- | ------------ | ------------ | ---------------- | --------- | ---------------- |
|               0.75 |            15 |       0.524… |       0.105… |           0.780… |    1.960… |           0.671… |