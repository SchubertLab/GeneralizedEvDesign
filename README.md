# Generalized Epitope-based Vaccine Design
General graph-based framework to design epitope vaccines.

# Setup
This framework makes heavy use of [Fred2](https://github.com/SchubertLab/Fred2) and one of the persistent solvers supported by Pyomo (currently either [Gurobi](https://www.gurobi.com/), which is the default we use, or [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)).

# Quick Start
For your convenience, we provide a makefile that invokes the necessary commands in the right order and with the right arguments, using the sample files in the `resources` folder:

```
make
```

This will process the data, produce the vaccine using all the supported methods, and save their evaluation in the `dev` folder. The vaccines are files named `made-<method>-vaccine.csv` and their evaluations are named `made-<method>-evaluation.csv.`

For a little more control, the following command allows you to specify the required input files and the output folder:

```
make all alleles=resources/alleles-small.csv proteins=resources/hiv1-bc-env-small.fasta BASE_DIR=dev
```

You are encouraged to read the next sections to learn the details, or jump straight to the Makefile.

# Command Description
Because of the long processing time required to design a vaccine, this process is split several stages, with intermediate results saved in CSV files. In the following sections we give concrete examples that work out of the box on how to run the commands; of course, you are free to modify the parameters to suit your needs.

By convention, the last argument of each command is the output file, and usage help can be printed by using the `--help` option.

## Data Preparation

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

Following are the commands required to prepare the data:

0. (Only if using the Makefile) start by copying the input files in the work directory:

   ```
   make init alleles=resources/alleles-small.csv proteins=resources/hiv1-bc-env-small.fasta
   ```

   By default, the work directory is `./dev/`, but it can be customized by using the argument `BASE_DIR=<dir>`. Note that if you do this now, you must include this argument for all subsequent `make` usages, too.

1. Extract the peptides and their coverage from the proteins you want to target:

   ```
   make coverage proteins=resources/hiv1-bc-env-small.fasta  # customize options with COVERAGE_OPTS="..."
   ```

   Or:

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
   make affinities alleles=resources/alleles-small.csv  # customize options with AFFINITIES_OPTS="..."
   ```
   
   Or:

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
   make epitopes   # customize options with EPITOPES_OPTS="..."
   ```

   Note that make will automatically use the intermediate files produced previously.

   Manual invocation via:

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

4. Compute the cleavage scores between all pairs of epitopes:

   ```
   make cleavages  # customize options with CLEAVAGE_OPTS="..."
   ```

   Or:

   ```
   python data_preparation.py -v compute-cleavages dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-cleavages.csv
   ```
   
   Sample output:

   | from      | to        |   score |
   | --------- | --------- | ------- |
   | FLGAAGSTM | FLGAAGSTM | -0.001… |
   | FLGAAGSTM | NVWATHACV |  0.798… |
   | FLGAAGSTM | VTVYYGVPV |  0.673… |
   | FLGAAGSTM | FIMIVGGLI |  1.222… |
   | FLGAAGSTM | KLTPLCVTL |  1.310… |
   | ...       | ...       |     ... |

## Vaccine Design
 - Mosaic: use this generalized framework to design a mosaic vaccine
   
   ```
   make mosaic-vaccine  # customize options with MOSAIC_OPTS="..."
   ```
   
   Or:

   ```
   python design.py -v mosaic dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-vaccine-mosaic.csv
   ```

 - String of Beads: use this generalized framework to design a string-of-beads vaccine
   
   ```
   make string-of-beads-vaccine  # customize options with STRING_OF_BEADS_OPTS="..."
   ```

   Or:

   ```
   python design.py -v string-of-beads dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-cleavages dev/hiv1-bc-env-small-vaccine-string-of-beads.csv
   ```

 - OptiTope: based on [1] and [2]

   ```
   make optitope-vaccine  # customize options with OPTITOPE_OPTS="..."
   ```

   Or:

   ```
   python design.py -v optitope dev/hiv1-bc-env-small-affinities.csv resources/alleles-small.csv dev/hiv1-bc-env-small-vaccine-optitope.csv
   ```

- PopCover: based on [3]

   ```
   make popcover-vaccine  # customize options with POPCOVER_OPTS="..."
   ```

   Or:

   ```
   python design.py -v popcover dev/hiv1-bc-env-small-coverage.csv dev/hiv1-bc-env-small-affinities.csv resources/alleles-small.csv dev/hiv1-bc-env-small-vaccine-popcover.csv
   ```

Sample output (same format for all methods):

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
make mosaic-evaluation  # make mosaic is also valid
```

Make can of course evaluate the vaccines produced by all methods via `make <method>-evaluation`, or simply `make <method>`.

The full command is as follows, make sure to use the correct input file for the vaccine:

```
python evaluation.py -v resources/hiv1-bc-env-small.fasta dev/hiv1-bc-env-small-coverage.csv resources/alleles-small.csv dev/hiv1-bc-env-small-epitopes.csv dev/hiv1-bc-env-small-vaccine-mosaic.csv dev/hiv1-bc-env-small-evaluation-mosaic.csv
```

Sample output:

| norm_prot_coverage | prot_coverage | pop_coverage | conservation | rel_pop_coverage | immunogen | max_pop_coverage |
| ------------------ | ------------- | ------------ | ------------ | ---------------- | --------- | ---------------- |
|               0.75 |            15 |       0.524… |       0.105… |           0.780… |    1.960… |           0.671… |

# References
[1] Toussaint NC, D ̈onnes P, Kohlbacher O. A mathematical framework for the selection of an optimal set of peptides forepitope-based vaccines.PLoS Comput Biol2008:4: e1000246.17.

[2] Toussaint NC, Kohlbacher O. OptiTope – a web server for theselection of an optimal set of peptides for epitope-basedvaccines.Nucleic Acids Res2009:37(suppl 2): W617–W622

[3] Lundegaard C, Buggert M, Karlsson A, Lund O, Perez C,Nielsen M. PopCover: a method for selecting of peptides withoptimal population and pathogen coverage. In:Proceedings ofthe First ACM International Conference on Bioinformatics andComputational Biology,2010.ACM