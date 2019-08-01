.PHONY: all init alleles proteins coverage affinities epitopes cleavages mosaic-vaccine mosaic-eval mosaic optitope-vaccine optitope-eval optitope popcover-vaccine popcover-eval popcover clean

BASE_DIR=./dev
alleles:=./resources/alleles-small.csv
ALLELES:=$(BASE_DIR)/made-alleles.csv
proteins:=./resources/hiv1-bc-env-small.fasta 
PROTEINS:=$(BASE_DIR)/made-proteins.csv
COVERAGE:=$(BASE_DIR)/made-coverage.csv
COVERAGE_OPTS:=
AFFINITIES:=$(BASE_DIR)/made-affinities.csv
AFFINITIES_OPTS:=
EPITOPES:=$(BASE_DIR)/made-epitopes.csv
EPITOPES_OPTS:=
CLEAVAGES:=$(BASE_DIR)/made-cleavages.csv
CLEAVAGES_OPTS:=
MOSAIC_VACCINE:=$(BASE_DIR)/made-mosaic-vaccine.csv
MOSAIC_OPTS:=
MOSAIC_EVAL:=$(BASE_DIR)/made-mosaic-evaluation.csv
POPCOVER_VACCINE:=$(BASE_DIR)/made-popcover-vaccine.csv
POPCOVER_OPTS:=
POPCOVER_EVAL:=$(BASE_DIR)/made-popcover-evaluation.csv
OPTITOPE_VACCINE:=$(BASE_DIR)/made-optitope-vaccine.csv
OPTITOPE_OPTS:=
OPTITOPE_EVAL:=$(BASE_DIR)/made-optitope-evaluation.csv

all: mosaic optitope popcover
init: $(ALLELES) $(PROTEINS)
alleles: $(ALLELES)
proteins: $(PROTEINS)
coverage: $(COVERAGE)
affinities: $(AFFINITIES)
epitopes: $(EPITOPES)
cleavages: $(CLEAVAGES)
mosaic-vaccine: $(MOSAIC_VACCINE)
mosaic-eval: $(MOSAIC_EVAL)
mosaic: mosaic-eval
optitope-vaccine: $(OPTITOPE_VACCINE)
optitope-eval: $(OPTITOPE_EVAL)
optitope: optitope-eval
popcover-vaccine: $(POPCOVER_VACCINE)
popcover-eval: $(POPCOVER_EVAL)
popcover: popcover-eval

$(ALLELES):
	mkdir -p $(BASE_DIR)
	cp $(alleles) $(ALLELES)

$(PROTEINS):
	mkdir -p $(BASE_DIR)
	cp $(proteins) $(PROTEINS)

$(COVERAGE): $(PROTEINS)
	python data_preparation.py -v extract-peptides $(PROTEINS) $(COVERAGE) $(COVERAGE_OPTS)

$(AFFINITIES): $(COVERAGE) $(ALLELES)
	python data_preparation.py -v compute-affinities $(ALLELES) $(COVERAGE) $(AFFINITIES) $(AFFINITIES_OPTS)

$(EPITOPES): $(ALLELES) $(COVERAGE) $(AFFINITIES)
	python data_preparation.py -v extract-epitopes $(ALLELES) $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(EPITOPES_OPTS)

$(CLEAVAGES): $(EPITOPES)
	python data_preparation.py -v compute-cleavages $(EPITOPES) $(CLEAVAGES) $(CLEAVAGES_OPTS)

$(MOSAIC_VACCINE): $(EPITOPES)
	python design.py -v mosaic $(EPITOPES) $(MOSAIC_VACCINE) $(MOSAIC_OPTS)

$(MOSAIC_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL)

$(OPTITOPE_VACCINE): $(AFFINITIES) $(ALLELES)
	python design.py -v optitope $(AFFINITIES) $(ALLELES) $(OPTITOPE_VACCINE) $(OPTITOPE_OPTS)

$(OPTITOPE_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(EPITOPES)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)

$(POPCOVER_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v popcover $(COVERAGE) $(AFFINITIES) $(ALLELES) $(POPCOVER_VACCINE) $(OPTITOPE_OPTS)

$(POPCOVER_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE) $(POPCOVER_EVAL)

clean:
	rm $(ALLELES) $(PROTEINS) $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(CLEAVAGES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL) $(POPCOVER_VACCINE) $(POPCOVER_EVAL) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)