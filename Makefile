.PHONY: all init alleles proteins coverage affinities epitopes cleavages preparation mosaic-vaccine mosaic-eval mosaic string-of-beads-vaccine string-of-beads-eval string-of-beads optitope-vaccine optitope-eval optitope popcover-vaccine popcover-eval popcover clean

# include custom configuration (variables from command line have precedence)
BASE_DIR?=./dev
CONFIG?="$(BASE_DIR)/config.mak"
-include $(CONFIG) 

# set default values for all variables that were not defined in the config or from the command line
alleles							?= ./resources/alleles-small.csv
ALLELES_NAME					?= made-alleles.csv
ALLELES							?= $(BASE_DIR)/$(ALLELES_NAME)
proteins						?= ./resources/hiv1-bc-env-small.fasta 
PROTEINS_NAME					?= made-proteins.fasta
PROTEINS						?= $(BASE_DIR)/$(PROTEINS_NAME)
COVERAGE_NAME					?= made-coverage.csv
COVERAGE						?= $(BASE_DIR)/$(COVERAGE_NAME)
COVERAGE_OPT					?= 
AFFINITIES_NAME					?= made-affinities.csv
AFFINITIES						?= $(BASE_DIR)/$(AFFINITIES_NAME)
AFFINITIES_OPTS					?= 
EPITOPES_NAME					?= made-epitopes.csv
EPITOPES						?= $(BASE_DIR)/$(EPITOPES_NAME)
EPITOPES_OPTS					?= 
CLEAVAGES_NAME					?= made-cleavages.csv
CLEAVAGES						?= $(BASE_DIR)/$(CLEAVAGES_NAME)
CLEAVAGES_OPTS					?= 
MOSAIC_VACCINE_NAME				?= made-mosaic-vaccine.csv
MOSAIC_VACCINE					?= $(BASE_DIR)/$(MOSAIC_VACCINE_NAME)
MOSAIC_OPTS						?= 
MOSAIC_EVAL_NAME				?= made-mosaic-evaluation.csv
MOSAIC_EVAL						?= $(BASE_DIR)/$(MOSAIC_EVAL_NAME)
STRING_OF_BEADS_VACCINE_NAME 	?= made-string-of-beads-vaccine.csv
STRING_OF_BEADS_VACCINE			?= $(BASE_DIR)/$(STRING_OF_BEADS_VACCINE_NAME)
STRING_OF_BEADS_OPTS			?= 
STRING_OF_BEADS_EVAL_NAME		?= made-string-of-beads-evaluation.csv
STRING_OF_BEADS_EVAL			?= $(BASE_DIR)/$(STRING_OF_BEADS_EVAL_NAME)
POPCOVER_VACCINE_NAME			?= made-popcover-vaccine.csv
POPCOVER_VACCINE				?= $(BASE_DIR)/$(POPCOVER_VACCINE_NAME)
POPCOVER_OPTS					?= 
POPCOVER_EVAL_NAME				?= made-popcover-evaluation.csv
POPCOVER_EVAL					?= $(BASE_DIR)/$(POPCOVER_EVAL_NAME)
OPTITOPE_VACCINE_NAME			?= made-optitope-vaccine.csv
OPTITOPE_VACCINE				?= $(BASE_DIR)/$(OPTITOPE_VACCINE_NAME)
OPTITOPE_OPTS					?= 
OPTITOPE_EVAL_NAME				?= made-optitope-evaluation.csv
OPTITOPE_EVAL					?= $(BASE_DIR)/$(OPTITOPE_EVAL_NAME)

# define aliases
all: mosaic optitope popcover
init: $(ALLELES) $(PROTEINS)
alleles: $(ALLELES)
proteins: $(PROTEINS)
coverage: $(COVERAGE)
affinities: $(AFFINITIES)
epitopes: $(EPITOPES)
cleavages: $(CLEAVAGES)
preparation: cleavages epitopes affinities coverage
mosaic-vaccine: $(MOSAIC_VACCINE)
mosaic-eval: $(MOSAIC_EVAL)
mosaic: mosaic-eval
string-of-beads-vaccine: $(STRING_OF_BEADS_VACCINE)
string-of-beads-eval: $(STRING_OF_BEADS_EVAL)
string-of-beads: string-of-beads-eval
optitope-vaccine: $(OPTITOPE_VACCINE)
optitope-eval: $(OPTITOPE_EVAL)
optitope: optitope-eval
popcover-vaccine: $(POPCOVER_VACCINE)
popcover-eval: $(POPCOVER_EVAL)
popcover: popcover-eval

# we only copy the input files when the targets are missing
# this is to prevent somebody from accedentally overwriting them with the samples in the resources folder
# when they run make -B <something>
# NB: make -B will redo *everything*, in order to only make the target you specify, you should touch one of its dependencies
$(ALLELES):
	mkdir -p $(BASE_DIR)
	[ ! -f $(ALLELES) ] && cp $(alleles) $(ALLELES) || echo "Not overwriting allele file, remove it manually first!"

$(PROTEINS):
	mkdir -p $(BASE_DIR)
	[ ! -f $(PROTEINS) ] && cp $(proteins) $(PROTEINS) || echo "Not overwriting proteins file, remove it manually first!"

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

$(STRING_OF_BEADS_VACCINE): $(EPITOPES) $(CLEAVAGES)
	python design.py -v string-of-beads $(EPITOPES) $(CLEAVAGES) $(STRING_OF_BEADS_VACCINE)

$(STRING_OF_BEADS_EVAL): $(STRING_OF_BEADS_VACCINE)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(STRING_OF_BEADS_VACCINE) $(STRING_OF_BEADS_EVAL)

$(OPTITOPE_VACCINE): $(AFFINITIES) $(ALLELES)
	python design.py -v optitope $(AFFINITIES) $(ALLELES) $(OPTITOPE_VACCINE) $(OPTITOPE_OPTS)

$(OPTITOPE_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(EPITOPES)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)

$(POPCOVER_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v popcover $(COVERAGE) $(AFFINITIES) $(ALLELES) $(POPCOVER_VACCINE) $(OPTITOPE_OPTS)

$(POPCOVER_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE)
	python evaluation.py -v $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE) $(POPCOVER_EVAL)

clean:
	rm -f $(ALLELES) $(PROTEINS) $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(CLEAVAGES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL) $(POPCOVER_VACCINE) $(POPCOVER_EVAL) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)