.PHONY: all init alleles proteins coverage affinities epitopes cleavages preparation mosaic-vaccine mosaic-eval mosaic string-of-beads-vaccine string-of-beads-eval string-of-beads optitope-vaccine optitope-eval optitope popcover-vaccine popcover-eval popcover clean

# include custom configuration (variables from command line have precedence)
BASE_DIR?=./dev
CONFIG?="$(BASE_DIR)/config.mak"
-include $(CONFIG) 

# set default values for all variables that were not defined in the config or from the command line
alleles					?= ./resources/alleles-small.csv
ALLELES_NAME			?= made-alleles.csv
ALLELES					?= $(BASE_DIR)/$(ALLELES_NAME)
proteins				?= ./resources/hiv1-bc-env-small.fasta 
PROTEINS_NAME			?= made-proteins.fasta
PROTEINS				?= $(BASE_DIR)/$(PROTEINS_NAME)

COVERAGE_NAME			?= made-coverage.csv
COVERAGE_LOG_NAME		?= made-coverage.log
COVERAGE				?= $(BASE_DIR)/$(COVERAGE_NAME)
COVERAGE_LOG			?= $(BASE_DIR)/$(COVERAGE_LOG_NAME)
COVERAGE_OPT			?= 

AFFINITIES_NAME			?= made-affinities.csv
AFFINITIES_LOG_NAME		?= made-affinities.log
AFFINITIES				?= $(BASE_DIR)/$(AFFINITIES_NAME)
AFFINITIES_LOG			?= $(BASE_DIR)/$(AFFINITIES_LOG_NAME)
AFFINITIES_OPTS			?= 

EPITOPES_NAME			?= made-epitopes.csv
EPITOPES_LOG_NAME		?= made-epitopes.log
EPITOPES				?= $(BASE_DIR)/$(EPITOPES_NAME)
EPITOPES_LOG			?= $(BASE_DIR)/$(EPITOPES_LOG_NAME)
EPITOPES_OPTS			?= 

CLEAVAGES_NAME			?= made-cleavages.csv
CLEAVAGES_LOG_NAME		?= made-cleavages.log
CLEAVAGES				?= $(BASE_DIR)/$(CLEAVAGES_NAME)
CLEAVAGES_LOG			?= $(BASE_DIR)/$(CLEAVAGES_LOG_NAME)
CLEAVAGES_OPTS			?= 

MOSAIC_VACCINE_NAME		?= made-mosaic-vaccine.csv
MOSAIC_LOG_NAME			?= made-mosaic-vaccine.log
MOSAIC_VACCINE			?= $(BASE_DIR)/$(MOSAIC_VACCINE_NAME)
MOSAIC_LOG				?= $(BASE_DIR)/$(MOSAIC_LOG_NAME)
MOSAIC_OPTS				?= 

MOSAIC_EVAL_NAME		?= made-mosaic-evaluation.csv
MOSAIC_EVAL_LOG_NAME	?= made-mosaic-evaluation.log
MOSAIC_EVAL				?= $(BASE_DIR)/$(MOSAIC_EVAL_NAME)
MOSAIC_EVAL_LOG			?= $(BASE_DIR)/$(MOSAIC_EVAL_LOG_NAME)

STROBE_VACCINE_NAME 	?= made-string-of-beads-vaccine.csv
STROBE_LOG_NAME 		?= made-string-of-beads-vaccine.log
STROBE_VACCINE			?= $(BASE_DIR)/$(STROBE_VACCINE_NAME)
STROBE_LOG				?= $(BASE_DIR)/$(STROBE_LOG_NAME)
STROBE_OPTS				?= 

STROBE_EVAL_NAME		?= made-string-of-beads-evaluation.csv
STROBE_EVAL_LOG_NAME	?= made-string-of-beads-evaluation.log
STROBE_EVAL				?= $(BASE_DIR)/$(STROBE_EVAL_NAME)
STROBE_EVAL_LOG			?= $(BASE_DIR)/$(STROBE_EVAL_LOG_NAME)

POPCOVER_VACCINE_NAME	?= made-popcover-vaccine.csv
POPCOVER_LOG_NAME		?= made-popcover-vaccine.log
POPCOVER_VACCINE		?= $(BASE_DIR)/$(POPCOVER_VACCINE_NAME)
POPCOVER_LOG			?= $(BASE_DIR)/$(POPCOVER_LOG_NAME)
POPCOVER_OPTS			?= 

POPCOVER_EVAL_NAME		?= made-popcover-evaluation.csv
POPCOVER_EVAL_LOG_NAME	?= made-popcover-evaluation.log
POPCOVER_EVAL			?= $(BASE_DIR)/$(POPCOVER_EVAL_NAME)
POPCOVER_EVAL_LOG		?= $(BASE_DIR)/$(POPCOVER_EVAL_LOG_NAME)

OPTITOPE_VACCINE_NAME	?= made-optitope-vaccine.csv
OPTITOPE_LOG_NAME		?= made-optitope-vaccine.log
OPTITOPE_VACCINE		?= $(BASE_DIR)/$(OPTITOPE_VACCINE_NAME)
OPTITOPE_LOG			?= $(BASE_DIR)/$(OPTITOPE_LOG_NAME)
OPTITOPE_OPTS			?= 

OPTITOPE_EVAL_NAME		?= made-optitope-evaluation.csv
OPTITOPE_EVAL_LOG_NAME	?= made-optitope-evaluation.log
OPTITOPE_EVAL			?= $(BASE_DIR)/$(OPTITOPE_EVAL_NAME)
OPTITOPE_EVAL_LOG		?= $(BASE_DIR)/$(OPTITOPE_EVAL_LOG_NAME

AGGREGATE_EVAL_NAME		?= made-aggregate-evaluation.csv
AGGREGATE_EVAL_LOG_NAME	?= made-aggregate-evaluation.log
AGGREGATE_EVAL			?= $(BASE_DIR)/$(AGGREGATE_EVAL_NAME)
AGGREGATE_EVAL_LOG		?= $(BASE_DIR)/$(AGGREGATE_EVAL_LOG_NAME)
AGGREGATE_PATH_SPEC		?= $(BASE_DIR)/made-*-evaluation.csv

# define aliases
all: mosaic string-of-beads optitope popcover
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
string-of-beads-vaccine: $(STROBE_VACCINE)
string-of-beads-eval: $(STROBE_EVAL)
string-of-beads: string-of-beads-eval
optitope-vaccine: $(OPTITOPE_VACCINE)
optitope-eval: $(OPTITOPE_EVAL)
optitope: optitope-eval
popcover-vaccine: $(POPCOVER_VACCINE)
popcover-eval: $(POPCOVER_EVAL)
popcover: popcover-eval
aggregate-eval: $(AGGREGATE_EVAL)

# we only copy the input files when the targets are missing
# this is to prevent somebody from accedentally overwriting them with the samples in the resources folder
# when they run make -B <something>
# NB: make -B will redo *everything*, in order to only make the target you specify, you should remove its output 
$(ALLELES):
	mkdir -p $(BASE_DIR)
	[ ! -f $(ALLELES) ] && cp $(alleles) $(ALLELES) || echo "Not overwriting allele file, remove it manually first!"

$(PROTEINS):
	mkdir -p $(BASE_DIR)
	[ ! -f $(PROTEINS) ] && cp $(proteins) $(PROTEINS) || echo "Not overwriting proteins file, remove it manually first!"

$(COVERAGE): $(PROTEINS)
	python data_preparation.py -v -l $(COVERAGE_LOG) extract-peptides $(PROTEINS) $(COVERAGE) $(COVERAGE_OPTS)

$(AFFINITIES): $(COVERAGE) $(ALLELES)
	python data_preparation.py -v -l $(AFFINITIES_LOG) compute-affinities $(ALLELES) $(COVERAGE) $(AFFINITIES) $(AFFINITIES_OPTS)

$(EPITOPES): $(ALLELES) $(COVERAGE) $(AFFINITIES)
	python data_preparation.py -v -l $(EPITOPES_LOG) extract-epitopes $(ALLELES) $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(EPITOPES_OPTS)

$(CLEAVAGES): $(EPITOPES)
	python data_preparation.py -v -l $(CLEAVAGES_LOG) compute-cleavages $(EPITOPES) $(CLEAVAGES) $(CLEAVAGES_OPTS)

$(MOSAIC_VACCINE): $(EPITOPES)
	python design.py -v -l $(MOSAIC_LOG) mosaic $(EPITOPES) $(MOSAIC_VACCINE) $(MOSAIC_OPTS)

$(MOSAIC_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE)
	python evaluation.py -v vaccine -l $(MOSAIC_EVAL_LOG) $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL)

$(STROBE_VACCINE): $(EPITOPES) $(CLEAVAGES)
	python design.py -v string-of-beads -l $(STROBE_LOG) $(EPITOPES) $(CLEAVAGES) $(STROBE_VACCINE) $(STROBE_OPTS)

$(STROBE_EVAL): $(STROBE_VACCINE)
	python evaluation.py -v vaccine -l $(STROBE_EVAL_LOG) $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(STROBE_VACCINE) $(STROBE_EVAL)

$(OPTITOPE_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v optitope -l $(OPTITOPE_LOG) $(COVERAGE) $(AFFINITIES) $(ALLELES) $(OPTITOPE_VACCINE) $(OPTITOPE_OPTS)

$(OPTITOPE_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(EPITOPES)
	python evaluation.py -v vaccine -l $(OPTITOPE_EVAL_LOG) $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)

$(POPCOVER_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v popcover -l $(POPCOVER_LOG) $(COVERAGE) $(AFFINITIES) $(ALLELES) $(POPCOVER_VACCINE) $(OPTITOPE_OPTS)

$(POPCOVER_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE)
	python evaluation.py -v vaccine -l $(POPCOVER_EVAL_LOG) $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE) $(POPCOVER_EVAL)

$(AGGREGATE_EVAL): $(MOSAIC_EVAL) $(STROBE_EVAL) $(OPTITOPE_EVAL) $(POPCOVER_EVAL)
	python evaluation.py -v aggregate -l $(AGGREGATE_EVAL_LOG) $(AGGREGATE_PATH_SPEC) $(AGGREGATE_EVAL)

clean:
	rm -f $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(CLEAVAGES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL) $(POPCOVER_VACCINE) $(POPCOVER_EVAL) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)