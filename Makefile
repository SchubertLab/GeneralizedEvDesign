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

COVERAGE_NAME			?= made-coverage
COVERAGE				?= $(BASE_DIR)/$(COVERAGE_NAME).csv
COVERAGE_LOG			?= $(BASE_DIR)/$(COVERAGE_NAME).log
COVERAGE_OPT			?=

AFFINITIES_NAME			?= made-affinities
AFFINITIES				?= $(BASE_DIR)/$(AFFINITIES_NAME).csv
AFFINITIES_LOG			?= $(BASE_DIR)/$(AFFINITIES_NAME).log
AFFINITIES_OPTS			?=

EPITOPES_NAME			?= made-epitopes
EPITOPES				?= $(BASE_DIR)/$(EPITOPES_NAME).csv
EPITOPES_LOG			?= $(BASE_DIR)/$(EPITOPES_NAME).log
EPITOPES_OPTS			?=

CLEAVAGES_NAME			?= made-cleavages
CLEAVAGES				?= $(BASE_DIR)/$(CLEAVAGES_NAME).csv
CLEAVAGES_LOG			?= $(BASE_DIR)/$(CLEAVAGES_NAME).log
CLEAVAGES_OPTS			?=

OVERLAPS_NAME			?= made-overlaps
OVERLAPS				?= $(BASE_DIR)/$(OVERLAPS_NAME).csv
OVERLAPS_LOG			?= $(BASE_DIR)/$(OVERLAPS_NAME).log
OVERLAPS_OPTS			?=

MOSAIC_NAME				?= made-mosaic-vaccine
MOSAIC_VACCINE			?= $(BASE_DIR)/$(MOSAIC_NAME).csv
MOSAIC_LOG				?= $(BASE_DIR)/$(MOSAIC_NAME).log
MOSAIC_OPTS				?=

MOSAIC_EVAL_NAME		?= $(MOSAIC_NAME)-evaluation
MOSAIC_EVAL				?= $(BASE_DIR)/$(MOSAIC_EVAL_NAME).csv
MOSAIC_EVAL_LOG			?= $(BASE_DIR)/$(MOSAIC_EVAL_LOG_NAME).csv

STROBE_NAME				?= made-string-of-beads-vaccine
STROBE_VACCINE			?= $(BASE_DIR)/$(STROBE_NAME).csv
STROBE_LOG				?= $(BASE_DIR)/$(STROBE_NAME).log
STROBE_OPTS				?=

STROBE_EVAL_NAME		?= $(STROBE_NAME)-evaluation
STROBE_EVAL				?= $(BASE_DIR)/$(STROBE_EVAL_NAME).csv
STROBE_EVAL_LOG			?= $(BASE_DIR)/$(STROBE_EVAL_NAME).log

POPCOVER_NAME			?= made-popcover-vaccine
POPCOVER_VACCINE		?= $(BASE_DIR)/$(POPCOVER_NAME).csv
POPCOVER_LOG			?= $(BASE_DIR)/$(POPCOVER_NAME).log
POPCOVER_OPTS			?= 

POPCOVER_EVAL_NAME		?= $(POPCOVER_NAME)-evaluation
POPCOVER_EVAL			?= $(BASE_DIR)/$(POPCOVER_EVAL_NAME).csv
POPCOVER_EVAL_LOG		?= $(BASE_DIR)/$(POPCOVER_EVAL_NAME).log

OPTITOPE_NAME			?= made-optitope-vaccine
OPTITOPE_VACCINE		?= $(BASE_DIR)/$(OPTITOPE_NAME).csv
OPTITOPE_LOG			?= $(BASE_DIR)/$(OPTITOPE_NAME).log
OPTITOPE_OPTS			?= 

OPTITOPE_EVAL_NAME		?= $(OPTITOPE_NAME)-evaluation
OPTITOPE_EVAL			?= $(BASE_DIR)/$(OPTITOPE_EVAL_NAME).csv
OPTITOPE_EVAL_LOG		?= $(BASE_DIR)/$(OPTITOPE_EVAL_NAME).log

AGGREGATE_EVAL_NAME		?= made-aggregate-evaluation
AGGREGATE_EVAL			?= $(BASE_DIR)/$(AGGREGATE_EVAL_NAME).csv
AGGREGATE_EVAL_LOG		?= $(BASE_DIR)/$(AGGREGATE_EVAL_NAME).log
AGGREGATE_PATH_SPEC		?= $(BASE_DIR)/made-*-evaluation.csv
AGGREGATE_OPTS			?=

# define aliases
all: mosaic string-of-beads optitope popcover
init: $(ALLELES) $(PROTEINS)
alleles: $(ALLELES)
proteins: $(PROTEINS)
coverage: $(COVERAGE)
affinities: $(AFFINITIES)
epitopes: $(EPITOPES)
overlaps: $(OVERLAPS)
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

$(OVERLAPS): $(EPITOPES)
	python data_preparation.py -v -l $(OVERLAPS_LOG) compute-overlaps $(EPITOPES) $(OVERLAPS) $(OVERLAPS_OPTS)

$(MOSAIC_VACCINE): $(EPITOPES) $(OVERLAPS)
	python design.py -v -l $(MOSAIC_LOG) mosaic $(EPITOPES) $(OVERLAPS) $(MOSAIC_VACCINE) $(MOSAIC_OPTS)

$(MOSAIC_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE)
	python evaluation.py -v -l $(MOSAIC_EVAL_LOG) vaccine $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL)

$(STROBE_VACCINE): $(EPITOPES) $(CLEAVAGES)
	python design.py -v -l $(STROBE_LOG) string-of-beads  $(EPITOPES) $(CLEAVAGES) $(STROBE_VACCINE) $(STROBE_OPTS)

$(STROBE_EVAL): $(STROBE_VACCINE)
	python evaluation.py -v -l $(STROBE_EVAL_LOG) vaccine $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(STROBE_VACCINE) $(STROBE_EVAL)

$(OPTITOPE_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v -l $(OPTITOPE_LOG) optitope  $(COVERAGE) $(AFFINITIES) $(ALLELES) $(OPTITOPE_VACCINE) $(OPTITOPE_OPTS)

$(OPTITOPE_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(EPITOPES)
	python evaluation.py -v -l $(OPTITOPE_EVAL_LOG) vaccine $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)

$(POPCOVER_VACCINE): $(COVERAGE) $(AFFINITIES) $(ALLELES)
	python design.py -v -l $(POPCOVER_LOG) popcover  $(COVERAGE) $(AFFINITIES) $(ALLELES) $(POPCOVER_VACCINE) $(POPCOVER_OPTS)

$(POPCOVER_EVAL): $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE)
	python evaluation.py -v -l $(POPCOVER_EVAL_LOG) vaccine $(PROTEINS) $(COVERAGE) $(ALLELES) $(EPITOPES) $(POPCOVER_VACCINE) $(POPCOVER_EVAL)

$(AGGREGATE_EVAL): $(MOSAIC_EVAL) $(STROBE_EVAL) $(OPTITOPE_EVAL) $(POPCOVER_EVAL)
	python evaluation.py -v -l $(AGGREGATE_EVAL_LOG) aggregate $(AGGREGATE_PATH_SPEC) $(AGGREGATE_EVAL) $(AGGREGATE_OPTS)

clean:
	rm -f $(COVERAGE) $(AFFINITIES) $(EPITOPES) $(CLEAVAGES) $(MOSAIC_VACCINE) $(MOSAIC_EVAL) $(POPCOVER_VACCINE) $(POPCOVER_EVAL) $(OPTITOPE_VACCINE) $(OPTITOPE_EVAL)