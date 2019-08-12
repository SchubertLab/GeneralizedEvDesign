#!/bin/bash
set -ex

BASEDIR="experiments/pareto"

mkdir -p "$BASEDIR/results"
if [ ! -f "$BASEDIR/results/made-proteins.fasta" ]; then
    python data_preparation.py -v random-sequences "experiments/resources/hiv1-bc-env.fasta" 300 "$BASEDIR/made-proteins.fasta"
else
    echo "Not overwriting proteins!"
fi;

if [ ! -f "$BASEDIR/results/made-alleles.csv" ]; then
    cp "experiments/resources/alleles.csv" "$BASEDIR/results/made-alleles.csv"
else
    echo "Not overwriting alleles!"
fi;

make epitopes cleavages BASE_DIR="$BASEDIR/results" COVERAGE_OPTS="--max-edits 0 --top-n -1"
python tradeoff.py \
    "$BASEDIR/results/made-epitopes.csv" \
    "$BASEDIR/results/made-cleavages.csv" \
    "$BASEDIR/results/made-tradeoff.csv" \
    --pareto-steps 11 --cocktail 1 \
    --max-aminoacids 0 --max-epitopes 10 \
    --min-alleles 0 --min-proteins 0