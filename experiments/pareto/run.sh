#!/bin/bash
set -ex

BASEDIR="experiments/pareto"

mkdir -p "$BASEDIR/results"
if [ ! -f "$BASEDIR/results/made-proteins.fasta" ]; then
    python data_preparation.py -v random-sequences "experiments/resources/hiv1-bc-env.fasta" 300 "$BASEDIR/results/made-proteins.fasta"
else
    echo "Not overwriting proteins!"
fi;

if [ ! -f "$BASEDIR/results/made-alleles.csv" ]; then
    cp "experiments/resources/alleles.csv" "$BASEDIR/results/made-alleles.csv"
else
    echo "Not overwriting alleles!"
fi;

make epitopes cleavages BASE_DIR="$BASEDIR/results" CONFIG="$BASEDIR/config.mak"
python tradeoff.py \
    "$BASEDIR/results/made-epitopes.csv" \
    "$BASEDIR/results/made-cleavages.csv" \
    "$BASEDIR/results/made-tradeoff.csv" \
    --log-file "$BASEDIR/results/made-tradeoff.log" \
    --pareto-steps 11 \
    --max-aminoacids 0 --max-epitopes 10 \
    --min-alleles 0 --min-proteins 0 \
    --cocktail 1 --verbose