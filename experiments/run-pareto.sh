#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make epitopes cleavages \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGEs_OPTS="--top-immunogen 1500" \
        BASE_DIR="${dir%/}" \

    python tradeoff.py \
        "${dir}/made-epitopes.csv" \
        "${dir}/made-cleavages.csv" \
        "${dir}/made-tradeoff.csv" \
        --log-file "${dir}/made-tradeoff.log" \
        --pareto-steps 11 \
        --max-aminoacids 0 --max-epitopes 10 \
        --min-alleles 0 --min-proteins 0 \
        --cocktail 1 --verbose
done