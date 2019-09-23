#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make epitopes cleavages \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGES_OPTS="--top-immunogen 1500" \
        CLEAVAGES_NAME="cleavages-1500ig" \
        BASE_DIR="${dir%/}" \

    if [ ! -f "${dir}/made-tradeoff.csv" ]; then
    python tradeoff.py \
        "${dir}/made-epitopes.csv" \
        "${dir}/cleavages-1500ig.csv" \
        "${dir}/made-tradeoff.csv" \
        --log-file "${dir}/made-tradeoff.log" \
        --pareto-steps 11 \
        --max-aminoacids 0 --max-epitopes 10 \
        --min-alleles 0 --min-proteins 0 \
        --cocktail 1 --verbose
    fi
done
