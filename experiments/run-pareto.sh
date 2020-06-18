#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make pareto \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGES_NAME="cleavages-1500ig" \
        CLEAVAGES_OPTS="--top-immunogen 1500" \
        PARETO_NAME="made-tradeoff" \
        PARETO_OPTS="--pareto-steps 11 --max-aminoacids 0 --max-epitopes 10" \
        BASE_DIR="${dir%/}"
done
