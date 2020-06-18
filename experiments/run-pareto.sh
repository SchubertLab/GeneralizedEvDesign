#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    # pareto without spacers
    make pareto \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGES_NAME="cleavages-1500ig" \
        CLEAVAGES_OPTS="--top-immunogen 1500" \
        PARETO_NAME="made-tradeoff" \
        PARETO_OPTS="--pareto-steps 11 --max-aminoacids 0 --max-epitopes 10" \
        BASE_DIR="${dir%/}"

    # pareto with spacers
    make spacers \
         SPACERS_NAME="spacers-1500ig" \
         SPACERS_OPTS="-l 0 -l 4 --top-immunogen 1500" \
         BASE_DIR="${dir%/}"

    make pareto \
         CONFIG="experiments/base-config.mak" \
         CLEAVAGES_NAME="spacers-1500ig" \
         PARETO_NAME="made-spacers-tradeoff" \
         PARETO_OPTS="--pareto-steps 11 --max-aminoacids 0 --max-epitopes 10" \
         BASE_DIR="${dir%/}"
done
