#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    # ignore errors because optitope fails since the problem is infeasible
    make mosaic optitope --ignore-errors \
        BASE_DIR="${dir%/}" \
        CONFIG="experiments/base-config.mak" \
        MOSAIC_OPTS="--max-epitopes 0 --max-aminoacids 135 --top-proteins 0 --top-alleles 0 --top-immunogen 0 --cocktail 1 --min-alleles 0.999 --min-proteins 0.99 --min-overlap 1 --greedy-subtour" \
        MOSAIC_NAME="mosaic-highest-coverage" \
        OPTITOPE_OPTS="--epitopes 15 --min-alleles 0.999 --min-proteins 0.99" \
        OPTITOPE_NAME="optitope-highest-coverage"
done

python evaluation.py -v aggregate \
    'experiments/results/nef-300-*/*-highest-coverage-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/(P<vax>\w+)-highest-coverage-evaluation.csv' \
    --summary-by vax \
    --output-summary 'experiments/results/mosaic-highest-coverage-evaluation-aggregate.csv'