#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make string-of-beads optitope \
        BASE_DIR="${dir%/}" \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGES_OPTS="--top-proteins 500 --top-alleles 500 --top-immunogen 500 --penalty 0.1 --cleavage-window 5 --cleavage-model PCM" \
        STROBE_OPTS="--max-epitopes 10 --max-aminoacids 0 --cocktail 1 --min-alleles 0.8 --min-proteins 0.95 --greedy-subtour" \
        STROBE_NAME="strobe-high-coverage" \
        OPTITOPE_OPTS="--epitopes 10 --min-alleles 0.8 --min-proteins 0.95" \
        OPTITOPE_NAME="optitope-high-coverage"

done

python evaluation.py -v aggregate \
    'experiments/results/nef-*/*-high-coverage-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/(?P<vax>\w+)-high-coverage-evaluation.csv' \
    --summary-by vax \
    --output-summary 'experiments/results/high-coverage-evaluation-aggregate.csv'