#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make CONFIG="experiments/base-config.mak" BASE_DIR="${dir%/}" POPCOVER_OPTS="--epitopes 10" popcover
done

python evaluation.py -v aggregate \
    'experiments/results/nef-300-*/made-popcover-vaccine-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/made-popcover-vaccine-evaluation.csv' \
    --output-summary 'experiments/results/popcover-evaluation-aggregate.csv'