#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make popcover \
        CONFIG="experiments/base-config.mak" \
        BASE_DIR="${dir%/}" \
        POPCOVER_OPTS="--epitopes 10" \
        POPCOVER_NAME="popcover-high-coverage"

    conservation=$(csvcut ${dir}/made-popcover-vaccine-evaluation.csv -c conservation | tail -n 1)
    min_proteins=$(csvcut ${dir}/made-popcover-vaccine-evaluation.csv -c prot_coverage | tail -n 1)
    alleles=$(csvcut ${dir}/made-popcover-vaccine-evaluation.csv -c alleles | tail -n 1)

    make string-of-beads \
        BASE_DIR="${dir%/}" \
        CONFIG="experiments/base-config.mak" \
        CLEAVAGES_OPTS="--top-proteins 500 --top-alleles 500 --top-immunogen 500 --penalty 0.1 --cleavage-window 5 --cleavage-model PCM" \
        STROBE_OPTS="--max-epitopes 10 --max-aminoacids 0 --cocktail 1 --min-alleles $alleles --min-proteins $min_proteins --min-avg-prot-conservation $conservation" \
        STROBE_NAME="strobe-high-coverage"
done

python evaluation.py -v aggregate \
    'experiments/results/nef-*/*-high-coverage-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/(?P<vax>\w+)-high-coverage-evaluation.csv' \
    --summary-by vax \
    --output-summary 'experiments/results/high-coverage-evaluation-aggregate.csv'