#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    for overlaps in 0 1 2 3 4 5 6 7 8 9
    do
        make mosaic \
            BASE_DIR="${dir%/}" \
            MOSAIC_NAME="made-mosaic-o$overlaps" \
            MOSAIC_OPTS="--top-immunogen -1 --top-proteins -1 --top-alleles -1 --max-epitopes -1 --max-aminoacids 90 --min-alleles 0 --min-proteins 0 --cocktail 1 --min-overlap $overlaps" \
            COVERAGE_OPTS="--max-edits 0 --top-n -1"
    done
done

python evaluation.py -v aggregate "experiments/results/nef-300-1/made-*-o*-evaluation.csv" \
    --path-format "experiments/results/nef-300-1/made-(?P<vax>\w+)-o(?P<overlap>\d+)-evaluation.csv" \
    --summary-by size --summary-by vax --output-summary "experiments/results/overlaps-evaluation-summary.csv"