#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    for aminoacids in 45 90 180 360 720
    do
        let epitopes=$aminoacids/9
        make mosaic optitope \
            BASE_DIR="${dir%/}" \
            MOSAIC_NAME="made-mosaic-a$aminoacids" \
            MOSAIC_OPTS="--top-immunogen -1 --top-proteins -1 --top-alleles -1 --max-epitopes -1 --max-aminoacids $aminoacids --min-alleles 0 --min-proteins 0 --cocktail 1 --min-overlap 2" \
            OPTITOPE_NAME="made-optitope-a$aminoacids" \
            OPTITOPE_OPTS="--epitopes $epitopes --min-alleles 0 --min-proteins 0" \
            STROBE_NAME="made-string-of-beads-a$aminoacids" \
            STROBE_OPTS="--max-aminoacids $aminoacids --max-epitopes 0 --min-alleles 0 --min-proteins 0 --cocktail 1" \
            CLEAVAGES_OPTS="--top-immunogen 1000" \
            COVERAGE_OPTS="--max-edits 0 --top-n -1"
    done
done

python evaluation.py -v aggregate "experiments/results/nef-300-*/made-*-a*-evaluation.csv" \
    --path-format "experiments/results/nef-300-\d+/made-(?P<vax>\w+)-a(?P<size>\d+)-evaluation.csv" \
    --summary-by size --summary-by vax --output-summary "experiments/results/advantage-evaluation-summary.csv"