#!/bin/bash
set -ex

# PREREQUISITE: run their online tool on each sample and put the result in fischer.fasta

for dir in experiments/results/nef-300-*; do
    # parse fischer's vaccine into our csv format
    python parse_fischer.py -v -l "${dir}/mosaic-fischer-vaccine.log" \
        "${dir}/made-coverage.csv" "${dir}/fischer.fasta" "${dir}/mosaic-fischer-vaccine.csv"

    # evaluate fischer's vaccine
    python evaluation.py vaccine \
        "${dir}/made-proteins.fasta" \
        "${dir}/made-coverage.csv" \
        "${dir}/made-alleles.csv" \
        "${dir}/made-epitopes.csv" \
        "${dir}/mosaic-fischer-vaccine.csv" \
        "${dir}/mosaic-fischer-vaccine-evaluation.csv"

    # also evaluate on the full set
    python evaluation.py vaccine \
        "experiments/results/hiv1bc-full/made-proteins.fasta" \
        "experiments/results/hiv1bc-full/made-coverage.csv" \
        "experiments/results/hiv1bc-full/made-alleles.csv" \
        "experiments/results/hiv1bc-full/made-epitopes.csv" \
        "${dir}/mosaic-fischer-vaccine.csv" \
        "${dir}/mosaic-fischer-val-vaccine-evaluation.csv"

    # extract properties
    let aminoacids=$(wc -l ${dir}/mosaic-fischer-vaccine.csv | cut -f 1 -d' ')-1+8
    conservation=$(csvcut ${dir}/mosaic-fischer-vaccine-evaluation.csv -c conservation | tail -n 1)
    min_proteins=$(csvcut ${dir}/mosaic-fischer-vaccine-evaluation.csv -c prot_coverage | tail -n 1)
    alleles=$(csvcut ${dir}/mosaic-fischer-vaccine-evaluation.csv -c alleles | tail -n 1)

    # run mosaic to match
    make mosaic \
        BASE_DIR="${dir%/}" \
        MOSAIC_OPTS="--min-overlap 8 --max-epitopes 0 --max-aminoacids $aminoacids --greedy-subtour --min-proteins $min_proteins --min-avg-prot-conservation $conservation --min-alleles $alleles" \
        MOSAIC_NAME="mosaic-genev-vaccine"

    # also evaluate on the full set
    make mosaic \
        BASE_DIR="experiments/results/hiv1bc-full" \
        MOSAIC_VACCINE="${dir}/mosaic-genev-vaccine.csv" \
        MOSAIC_EVAL="${dir}/mosaic-genev-val-vaccine-evaluation.csv"

done

python evaluation.py -v aggregate \
    'experiments/results/nef-*/mosaic-*-vaccine-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/mosaic-(?P<vax>[a-z\-]+)-vaccine-evaluation.csv' \
    --summary-by vax \
    --output-summary 'experiments/results/mosaic-fischer-evaluation-aggregate.csv'
