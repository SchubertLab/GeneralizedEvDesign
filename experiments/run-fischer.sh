#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do

    # run their online tool on each sample and put the result in fischer.fasta
    python parse_fischer.py -v -l "${dir}/mosaic-fischer-vaccine.log" \
        "${dir}/made-coverage.csv" "${dir}/fischer.fasta" "${dir}/mosaic-fischer-vaccine.csv"
    
    python evaluation.py vaccine \
        "${dir}/made-proteins.fasta" \
        "${dir}/made-coverage.csv" \
        "${dir}/made-alleles.csv" \
        "${dir}/made-epitopes.csv" \
        "${dir}/mosaic-fischer-vaccine.csv" \
        "${dir}/mosaic-fischer-vaccine-evaluation.csv"

    let aminoacids=$(wc -l ${dir}/mosaic-fischer-vaccine.csv | cut -f 1 -d' ')-1+8
    conservation=$(csvcut ${dir}/mosaic-fischer-vaccine-evaluation.csv -c conservation | tail -n 1)
    min_proteins=$(csvcut ${dir}/mosaic-fischer-vaccine-evaluation.csv -c prot_coverage | tail -n 1)

    make mosaic \
        BASE_DIR="${dir%/}" \
        MOSAIC_OPTS="--min-overlap 8 --max-epitopes 0 --max-aminoacids $aminoacids --greedy-subtour --min-proteins $min_proteins --min-avg-prot-conservation $conservation" \
        MOSAIC_NAME="mosaic-genev-vaccine"

done

python evaluation.py -v aggregate \
    'experiments/results/nef-*/mosaic-*-vaccine-evaluation.csv' \
    --path-format 'experiments/results/nef-300-(?P<rep>\d+)/mosaic-(?P<vax>\w+)-vaccine-evaluation.csv' \
    --summary-by vax \
    --output-summary 'experiments/results/mosaic-fischer-evaluation-aggregate.csv'