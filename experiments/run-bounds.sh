#!/bin/bash
set -ex

# run this after mosaic with make, so the inputs are already prepared

function run_bound() {
    if [ ! -f "${1}/bound-$2-evaluation.csv" ]; then
        python mosaic_bounds.py $2 \
            "${1}/made-proteins.fasta" \
            "${1}/made-alleles.csv" \
            "${1}/made-epitopes.csv" \
            "${1}/made-overlaps.csv" \
            "${1}/bound-$2.csv" --min-overlap 8 -e 0 -a 90 --greedy-subtour \
            --verbose --log-file "${1}/bound-$2.log"

        python evaluation.py \
            --verbose --log-file "${1}/bound-$2-evaluation.log" \
            vaccine \
            "${1}/made-proteins.fasta" \
            "${1}/made-coverage.csv" \
            "${1}/made-alleles.csv" \
            "${1}/made-epitopes.csv" \
            "${1}/bound-$2.csv" \
            "${1}/bound-$2-evaluation.csv"
    fi
}

for dir in experiments/results/nef-300-*; do
    run_bound "${dir}" "immunogen"
    run_bound "${dir}" "coverage"
    run_bound "${dir}" "conservation"
done