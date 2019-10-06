#!/bin/bash
set -ex

# run this after mosaic with make, so the inputs are already prepared

function run_bound() {
    basename="ig-bound-$2-a$3"
    if [ ! -f "${1}/$basename.csv" ]; then
        python mosaic_bounds.py $2 \
            "${1}/made-proteins.fasta" \
            "${1}/made-alleles.csv" \
            "${1}/made-epitopes.csv" \
            "${1}/made-overlaps.csv" \
            "${1}/$basename.csv" --min-overlap 8 -e 0 -a $3 --greedy-subtour \
            --verbose --log-file "${1}/$basename.log"

        python evaluation.py \
            --verbose --log-file "${1}/$basename-evaluation.log" \
            vaccine \
            "${1}/made-proteins.fasta" \
            "${1}/made-coverage.csv" \
            "${1}/made-alleles.csv" \
            "${1}/made-epitopes.csv" \
            "${1}/$basename.csv" \
            "${1}/$basename-evaluation.csv"
    fi
}

for dir in experiments/results/nef-300-*; do
    for aminoacids in 45 72 90 135
    do
        run_bound "${dir}" "immunogen" $aminoacids
        run_bound "${dir}" "coverage" $aminoacids
        run_bound "${dir}" "conservation" $aminoacids
    done
done

python evaluation.py aggregate \
    "experiments/results/nef-300-*/ig-bound-*-evaluation.csv" \
    --path-format "experiments/results/nef-300-(?P<rep>\d+)/ig-bound-(?P<vax>\w+)-a(?P<aa>\d+)-evaluation.csv" \
    -s aa -s vax -S experiments/results/ig-bound-evaluation-aggregate.csv
