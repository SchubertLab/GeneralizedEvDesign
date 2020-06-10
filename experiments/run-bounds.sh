#!/bin/bash
set -ex

function run_bound() {
    # arguments: <base-dir> <optimize-for> <length-aminoacids> <immunogen-predictor>

    # create necessary inputs
    make overlaps affinities epitopes \
         CONFIG="experiments/base-config.mak" \
         AFFINITIES_OPTS="--predictor $4" \
         AFFINITIES_NAME="made-$4-affinities" \
         EPITOPES_NAME="made-$4-epitopes" \
         BASE_DIR="${1%/}"

    basename="ig-bound-$2-a$3-$4"
    if [ ! -f "${1}/$basename.csv" ]; then
        python mosaic_bounds.py $2 \
            "${1}/made-proteins.fasta" \
            "${1}/made-alleles.csv" \
            "${1}/made-$4-epitopes.csv" \
            "${1}/made-overlaps.csv" \
            "${1}/$basename.csv" \
            --min-overlap 8 -e 0 -a $3 --greedy-subtour \
            --verbose --log-file "${1}/$basename.log"
    fi

    if [ ! -f "${1}/$basename-evaluation.csv" ]; then
        python evaluation.py \
            --verbose --log-file "${1}/$basename-evaluation.log" \
            vaccine \
            "${1}/made-proteins.fasta" \
            "${1}/made-coverage.csv" \
            "${1}/made-alleles.csv" \
            "${1}/made-$4-epitopes.csv" \
            "${1}/$basename.csv" \
            "${1}/$basename-evaluation.csv"
    fi
}

for ig_pred in "netmhcpan-rank" "pickpocket" "mhcflurry" "netmhcpan"
do
    for dir in experiments/results/nef-300-*
    do
        for aminoacids in 45 72 90 135
        do
            run_bound "${dir}" "immunogen" $aminoacids $ig_pred
            run_bound "${dir}" "coverage" $aminoacids $ig_pred
            run_bound "${dir}" "conservation" $aminoacids $ig_pred
        done
    done

    python evaluation.py aggregate \
           "experiments/results/nef-300-*/ig-bound-*-evaluation.csv" \
           --path-format "experiments/results/nef-300-(?P<rep>\d+)/ig-bound-(?P<vax>\w+)-a(?P<aa>\d+)-$ig_pred-evaluation.csv" \
           -s aa -s vax -S experiments/results/ig-bound-evaluation-aggregate.csv
done

