#!/bin/bash
set -ex

for immunogen in "netmhcpan-rank" "pickpocket" "mhcflurry"; do
    for dir in experiments/results/nef-300-*; do
        make cleavages affinities epitopes \
            CONFIG="experiments/base-config.mak" \
            AFFINITIES_OPTS="--predictor $immunogen" \
            AFFINITIES_NAME="made-$immunogen-affinities" \
            CLEAVAGES_OPTS="--top-immunogen 1500" \
            CLEAVAGES_NAME="cleavages-$immunogen-1500ig" \
            EPITOPES_NAME="made-$immunogen-epitopes" \
            BASE_DIR="${dir%/}"

        if [ ! -f "${dir}/made-$immunogen-tradeoff.csv" ]; then
            python tradeoff.py \
                "${dir}/made-$immunogen-epitopes.csv" \
                "${dir}/cleavages-$immunogen-1500ig.csv" \
                "${dir}/made-$immunogen-tradeoff.csv" \
                --log-file "${dir}/made-$immunogen-tradeoff.log" \
                --pareto-steps 11 \
                --max-aminoacids 0 --max-epitopes 10 \
                --min-alleles 0 --min-proteins 0 \
                --cocktail 1 --verbose
        fi
    done
done
