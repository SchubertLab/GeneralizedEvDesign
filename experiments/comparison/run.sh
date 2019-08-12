#!/bin/bash
set -ex

BASEDIR="experiments/comparison"

for i in {1..10}
do
    echo "Trial $i"
    RESDIR="$BASEDIR/results/nef-300-$i"
    mkdir -p $RESDIR
    if [ ! -f "$RESDIR/made-proteins.fasta" ]; then
        python data_preparation.py -v random-sequences "experiments/resources/hiv1-bc-nef.fasta" 300 "$RESDIR/made-proteins.fasta"
    else
        echo "Not overwriting proteins!"
    fi;

    if [ ! -f "$RESDIR/made-alleles.csv" ]; then
        cp "experiments/resources/alleles.csv" "$RESDIR/made-alleles.csv"
    else
        echo "Not overwriting alleles!"
    fi;

    make CONFIG="$BASEDIR/config.mak" BASE_DIR=$RESDIR mosaic optitope popcover
done

python evaluation.py -v aggregate "$BASEDIR/results/nef-300-*/*-evaluation.csv" \
    "$BASEDIR/results/aggregate-evaluation.csv" \
    --path-format "$BASEDIR/results/nef-300-(?P<rep>\d+)/made-(?P<vax>[a-z0-9\-]+)-evaluation.csv" \
    --summary-by vax --output-summary "$BASEDIR/results/evaluation-summary.csv"