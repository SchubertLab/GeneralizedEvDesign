#!/bin/bash
set -ex

BASEDIR="experiments"

for i in {1..5}
do
    echo "Bootstrap $i"
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
done