#!/bin/bash
set -ex

BASEDIR="experiments/advantage"

mkdir -p "$BASEDIR/results"
if [ ! -f "$BASEDIR/results/made-proteins.fasta" ]; then
    python data_preparation.py -v random-sequences "experiments/resources/hiv1-bc-env.fasta" 300 "$BASEDIR/results/made-proteins.fasta"
else
    echo "Not overwriting proteins!"
fi;

if [ ! -f "$BASEDIR/results/made-alleles.csv" ]; then
    cp "experiments/resources/alleles.csv" "$BASEDIR/results/made-alleles.csv"
else
    echo "Not overwriting alleles!"
fi;

for aminoacids in 90 180 360 720
do
    let epitopes=$aminoacids/9
    make mosaic optitope \
        BASE_DIR="$BASEDIR/results" \
        MOSAIC_VACCINE_NAME="made-mosaic-a$aminoacids.csv" \
        MOSAIC_LOG_NAME="made-mosaic-a$aminoacids.log" \
        MOSAIC_EVAL_NAME="made-mosaic-a$aminoacids-eval.csv" \
        MOSAIC_OPTS="--top-immunogen 1000 --top-proteins -1 --top-alleles -1 --max-epitopes -1 --max-aminoacids $aminoacids --min-alleles 0 --min-proteins 0 --cocktail 1" \
        OPTITOPE_VACCINE_NAME="made-optitope-a$aminoacids.csv" \
        OPTITOPE_LOG_NAME="made-optitope-a$aminoacids.log" \
        OPTITOPE_EVAL_NAME="made-optitope-a$aminoacids-eval.csv" \
        OPTITOPE_OPTS="--epitopes $epitopes --min-alleles 0 --min-proteins 0" \
        COVERAGE_OPTS="--max-edits 0 --top-n -1"
done

python evaluation.py -v aggregate "$BASEDIR/results/*-eval.csv" \
    "$BASEDIR/results/aggregate-evaluation.csv" \
    --path-format "$BASEDIR/results/made-(?P<vax>\w+)-a(?P<size>\d+)-eval.csv" \
    --summary-by size --summary-by vax --output-summary "$BASEDIR/results/evaluation-summary.csv"