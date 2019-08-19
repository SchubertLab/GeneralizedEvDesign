#!/bin/bash
set -ex

for dir in experiments/results/nef-300-*; do
    make CONFIG="experiments/base-config.mak" BASE_DIR="${dir%/}" POPCOVER_OPTS="--epitopes 10" popcover
done