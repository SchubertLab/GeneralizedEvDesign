make CONFIG="experiments/base-config.mak" \
    BASE_DIR="experiments/results/hiv1bc-full" \
    proteins="experiments/resources/hiv1-bc-nef.fasta" \
    alleles="experiments/resources/alleles.csv" \
    OPTITOPE_OPTS="--min-alleles 0.99 --min-proteins 0.99 --epitopes 24" \
    MOSAIC_OPTS="--max-aminoacids 45 --max-epitopes 0 --cocktail 6 --top-immunogen 500 --top-proteins 500 --top-alleles 0 --min-proteins 0.99 --min-overlap 8 --greedy-subtour" \
    MOSAIC_NAME="mosaic-6cocktail" \
    optitope mosaic