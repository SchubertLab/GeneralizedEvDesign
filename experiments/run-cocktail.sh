mkdir -p experiments/results/hiv1bc-full

make CONFIG="experiments/base-config.mak" \
    BASE_DIR="experiments/results/hiv1bc-full" \
    proteins="experiments/resources/hiv1-bc-nef.fasta" \
    alleles="experiments/resources/alleles.csv" \
    OPTITOPE_OPTS="--min-alleles 0.99 --min-proteins 0.99 --epitopes 24" \
    MOSAIC_OPTS="--max-aminoacids 54 --max-epitopes 0 --cocktail 4 --top-immunogen 0 --top-proteins 25000 --top-alleles 0 --min-proteins 0.99 --min-overlap 8 --greedy-subtour" \
    optitope mosaic