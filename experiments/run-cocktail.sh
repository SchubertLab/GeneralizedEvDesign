make CONFIG="experiments/base-config.mak" \
    BASE_DIR="experiments/results/hiv1bc-full" \
    proteins="experiments/resources/hiv1-bc-nef.fasta" \
    alleles="experiments/resources/alleles.csv" \
    OPTITOPE_OPTS="--min-alleles 0.99 --min-proteins 0.99 --epitopes 24" \
    OPTITOPE_NAME="optitope-p99" \
    MOSAIC_OPTS="--max-aminoacids 54 --max-epitopes 0 --cocktail 4 --top-immunogen 1000 --top-proteins 1000 --top-alleles 0 --min-proteins 0.99 --min-overlap 8 --greedy-subtour" \
    MOSAIC_NAME="mosaic-4cocktail" \
    MOSAIC_EVAL_LOG_NAME="mosaic-4cocktail-evaluation" \
    optitope mosaic
