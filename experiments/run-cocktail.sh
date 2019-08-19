mkdir -p experiments/results/hiv1bc-full

make CONFIG="experiments/base-config.mak" \
    BASE_DIR="experiments/results/hiv1bc-full" \
    proteins="experiments/resources/hiv1-bc-nef.fasta" \
    alleles="experiments/resources/alleles.csv" \
    OPTITOPE_OPTS="--min-alleles 0.99 --min-proteins 0.99 --epitopes 24" \
    MOSAIC_OPTS="--max-aminoacids 36 --max-epitopes 0 --cocktail 6 --top-immunogen 400 --top-proteins 500 --top-alleles 300 --min-alleles 0.99 --min-proteins 0.99" \
    optitope mosaic
