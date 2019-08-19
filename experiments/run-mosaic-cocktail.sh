make mosaic \
    BASE_DIR="experiments/results/hiv1bc-full" \
    CONFIG="experiments/base-config.mak" \
    MOSAIC_OPTS="--max-epitopes 0 --max-aminoacids 36 --top-proteins 500 --top-alleles 500 --top-immunogen 500 --cocktail 6 --min-alleles 0.8 --min-proteins 0.95" \
    MOSAIC_NAME="mosaic-high-coverage"