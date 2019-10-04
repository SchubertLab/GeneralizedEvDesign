make mosaic \
    BASE_DIR="experiments/results/hiv1bc-full"\
    MOSAIC_NAME="mosaic-entropy-mid" \
    MOSAIC_OPTS="--min-overlap 8 --greedy-subtour --max-epitopes 0 --max-aminoacids 90 --top-immunogen 5000 --top-proteins 5000"

make mosaic \
    BASE_DIR="experiments/results/hiv1bc-full"\
    MOSAIC_NAME="mosaic-entropy-long" \
    MOSAIC_OPTS="--min-overlap 8 --greedy-subtour --max-epitopes 0 --max-aminoacids 180 --top-immunogen 5000 --top-proteins 5000"

make mosaic \
    BASE_DIR="experiments/results/hiv1bc-full"\
    MOSAIC_NAME="mosaic-entropy-short" \
    MOSAIC_OPTS="--min-overlap 8 --greedy-subtour --max-epitopes 0 --max-aminoacids 28 --top-immunogen 5000 --top-proteins 5000"

make optitope \
    BASE_DIR="experiments/results/hiv1bc-full"\
    OPTITOPE_NAME="optitope-entropy" \
    OPTITOPE_OPTS="--epitopes 20"
