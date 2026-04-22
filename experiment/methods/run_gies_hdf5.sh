#!/bin/bash

INPUT_DIR="$HOME/input"
OUTPUT_DIR="$HOME/output/gies"

for f in train train_fold1 train_fold2 train_fold3 train_fold4 train_fold5; do
    Y_FILE="$INPUT_DIR/Y_matrix_essential_${f}.csv"
    TARGETS_FILE="$INPUT_DIR/targets_essential_${f}.txt"
    OUT_PREFIX="$OUTPUT_DIR/G_gies_${f}"

    Rscript run_gies.R \
        "$Y_FILE" \
        "$TARGETS_FILE" \
        "$OUT_PREFIX"
done

echo "All runs completed"