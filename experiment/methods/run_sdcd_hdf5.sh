#!/bin/bash

INPUT_DIR="$HOME/input"
OUTPUT_DIR="$HOME/output/sdcd"

for f in train train_fold1 train_fold2 train_fold3 train_fold4 train_fold5; do
    Y_FILE="$INPUT_DIR/Y_matrix_essential_${f}.csv"
    TARGETS_FILE="$INPUT_DIR/targets_essential_${f}.txt"
    OUT_FILE="$OUTPUT_DIR/G_sdcd_${f}.csv"

    python run_sdcd.py \
        --y_path "$Y_FILE" \
        --targets_path "$TARGETS_FILE" \
        --out_path "$OUT_FILE"
done

echo "All runs completed"