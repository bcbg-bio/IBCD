#!/bin/bash

INPUT_DIR="$HOME/input"
OUTPUT_DIR="$HOME/output/igsp"

for f in train train_fold1 train_fold2 train_fold3 train_fold4 train_fold5; do
    Y_FILE="$INPUT_DIR/Y_matrix_essential_${f}.csv"
    TARGETS_FILE="$INPUT_DIR/targets_essential_${f}.txt"
    OUT_FILE="$OUTPUT_DIR/G_igsp_${f}.csv"

    python run_igsp.py \
        --y_path "$Y_FILE"  \
        --t_path "$T_FILE" \
        --out_path "$OUT_FILE"
done

echo "All runs completed"