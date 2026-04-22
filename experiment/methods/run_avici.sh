#!/bin/bash

INPUT_DIR="$HOME/input"

for cov in 5 15 25 50; do
    for graph in er sf; do
        for seed in {42..51}; do
            RUN_DIR="$INPUT_DIR/50d/$cov/$graph/$seed"
            Y_FILE="$RUN_DIR/Y_matrix.csv"
            OUT_FILE="$RUN_DIR/G_avici.csv"

            python run_avici.py \
                --y_path "$Y_FILE" \
                --out_path "$OUT_FILE"
        done
    done
done

echo "All runs completed"