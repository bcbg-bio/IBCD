#!/bin/bash

INPUT_DIR="$HOME/input"

for cov in 5 15 25 50 75 100; do
    for graph in er sf; do
        for seed in {42..51}; do
            RUN_DIR="$INPUT_DIR/50d/$cov/$graph/$seed"
            Y_FILE="$RUN_DIR/Y_matrix.csv"

            julia run_llcb.jl "$Y_FILE"
        done
    done
done

echo "All runs completed"