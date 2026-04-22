#!/bin/bash

#INPUT_DIR="$HOME/input"

for dim in 50d 150d 250d 500d; do
    if [ "$dim" = "50d" ]; then
        cov_list="5 15 25 50 75 100"
    else
        cov_list="none"
    fi

    for cov in $cov_list; do
        for graph in er sf; do
            for seed in {42..51}; do

                if [ "$dim" = "50d" ]; then
                    RUN_DIR="$INPUT_DIR/$dim/$cov/$graph/$seed"
                else
                    RUN_DIR="$INPUT_DIR/$dim/$graph/$seed"
                fi

                Y_FILE="$RUN_DIR/Y_matrix.csv"
                T_FILE="$RUN_DIR/targets.txt"
                OUT_FILE="$RUN_DIR/G_igsp.csv"

                python run_igsp.py \
                    --y_path "$Y_FILE" \
                    --t_path "$T_FILE" \
                    --out_path "$OUT_FILE"
            done
        done
    done
done

echo "All runs completed"