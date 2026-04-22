#!/bin/bash

#BSUB -J IBCD
#BSUB -q dbeigpu
#BSUB -gpu "num=3"
#BSUB -n 3
#BSUB -o ibcd_output.txt
#BSUB -e ibcd_error.txt

INPUT_DIR="$HOME/input"

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

                DATA_FILE="$RUN_DIR/Y_matrix.csv"

                python ibcd.py \
                    --data "$DATA_FILE" \
                    --prior "$graph" \
                    --causal_order true \
                    --output_dir "$RUN_DIR"
            done
        done
    done
done

echo "All runs completed"