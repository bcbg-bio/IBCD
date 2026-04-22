#!/bin/bash

#BSUB -J IBCD
#BSUB -q dbeigpu
#BSUB -gpu "num=3"
#BSUB -n 3
#BSUB -o ibcd_output.txt
#BSUB -e ibcd_error.txt

INPUT_DIR="$HOME/input"
OUTPUT_DIR="$HOME/output/ibcd"

for f in train train_fold1 train_fold2 train_fold3 train_fold4 train_fold5; do
    DATA_FILE="$INPUT_DIR/Y_matrix_essential_${f}.csv"
    OUT_DIR="$OUTPUT_DIR/${f}"

    python ibcd.py \
        --data "$DATA_FILE" \
        --prior sf \
        --causal_order false \
        --output_dir "$OUT_DIR"
done

echo "All runs completed"