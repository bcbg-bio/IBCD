#!/bin/bash
# SF sparsity anchor x target acceptance, D = 150, with the optimised start.
#
#   anchor  em    global spike weight from the EM of eq 12 (the default)
#           none  the published max-normalisation, --sf_anchor none
#   accept  0.9   the default
#           0.8   numpyro's own default; larger steps, so shorter trajectories,
#                 at the cost of more divergences
#
# The em / 0.9 cell, for SF and for ER, is init_optimized_none in the
# init-strategy comparison (run_ibcd_init_compare.sh), which these defaults
# reproduce exactly, so it is not rerun here. The anchor only exists on the
# SF path, so ER only needs the acceptance arm.
#
#   tasks  1-15  SF: none/0.9, em/0.8, none/0.8  x 5 seeds
#   tasks 16-20  ER: 0.8                          x 5 seeds
#
# The start is left to the default (optimized), as are the other sampler
# settings.
#
#   mkdir -p logs && bsub < run_ibcd_anchor_accept.sh

#BSUB -J "IBCDanch[1-20]%20"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 06:00
#BSUB -o logs/ibcd_anchacc_%I.out
#BSUB -e logs/ibcd_anchacc_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=150
N_INT=100
SEEDS=(42 43 44 45 46)
SF_ANCHORS=(none em none)
SF_ACCEPTS=(0.9 0.8 0.8)

# 1..20 -> (graph, anchor, accept, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 5 ))]}
if (( i < 15 )); then
    GRAPH=sf
    ANCHOR=${SF_ANCHORS[$(( i / 5 ))]}
    ACCEPT=${SF_ACCEPTS[$(( i / 5 ))]}
    ARM_FLAGS=(--sf_anchor "$ANCHOR" --target_accept_prob "$ACCEPT")
    OUT_DIR="$OUT_ROOT/anchor_${ANCHOR}_accept_${ACCEPT}/${D}d/${N_INT}/${GRAPH}/${SEED}"
else
    GRAPH=er
    ACCEPT=0.8
    ARM_FLAGS=(--target_accept_prob "$ACCEPT")
    OUT_DIR="$OUT_ROOT/accept_${ACCEPT}/${D}d/${N_INT}/${GRAPH}/${SEED}"
fi

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: graph=$GRAPH seed=$SEED flags=${ARM_FLAGS[*]}"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --save_diagnostics \
    "${ARM_FLAGS[@]}"

echo "task ${LSB_JOBINDEX} finished with status $?"
