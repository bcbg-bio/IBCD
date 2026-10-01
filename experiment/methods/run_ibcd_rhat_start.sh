#!/bin/bash
# The R_hat start: NUTS begun from the total effects soft-thresholded at 3 SE.
#
#   rhat_adam  --init_strategy rhat: sign(R_hat) max(|R_hat| - 3 SE, 0), then
#              2000 Adam steps on the log posterior, then NUTS
#   rhat_raw   the same start with --init_opt_steps 0: NUTS begins at the
#              thresholded R_hat itself (jittered)
#
# At D = 500 the empty start, even with the Adam descent, lands on a wrong
# graph 3,000-5,500 nats below the optimum near the true G; the R_hat start
# lands within ~650-900 nats of it with valid rho. D = 150 SF checks the new
# start against the current default (optimized), and D = 150 ER that it does
# not disturb a family that already works; their optimized-start reference
# runs are init_optimized_none in run_ibcd_init_compare.sh.
#
#   tasks  1-10  SF D = 500   rhat_adam, rhat_raw  x 5 seeds
#   tasks 11-20  SF D = 150   rhat_adam, rhat_raw  x 5 seeds
#   tasks 21-30  ER D = 150   rhat_adam, rhat_raw  x 5 seeds
#
# D = 500 first, since those take ~8 h. One array has to request the D = 500
# resources for every task (40 GB, 24 h); submit [11-30] separately with
# -gpu "num=1:gmem=20G" -W 06:00 if the larger request slows scheduling.
#
#   mkdir -p logs && bsub < run_ibcd_rhat_start.sh

#BSUB -J "IBCDrhat[1-30]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_rhat_%I.out
#BSUB -e logs/ibcd_rhat_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

N_INT=100
SEEDS=(42 43 44 45 46)
ARMS=(rhat_adam rhat_raw)
BLOCK_D=(500 150 150)
BLOCK_GRAPH=(sf sf er)

# 1..30 -> (block, arm, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 5 ))]}
ARM=${ARMS[$(( (i / 5) % 2 ))]}
D=${BLOCK_D[$(( i / 10 ))]}
GRAPH=${BLOCK_GRAPH[$(( i / 10 ))]}

case "$ARM" in
    rhat_adam) ARM_FLAGS=(--init_strategy rhat) ;;
    rhat_raw)  ARM_FLAGS=(--init_strategy rhat --init_opt_steps 0) ;;
    *)         echo "unknown arm '$ARM'" >&2; exit 1 ;;
esac

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/start_${ARM}/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: arm=$ARM D=$D graph=$GRAPH seed=$SEED"
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
