#!/bin/bash
# Starting point x soft rho constraint, D = 150, SF and ER.
#
#   init  median     init_to_median(num_samples=50): an essentially empty G,
#                    the default. Every chain started here at D = 150 SF left
#                    rho(G) < 1 within 20 warmup iterations.
#         optimized  the same start, then 2000 Adam steps on the log posterior
#                    before NUTS (posterior unchanged), each chain from its own
#                    start and jittered afterwards
#   pen   none       the model as published
#         soft       barrier on the Gelfand bound, zero up to rho = 1.0, one nat
#                    at 1.2 (the defaults of --rho_barrier_start/_width)
#
# 2 x 2 arms x 5 seeds x 2 graph families = 40 array elements. SF is 1-20,
# the family that fails from the default start; ER is 21-40, the control for
# whether either change disturbs a posterior that already works. The SF block
# keeps the indices it had before ER was added, so an SF-only submission of
# [1-20] is unchanged; submit [21-40] to add ER alone.
#
# median/none (SF 1-5, ER 21-25) is the current default, so it should
# reproduce the Figure 2 sweep's inverse/150d/100/<graph> runs byte-for-byte;
# SF seed 44 is the one that sweep lost to a cluster abort, and all ten ER
# runs there completed, so ER 21-25 can be skipped. Sampler settings are the
# defaults, as in the other sweeps.
#
#   mkdir -p logs && bsub < run_ibcd_init_compare.sh

#BSUB -J "IBCDinit2[1-40]%20"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 06:00
#BSUB -o logs/ibcd_initcmp_%I.out
#BSUB -e logs/ibcd_initcmp_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=150
N_INT=100
SEEDS=(42 43 44 45 46)
ARMS=(median_none median_soft optimized_none optimized_soft)
GRAPHS=(sf er)

# 1..40 -> (graph, arm, seed); seed varies fastest, graph slowest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 5 ))]}
ARM=${ARMS[$(( (i / 5) % 4 ))]}
GRAPH=${GRAPHS[$(( i / 20 ))]}

case "$ARM" in
    median_none)    ARM_FLAGS=(--init_strategy median) ;;
    median_soft)    ARM_FLAGS=(--init_strategy median    --rho_penalty barrier --rho_estimator gelfand) ;;
    optimized_none) ARM_FLAGS=(--init_strategy optimized) ;;
    optimized_soft) ARM_FLAGS=(--init_strategy optimized --rho_penalty barrier --rho_estimator gelfand) ;;
    *)              echo "unknown arm '$ARM'" >&2; exit 1 ;;
esac

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/init_${ARM}/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: arm=$ARM graph=$GRAPH seed=$SEED"
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
