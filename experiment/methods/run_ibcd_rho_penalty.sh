#!/bin/bash
# Constraining rho(G) at D = 150: three versions compared.
#
#   appH             Appendix H as written: Gaussian prior N(0, 0.5^2) on rho,
#                    rho from 50 steps of power iteration (port of 42660bd)
#   gelfand_gauss    the same Gaussian prior, rho from the Gelfand upper bound
#   gelfand_barrier  a barrier that is zero below rho = 0.9, Gelfand bound
#
# 3 arms x 2 graph families x 5 seeds = 30 array elements. SF comes first
# (1-15), since it is the family that fails without a constraint; ER (16-30)
# is the control for whether the constraint disturbs a posterior that already
# works. Submit [1-15] instead of [1-30] to run SF only.
#
# The unconstrained baseline for the same datasets is the Figure 2 sweep's
# inverse/150d/100/<graph>/<seed>, so it is not rerun here. Everything else
# is left at the defaults, as in the other sweeps.
#
#   mkdir -p logs && bsub < run_ibcd_rho_penalty.sh

#BSUB -J "IBCDrho[1-30]%15"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 06:00
#BSUB -o logs/ibcd_rho_%I.out
#BSUB -e logs/ibcd_rho_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=150
N_INT=100
SEEDS=(42 43 44 45 46)
ARMS=(appH gelfand_gauss gelfand_barrier)
GRAPHS=(sf er)

# 1..30 -> (graph, arm, seed); seed varies fastest, graph slowest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 5 ))]}
ARM=${ARMS[$(( (i / 5) % 3 ))]}
GRAPH=${GRAPHS[$(( i / 15 ))]}

case "$ARM" in
    appH)            ARM_FLAGS=(--rho_penalty gaussian --rho_estimator power) ;;
    gelfand_gauss)   ARM_FLAGS=(--rho_penalty gaussian --rho_estimator gelfand) ;;
    gelfand_barrier) ARM_FLAGS=(--rho_penalty barrier  --rho_estimator gelfand) ;;
    *)               echo "unknown arm '$ARM'" >&2; exit 1 ;;
esac

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/rho_${ARM}/${D}d/${N_INT}/${GRAPH}/${SEED}"

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
