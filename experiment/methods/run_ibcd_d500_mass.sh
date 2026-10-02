#!/bin/bash
# Where is the posterior mass at D = 500? SF, 3 seeds, Adam -> NUTS throughout.
#
# Every D = 500 chain started from the R_hat start began inside rho(G) < 1
# (rho 0.52-0.77) and left it during warmup; all sampling draws had
# rho >= 143. That points at the posterior itself rather than the start.
#
#   truth         --init_strategy file --init_G <true G>, current prior.
#                 If these chains also leave rho < 1, no start will fix it
#                 and the model has to change. If they stay, the problem is
#                 in warmup.
#   rhat_barrier  the R_hat start + the soft barrier on rho (zero below 1.0,
#                 one nat at 1.2), Gelfand estimate
#   rhat_slab     the R_hat start + the regularised horseshoe, each entry's
#                 slab scale capped at the EM slab scale (--slab_width em).
#                 Removes the Cauchy tail (a prior draw at D = 500 has ~350
#                 entries with |G| > 1 under the horseshoe, none here) while
#                 keeping shrinkage proportional to 1 - pi0 below the cap.
#
# All arms use the default Adam phase (lr 0.01, 2000 steps) before NUTS.
#
# 3 arms x 3 seeds = 9 array elements. Shortened to 700 warmup / 300
# samples, enough to tell "every draw rejected" from "it stayed"; a chain
# leaving only after iteration 1000 would be missed, but every escape seen
# so far happened early in warmup.
#
#   mkdir -p logs && bsub < run_ibcd_d500_mass.sh

#BSUB -J "IBCDmass[1-9]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 12:00
#BSUB -o logs/ibcd_mass_%I.out
#BSUB -e logs/ibcd_mass_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=500
N_INT=100
GRAPH=sf
SEEDS=(42 43 44)
ARMS=(truth rhat_barrier rhat_slab)

# 1..9 -> (arm, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 3 ))]}
ARM=${ARMS[$(( i / 3 ))]}

SIM_DIR="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}"
case "$ARM" in
    truth)        ARM_FLAGS=(--init_strategy file --init_G "$SIM_DIR/G_matrix.csv") ;;
    rhat_barrier) ARM_FLAGS=(--init_strategy rhat --rho_penalty barrier --rho_estimator gelfand) ;;
    rhat_slab)    ARM_FLAGS=(--init_strategy rhat --slab_width em) ;;
    *)            echo "unknown arm '$ARM'" >&2; exit 1 ;;
esac

DATA_FILE="$SIM_DIR/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/d500_${ARM}/${D}d/${N_INT}/${GRAPH}/${SEED}"

for f in "$DATA_FILE" "$SIM_DIR/G_matrix.csv"; do
    if [[ ! -f "$f" ]]; then echo "missing input: $f" >&2; exit 1; fi
done
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: arm=$ARM seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --num_warmup 700 \
    --num_samples 300 \
    --save_diagnostics \
    "${ARM_FLAGS[@]}"

echo "task ${LSB_JOBINDEX} finished with status $?"
