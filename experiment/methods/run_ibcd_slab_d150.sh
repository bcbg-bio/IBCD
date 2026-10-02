#!/bin/bash
# The slab width at D = 150, at target acceptance 0.9 and 0.8, SF and ER.
#
# Checks the candidate method -- the R_hat start plus the regularised
# horseshoe with its slab width from the EM (--slab_width em) -- where the
# published model already worked, and whether a target acceptance of 0.8 is
# safe with it. At D = 500 the slab width gave zero divergences at 0.9 with
# step sizes 4-5x larger than under the horseshoe, whose funnel is what made
# 0.8 fail before (stuck chains on ER).
#
#   accept  0.9 (the default) and 0.8
#   graph   SF and ER
#
# 2 x 2 x 5 seeds = 20 array elements, full 1000 warmup / 1000 samples. The
# comparisons are the D = 150 runs without the slab width: init_optimized_none
# (run_ibcd_init_compare.sh) and start_rhat_adam (run_ibcd_rhat_start.sh).
# Flags are passed explicitly rather than left to the defaults.
#
#   mkdir -p logs && bsub < run_ibcd_slab_d150.sh

#BSUB -J "IBCDslab150[1-20]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 06:00
#BSUB -o logs/ibcd_slab150_%I.out
#BSUB -e logs/ibcd_slab150_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=150
N_INT=100
SEEDS=(42 43 44 45 46)
ACCEPTS=(0.9 0.8)
GRAPHS=(sf er)

# 1..20 -> (graph, accept, seed); seed varies fastest, graph slowest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 5 ))]}
ACCEPT=${ACCEPTS[$(( (i / 5) % 2 ))]}
GRAPH=${GRAPHS[$(( i / 10 ))]}

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/slab_accept_${ACCEPT}/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: graph=$GRAPH accept=$ACCEPT seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --init_strategy rhat \
    --slab_width em \
    --target_accept_prob "$ACCEPT" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
