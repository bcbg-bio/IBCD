#!/bin/bash
# D=50 sweep: exact inverse vs truncated series at K=6 and K=24.
#
# 3 arms x 2 graph families x 10 seeds = 60 array elements, one dataset each.
# Chains are drawn vectorized within a job: CUDA binds a process to a single
# MIG instance, so reserving more slices does not let numpyro draw chains in
# parallel. See run_ibcd_chains.sh for the one-chain-per-job layout, which is
# what the larger runs will need.
#
# Generate the data first (CPU only, a few minutes at D=50):
#
#   Rscript $HOME/projects/IBCD/experiment/simulation/generate_sim_data.R \
#       --out_dir $PROJECT/IBCD_results/sim_data \
#       --dims 50 --graphs er,sf --n_int 100 --seeds 42:51
#
# Then:  bsub < run_ibcd_d50_sweep.sh

#BSUB -J "IBCD50[1-60]%40"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 08:00
#BSUB -o logs/ibcd50_%I.out
#BSUB -e logs/ibcd50_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=50
N_INT=100
SEEDS=(42 43 44 45 46 47 48 49 50 51)
GRAPHS=(er sf)
ARMS=(inverse k6 k24)

# 1..60 -> (arm, graph, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
GRAPH=${GRAPHS[$(( (i / 10) % 2 ))]}
ARM=${ARMS[$(( i / 20 ))]}

case "$ARM" in
    inverse) ARM_FLAGS=() ;;                                     # exact (I-G)^-1
    k6)      ARM_FLAGS=(--truncated_series --series_order 6) ;;
    k24)     ARM_FLAGS=(--truncated_series --series_order 24) ;;
    *)       echo "unknown arm '$ARM'" >&2; exit 1 ;;
esac

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/${ARM}/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run generate_sim_data.R first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: arm=$ARM graph=$GRAPH seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"

# Confirm the MIG profile and memory ceiling actually handed to this process.
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --save_diagnostics \
    "${ARM_FLAGS[@]}"

echo "task ${LSB_JOBINDEX} finished with status $?"
