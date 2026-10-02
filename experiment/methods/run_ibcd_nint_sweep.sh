#!/bin/bash
# Figure 3: intervention sample size at D = 50, exact inverse.
#
# n_int in {5, 15, 25, 50, 75} x 2 graph families x 10 seeds = 100 array
# elements. Control sample size tracks the interventional one: the generator
# sets N_cont = D * n_int, so each dataset has 2 * D * n_int rows in total.
#
# Sampler settings are the defaults as of this commit: target_accept_prob 0.9,
# num_warmup 1000, max_tree_depth 12, 3 vectorized chains. Left implicit so
# every sweep moves together if a default changes; each run records what it
# actually used in diagnostics.json.
#
#
# Generate the data first:
#
#   Rscript $HOME/projects/IBCD/experiment/simulation/generate_sim_data.R \
#       --out_dir $PROJECT/IBCD_results/sim_data \
#       --dims 50 --graphs er,sf --n_int 5,15,25,50,75 --seeds 42:51 \
#       --write_iv false
#
#   mkdir -p logs && bsub < run_ibcd_nint_sweep.sh

#BSUB -J "IBCDnint[1-100]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 08:00
#BSUB -o logs/ibcd_nint_%I.out
#BSUB -e logs/ibcd_nint_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

D=50
SEEDS=(42 43 44 45 46 47 48 49 50 51)
GRAPHS=(er sf)
N_INTS=(5 15 25 50 75)

# 1..100 -> (n_int, graph, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
GRAPH=${GRAPHS[$(( (i / 10) % 2 ))]}
N_INT=${N_INTS[$(( i / 20 ))]}

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/inverse/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run generate_sim_data.R first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: n_int=$N_INT graph=$GRAPH seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
