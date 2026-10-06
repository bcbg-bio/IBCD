#!/bin/bash
# Figure 2: D in {50, 150, 250, 500}, exact inverse.
#
# 4 dims x 2 graph families x 10 seeds = 80 array elements.
#
# Sampler settings are the defaults as of this commit: target_accept_prob 0.9,
# num_warmup 1000, max_tree_depth 12, 3 vectorized chains. They are left
# implicit so one change of default moves every sweep together; each run
# records what it actually used in diagnostics.json.
#
# Generate the data first. --write_iv false if not running a method that
# consumes the IVs directly.
#
#   Rscript $HOME/projects/IBCD/experiment/simulation/generate_sim_data.R \
#       --out_dir $PROJECT/IBCD_results/sim_data \
#       --dims 50,150,250,500 --graphs er,sf --n_int 100 --seeds 42:51 \
#       --write_iv false
#
# Disk: at D=500 each Y_matrix.csv and Y_with_targets.csv is about 1 GB and
# each G_draws.npy about 3 GB, so that tier alone needs roughly 100 GB.
#
# SF_ANCHOR selects the SF prior's sparsity level (ibcd.py --sf_anchor): em,
# the default, or none, the published max-normalisation. none writes to a
# separate directory (inverse_sfnone/). It has no effect on ER tasks, so
# submit only the SF ones, e.g.
#
#   SF_ANCHOR=none bsub -J "IBCDdim[11-20,31-40,51-60,71-80]" < run_ibcd_dims_sweep.sh
#
#   mkdir -p logs && bsub < run_ibcd_dims_sweep.sh

#BSUB -J "IBCDdim[1-80]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_dim_%I.out
#BSUB -e logs/ibcd_dim_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs"

N_INT=100
SEEDS=(42 43 44 45 46 47 48 49 50 51)
GRAPHS=(er sf)
DIMS=(50 150 250 500)

# 1..80 -> (dim, graph, seed); seed varies fastest
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
GRAPH=${GRAPHS[$(( (i / 10) % 2 ))]}
D=${DIMS[$(( i / 20 ))]}

DATA_FILE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/Y_with_targets.csv"
SF_ANCHOR="${SF_ANCHOR:-em}"
RUN_DIR=inverse
[[ $SF_ANCHOR == none ]] && RUN_DIR=inverse_sfnone
OUT_DIR="$OUT_ROOT/$RUN_DIR/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run generate_sim_data.R first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX} (sf_anchor $SF_ANCHOR): D=$D graph=$GRAPH seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --sf_anchor "$SF_ANCHOR" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
