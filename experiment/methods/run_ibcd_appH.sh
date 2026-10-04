#!/bin/bash
# Appendix H / Table 6: the Figure 2 sweep again with the soft acyclicity
# constraint, a Gaussian prior N(0, 0.5^2) on the spectral radius of G
# computed by power iteration, as on main. Everything else is the default, so
# each run differs from its Figure 2 counterpart only by the constraint.
#
# 4 dims x 2 graph families x 10 seeds = 80 array elements, same mapping and
# data as run_ibcd_dims_sweep.sh.
#
# Outputs: $PROJECT/IBCD_results/ibcd_runs/appH/<D>d/100/<graph>/<seed>
#
#   mkdir -p logs && bsub < run_ibcd_appH.sh

#BSUB -J "IBCDappH[1-80]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_appH_%I.out
#BSUB -e logs/ibcd_appH_%I.err

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
OUT_DIR="$OUT_ROOT/appH/${D}d/${N_INT}/${GRAPH}/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: D=$D graph=$GRAPH seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$GRAPH" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --rho_penalty gaussian \
    --rho_estimator power \
    --rho_sigma 0.5 \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
