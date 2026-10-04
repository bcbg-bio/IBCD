#!/bin/bash
# SF graphs with reduced R^2-sortability (Appendix F.2): Table 3 ('full') and
# Table 4 ('thresholded', edges below 0.125 removed), on the data from
# experiment/simulation/run_generate_extra.sh. SF prior, default settings.
#
#   2 versions x 2 dims (50, 250) x 10 seeds = 40 array elements;
#   seed fastest, then D, then version.
#
# Outputs: $PROJECT/IBCD_results/ibcd_runs/r2/<version>/<D>d/sf/<seed>
#
#   mkdir -p logs && bsub < run_ibcd_r2.sh

#BSUB -J "IBCDr2[1-40]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_r2_%I.out
#BSUB -e logs/ibcd_r2_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data_r2"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs/r2"

SEEDS=(42 43 44 45 46 47 48 49 50 51)
DIMS=(50 250)
VERSIONS=(full thresholded)

i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
D=${DIMS[$(( (i / 10) % 2 ))]}
VERSION=${VERSIONS[$(( i / 20 ))]}

DATA_FILE="$DATA_ROOT/$VERSION/${D}d/100/sf/${SEED}/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/$VERSION/${D}d/sf/${SEED}"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run experiment/simulation/run_generate_extra.sh first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: version=$VERSION D=$D seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior sf \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
