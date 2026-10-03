#!/bin/bash
# IBCD on the K562 Perturb-seq data (paper Section 4.6): Figure 5, Figure 8 and
# IBCD's row of Table 2.
#
#   screen  essential, gwps
#   prior   sf, er
#   split   train        all cells           -> Figure 8 (essential vs GWPS)
#           train_fold1-5  80% CV training sets -> Figure 5 and Table 2 (fold pairs)
#
# 2 x 2 x 6 = 24 array elements, at D = 521 each about as costly as one D = 500
# run of the Figure 2 sweep. Inputs are the files run_preprocess.sh writes;
# sampler settings are the defaults, as in the sweeps, and are recorded in each
# run's diagnostics.json. Analyse with pip_agreement.py once all are done.
#
# Disk: each G_draws.npy is about 3.3 GB, so about 80 GB in all.
#
#   mkdir -p logs && bsub < run_ibcd_real.sh

#BSUB -J "IBCDreal[1-24]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_real_%I.out
#BSUB -e logs/ibcd_real_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
ROOT="$PROJECT/IBCD_results/real_data"

SPLITS=(train train_fold1 train_fold2 train_fold3 train_fold4 train_fold5)
PRIORS=(sf er)
SCREENS=(essential gwps)

# 1..24 -> (screen, prior, split); split varies fastest
i=$((LSB_JOBINDEX - 1))
SPLIT=${SPLITS[$(( i % 6 ))]}
PRIOR=${PRIORS[$(( (i / 6) % 2 ))]}
SCREEN=${SCREENS[$(( i / 12 ))]}

DATA_FILE="$ROOT/$SCREEN/input/Y_matrix_${SCREEN}_${SPLIT}.csv"
OUT_DIR="$ROOT/ibcd/$SCREEN/$PRIOR/$SPLIT"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run run_preprocess.sh first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: screen=$SCREEN prior=$PRIOR split=$SPLIT"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$PRIOR" \
    --output_dir "$OUT_DIR" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
