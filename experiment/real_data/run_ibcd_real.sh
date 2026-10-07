#!/bin/bash
# IBCD on the K562 Perturb-seq data (paper Section 4.6): Figure 5, Figure 8 and
# IBCD's row of Table 2.
#
#   screen  essential, gwps
#   prior   sf, er
#   split   train        all cells           -> Figure 8 (essential vs GWPS)
#           train_fold1-5  80% CV training sets -> Figure 5 and Table 2 (fold pairs)
#
# 2 x 2 x 6 = 24 array elements. With the EM-anchored SF prior each takes about
# 14.5 h (26 s per iteration), about 5x a D = 500 run of the Figure 2 sweep:
# on these screens the EM puts about 2/3 of R_hat in the slab, so both priors
# are close to dense. Inputs are the files run_preprocess.sh writes;
# sampler settings are the defaults, as in the sweeps, and are recorded in each
# run's diagnostics.json. Analyse with pip_agreement.py once all are done.
#
# Disk: each G_draws.npy is about 3.3 GB, so about 80 GB in all.
#
# SF_ANCHOR selects the SF prior's sparsity level (ibcd.py --sf_anchor): em,
# the default, or none, the published max-normalisation. none writes to a
# separate directory (ibcd_sfnone/). It has no effect on ER tasks, so
# submit only the SF ones, e.g.
#
#   SF_ANCHOR=none bsub -J "IBCDreal[1-6,13-18]" < run_ibcd_real.sh
#
# INIT_STRATEGY sets where NUTS starts (ibcd.py --init_strategy): rhat, the
# default, or optimized, the empty graph plus the same Adam descent. Anything
# but rhat appends _init<strategy> to the output directory, e.g.
#
#   INIT_STRATEGY=optimized bsub -J "IBCDreal[1-24]" < run_ibcd_real.sh
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
SF_ANCHOR="${SF_ANCHOR:-em}"
RUN_DIR=ibcd
[[ $SF_ANCHOR == none ]] && RUN_DIR=ibcd_sfnone
INIT_STRATEGY="${INIT_STRATEGY:-rhat}"
[[ $INIT_STRATEGY != rhat ]] && RUN_DIR=${RUN_DIR}_init${INIT_STRATEGY}
OUT_DIR="$ROOT/$RUN_DIR/$SCREEN/$PRIOR/$SPLIT"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run run_preprocess.sh first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX} (sf_anchor $SF_ANCHOR, init $INIT_STRATEGY): screen=$SCREEN prior=$PRIOR split=$SPLIT"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$PRIOR" \
    --output_dir "$OUT_DIR" \
    --sf_anchor "$SF_ANCHOR" \
    --init_strategy "$INIT_STRATEGY" \
    --save_diagnostics

echo "task ${LSB_JOBINDEX} finished with status $?"
