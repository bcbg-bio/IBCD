#!/bin/bash
# Ablations (Tables 7 and 8) on the data from run_generate_extra.sh
# (int_beta -1.5, v 0.2). Sampler settings are the defaults, as in the sweeps.
#
#   tasks   1-70   Table 7, D = 50
#                    1-10   ER data, ER prior         (ER baseline)
#                   11-20   ER data, ER prior, oracle spike weight (--oracle_G)
#                   21-30   ER data, ER prior, MVN likelihood (--likelihood mvn)
#                   31-40   ER data, ER prior, one global spike weight (--global_prior)
#                   41-50   ER data, SF prior
#                   51-60   SF data, SF prior         (SF baseline)
#                   61-70   SF data, ER prior
#   tasks  71-130  Table 8, ER data at D = 150, 250, 500 under the ER and the
#                  SF prior (D slowest, then prior, then seed)
#   tasks 131-190  the same for SF data under the SF and the ER prior, which
#                  the paper did not have; submit [1-130] to leave it out
#
# Outputs: $PROJECT/IBCD_results/ibcd_runs/ablation/<D>d/<graph>_data/<arm>/<seed>
# with arm one of er_prior, sf_prior, oracle, mvn, global.
#
#   mkdir -p logs && bsub < run_ibcd_ablation.sh

#BSUB -J "IBCDabl[1-190]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o logs/ibcd_abl_%I.out
#BSUB -e logs/ibcd_abl_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

IBCD="$HOME/projects/IBCD/src/ibcd.py"
DATA_ROOT="$PROJECT/IBCD_results/sim_data_ablation"
OUT_ROOT="$PROJECT/IBCD_results/ibcd_runs/ablation"

SEEDS=(42 43 44 45 46 47 48 49 50 51)
i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
block=$(( i / 10 ))

# block -> D, data graph, prior, arm
if (( block < 7 )); then
    D=50
    GRAPHS=(er er er er er sf sf)
    PRIORS=(er er er er sf sf er)
    ARMS=(er_prior oracle mvn global sf_prior sf_prior er_prior)
    GRAPH=${GRAPHS[$block]}; PRIOR=${PRIORS[$block]}; ARM=${ARMS[$block]}
else
    b=$(( block - 7 ))                       # 0..11
    DIMS=(150 250 500)
    GRAPH=$( (( b < 6 )) && echo er || echo sf )
    b=$(( b % 6 ))
    D=${DIMS[$(( b / 2 ))]}
    if [[ $GRAPH == er ]]; then PRIORS=(er sf); else PRIORS=(sf er); fi
    PRIOR=${PRIORS[$(( b % 2 ))]}
    ARM=${PRIOR}_prior
fi

DATA_DIR="$DATA_ROOT/${D}d/100/${GRAPH}/${SEED}"
DATA_FILE="$DATA_DIR/Y_with_targets.csv"
OUT_DIR="$OUT_ROOT/${D}d/${GRAPH}_data/${ARM}/${SEED}"

EXTRA=()
case $ARM in
    oracle) EXTRA=(--oracle_G "$DATA_DIR/G_matrix.csv") ;;
    mvn)    EXTRA=(--likelihood mvn) ;;
    global) EXTRA=(--global_prior) ;;
esac

if [[ ! -f "$DATA_FILE" ]]; then
    echo "missing input: $DATA_FILE" >&2
    echo "run experiment/simulation/run_generate_extra.sh first" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

echo "task ${LSB_JOBINDEX}: D=$D data=$GRAPH prior=$PRIOR arm=$ARM seed=$SEED"
echo "  in  $DATA_FILE"
echo "  out $OUT_DIR"
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; d=jax.devices(); print('devices', d, 'bytes_limit', d[0].memory_stats()['bytes_limit'])"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior "$PRIOR" \
    --output_dir "$OUT_DIR" \
    --seed "$SEED" \
    --save_diagnostics \
    "${EXTRA[@]}"

echo "task ${LSB_JOBINDEX} finished with status $?"
