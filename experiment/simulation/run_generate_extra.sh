#!/bin/bash
# Simulated data for the ablations (Tables 7 and 8) and the reduced
# R^2-sortability SF graphs (Appendix F.2, Tables 3 and 4), one
# (setting, D, graph) per task, 10 seeds each, n_int = 100.
#
#   tasks 1-8    ablations: int_beta -1.5, v 0.2 (Appendix D.1), ER and SF at
#                D = 50, 150, 250, 500 -> sim_data_ablation/
#   tasks 9-12   R^2: SF, int_beta -2, v 0.25, net_vars 0.4 (D = 50) or
#                0.5 (D = 250); 'full' without and 'thresholded' with
#                min_v = 0.125 -> sim_data_r2/{full,thresholded}/
#
# Same generator and seeds as the Figure 2 sweep (generate_sim_data.R), so the
# settings above are the only differences. IV summaries are not written:
# ibcd.py computes its own.
#
# Needs inspre and mc2d in your R library, as for the sweep data.
#
#   mkdir -p logs && bsub < run_generate_extra.sh

#BSUB -J "IBCDgen[1-12]"
#BSUB -q i2c2_normal
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -W 8:00
#BSUB -o logs/ibcd_gen_%I.out
#BSUB -e logs/ibcd_gen_%I.err

module load R/4.5

GEN="$HOME/projects/IBCD/experiment/simulation/generate_sim_data.R"
ROOT="$PROJECT/IBCD_results"
SEEDS="42:51"

i=$((LSB_JOBINDEX - 1))
if (( LSB_JOBINDEX <= 8 )); then
    DIMS=(50 150 250 500)
    GRAPHS=(er sf)
    D=${DIMS[$(( i % 4 ))]}
    GRAPH=${GRAPHS[$(( i / 4 ))]}
    echo "task ${LSB_JOBINDEX}: ablation data, D=$D graph=$GRAPH"
    Rscript "$GEN" --out_dir "$ROOT/sim_data_ablation" --dims "$D" --graphs "$GRAPH" \
        --n_int 100 --seeds "$SEEDS" --write_iv false --int_beta -1.5 --v 0.2
else
    j=$((LSB_JOBINDEX - 9))
    DIMS=(50 250)
    NET_VARS=(0.4 0.5)
    VERSIONS=(full thresholded)
    MIN_V=(none 0.125)
    D=${DIMS[$(( j % 2 ))]}
    NV=${NET_VARS[$(( j % 2 ))]}
    VERSION=${VERSIONS[$(( j / 2 ))]}
    MV=${MIN_V[$(( j / 2 ))]}
    echo "task ${LSB_JOBINDEX}: R^2 data ($VERSION), D=$D net_vars=$NV min_v=$MV"
    Rscript "$GEN" --out_dir "$ROOT/sim_data_r2/$VERSION" --dims "$D" --graphs sf \
        --n_int 100 --seeds "$SEEDS" --write_iv false --net_vars "$NV" --min_v "$MV"
fi

echo "task ${LSB_JOBINDEX} finished with status $?"
