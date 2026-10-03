#!/bin/bash
# SID for the IBCD Figure 2 and 3 sweeps, one graph per task.
#
# Uses experiment/simulation/sid_one.R, which computes SID exactly as
# metrics.R does (both matrices thresholded at 0.075, no cycle resolution).
# IBCD's estimates are no longer DAGs, and SID treats their two-cycles as
# undirected CPDAG edges, which can take hours at D >= 250. As in Table 12, a
# graph whose SID does not finish within 12 hours is recorded as NA.
#
#   tasks   1-80   the dimension sweep  (same mapping as run_ibcd_dims_sweep.sh)
#   tasks 81-180   the n_int sweep      (same mapping as run_ibcd_nint_sweep.sh)
#
# SID needs CPU only, so this runs on i2c2_normal. R comes from the R/4.5
# module; SID and its Bioconductor dependencies have to be installed in your
# R library first (see the install line below). Results go to
# $PROJECT/IBCD_results/sid/, one CSV per graph.
#
#   module load R/4.5 && Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("graph", "RBGL", "SID"), ask = FALSE, update = FALSE)'
#
#   mkdir -p logs && bsub < run_ibcd_sid.sh

#BSUB -J "IBCDsid[1-180]"
#BSUB -q i2c2_normal
#BSUB -n 1
#BSUB -W 13:00
#BSUB -o logs/ibcd_sid_%I.out
#BSUB -e logs/ibcd_sid_%I.err

module load R/4.5

SID_R="$HOME/projects/IBCD/experiment/simulation/sid_one.R"
DATA_ROOT="$PROJECT/IBCD_results/sim_data"
RUN_ROOT="$PROJECT/IBCD_results/ibcd_runs/inverse"
OUT_ROOT="$PROJECT/IBCD_results/sid"
mkdir -p "$OUT_ROOT"

if ! Rscript -e 'suppressMessages(library(SID))' >/dev/null 2>&1; then
    echo "R package SID is not available to Rscript after 'module load R/4.5';" \
         "install it with the line at the top of this script" >&2
    exit 1
fi

SEEDS=(42 43 44 45 46 47 48 49 50 51)
GRAPHS=(er sf)
DIMS=(50 150 250 500)
N_INTS=(5 15 25 50 75)

if (( LSB_JOBINDEX <= 80 )); then
    i=$((LSB_JOBINDEX - 1))
    D=${DIMS[$(( i / 20 ))]}; N_INT=100
else
    i=$((LSB_JOBINDEX - 81))
    D=50; N_INT=${N_INTS[$(( i / 20 ))]}
fi
SEED=${SEEDS[$(( i % 10 ))]}
GRAPH=${GRAPHS[$(( (i / 10) % 2 ))]}

G_TRUE="$DATA_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/G_matrix.csv"
G_EST="$RUN_ROOT/${D}d/${N_INT}/${GRAPH}/${SEED}/G.csv"
LABEL="D=${D},n_int=${N_INT},graph=${GRAPH},seed=${SEED}"
OUT="$OUT_ROOT/${D}d_${N_INT}_${GRAPH}_${SEED}.csv"

echo "task ${LSB_JOBINDEX}: $LABEL"
for f in "$G_TRUE" "$G_EST"; do
    if [[ ! -f "$f" ]]; then echo "missing input: $f" >&2; exit 1; fi
done

timeout 12h Rscript "$SID_R" --g_true "$G_TRUE" --g_est "$G_EST" --out "$OUT" --label "$LABEL"
status=$?
if (( status == 124 )); then
    echo "label,sid,seconds,status" > "$OUT"
    echo "\"$LABEL\",NA,43200,timeout" >> "$OUT"
    echo "SID timed out after 12 h; recorded as NA"
fi
echo "task ${LSB_JOBINDEX} finished with status $status"
