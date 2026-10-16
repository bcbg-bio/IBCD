#!/bin/bash
# SID for the reduced-R^2-sortability SF graphs (Tables 3 and 4), one graph
# per task, from the runs of run_ibcd_r2.sh.
#
# Uses experiment/simulation/sid_one.R, which thresholds both the estimate
# and the truth at --eps and passes them to SID::structIntervDist with no
# cycle resolution, as for Tables 12 and 16. The thresholds are the paper's:
#
#   full         (Table 3)  0.07 at D = 50, 0.01 at D = 250
#   thresholded  (Table 4)  0.075
#
# The full version at D = 250 keeps many small entries, so its cyclic
# estimates can need a lot of memory and time; as in Table 12, a graph whose
# SID does not finish within 12 hours is recorded as NA.
#
#   2 versions x 2 dims x 10 seeds = 40 array elements, same mapping as
#   run_ibcd_r2.sh (seed fastest, then D, then version).
#
# Results go to $PROJECT/IBCD_results/sid_r2/, one CSV per graph. SID and its
# Bioconductor dependencies must be installed (see run_ibcd_sid.sh).
#
#   mkdir -p logs && bsub < run_ibcd_sid_r2.sh

#BSUB -J "IBCDsidr2[1-40]"
#BSUB -q i2c2_normal
#BSUB -n 1
#BSUB -M 64000
#BSUB -R "rusage[mem=64000]"
#BSUB -W 13:00
#BSUB -o logs/ibcd_sid_r2_%I.out
#BSUB -e logs/ibcd_sid_r2_%I.err

module load R/4.5

SID_R="$HOME/projects/IBCD/experiment/simulation/sid_one.R"
DATA_ROOT="$PROJECT/IBCD_results/sim_data_r2"
RUN_ROOT="$PROJECT/IBCD_results/ibcd_runs/r2"
OUT_ROOT="$PROJECT/IBCD_results/sid_r2"
mkdir -p "$OUT_ROOT"

if ! Rscript -e 'suppressMessages(library(SID))' >/dev/null 2>&1; then
    echo "R package SID is not available to Rscript after 'module load R/4.5';" \
         "install it with the line at the top of run_ibcd_sid.sh" >&2
    exit 1
fi

SEEDS=(42 43 44 45 46 47 48 49 50 51)
DIMS=(50 250)
VERSIONS=(full thresholded)

i=$((LSB_JOBINDEX - 1))
SEED=${SEEDS[$(( i % 10 ))]}
D=${DIMS[$(( (i / 10) % 2 ))]}
VERSION=${VERSIONS[$(( i / 20 ))]}

if [[ $VERSION == thresholded ]]; then
    EPS=0.075
elif (( D == 50 )); then
    EPS=0.07
else
    EPS=0.01
fi

G_TRUE="$DATA_ROOT/$VERSION/${D}d/100/sf/${SEED}/G_matrix.csv"
G_EST="$RUN_ROOT/$VERSION/${D}d/sf/${SEED}/G.csv"
LABEL="version=${VERSION},D=${D},seed=${SEED},eps=${EPS}"
OUT="$OUT_ROOT/${VERSION}_${D}d_${SEED}.csv"

echo "task ${LSB_JOBINDEX}: $LABEL"
for f in "$G_TRUE" "$G_EST"; do
    if [[ ! -f "$f" ]]; then echo "missing input: $f" >&2; exit 1; fi
done

timeout 12h Rscript "$SID_R" --g_true "$G_TRUE" --g_est "$G_EST" --eps "$EPS" --out "$OUT" --label "$LABEL"
status=$?
if (( status == 124 )); then
    echo "label,sid,seconds,status" > "$OUT"
    echo "\"$LABEL\",NA,43200,timeout" >> "$OUT"
    echo "SID timed out after 12 h; recorded as NA"
fi
echo "task ${LSB_JOBINDEX} finished with status $status"
