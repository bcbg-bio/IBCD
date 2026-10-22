#!/bin/bash
# Preprocess the K562 Perturb-seq screens into IBCD inputs (see
# preprocess_screen.R): cell-cycle regression, per-GEM normalisation against
# the non-targeting cells, the full data and 5 CV folds.
#
#   task 1   essential   K562_essential_raw_singlecell_01.h5ad  (~10 GB)
#   task 2   GWPS        K562_gwps_raw_singlecell_01.h5ad       (~66 GB)
#
# The gene list defaults to experiment/real_data/common_genes.csv, the 521
# genes used in the paper; set GENES to use another (e.g. the common_genes.csv
# that run_select_genes.sh writes to its output directory):
#
#   GENES=$PROJECT/IBCD_results/real_data/genes/common_genes.csv bsub < run_preprocess.sh
#
# Needs hdf5r, tidyverse, matrixStats and caret in your R library:
#
#   module load R/4.5 hdf5/1.12.1 && Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages(c("hdf5r", "tidyverse", "matrixStats", "caret"))'
#
# The script streams the raw counts in blocks of cells, so memory stays at a
# few GB; it reads the whole count matrix about four times, which for GWPS is
# the bulk of the run time.
#
#   mkdir -p logs && bsub < run_preprocess.sh

#BSUB -J "IBCDprep[1-2]"
#BSUB -q i2c2_normal
#BSUB -n 1
#BSUB -M 64000
#BSUB -R "rusage[mem=64000]"
#BSUB -W 72:00
#BSUB -o logs/ibcd_preprocess_%I.out
#BSUB -e logs/ibcd_preprocess_%I.err

module load R/4.5
module load hdf5/1.12.1

REPO="$HOME/projects/IBCD"
DATA="$PROJECT/datasets/perturb_seq"
OUT_ROOT="$PROJECT/IBCD_results/real_data"
GENES="${GENES:-$REPO/experiment/real_data/common_genes.csv}"

SCREENS=(essential gwps)
SCREEN=${SCREENS[$(( LSB_JOBINDEX - 1 ))]}
RAW="$DATA/K562_${SCREEN}_raw_singlecell_01.h5ad"
OUT="$OUT_ROOT/$SCREEN"

echo "task ${LSB_JOBINDEX}: $SCREEN"
echo "  raw   $RAW"
echo "  genes $GENES"
echo "  out   $OUT"
for f in "$RAW" "$GENES"; do
    if [[ ! -f "$f" ]]; then echo "missing input: $f" >&2; exit 1; fi
done
mkdir -p "$OUT"

Rscript "$REPO/experiment/real_data/preprocess_screen.R" \
    --raw_h5ad "$RAW" \
    --genes "$GENES" \
    --screen "$SCREEN" \
    --out_dir "$OUT"

echo "task ${LSB_JOBINDEX} finished with status $?"
