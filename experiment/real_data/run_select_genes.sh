#!/bin/bash
# Gene selection for the K562 Perturb-seq analysis (paper Appendix D.3): genes
# whose best guide knocks them down by more than 0.75 SD with more than 50
# cells, in both the essential and the genome-wide screens. The paper reports
# 521. See select_genes.R.
#
# Needs hdf5r and dplyr in your R library:
#
#   module load R/4.5 && Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages(c("hdf5r", "dplyr"))'
#
# hdf5r compiles against the system HDF5 library; if its install fails, load
# an hdf5 module first.
#
# Memory: the GWPS non-targeting cells are held in memory (all genes x control
# cells), a few GB; the request below leaves room. -M/rusage units depend on
# the cluster's LSF configuration (usually MB).
#
#   mkdir -p logs && bsub < run_select_genes.sh

#BSUB -J IBCDgenes
#BSUB -q i2c2_normal
#BSUB -n 1
#BSUB -M 64000
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/ibcd_select_genes.out
#BSUB -e logs/ibcd_select_genes.err

module load R/4.5

DATA="$PROJECT/datasets/perturb_seq"
OUT="$PROJECT/IBCD_results/real_data/genes"
mkdir -p "$OUT"

Rscript "$HOME/projects/IBCD/experiment/real_data/select_genes.R" \
    --essential "$DATA/K562_essential_normalized_singlecell_01.h5ad" \
    --gwps "$DATA/K562_gwps_normalized_singlecell_01.h5ad" \
    --out_dir "$OUT"

echo "finished with status $?"
