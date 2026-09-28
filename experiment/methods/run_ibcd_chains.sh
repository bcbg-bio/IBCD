#!/bin/bash
# One MCMC chain per LSF array job.
#
# CUDA exposes a single MIG device per process on this node, so num_chains > 1
# inside one job is drawn sequentially however many slices are reserved. Each
# array element instead draws one chain on its own slice, with the array index
# as the seed, and combine_chains.py merges them afterwards.
#
#   bsub < run_ibcd_chains.sh          # draws the chains
#   # then, once all three have finished:
#   python combine_chains.py --chain_dirs "$OUTPUT_DIR"/chain_{1,2,3} \
#                            --output_dir "$OUTPUT_DIR/combined"

#BSUB -J "IBCD[1-3]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 24:00
#BSUB -o ibcd_chain_%I.out
#BSUB -e ibcd_chain_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

DATA_FILE="/path/to/X.csv"
OUTPUT_DIR="/path/to/ibcd_out"
IBCD="$HOME/projects/IBCD/src/ibcd.py"

# Keep these: they are how the MIG profile and memory ceiling get confirmed
# on every run.
echo "CUDA_VISIBLE_DEVICES = ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import jax; print('bytes_limit', jax.devices()[0].memory_stats()['bytes_limit'])"

CHAIN_DIR="$OUTPUT_DIR/chain_${LSB_JOBINDEX}"
mkdir -p "$CHAIN_DIR"

python "$IBCD" \
    --data "$DATA_FILE" \
    --prior sf \
    --output_dir "$CHAIN_DIR" \
    --num_chains 1 \
    --seed "$LSB_JOBINDEX" \
    --save_diagnostics

echo "chain ${LSB_JOBINDEX} finished"
