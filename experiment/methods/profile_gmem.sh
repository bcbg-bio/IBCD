#!/bin/bash
# Peak GPU memory of an IBCD run, to size the gmem request of the large runs.
#
# JAX preallocates most of the device by default, so nvidia-smi shows the
# pool, not what a run needs. This records JAX's own allocator peak
# (memory_stats()["peak_bytes_in_use"]) after a short run. The stored draws
# grow linearly with --num_samples, so two sample counts give the per-draw
# cost and an extrapolation to the production 1000:
#
#   tasks 1-2    simulated D = 50 SF, seed 42       num_samples 100, 200
#   tasks 3-4    simulated D = 150 SF, seed 42      num_samples 100, 200
#   tasks 5-6    simulated D = 250 SF, seed 42      num_samples 100, 200
#   tasks 7-8    simulated D = 500 SF, seed 42      num_samples 100, 200
#   tasks 9-10   real essential screen, SF, train   num_samples 100, 200
#
# 20 warmup iterations each; other settings are the defaults (3 chains). The
# stored draws do not depend on the prior, so SF stands in for ER too. Writes
# gmem.json to each output directory and prints the extrapolation. Expect
# minutes at D <= 250, under an hour at D = 500 and up to about 2 hours for
# the real data.
#
#   mkdir -p logs && bsub < profile_gmem.sh

#BSUB -J "IBCDgmem[1-10]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=40G"
#BSUB -W 04:00
#BSUB -o logs/ibcd_gmem_%I.out
#BSUB -e logs/ibcd_gmem_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

SRC="$HOME/projects/IBCD/src"
ROOT="$PROJECT/IBCD_results"

i=$((LSB_JOBINDEX - 1))
SAMPLES=(100 200)
NS=${SAMPLES[$(( i % 2 ))]}
DIMS=(50 150 250 500)
k=$(( i / 2 ))
if (( k < 4 )); then
    D=${DIMS[$k]}
    NAME=sim${D}_sf; DATA="$ROOT/sim_data/${D}d/100/sf/42/Y_with_targets.csv"
else
    NAME=real_essential_sf; DATA="$ROOT/real_data/essential/input/Y_matrix_essential_train.csv"
fi
OUT="$ROOT/gmem_profile/${NAME}_ns${NS}"
mkdir -p "$OUT"
echo "task ${LSB_JOBINDEX}: $NAME, num_samples $NS -> $OUT"

python - "$SRC" "$DATA" "$OUT" "$NS" <<'PY'
import json, runpy, sys
src, data, out, ns = sys.argv[1:5]
sys.path.insert(0, src)
sys.argv = [f"{src}/ibcd.py", "--data", data, "--prior", "sf", "--output_dir", out,
            "--num_warmup", "20", "--num_samples", ns, "--save_diagnostics"]
runpy.run_path(f"{src}/ibcd.py", run_name="__main__")
import jax
st = jax.devices()[0].memory_stats()
if st is None:
    sys.exit(f"no memory_stats on {jax.devices()[0]}; this needs a GPU")
rec = {k: st.get(k) for k in ("peak_bytes_in_use", "bytes_limit", "bytes_in_use")}
rec["num_samples"] = int(ns)
json.dump(rec, open(f"{out}/gmem.json", "w"), indent=2)
print(f"peak {rec['peak_bytes_in_use'] / 1e9:.2f} GB of {rec['bytes_limit'] / 1e9:.2f} GB at num_samples {ns}")
PY

# once both sample counts are done, extrapolate to num_samples 1000
python - "$ROOT/gmem_profile/${NAME}" <<'PY'
import json, os, sys
base = sys.argv[1]
recs = {}
for ns in (100, 200):
    p = f"{base}_ns{ns}/gmem.json"
    if os.path.exists(p):
        recs[ns] = json.load(open(p))["peak_bytes_in_use"] / 1e9
if len(recs) == 2:
    per = (recs[200] - recs[100]) / 100
    print(f"{os.path.basename(base)}: {per * 1000:.2f} GB per 1000 draws; projected peak at 1000 draws "
          f"{recs[100] + 900 * per:.1f} GB")
PY
