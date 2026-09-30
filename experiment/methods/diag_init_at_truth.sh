#!/bin/bash
# Diagnostic: is the SF runaway (rho(G) in the hundreds at D=150) where the
# posterior actually is, or a basin the sampler falls into from its default
# start?
#
# Runs NUTS on an existing SF D=150 run's inputs with the production sampler
# settings, from two starting points:
#   truth   - initialised at the true G
#   median  - the pipeline's default, init_to_median(num_samples=50)
# and records rho(G) every 10 iterations through warmup and sampling.
#
# If the truth-initialised chains stay at rho < 1 the runaway is a trap and the
# fix is initialisation / geometry. If they also run away, the posterior puts
# mass there and the fix has to be in the prior.
#
# Reads the R.csv, U.csv, V.csv, xi.csv, SE_hat.csv that ibcd.py already wrote
# into each run directory, plus the true G_matrix.csv. Writes one json per task.
#
#   mkdir -p logs && bsub < diag_init_at_truth.sh

#BSUB -J "IBCDinit[1-6]"
#BSUB -q dbeigpu
#BSUB -n 1
#BSUB -gpu "num=1:gmem=20G"
#BSUB -W 06:00
#BSUB -o logs/diag_init_%I.out
#BSUB -e logs/diag_init_%I.err

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibcd

SEEDS=(42 43 44)
ARMS=(truth median)
i=$((LSB_JOBINDEX - 1))
export SEED=${SEEDS[$(( i % 3 ))]}
export ARM=${ARMS[$(( i / 3 ))]}
export RUN_DIR="$PROJECT/IBCD_results/ibcd_runs/inverse/150d/100/sf/$SEED"
export TRUE_G="$PROJECT/IBCD_results/sim_data/150d/100/sf/$SEED/G_matrix.csv"
export OUT="$PROJECT/IBCD_results/diag_init/${ARM}_${SEED}.json"
export IBCD_SRC="$HOME/projects/IBCD/src"
mkdir -p "$(dirname "$OUT")"
echo "task $LSB_JOBINDEX: arm=$ARM seed=$SEED"

python - <<'PY'
import os, sys, json, time
import numpy as np, pandas as pd, jax, jax.numpy as jnp
sys.path.insert(0, os.environ["IBCD_SRC"])
from numpyro import infer
from numpyro.infer import MCMC, NUTS
from empirical_prior import (scale_free_degree, solve_edge_weights_rowwise,
                             load_R_and_SE_hat, empirical_bayes_em)
from model import matrix_model_spike_horseshoe

arm, seed, d = os.environ["ARM"], int(os.environ["SEED"]), os.environ["RUN_DIR"]
Rh = pd.read_csv(f"{d}/R.csv").values
xi = pd.read_csv(f"{d}/xi.csv").values
D = xi.shape[0]

# the SF prior exactly as ibcd.py builds it (defaults alpha_er=2, pi0_floor=0.05)
w, se = load_R_and_SE_hat(f"{d}/R.csv", f"{d}/SE_hat.csv")
pi0_global, _, _ = empirical_bayes_em(w, se, alpha_er=2.0)
pi0_ij, _ = solve_edge_weights_rowwise(
    xi, scale_free_degree(Rh, pi0_global=float(pi0_global), pi0_floor=0.05))
UL = jnp.linalg.cholesky(jnp.array(pd.read_csv(f"{d}/U.csv").values))
VL = jnp.linalg.cholesky(jnp.array(pd.read_csv(f"{d}/V.csv").values))

if arm == "truth":
    # G = pi0*spike + (1-pi0)*tau*lam*eps: put the true edges in the slab
    Gt = pd.read_csv(os.environ["TRUE_G"]).values.astype(float)
    Gt[np.abs(Gt) < 1e-8] = 0.0                       # generator leaves denormal dust
    edge = Gt != 0
    mult = (1.0 - pi0_ij) * 0.1
    lam = np.where(edge, np.abs(Gt) / np.maximum(mult, 1e-12), 1.0)
    eps = np.where(edge, np.sign(Gt), 0.0)
    init = infer.init_to_value(values={"lam": jnp.array(lam), "eps": jnp.array(eps),
                                       "spike": jnp.zeros((D, D))})
else:
    init = infer.init_to_median(num_samples=50)

kernel = NUTS(matrix_model_spike_horseshoe, target_accept_prob=0.9,
              max_tree_depth=12, init_strategy=init)
mcmc = MCMC(kernel, num_warmup=1000, num_samples=1000, num_chains=3,
            chain_method="vectorized", progress_bar=False)
kw = dict(obs_data=Rh, pi0_ij=pi0_ij, U_lower=UL, V_lower=VL, D=D)

def rho_trace(G):                                     # G: (chains, iters, D, D)
    G = np.asarray(G)[:, ::10]
    return np.abs(np.linalg.eigvals(G.reshape(-1, D, D))).max(1).reshape(G.shape[0], -1)

t0 = time.time()
key = jax.random.PRNGKey(seed)
mcmc.warmup(key, collect_warmup=True, extra_fields=("diverging",), **kw)
rho_w = rho_trace(mcmc.get_samples(group_by_chain=True)["G"])
div_w = float(np.mean(mcmc.get_extra_fields()["diverging"]))
mcmc.post_warmup_state = mcmc.last_state
mcmc.run(mcmc.post_warmup_state.rng_key, extra_fields=("diverging",), **kw)
rho_s = rho_trace(mcmc.get_samples(group_by_chain=True)["G"])
div_s = float(np.mean(mcmc.get_extra_fields()["diverging"]))

def first_escape(r):                                  # iteration index, or None
    k = np.flatnonzero(r >= 1.0)
    return int(k[0]) * 10 if k.size else None

out = {
    "arm": arm, "seed": seed, "D": int(D), "pi0_global": float(pi0_global),
    "runtime_s": time.time() - t0,
    "divergence_pct": {"warmup": 100 * div_w, "sampling": 100 * div_s},
    "sampling_frac_rho_below_1": float(np.mean(rho_s < 1.0)),
    "sampling_rho_median_per_chain": [float(np.median(r)) for r in rho_s],
    "warmup_first_iter_rho_ge_1_per_chain": [first_escape(r) for r in rho_w],
    "rho_every_10_iters": {"warmup": rho_w.round(4).tolist(),
                           "sampling": rho_s.round(4).tolist()},
}
json.dump(out, open(os.environ["OUT"], "w"))
print(f"{arm} seed {seed}: frac rho<1 during sampling {out['sampling_frac_rho_below_1']:.3f}; "
      f"median rho per chain {[round(x, 3) for x in out['sampling_rho_median_per_chain']]}; "
      f"first warmup iter with rho>=1 {out['warmup_first_iter_rho_ge_1_per_chain']}; "
      f"{out['runtime_s'] / 3600:.2f} h")
PY
echo "task ${LSB_JOBINDEX} finished with status $?"
