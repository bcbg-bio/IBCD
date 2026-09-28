"""Combine per-chain IBCD runs into one set of posterior summaries.

On a MIG-partitioned node CUDA exposes a single device per process, so
`--num_chains 3` cannot run in parallel inside one job. The alternative is one
chain per job, each with its own `--seed` and `--output_dir`, combined here.

Each input directory must contain the `G_draws.npy` written by `ibcd.py`, of
shape (chains, draws, D, D); chains are concatenated along the first axis.
Column names are taken from the first directory's `G.csv`.

    python combine_chains.py --chain_dirs out/chain1 out/chain2 out/chain3 \\
                             --output_dir out/combined

Draws whose spectral radius reaches 1 are excluded from every summary, exactly
as in a single run, and R-hat and ESS are computed across the combined chains.
"""

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd

from model import compute_lfsr, convergent_draws, posterior_diagnostics


def load_chains(chain_dirs):
    """Stack per-chain draws and recover the variable names.

    Returns:
        draws (np.ndarray): shape (n_chains, n_draws, D, D).
        colnames (list): variable names, or None if no G.csv was found.
    """
    blocks, colnames = [], None
    for d in chain_dirs:
        path = os.path.join(d, "G_draws.npy")
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found")
        g = np.load(path)
        if g.ndim != 4:
            raise ValueError(f"{path} has shape {g.shape}, expected 4 dimensions")
        blocks.append(g)
        if colnames is None:
            g_csv = os.path.join(d, "G.csv")
            if os.path.exists(g_csv):
                colnames = pd.read_csv(g_csv, index_col=0).columns.tolist()

    shapes = {b.shape[1:] for b in blocks}
    if len(shapes) > 1:
        raise ValueError(f"chains disagree on draw shape: {sorted(shapes)}")
    return np.concatenate(blocks, axis=0), colnames


def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    posterior, colnames = load_chains(args.chain_dirs)
    n_chains, n_draws, D, _ = posterior.shape
    if colnames is None:
        colnames = [f"V{i + 1}" for i in range(D)]
    print(f"Combined {n_chains} chains x {n_draws} draws at D={D}")

    flat = posterior.reshape(-1, D, D)
    keep, rho = convergent_draws(flat)
    n_drop = int((~keep).sum())
    if n_drop:
        pct = 100.0 * n_drop / keep.size
        per_chain = (~keep).reshape(n_chains, n_draws).sum(axis=1)
        detail = ", ".join(f"chain {c}: {int(k)}/{n_draws}"
                           for c, k in enumerate(per_chain))
        message = (
            f"{n_drop} of {keep.size} draws ({pct:.1f}%) have spectral radius "
            f">= 1 and were excluded ({detail}); max rho = {rho.max():.3g}"
        )
        if pct >= 10.0:
            warnings.warn(message + ". Treat these results with caution.",
                          RuntimeWarning)
        else:
            print(message)
    if not keep.any():
        raise RuntimeError("Every draw has spectral radius >= 1.")

    diagnostics = posterior_diagnostics(posterior, rho, keep, seed=args.seed)
    diagnostics["config"] = {"chain_dirs": list(args.chain_dirs),
                             "epsilon": args.epsilon}
    with open(os.path.join(args.output_dir, "diagnostics.json"), "w") as fh:
        json.dump(diagnostics, fh, indent=2)

    issues = []
    if diagnostics.get("r_hat", {}).get("max", 0.0) > 1.05:
        issues.append(f"max r_hat {diagnostics['r_hat']['max']:.3f} > 1.05")
    if diagnostics.get("ess", {}).get("n_nonpositive_raw", 0) > 0:
        issues.append(
            f"{diagnostics['ess']['n_nonpositive_raw']} entries had a "
            "non-positive ESS estimate"
        )
    if issues:
        warnings.warn("Chains did not converge cleanly: " + "; ".join(issues)
                      + ". See diagnostics.json.", RuntimeWarning)

    kept = flat[keep]
    frame = lambda m: pd.DataFrame(m, columns=colnames, index=colnames)  # noqa: E731
    frame(kept.mean(axis=0)).to_csv(f"{args.output_dir}/G.csv", index=True)
    frame((np.abs(kept) > args.epsilon).mean(axis=0)).to_csv(
        f"{args.output_dir}/pip.csv", index=True)
    frame(compute_lfsr(kept)).to_csv(f"{args.output_dir}/lfsr.csv", index=True)

    print(
        "Diagnostics: "
        f"max r_hat {diagnostics.get('r_hat', {}).get('max', float('nan')):.4f}, "
        f"min ESS {diagnostics.get('ess', {}).get('min', float('nan')):.0f}, "
        f"rho median {diagnostics['spectral_radius']['median']:.3f} "
        f"max {diagnostics['spectral_radius']['max']:.3g}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Combine per-chain IBCD runs into one posterior mean, PIP and "
            "LFSR, with cross-chain R-hat and ESS."
        )
    )
    parser.add_argument("--chain_dirs", required=True, nargs="+",
                        help="Output directories of the individual chain runs.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write the combined summaries to.")
    parser.add_argument("--epsilon", type=float, default=0.05,
                        help="PIP threshold: edges with |G| > epsilon are "
                             "counted as active. Default = 0.05.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Seed for the ESS subsample. Default = 42.")
    main(parser.parse_args())
