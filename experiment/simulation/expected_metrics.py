"""Expected F1 and SHD over the posterior draws of one IBCD run (Table 9).

Each draw of G is thresholded at |G| > eps and scored against the true graph
thresholded the same way (off-diagonal entries), as experiment/simulation/
metrics.R scores every method; F1 and SHD are then averaged over the draws.
Run on the cluster after ibcd.py, so the draws never have to be copied.

Usage:

    python expected_metrics.py --run_dir <ibcd output dir> --g_true G_matrix.csv \
        [--eps 0.075]

Writes <run_dir>/expected_metrics.json with E_F1, E_SHD, eps and n_draws.
"""

import argparse
import json
import os

import numpy as np
import pandas as pd


def expected_metrics(draws, G_true, eps=0.075, chunk=200):
    """E-F1 and E-SHD of thresholded posterior draws against a true graph.

    Args:
        draws (array): Posterior draws of G, shape (..., D, D).
        G_true (np.ndarray): D x D true graph.
        eps (float): Threshold on |G| for an edge, in each draw and in the truth.
        chunk (int): Draws processed at a time, to bound memory.

    Returns:
        dict: E_F1, E_SHD (averages over draws) and n_draws.
    """
    D = G_true.shape[0]
    flat = draws.reshape(-1, D, D)
    off = ~np.eye(D, dtype=bool)
    truth = np.abs(np.asarray(G_true, dtype=float))[off] > eps
    f1, shd = [], []
    for i in range(0, flat.shape[0], chunk):
        a = np.abs(np.asarray(flat[i:i + chunk]))[:, off] > eps
        tp = (a & truth).sum(axis=1)
        fp = (a & ~truth).sum(axis=1)
        fn = (~a & truth).sum(axis=1)
        denom = 2 * tp + fp + fn
        f1.append(np.where(denom > 0, 2 * tp / np.maximum(denom, 1), 1.0))
        shd.append(fp + fn)
    f1, shd = np.concatenate(f1), np.concatenate(shd)
    return {"E_F1": float(f1.mean()), "E_SHD": float(shd.mean()), "n_draws": int(f1.size)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_dir", required=True, help="ibcd.py output directory holding G_draws.npy.")
    ap.add_argument("--g_true", required=True, help="CSV of the true G (no index column).")
    ap.add_argument("--eps", type=float, default=0.075, help="Edge threshold on |G|. Default 0.075.")
    args = ap.parse_args()

    draws = np.load(os.path.join(args.run_dir, "G_draws.npy"), mmap_mode="r")
    G_true = pd.read_csv(args.g_true).values
    out = dict(expected_metrics(draws, G_true, args.eps), eps=args.eps)
    with open(os.path.join(args.run_dir, "expected_metrics.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"E-F1 {out['E_F1']:.4f}, E-SHD {out['E_SHD']:.1f} over {out['n_draws']} draws (eps {args.eps})")


if __name__ == "__main__":
    main()
