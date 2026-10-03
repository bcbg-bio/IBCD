"""
PIP agreement on the K562 Perturb-seq data (paper Section 4.6), from the runs
of run_ibcd_real.sh:

  Figure 5   cross-validation fold pairs: PIP(i) against PIP(j), i < j, all ten
             pairs overlaid, ER and SF side by side, one figure per screen.
  Figure 8   full essential screen against full GWPS screen, ER and SF.
  Table 2    IBCD's row: F1 and SHD between the fold pairs' posterior means
             thresholded at 0.075 (as experiment/methods/cv_comparison.R on
             main computes them), and the pooled fold-pair PIP correlation.

Only off-diagonal entries are used. Pearson r is reported both on PIP and on
log10 PIP; a PIP of 0 cannot go on a log axis, so for the log scale (and the
plots) zeros are set to half the smallest nonzero PIP, 1 / (2 * draws).

PIPs come from each run's pip.csv (the epsilon it was run with, 0.05 by
default) unless --epsilon is given, in which case they are recomputed from
G_draws.npy, keeping every draw as ibcd.py does by default.

Usage:

  python pip_agreement.py --root $PROJECT/IBCD_results/real_data/ibcd \
      --out_dir $PROJECT/IBCD_results/real_data/figures [--epsilon 0.05]

Expects <root>/<screen>/<prior>/<split>/ with split in train, train_fold1-5.
Missing runs are skipped with a note.
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SCREENS = ("essential", "gwps")
SCREEN_LABEL = {"essential": "Essential", "gwps": "GWPS"}
PRIORS = ("er", "sf")
FOLDS = tuple(f"train_fold{k}" for k in range(1, 6))


def n_draws(run_dir):
    draws = np.load(os.path.join(run_dir, "G_draws.npy"), mmap_mode="r")
    return int(np.prod(draws.shape[:-2]))


def load_run(run_dir, genes, epsilon=None):
    """PIP, posterior mean and draw count of one run, rows and columns in `genes` order."""
    G = pd.read_csv(os.path.join(run_dir, "G.csv"), index_col=0)
    if epsilon is None:
        pip = pd.read_csv(os.path.join(run_dir, "pip.csv"), index_col=0)
    else:
        draws = np.load(os.path.join(run_dir, "G_draws.npy"), mmap_mode="r")
        D = draws.shape[-1]
        flat = draws.reshape(-1, D, D)
        count = np.zeros((D, D))
        for s in range(0, flat.shape[0], 250):
            count += (np.abs(np.asarray(flat[s:s + 250])) > epsilon).sum(axis=0)
        pip = pd.DataFrame(count / flat.shape[0], index=G.index, columns=G.columns)
    if genes is None:
        genes = list(G.columns)
    if set(genes) != set(G.columns):
        raise ValueError(f"{run_dir}: genes differ from the first run loaded")
    return (pip.loc[genes, genes].to_numpy(), G.loc[genes, genes].to_numpy(),
            n_draws(run_dir), genes)


def offdiag(M):
    return M[~np.eye(M.shape[0], dtype=bool)]


def fold_pair_metrics(G1, G2, eps=0.075):
    """F1 and SHD of two thresholded graphs, as cv_comparison.R's calc_metrics(abs(G1), abs(G2), eps)."""
    a = offdiag(np.abs(G1)) >= eps
    b = offdiag(np.abs(G2)) >= eps
    tp = np.sum(a & b)
    fp = np.sum(a & ~b)
    fn = np.sum(~a & b)
    return 2 * tp / (2 * tp + fp + fn), int(fp + fn)


def log_pip(p, floor):
    return np.log10(np.maximum(p, floor))


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def scatter_panel(fig, spec, x, y, floor, title, xlabel, ylabel, side_hist):
    """Log-log scatter with a marginal histogram on top (and on the right if side_hist)."""
    lx, ly = log_pip(x, floor), log_pip(y, floor)
    lo = np.log10(floor) - 0.1
    lim = (lo, 0.05)
    bins = np.linspace(lo, 0.0, 41)
    if side_hist:
        gs = spec.subgridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4), wspace=0.03, hspace=0.03)
    else:
        gs = spec.subgridspec(2, 1, height_ratios=(1, 4), hspace=0.03)
    ax = fig.add_subplot(gs[1, 0])
    top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax.scatter(lx, ly, s=1, c="0.2", alpha=0.03, linewidths=0, rasterized=True)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ticks = np.arange(np.ceil(lo), 1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
    ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(0.03, 0.97, f"Pearson r = {pearson(x, y):.3f} (PIP)\n"
                        f"Pearson r = {pearson(lx, ly):.3f} (log PIP)",
            transform=ax.transAxes, va="top", fontsize=8)
    top.hist(lx, bins=bins, color="0.5", edgecolor="0.3", linewidth=0.3)
    top.axis("off")
    top.set_title(title)
    if side_hist:
        right = fig.add_subplot(gs[1, 1], sharey=ax)
        right.hist(ly, bins=bins, orientation="horizontal", color="0.5", edgecolor="0.3", linewidth=0.3)
        right.axis("off")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", required=True, help="Directory holding <screen>/<prior>/<split> runs.")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--epsilon", type=float, default=None,
                        help="Recompute PIP from G_draws.npy at this threshold. Default: use pip.csv.")
    parser.add_argument("--graph_eps", type=float, default=0.075,
                        help="Threshold on the posterior mean for the Table 2 F1 and SHD. Default 0.075.")
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    runs, genes = {}, None
    for screen, prior, split in itertools.product(SCREENS, PRIORS, ("train",) + FOLDS):
        run_dir = os.path.join(args.root, screen, prior, split)
        if not os.path.exists(os.path.join(run_dir, "G.csv")):
            print(f"missing: {run_dir}")
            continue
        pip, G, n, genes = load_run(run_dir, genes, args.epsilon)
        runs[screen, prior, split] = (offdiag(pip), G, n)
        print(f"loaded {screen} {prior} {split}: {n} draws, mean PIP {offdiag(pip).mean():.3f}")
    if not runs:
        raise SystemExit("no runs found")
    floor = 1.0 / (2 * max(n for _, _, n in runs.values()))
    eps_note = f"epsilon {args.epsilon}" if args.epsilon is not None else "pip.csv"

    # ---- Figure 5 and Table 2: fold pairs ----
    pair_rows, summary_rows = [], []
    for screen in SCREENS:
        panels = []
        for prior in PRIORS:
            folds = [f for f in FOLDS if (screen, prior, f) in runs]
            if len(folds) < 2:
                continue
            xs, ys = [], []
            for fi, fj in itertools.combinations(folds, 2):
                pi, Gi, _ = runs[screen, prior, fi]
                pj, Gj, _ = runs[screen, prior, fj]
                f1, shd = fold_pair_metrics(Gi, Gj, args.graph_eps)
                pair_rows.append(dict(screen=screen, prior=prior, fold_i=fi, fold_j=fj, F1=f1, SHD=shd,
                                      r_pip=pearson(pi, pj),
                                      r_log_pip=pearson(log_pip(pi, floor), log_pip(pj, floor))))
                xs.append(pi)
                ys.append(pj)
            x, y = np.concatenate(xs), np.concatenate(ys)
            pairs = pd.DataFrame([r for r in pair_rows if r["screen"] == screen and r["prior"] == prior])
            summary_rows.append(dict(
                comparison="fold pairs", screen=screen, prior=prior, n_pairs=len(pairs),
                F1_mean=pairs.F1.mean(), F1_sd=pairs.F1.std(), SHD_mean=pairs.SHD.mean(), SHD_sd=pairs.SHD.std(),
                r_pip_pooled=pearson(x, y), r_log_pip_pooled=pearson(log_pip(x, floor), log_pip(y, floor))))
            panels.append((prior, x, y))
        if not panels:
            continue
        fig = plt.figure(figsize=(4.2 * len(panels), 4.2))
        outer = fig.add_gridspec(1, len(panels), wspace=0.35)
        for k, (prior, x, y) in enumerate(panels):
            scatter_panel(fig, outer[k], x, y, floor, prior.upper(),
                          "PIP fold $i$", "PIP fold $j$", side_hist=False)
        fig.suptitle(f"K562 {SCREEN_LABEL[screen]} screen, all fold pairs $i < j$ ({eps_note})", fontsize=9, y=1.0)
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(args.out_dir, f"figure5_{screen}.{ext}"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ---- Figure 8: full essential against full GWPS ----
    panels = []
    for prior in PRIORS:
        if ("essential", prior, "train") in runs and ("gwps", prior, "train") in runs:
            x = runs["gwps", prior, "train"][0]
            y = runs["essential", prior, "train"][0]
            summary_rows.append(dict(comparison="essential vs gwps", screen="both", prior=prior,
                                     r_pip_pooled=pearson(x, y),
                                     r_log_pip_pooled=pearson(log_pip(x, floor), log_pip(y, floor))))
            panels.append((prior, x, y))
    if panels:
        fig = plt.figure(figsize=(4.6 * len(panels), 4.6))
        outer = fig.add_gridspec(1, len(panels), wspace=0.3)
        for k, (prior, x, y) in enumerate(panels):
            scatter_panel(fig, outer[k], x, y, floor, prior.upper(), "PIP GWPS", "PIP Essential", side_hist=True)
        fig.suptitle(f"Full screens ({eps_note})", fontsize=9, y=1.0)
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(args.out_dir, f"figure8.{ext}"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    pd.DataFrame(pair_rows).to_csv(os.path.join(args.out_dir, "fold_pairs.csv"), index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(args.out_dir, "pip_agreement_summary.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_columns", None, "display.precision", 3):
        print(summary.to_string(index=False))
    print(f"zeros floored at {floor:.2e} for log PIP; figures and tables in {args.out_dir}")


if __name__ == "__main__":
    main()
