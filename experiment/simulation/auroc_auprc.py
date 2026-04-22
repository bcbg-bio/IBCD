import os
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

def compute_auroc_auprc_from_mats(X_est, X_true, tol=1e-8):
    X_est  = np.asarray(X_est, float)
    X_true = np.asarray(X_true, float)
    assert X_est.shape == X_true.shape

    D = X_true.shape[0]
    mask = ~np.eye(D, dtype=bool)
    scores = np.abs(X_est)[mask].ravel()
    y_true = (np.abs(X_true)[mask] > tol).astype(int)

    if y_true.min() == y_true.max():
        return np.nan, np.nan

    auroc = roc_auc_score(y_true, scores)
    auprc = average_precision_score(y_true, scores)
    return auroc, auprc

def parse_seeds(s):
    s = s.strip()
    if ":" in s:
        a, b = s.split(":")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True
                        , help="Method name")
    parser.add_argument("--g_true", required=True,
                        help="Path to true G, use {seed}")
    parser.add_argument("--g_est", required=True,
                        help="Path to estimated G, use {seed}")
    parser.add_argument("--seeds", default="42:51",
                        help="Seed range, e.g. 42:51")
    parser.add_argument("--out_csv", default="summary.csv",
                        help="Output csv file")

    args = parser.parse_args()
    seeds = parse_seeds(args.seeds)

    aurocs, auprcs = [], []

    for seed in seeds:
        true_path = args.g_true.replace("{seed}", str(seed))
        est_path  = args.g_est.replace("{seed}", str(seed))

        X_true = pd.read_csv(true_path).values
        X_est  = pd.read_csv(est_path).values

        auroc, auprc = compute_auroc_auprc_from_mats(X_est, X_true)
        aurocs.append(auroc)
        auprcs.append(auprc)

    row = pd.DataFrame([{
        "method": args.method,
        "AUROC_mean": np.mean(aurocs),
        "AUROC_sd": np.std(aurocs),
        "AUPRC_mean": np.mean(auprcs),
        "AUPRC_sd": np.std(auprcs),
        "n_seeds": len(seeds)
    }])

    if os.path.exists(args.out_csv):
        row.to_csv(args.out_csv, mode="a", header=False, index=False)
    else:
        row.to_csv(args.out_csv, index=False)

    print(f"Updated summary: {args.out_csv}")

if __name__ == "__main__":
    main()
