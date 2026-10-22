# generate_from_csv.py
# -----------------------------------------------------------------------------
# Utility script to convert a *single* tabular CSV file into the same directory
# structure expected by BayesDAG.
#
# It creates two sibling folders in the *current* directory:
#   ├── <dataset_name>            (mean‑centered, **not** standardised)
#   └── <dataset_name>_std        (mean‑centred **and** standardised)
# and writes the following artefacts into each one:
#   ├── adj_matrix.csv            (ground‑truth DAG, here a zero matrix placeholder)
#   ├── all.csv                   (all rows from your CSV)
#   ├── train.csv                 (first N rows – see --train‑ratio)
#   ├── test.csv                  (remaining rows)
#   ├── held_out_interventions.pkl  -> `None` (kept for API compatibility)
#   ├── true_posterior.pkl          -> `None`
#   └── graph_args.pkl              -> metadata dictionary
# -----------------------------------------------------------------------------
# Example
# -------
#   python generate_from_csv.py --csv Y_matrix.csv --train-ratio 0.8 --seed 42
# -----------------------------------------------------------------------------

import argparse
import json
import pickle as pkl
from pathlib import Path

import numpy as np
import pandas as pd


def save_data(savedir: Path, adj_matrix: np.ndarray, X: np.ndarray, num_samples_train: int, graph_args: dict) -> None:
    """Write the split datasets and metadata expected by downstream pipelines."""
    savedir.mkdir(parents=True, exist_ok=True)

    # Core tables ----------------------------------------------------------------
    np.savetxt(savedir / "adj_matrix.csv", adj_matrix, delimiter=",", fmt="%i")
    np.savetxt(savedir / "all.csv", X, delimiter=",")
    np.savetxt(savedir / "train.csv", X[:num_samples_train, :], delimiter=",")
    np.savetxt(savedir / "test.csv", X[num_samples_train:, :], delimiter=",")

    # Pickled placeholders -------------------------------------------------------
    (savedir / "held_out_interventions.pkl").write_bytes(pkl.dumps(None))
    (savedir / "true_posterior.pkl").write_bytes(pkl.dumps(None))
    (savedir / "graph_args.pkl").write_bytes(pkl.dumps(graph_args))


def main():
    parser = argparse.ArgumentParser(description="Convert a CSV file into GraN‑DAG format.")
    parser.add_argument("--csv", type=Path, required=True, help="Path to the source CSV file.")
    parser.add_argument(
        "--train-ratio",
        type=float,
        default=0.8,
        help="Fraction of rows to use for the training split (0 < r < 1).",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed for reproducibility.")
    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # 1. Load the data ---------------------------------------------------------
    # -------------------------------------------------------------------------
    df = pd.read_csv(args.csv)
    X_raw = df.values.astype(np.float32)

    num_samples = X_raw.shape[0]
    num_train = int(round(args.train_ratio * num_samples))

    # -------------------------------------------------------------------------
    # 2. Create a *placeholder* adjacency matrix. -----------------------------
    # -------------------------------------------------------------------------
    # If you *know* the ground‑truth DAG, replace this with your own matrix.
    adj_matrix = np.zeros((X_raw.shape[1], X_raw.shape[1]), dtype=int)

    # -------------------------------------------------------------------------
    # 3. Prepare the mean‑centred and standardised versions -------------------
    # -------------------------------------------------------------------------
    mean = X_raw[:num_train].mean(axis=0, keepdims=True)
    std = X_raw[:num_train].std(axis=0, keepdims=True)
    std[std == 0] = 1.0  # Guard against constant columns.

    X_mean_centered = X_raw - mean
    X_standardised = (X_raw - mean) / std

    # -------------------------------------------------------------------------
    # 4. Metadata to mimic GraN‑DAG expectations ------------------------------
    # -------------------------------------------------------------------------
    graph_args = {
        "num_variables": X_raw.shape[1],
        "exp_edges": int(adj_matrix.sum()),
        "seed": args.seed,
        "graph_type": args.csv.stem,
        "exp_edges_per_node": float(adj_matrix.sum() / X_raw.shape[1]),
    }

    # -------------------------------------------------------------------------
    # 5. Persist to disk -------------------------------------------------------
    # -------------------------------------------------------------------------
    # Folder names: <dataset_name> and <dataset_name>_std
    out_root = Path.cwd() / args.csv.stem
    save_data(out_root, adj_matrix, X_mean_centered, num_train, graph_args)
    save_data(out_root.with_name(f"{args.csv.stem}_std"), adj_matrix, X_standardised, num_train, graph_args)

    # -------------------------------------------------------------------------
    # 6. Write the global dataset configuration -------------------------------
    # -------------------------------------------------------------------------
    config = {
        "dataset_format": "csv",
        "use_predefined_dataset": False,
        "test_fraction": 0.1,
        "val_fraction": 0.1,
        "random_seed": [0, 0],
        "negative_sample": False,
    }
    config_path = out_root / "dataset_config.json"
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    print(
        f"✓  Converted '{args.csv.name}' -> '{out_root.name}/' & '{out_root.name}_std/' \n"
        f"   Train samples : {num_train}\n"
        f"   Test  samples : {num_samples - num_train}\n"
        f"   Variables     : {X_raw.shape[1]}"
    )


if __name__ == "__main__":
    main()
