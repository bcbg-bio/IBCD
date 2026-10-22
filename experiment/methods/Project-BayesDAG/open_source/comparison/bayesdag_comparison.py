"""This script compares the results of the BayesDAG algorithm with the ground truth matrix."""

import pandas as pd
from calc_metrics import calc_metrics

# Read Ground Truth Matrix
df = pd.read_csv("G_matrix.csv")

# Rename columns
df.columns = range(df.shape[1])

# Convert to binary matrix
X_true = (df > 0).astype(int) + (df.T < 0).astype(int)

# Read BayesDAG results
df = pd.read_csv("best_edge_probabilities_r1.csv")

# First binary comparison
X = (df > 0.50).astype(int) # Convert to binary matrix

# Compare and save results
metrics = calc_metrics(X.values, X_true.values)
pd.DataFrame([metrics]).to_csv("bayesdag_results_r1!.csv", index=False)

# Compare probability results and save it
metrics = calc_metrics(df.values, X_true.values)
pd.DataFrame([metrics]).to_csv("bayesdag_probability_results_r1!.csv", index=False)
