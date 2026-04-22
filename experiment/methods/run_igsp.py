import argparse
import pandas as pd
import numpy as np
import helper_functions
import conditional_independence

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--y_path", required=True, help="Path to Y_matrix.csv")
    parser.add_argument("--t_path", required=True, help="Path to targets.txt")
    parser.add_argument("--out_path", required=True, help="Path to save output matrix")
    args = parser.parse_args()

    X_df = pd.read_csv(args.y_path)
    colnames = X_df.columns.to_list()
    X = X_df.values

    targets = np.squeeze(pd.read_csv(args.t_path, header=None).values)

    X_dt = helper_functions.make_dotears_data(X, targets, colnames)
    G = helper_functions.run_igsp(X_dt)

    pd.DataFrame(G).to_csv(args.out_path, index=False)
