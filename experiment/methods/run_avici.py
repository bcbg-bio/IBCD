import argparse
import pandas as pd
import avici

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--y_path", required=True, help="Path to Y_matrix.csv")
    parser.add_argument("--out_path", required=True, help="Path to save output matrix")
    args = parser.parse_args()

    X = pd.read_csv(args.y_path).values

    model = avici.load_pretrained(download="scm-v0")
    G_prob = model(x=X)

    pd.DataFrame(G_prob).to_csv(args.out_path, index=False)
