import os
import argparse
import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import cvxpy as cp
from numpyro import infer
from numpyro.infer import MCMC, NUTS

from empirical_prior import (
    scale_free_degree,
    solve_edge_weights_rowwise,
    load_R_and_SE_hat,
    empirical_bayes_em,
    solve_spike_slab_diagonal_spike,
)
from model import matrix_model_spike_horseshoe, compute_lfsr
from iv_regression import xi_norm, run_all_IV


def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.data)

    xi_norm(df, out_dir=args.output_dir)

    #print("2) Running 2SLS IV regressions...")
    run_all_IV(df, out_dir=args.output_dir)

    # 2SLS outputs
    Rhat_path = f"{args.output_dir}/R.csv"
    SE_hat_path = f"{args.output_dir}/SE_hat.csv"
    xi_path = f"{args.output_dir}/xi.csv"
    U_path = f"{args.output_dir}/U.csv"
    V_path = f"{args.output_dir}/V.csv"

    Rhat_df = pd.read_csv(Rhat_path)
    colnames = Rhat_df.columns.tolist()

    xi = pd.read_csv(xi_path).values
    U_mat = pd.read_csv(U_path)
    V_mat = pd.read_csv(V_path)

    D = xi.shape[0]

    #dag_flag = (args.dag.lower() == "true")
    causal_order_flag = (args.causal_order.lower() == "true")

    if args.prior.lower() == "sf":
        # -------- Scale-free (SF) prior --------
        print("2) Using SF prior (Scale-Free)...")
        R = Rhat_df.values
        pi0_mat = scale_free_degree(R)
        offdiag_mask = ~np.eye(D, dtype=bool)
        pi0_i = pi0_mat[offdiag_mask].reshape(D, D-1).mean(axis=1)

        print("Estimated spike weight:", pi0_i)
        print("3) Running edge specific weights for SF...")
        pi0_ij, pi_k_ij = solve_edge_weights_rowwise(
            xi,
            pi0_i,
            causal_order=causal_order_flag,
            alpha_sf=args.alpha_sf,
            solver=cp.ECOS,
        )

    elif args.prior.lower() == "er":
        # -------- Erdős–Rényi (ER) prior --------
        print("2) Using ER prior (Erdős–Rényi)...")
        # Load off-diagonal or upper-tri entries depending on DAG
        w, se_hat = load_R_and_SE_hat(Rhat_path, SE_hat_path, causal_order=causal_order_flag)
        #print(f"Shape of w: {w.shape}, Shape of se: {se_hat.shape}")
        pi0, pi_slabs, slab_scales = empirical_bayes_em(
            w,
            se_hat,
            alpha_er=args.alpha_er,
        )
        print("Estimated spike weight:", float(pi0))
        #print("Sum pi:", float(pi0 + pi_slabs.sum()))

        print("3) Running edge specific weights for ER...")
        pi0_ij, pi_k_ij, _ = solve_spike_slab_diagonal_spike(xi, pi0=pi0)

    else:
        raise ValueError("args.prior must be 'sf' or 'er'.")

    U_lower = jnp.linalg.cholesky(jnp.array(U_mat.values))
    V_lower = jnp.linalg.cholesky(jnp.array(V_mat.values))

    print("4) Running inference...")

    kernel = NUTS(
        matrix_model_spike_horseshoe,
        target_accept_prob=0.7,
        max_tree_depth=10,
        init_strategy=infer.init_to_median(num_samples=3),
    )

    mcmc = MCMC(
        kernel,
        num_warmup=args.num_warmup,
        num_samples=args.num_samples,
        num_chains=args.num_chains,
        progress_bar=True,
    )

    mcmc.run(
        jax.random.PRNGKey(42),
        obs_data=Rhat_df.values,
        pi0_ij=pi0_ij,
        U_lower=U_lower,
        V_lower=V_lower,
        D=D,
    )


    # ------------------------ outputs ---------------------------------
    posterior = np.asarray(
        jax.device_get(mcmc.get_samples(group_by_chain=True)["G"]),
        dtype=np.float32,
    )
    np.save(os.path.join(args.output_dir, "G_draws.npy"), posterior)

    flat = posterior.reshape(-1, D, D)

    posterior_mean = flat.mean(axis=0)
    pd.DataFrame(posterior_mean, columns=colnames, index=colnames).to_csv(
        f"{args.output_dir}/G.csv", index=True
    )    

    pip_matrix = (np.abs(flat) > args.epsilon).mean(axis=0)
    pd.DataFrame(pip_matrix, columns=colnames, index=colnames).to_csv(
        f"{args.output_dir}/pip.csv", index=True
    )

    lfsr_matrix = compute_lfsr(flat)
    pd.DataFrame(lfsr_matrix, columns=colnames, index=colnames).to_csv(
        f"{args.output_dir}/lfsr.csv", index=True
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "IBCD pipeline.\n"
            "1) Load data.csv (observation + intervention)\n"
            "2) Run 2SLS\n"
            "3) Choose SF (scale-free) or ER (Erdős–Rényi) empirical prior\n"
            "4) Fit empirical Bayesian spike-and-slab prior on matrix normal model\n"
            "5) Output G, PIP, and LFSR.\n"
        )
    )

    parser.add_argument(
        "--data",
        required=True,
        help="Path to input data CSV.",
    )

    parser.add_argument(
        "--prior",
        required=True,
        choices=["sf", "er"],
        help="Choice of empirical prior: 'sf' = scale-free, 'er' = Erdős–Rényi.",
    )

    parser.add_argument(
        "--causal_order",
        default="false",
        choices=["true", "false"],
        help="'false' = input variables are unordered, as in most real data (e.g. perturb-seq); 'true' = input variables are already in causal/topological order.",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save all outputs.",
    )

    parser.add_argument(
        "--alpha_sf",
        type=float,
        default=1.0,
        help="Penalty parameter for SF prior (scale-free row-wise optimization). Default=1.0.",
    )

    parser.add_argument(
        "--alpha_er",
        type=float,
        default=2.0,
        help="Alpha for EM in ER prior. Controls shrinkage strength. Default=2.0.",
    )

    parser.add_argument(
        "--num_warmup",
        type=int,
        default=300,
        help="Number of NUTS warm-up iterations. Default = 300.",
    )

    parser.add_argument(
        "--num_samples",
        type=int,
        default=1000,
        help="Number of posterior samples per chain after warm-up. Default = 1000.",
    )

    parser.add_argument(
        "--num_chains",
        type=int,
        default=3,
        help="Number of parallel MCMC chains. Default = 3.",
    )

    parser.add_argument(
        "--epsilon",
        type=float,
        default=0.05,
        help=(
            "Threshold for computing PIP: edges with |G| > epsilon are counted as active. "
            "Default = 0.05."
        ),
    )

    args = parser.parse_args()
    main(args)
