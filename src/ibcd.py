import os
import argparse
import json
import time
import warnings
import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
from numpyro import infer
from numpyro.infer import MCMC, NUTS

from empirical_prior import (
    scale_free_degree,
    solve_edge_weights_rowwise,
    load_R_and_SE_hat,
    empirical_bayes_em,
    em_slab_scale,
    solve_spike_slab_diagonal_spike,
)
from model import (matrix_model_spike_horseshoe, compute_lfsr, draw_spectral_radius,
                   posterior_diagnostics, optimized_init, report_spectral_radius)
from iv_regression import xi_norm, run_all_IV, compute_S_hat


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

    em_pi_k = em_sigma_k = None          # the EM's slab fit, if it is run
    if (args.oracle_G or args.global_prior) and args.prior.lower() != "er":
        raise ValueError("--oracle_G and --global_prior are ablations of the ER prior; use --prior er")
    if args.prior.lower() == "sf":
        # -------- Scale-free (SF) prior --------
        print("2) Using SF prior (Scale-Free)...")
        R = Rhat_df.values
        # theta = sum_j |R_hat_ij|^2 is in units of squared effect size, not
        # degree, so equation 13 fixes the shape of the degree profile but not
        # the overall sparsity. The EM of equation 12 supplies that level, the
        # same one equation 18 gives the ER prior; theta then only decides how
        # it is shared across rows.
        # --sf_anchor none restores the published max-normalisation, which
        # has no sparsity level of its own.
        if args.sf_anchor == "em":
            w, se_hat = load_R_and_SE_hat(Rhat_path, SE_hat_path)
            pi0_global, em_pi_k, em_sigma_k = empirical_bayes_em(w, se_hat, alpha_er=args.alpha_er)
            pi0_i = scale_free_degree(R, pi0_global=float(pi0_global),
                                      pi0_floor=args.pi0_floor)
            print("Estimated global spike weight:", float(pi0_global))
        else:
            pi0_i = scale_free_degree(R)
        print("Estimated spike weight:", pi0_i)
        print("3) Running edge specific weights for SF...")
        pi0_ij, pi_k_ij = solve_edge_weights_rowwise(
            xi,
            pi0_i,
        )

    elif args.prior.lower() == "er":
        # -------- Erdős–Rényi (ER) prior --------
        print("2) Using ER prior (Erdős–Rényi)...")
        # Load all off-diagonal entries
        w, se_hat = load_R_and_SE_hat(Rhat_path, SE_hat_path)
        #print(f"Shape of w: {w.shape}, Shape of se: {se_hat.shape}")
        pi0, pi_slabs, slab_scales = empirical_bayes_em(
            w,
            se_hat,
            alpha_er=args.alpha_er,
        )
        em_pi_k, em_sigma_k = pi_slabs, slab_scales
        print("Estimated spike weight:", float(pi0))
        #print("Sum pi:", float(pi0 + pi_slabs.sum()))
        if args.oracle_G:
            # Appendix E: the spike weight is the true graph's share of zeros
            G_true = pd.read_csv(args.oracle_G).values.astype(float)
            off = ~np.eye(G_true.shape[0], dtype=bool)
            pi0 = float(np.mean(np.abs(G_true[off]) < 1e-8))
            print("Oracle spike weight from the true graph:", pi0)

        if args.global_prior:
            # Table 7's global prior: one spike weight for every edge
            print("3) Global prior: the same spike weight for every edge...")
            pi0_ij = np.full((D, D), float(pi0))
        else:
            print("3) Running edge specific weights for ER...")
            pi0_ij, pi_k_ij = solve_spike_slab_diagonal_spike(xi, pi0=pi0)

    else:
        raise ValueError("args.prior must be 'sf' or 'er'.")

    slab_width = None
    if args.slab_width == "em":
        if em_pi_k is None:               # SF with --sf_anchor none skips the EM
            w, se_hat = load_R_and_SE_hat(Rhat_path, SE_hat_path)
            _, em_pi_k, em_sigma_k = empirical_bayes_em(w, se_hat, alpha_er=args.alpha_er)
        slab_width = em_slab_scale(em_pi_k, em_sigma_k)
        print(f"Regularised horseshoe, slab width {slab_width:.4f}")

    U_lower = jnp.linalg.cholesky(jnp.array(U_mat.values))
    V_lower = jnp.linalg.cholesky(jnp.array(V_mat.values))

    print("4) Running inference...")

    kernel = NUTS(
        matrix_model_spike_horseshoe,
        target_accept_prob=args.target_accept_prob,
        max_tree_depth=args.max_tree_depth,
        init_strategy=infer.init_to_median(num_samples=50),
    )

    mcmc = MCMC(
        kernel,
        num_warmup=args.num_warmup,
        num_samples=args.num_samples,
        num_chains=args.num_chains,
        chain_method=args.chain_method,
        progress_bar=True,
    )

    if args.num_chains > 1 and args.chain_method == "parallel":
        n_dev = jax.local_device_count()
        if n_dev < args.num_chains:
            warnings.warn(
                f"chain_method='parallel' with {args.num_chains} chains but only "
                f"{n_dev} visible device(s); numpyro will draw the chains "
                "sequentially. Use --chain_method vectorized to share one "
                "device, or run one chain per job with --num_chains 1 and a "
                "distinct --seed.",
                RuntimeWarning,
            )

    # one set of model arguments, shared by the optimised start and the sampler
    model_kwargs = dict(
        obs_data=Rhat_df.values,
        pi0_ij=pi0_ij,
        U_lower=U_lower,
        V_lower=V_lower,
        D=D,
        rho_penalty=None if args.rho_penalty == "none" else args.rho_penalty,
        rho_sigma=args.rho_sigma,
        slab_width=slab_width,
    )

    if args.likelihood == "mvn":
        # Table 7 ablation: the full covariance of vec(R_hat) (equation 9),
        # truncated to its top --mvn_rank eigenpairs plus a small diagonal
        S = compute_S_hat(df, Rhat_df.values)
        lam_S, Q_S = np.linalg.eigh(S)
        top = np.argsort(lam_S)[::-1][:args.mvn_rank]
        model_kwargs["mvn_factor"] = jnp.array(Q_S[:, top] * np.sqrt(np.maximum(lam_S[top], 0.0)))
        model_kwargs["mvn_diag"] = jnp.full(D * D, args.mvn_jitter)
        print(f"Multivariate-normal likelihood: rank {args.mvn_rank} of S "
              f"({lam_S[top].sum() / lam_S.clip(min=0).sum():.1%} of its trace) "
              f"plus {args.mvn_jitter:g} on the diagonal")

    if args.init_opt_steps > 0:
        print(f"4a) Optimising starting values ({args.init_opt_steps} Adam steps)...")
    t_init = time.perf_counter()
    init_params, init_info = optimized_init(
        matrix_model_spike_horseshoe, model_kwargs,
        jax.random.fold_in(jax.random.PRNGKey(args.seed), 1),
        num_chains=args.num_chains, steps=args.init_opt_steps,
        lr=args.init_opt_lr, jitter=args.init_jitter,
    )
    if args.num_chains == 1:
        init_params = jax.tree_util.tree_map(lambda x: x[0], init_params)
    init_info["seconds"] = time.perf_counter() - t_init
    print("    rho(G) at the start of each chain:",
          [round(r, 3) for r in init_info["rho_start_per_chain"]])

    t_start = time.perf_counter()
    mcmc.run(
        jax.random.PRNGKey(args.seed),
        init_params=init_params,
        extra_fields=("num_steps", "diverging", "accept_prob",
                      "adapt_state.step_size"),
        **model_kwargs,
    )


    # ------------------------ outputs ---------------------------------
    posterior = np.asarray(
        jax.device_get(mcmc.get_samples(group_by_chain=True)["G"]),
        dtype=np.float32,
    )
    elapsed = time.perf_counter() - t_start
    extra = jax.device_get(mcmc.get_extra_fields(group_by_chain=True))
    np.save(os.path.join(args.output_dir, "G_draws.npy"), posterior)

    flat = posterior.reshape(-1, D, D)

    # rho(G) for every draw; every draw is kept and rho is only reported
    rho = draw_spectral_radius(flat)
    report_spectral_radius(rho, posterior.shape[0])

    # Diagnostics are optional: split R-hat and ESS are O(D^2) over the entries
    # of G, which is the expensive part of this block at large D.
    if args.save_diagnostics:
        diagnostics = posterior_diagnostics(posterior, rho,
                                            extra_fields=extra, seed=args.seed,
                                            max_tree_depth=args.max_tree_depth)
        diagnostics["runtime_seconds"] = float(elapsed)
        diagnostics["config"] = {
            "data": args.data,
            "prior": args.prior,
            "seed": args.seed,
            "num_warmup": args.num_warmup,
            "num_samples": args.num_samples,
            "num_chains": args.num_chains,
            "target_accept_prob": args.target_accept_prob,
            "max_tree_depth": args.max_tree_depth,
            "epsilon": args.epsilon,
            "rho_penalty": args.rho_penalty,
            "rho_sigma": args.rho_sigma if args.rho_penalty == "gaussian" else None,
            "sf_anchor": args.sf_anchor if args.prior.lower() == "sf" else None,
            "init_strategy": "optimized",
            "init_opt_steps": args.init_opt_steps,
            "init_opt_lr": args.init_opt_lr,
            "init_jitter": args.init_jitter,
            "slab_width": args.slab_width,
            "slab_width_value": slab_width,
            "oracle_G": args.oracle_G,
            "global_prior": bool(args.global_prior),
            "likelihood": args.likelihood,
            "mvn_rank": args.mvn_rank if args.likelihood == "mvn" else None,
            "mvn_jitter": args.mvn_jitter if args.likelihood == "mvn" else None,
        }
        diagnostics["init"] = init_info
        if args.rho_penalty != "none":
            # what the penalty acted on, against the exact rho recorded above
            est = np.asarray(jax.device_get(mcmc.get_samples()["rho_estimate"])).ravel()
            diagnostics["rho_estimate"] = {
                k: float(v) for k, v in zip(
                    ["min", "median", "p95", "max"],
                    np.percentile(est, [0, 50, 95, 100]))
            }
            diagnostics["rho_estimate"]["median_ratio_to_exact"] = float(
                np.median(est / np.maximum(rho, 1e-12)))
        with open(os.path.join(args.output_dir, "diagnostics.json"), "w") as fh:
            json.dump(diagnostics, fh, indent=2)

        dg = diagnostics
        issues = []
        if dg.get("r_hat", {}).get("max", 0.0) > 1.05:
            issues.append(f"max r_hat {dg['r_hat']['max']:.3f} > 1.05")
        if dg.get("divergences", {}).get("pct", 0.0) > 5.0:
            issues.append(f"{dg['divergences']['pct']:.1f}% of transitions diverged")
        if dg.get("ess", {}).get("n_nonpositive_raw", 0) > 0:
            issues.append(
                f"{dg['ess']['n_nonpositive_raw']} entries had a non-positive ESS estimate"
            )
        if issues:
            warnings.warn(
                "Sampler did not converge cleanly: " + "; ".join(issues)
                + ". See diagnostics.json.",
                RuntimeWarning,
            )

        print(
            "5) Diagnostics: "
            f"divergences {dg.get('divergences', {}).get('n', 'NA')}"
            f" ({dg.get('divergences', {}).get('pct', float('nan')):.1f}%), "
            f"max r_hat {dg.get('r_hat', {}).get('max', float('nan')):.4f}, "
            f"min ESS {dg.get('ess', {}).get('min', float('nan')):.0f}, "
            f"step size {dg.get('step_size', {}).get('median', float('nan')):.2e}, "
            f"accept {dg.get('accept_prob', {}).get('mean', float('nan')):.2f}, "
            f"rho median {dg['spectral_radius']['median']:.3f} max {dg['spectral_radius']['max']:.3g}, "
            f"{elapsed:.1f}s"
        )

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
        "--output_dir",
        required=True,
        help="Directory to save all outputs.",
    )

    parser.add_argument(
        "--sf_anchor",
        choices=["em", "none"],
        default="em",
        help=(
            "SF prior only: where the overall sparsity level comes from. 'em' "
            "takes the global spike weight from the EM of equation 12 and lets "
            "the degree profile share it out across rows. 'none' is the "
            "published max-normalisation, pi0_i = 1 - theta_i / max(theta, "
            "phi), which has no level of its own. Default em."
        ),
    )
    parser.add_argument(
        "--pi0_floor",
        type=float,
        default=0.05,
        help=("SF prior only: smallest per-node spike proportion, so no row is "
              "left entirely unshrunk. The row budget is capped at "
              "(1 - pi0_floor) * (D - 1) and the excess redistributed."),
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
        default=1000,
        help="Number of NUTS warm-up iterations. Default = 1000.",
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
        "--target_accept_prob",
        type=float,
        default=0.9,
        help=(
            "NUTS target acceptance probability. Raising it shrinks the "
            "adapted step size, which reduces divergences at the cost of "
            "longer trajectories. Default = 0.9."
        ),
    )

    parser.add_argument(
        "--max_tree_depth",
        type=int,
        default=12,
        help=(
            "Maximum NUTS tree depth; each iteration costs at most "
            "2^depth - 1 leapfrog steps. Default = 12."
        ),
    )

    parser.add_argument(
        "--chain_method",
        choices=["parallel", "sequential", "vectorized"],
        default="vectorized",
        help=(
            "How to draw multiple chains. 'vectorized' maps them onto one "
            "device, which is the only form of within-process parallelism "
            "available when CUDA exposes a single device, and avoids the "
            "post-hoc stack that 'sequential' pays for. 'parallel' needs one "
            "visible device per chain and falls back to sequential otherwise. "
            "Default = vectorized."
        ),
    )

    parser.add_argument(
        "--save_diagnostics",
        action="store_true",
        help=(
            "Write diagnostics.json to the output directory: divergences, "
            "leapfrog steps, split R-hat, ESS, spectral-radius quantiles, "
            "rejected-draw counts, runtime and the run configuration. Off by "
            "default because R-hat and ESS are O(D^2) over the entries of G."
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="PRNG seed for MCMC. Default = 42.",
    )

    parser.add_argument(
        "--rho_penalty",
        choices=["none", "gaussian"],
        default="none",
        help=(
            "Soft acyclicity constraint on G: 'gaussian' is the Appendix H "
            "prior N(0, rho_sigma^2) on the spectral radius of G, estimated by "
            "50 steps of power iteration. Default none."
        ),
    )

    parser.add_argument(
        "--rho_sigma",
        type=float,
        default=0.5,
        help="Scale of the Gaussian rho penalty. Default 0.5, the Appendix H value.",
    )

    parser.add_argument(
        "--init_opt_steps",
        type=int,
        default=2000,
        help=("NUTS starts from init_to_median(num_samples=50), an essentially "
              "empty G, after this many Adam steps on the log posterior, which "
              "leave the posterior unchanged; 0 starts NUTS at the (jittered) "
              "init_to_median draws themselves. Default 2000."),
    )

    parser.add_argument(
        "--slab_width",
        choices=["none", "em"],
        default="em",
        help=(
            "'none' for the horseshoe as published; 'em' for the regularised "
            "horseshoe with its slab width set to the RMS scale of the EM's "
            "fitted slab. The regularised slab's tail beyond the width is "
            "Gaussian rather than Cauchy. Default em."
        ),
    )

    parser.add_argument(
        "--init_opt_lr",
        type=float,
        default=0.01,
        help="Adam learning rate for the descent before NUTS. Default 0.01.",
    )

    parser.add_argument(
        "--init_jitter",
        type=float,
        default=0.1,
        help=("Noise added to log(lam) and eps after optimising, so the "
              "chains do not start at one point. Default 0.1."),
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

    parser.add_argument(
        "--oracle_G",
        default=None,
        help=("Ablation (Table 7, Appendix E): CSV of the true G; the ER prior's "
              "spike weight is set to its share of zero off-diagonal entries "
              "instead of the EM estimate. ER prior only."),
    )

    parser.add_argument(
        "--global_prior",
        action="store_true",
        help=("Ablation (Table 7): give every edge the EM's global spike weight "
              "instead of the edge-specific weights. ER prior only."),
    )

    parser.add_argument(
        "--likelihood",
        choices=["mn", "mvn"],
        default="mn",
        help=("'mn' is the matrix-normal likelihood. 'mvn' is the Table 7 "
              "ablation: a multivariate normal on vec(R_hat) whose covariance is "
              "the full covariance S of equation 9 truncated to its top "
              "--mvn_rank components. S is D^2 x D^2, so small D only. Default mn."),
    )

    parser.add_argument(
        "--mvn_rank",
        type=int,
        default=10,
        help="Components of S kept for --likelihood mvn. Default 10, as in the paper.",
    )

    parser.add_argument(
        "--mvn_jitter",
        type=float,
        default=1e-5,
        help=("Diagonal added to the truncated S for --likelihood mvn, which is "
              "otherwise singular. Default 1e-5."),
    )

    args = parser.parse_args()
    main(args)
