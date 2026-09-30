import numpyro
import numpyro.distributions as dist
from numpyro.diagnostics import split_gelman_rubin, effective_sample_size
import jax.numpy as jnp
import arviz as az
import numpy as np

def power_iteration_radius(G, max_iter=50):
    """Spectral radius estimate by power iteration, as in the Appendix H penalty.

    A line-for-line port of `spectral_radius` in commit 42660bd on main, kept
    unchanged so the published penalty can be tested as written. It estimates
    |lambda_max| rather than bounding it, and when the dominant eigenvalues
    are a complex-conjugate pair the iterate rotates instead of converging.
    """
    v = jnp.ones(G.shape[0])
    for _ in range(max_iter):
        v = G @ v
        v = v / (jnp.linalg.norm(v) + 1e-8)
    return jnp.linalg.norm(G @ v)


def gelfand_radius(G, n_squarings=6):
    """Upper bound on the spectral radius from Gelfand's formula.

    rho(G) <= ||G^k||_F^(1/k), tightening as k grows; k = 2^n_squarings. G^k
    is formed by repeated squaring, renormalising at every step and summing
    the logs, so it neither overflows at large rho nor underflows for a
    near-nilpotent G. Uses only matmuls, so it runs on GPU and has a smooth
    gradient, unlike an eigendecomposition of a nonsymmetric matrix.
    """
    tiny = 1e-30
    A = G
    log_r = 0.0
    for i in range(n_squarings):
        s = jnp.sqrt(jnp.sum(A * A) + tiny)
        log_r = log_r + jnp.log(s) * 2.0 ** (-i)
        A = A / s
        A = A @ A
    log_r = log_r + jnp.log(jnp.sqrt(jnp.sum(A * A) + tiny)) * 2.0 ** (-n_squarings)
    return jnp.exp(log_r)


def rho_log_penalty(rho, kind, sigma=0.5, start=0.9, width=0.05):
    """Log-density term that discourages a large spectral radius.

    'gaussian' is the Appendix H prior, log N(rho; 0, sigma^2). It acts at
    every rho, so it also shrinks cycles well inside rho < 1. 'barrier' is
    exactly zero below `start` and grows as ((rho - start) / width)^2 above
    it, so it leaves the interior alone and enforces the rho < 1 that the
    path-sum definition of R requires. The hinge is squared so the gradient
    stays continuous, which NUTS needs.
    """
    if kind == "gaussian":
        return dist.Normal(0.0, sigma).log_prob(rho)
    if kind == "barrier":
        return -jnp.square(jnp.maximum(rho - start, 0.0) / width)
    raise ValueError(f"unknown rho penalty '{kind}'")


def matrix_model_spike_horseshoe(obs_data, pi0_ij, U_lower, V_lower, D, sigma0=0.001, tau=0.1,
                                 epsilon=1e-5, truncated_series=False, series_order=24,
                                 rho_penalty=None, rho_estimator="power", rho_sigma=0.5,
                                 rho_barrier_start=0.9, rho_barrier_width=0.05):
    """
    NumPyro model: Spike-and-horseshoe prior over matrix G, MatrixNormal likelihood.

    Parameters:
        obs_data (array): R_hat matrix.
        pi0_ij (array): Edge-specific spike weights.
        U_lower (array): Cholesky of row covariance.
        V_lower (array): Cholesky of column covariance.
        D (int): Dimension.
        sigma0 (float): Std for spike component.
        tau (float): Global scale for horseshoe slab.
        epsilon (float): Regularization for stability, used only for the
            explicit inverse.
        truncated_series (bool): Compute R as a truncated path sum instead of
            inverting (I - G). The sum is a polynomial in G, so it has no pole
            and bounded gradients, but it costs roughly an order of magnitude
            more per gradient. Equivalent to the inverse for a DAG once
            series_order reaches the longest directed path.
        series_order (int): Highest power retained in the truncated path sum.
            Ignored unless truncated_series is set.
        rho_penalty (str): None for no constraint on rho(G), 'gaussian' for
            the Appendix H prior N(0, rho_sigma^2), or 'barrier' for a term
            that is zero below rho_barrier_start. See `rho_log_penalty`.
        rho_estimator (str): 'power' (Appendix H's power iteration) or
            'gelfand' (an upper bound). Ignored unless rho_penalty is set.
        rho_sigma (float): Scale of the Gaussian penalty; 0.5 as on main.
        rho_barrier_start (float): Spectral radius where the barrier begins.
        rho_barrier_width (float): Distance over which the barrier costs 1 nat.
    """

    # Sample horseshoe local scales (HalfCauchy), shape (D, D)
    lam = numpyro.sample("lam", dist.HalfCauchy(1.).expand((D, D)), infer={"is_auxiliary": True})

    # Sample standard normal noise, for non-centered parameterization
    eps = numpyro.sample("eps", dist.Normal(0., 1.).expand((D, D)), infer={"is_auxiliary": True})

    # Slab values: non-centered horseshoe (continuous)
    slab_vals = tau * lam * eps

    # Spike values: zero-mean Normal with tiny variance
    spike_vals = numpyro.sample("spike", dist.Normal(0., sigma0).expand((D, D)), infer={"is_auxiliary": True})

    # Binary mixture via convex combination (no discrete Cat). The diagonal is
    # zeroed: the model has no self-loops, R_ii is fixed at 1 by definition,
    # and a nonzero diagonal would stop a DAG's G from being nilpotent, so the
    # path sum below would not terminate exactly.
    G = pi0_ij * spike_vals + (1. - pi0_ij) * slab_vals
    G = G * (1. - jnp.eye(D))
    numpyro.deterministic("G", G)           # keep only G

    if rho_penalty:
        if rho_estimator == "power":
            rho_est = power_iteration_radius(G)
        elif rho_estimator == "gelfand":
            rho_est = gelfand_radius(G)
        else:
            raise ValueError(f"unknown rho estimator '{rho_estimator}'")
        # recorded so each run can compare the estimate the penalty acted on
        # with the exact spectral radius computed afterwards
        numpyro.deterministic("rho_estimate", rho_est)
        numpyro.factor("spectral_radius_penalty",
                       rho_log_penalty(rho_est, rho_penalty, sigma=rho_sigma,
                                       start=rho_barrier_start,
                                       width=rho_barrier_width))

    # MatrixNormal mean. R is the sum over directed paths, sum_d G^d, which
    # equals (I - G)^-1 for a DAG. Accumulating the truncated sum by Horner
    # keeps R polynomial in G, so it has no pole and bounded gradients.
    I = jnp.eye(D)
    if truncated_series:
        R_mean = I
        for _ in range(series_order):
            R_mean = I + G @ R_mean
    else:
        I_minus_G = I - G + epsilon * I
        R_mean = jnp.linalg.solve(I_minus_G, I)

    # MatrixNormal likelihood
    numpyro.sample(
        "R_hat_obs",
        dist.MatrixNormal(
            loc=R_mean,
            scale_tril_row=U_lower[:D, :D],
            scale_tril_column=V_lower[:D, :D],
        ),
        obs=obs_data[:D, :D],
    )

def convergent_draws(G_draws, max_spectral_radius=1.0):
    """Flag posterior draws of G whose implied total-effect matrix exists.

    R is defined as the sum over directed paths, sum_d G^d, which converges
    only when the spectral radius of G is below one. Draws outside that region
    have left the part of the parameter space the model is defined on; they
    arise when the sampler crosses the singularity of (I - G).

    Parameters:
        G_draws (array): Posterior draws of G, shape (..., D, D).
        max_spectral_radius (float): Exclusive upper bound on rho(G).

    Returns:
        keep (np.ndarray): Boolean mask over flattened draws.
        rho (np.ndarray): Spectral radius of each draw.
    """
    g = np.asarray(G_draws)
    flat = g.reshape(-1, g.shape[-2], g.shape[-1])
    rho = np.abs(np.linalg.eigvals(flat)).max(axis=1)
    return rho < max_spectral_radius, rho


def posterior_diagnostics(G_draws, rho, keep, extra_fields=None,
                          max_ess_entries=5000, seed=0, max_tree_depth=12,
                          max_rhat_entries=20000):
    """Summarise sampler behaviour and draw validity for one run.

    Convergence statistics are computed over the entries of G. ESS is
    estimated on a random subsample of entries, since an autocorrelation
    estimate per entry is O(D^2) and D can be 500.

    Parameters:
        G_draws (array): Posterior draws, shape (chains, draws, D, D).
        rho (array): Spectral radius of each flattened draw.
        keep (array): Boolean mask over flattened draws.
        extra_fields (dict): numpyro extra fields grouped by chain; the keys
            'diverging' and 'num_steps' are used when present.
        max_ess_entries (int): Cap on how many entries of G enter the ESS
            estimate.
        seed (int): Seed for choosing that subsample.
        max_tree_depth (int): The sampler's tree-depth limit, used to work out
            how many leapfrog steps a saturated iteration takes.
        max_rhat_entries (int): Cap on how many entries of G enter the split
            R-hat estimate. Bounds peak memory, which is what this costs at
            large D; the estimate itself is cheap.

    Returns:
        dict: JSON-serialisable diagnostics.
    """
    g = np.asarray(G_draws)
    n_chains, n_draws, D, _ = g.shape
    keep = np.asarray(keep)
    rho = np.asarray(rho)

    out = {
        "n_chains": int(n_chains),
        "n_draws_per_chain": int(n_draws),
        "n_draws_total": int(keep.size),
        "D": int(D),
    }

    # draws outside the region where R = sum_d G^d converges
    per_chain = (~keep).reshape(n_chains, n_draws).sum(axis=1)
    out["nonconvergent"] = {
        "n": int((~keep).sum()),
        "pct": float(100.0 * (~keep).mean()),
        "per_chain": [int(x) for x in per_chain],
    }
    out["spectral_radius"] = {
        k: float(v) for k, v in zip(
            ["min", "median", "p95", "p99", "max"],
            np.percentile(rho, [0, 50, 95, 99, 100]),
        )
    }

    # convergence over the entries of G, on the retained draws only
    offdiag = ~np.eye(D, dtype=bool)
    kept = keep.reshape(n_chains, n_draws)
    usable = kept.all(axis=1)
    if usable.sum() >= 2:
        # Index the entries before the chains: g[usable][:, :, offdiag] would
        # materialise the whole off-diagonal block, which is 3 GB at D = 500 on
        # top of the draws themselves. Subsampling first keeps it to a few
        # hundred MB, and r_hat over a random subset of entries answers the
        # same question as r_hat over all of them.
        rng = np.random.default_rng(seed)
        off_idx = np.flatnonzero(offdiag)
        n_rhat = min(max_rhat_entries, off_idx.size)
        r_idx = np.sort(rng.choice(off_idx, size=n_rhat, replace=False))
        sub = g.reshape(n_chains, n_draws, D * D)[:, :, r_idx][usable]
        try:
            rhat = np.asarray(split_gelman_rubin(sub))
            out["r_hat"] = {
                "max": float(np.nanmax(rhat)),
                "median": float(np.nanmedian(rhat)),
                "frac_above_1_01": float(np.nanmean(rhat > 1.01)),
                "n_entries_used": int(n_rhat),
                "n_entries_total": int(off_idx.size),
            }
        except Exception as exc:          # noqa: BLE001
            out["r_hat"] = {"error": str(exc)}
        idx = rng.choice(sub.shape[-1], size=min(max_ess_entries, sub.shape[-1]),
                         replace=False)
        try:
            ess = np.asarray(effective_sample_size(sub[..., idx]))
            # The autocorrelation estimator can return non-positive values when
            # the chains are short or badly mixed. Those carry no information,
            # so floor them at zero and record how many there were.
            n_bad = int(np.sum(~(ess > 0)))
            out["ess"] = {
                "min": float(np.nanmin(np.maximum(ess, 0.0))),
                "median": float(np.nanmedian(np.maximum(ess, 0.0))),
                "n_entries_used": int(idx.size),
                "n_nonpositive_raw": n_bad,
            }
        except Exception as exc:          # noqa: BLE001
            out["ess"] = {"error": str(exc)}
    else:
        out["r_hat"] = {"note": "fewer than two chains free of non-convergent draws"}
        out["ess"] = {"note": "fewer than two chains free of non-convergent draws"}

    if extra_fields:
        # Step size and achieved acceptance are how a target_accept_prob change
        # shows up: raising the target shrinks the step, which trades
        # divergences for longer trajectories.
        if "adapt_state.step_size" in extra_fields:
            ss = np.asarray(extra_fields["adapt_state.step_size"]).reshape(n_chains, -1)
            out["step_size"] = {
                "per_chain": [float(x) for x in ss[:, -1]],
                "median": float(np.median(ss[:, -1])),
            }
        if "accept_prob" in extra_fields:
            ap = np.asarray(extra_fields["accept_prob"]).reshape(n_chains, -1)
            out["accept_prob"] = {
                "mean": float(np.nanmean(ap)),
                "per_chain": [float(x) for x in np.nanmean(ap, axis=1)],
            }
        if "diverging" in extra_fields:
            dv = np.asarray(extra_fields["diverging"])
            out["divergences"] = {
                "n": int(dv.sum()),
                "pct": float(100.0 * dv.mean()),
                "per_chain": [int(x) for x in dv.reshape(n_chains, -1).sum(axis=1)],
            }
        if "num_steps" in extra_fields:
            ns = np.asarray(extra_fields["num_steps"])
            # A saturated iteration takes 2^max_tree_depth - 1 steps, so the
            # cap depends on the sampler setting and cannot be hardcoded.
            cap = 2 ** int(max_tree_depth) - 1
            out["leapfrog"] = {
                "total": int(ns.sum()),
                "mean_per_iter": float(ns.mean()),
                "max_leapfrog_steps": int(cap),
                "pct_at_max_tree_depth": float(100.0 * (ns >= cap).mean()),
            }
    return out


def compute_lfsr(flat_samples):
    """
    Compute Local Falscale_free_degree Sign Rate for each edge from posterior samples.

    Parameters:
        flat_samples (array): Posterior samples of shape (samples, D, D)

    Returns:
        array: LFSR matrix of shape (D, D)
    """
    pos_probs = np.mean(flat_samples > 0, axis=0)
    neg_probs = np.mean(flat_samples < 0, axis=0)
    return np.minimum(pos_probs, neg_probs)

def plot_posterior_(idata, ax=None):
    """
    Plot posterior distribution of G from inference data.

    Parameters:
        idata (InferenceData): ArviZ-formatted object.
        ax (matplotlib axis): Optional axis to draw on.
    """
    if ax is None:
        az.plot_posterior(idata.posterior.G)
