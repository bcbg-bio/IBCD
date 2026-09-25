import numpyro
import numpyro.distributions as dist
from numpyro.diagnostics import split_gelman_rubin, effective_sample_size
import jax.numpy as jnp
import arviz as az
import numpy as np

def matrix_model_spike_horseshoe(obs_data, pi0_ij, U_lower, V_lower, D, sigma0=0.001, tau=0.1,
                                 epsilon=1e-5, truncated_series=False, series_order=24):
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
                          max_ess_entries=5000, seed=0):
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
        sub = g[usable][:, :, offdiag]
        try:
            rhat = np.asarray(split_gelman_rubin(sub))
            out["r_hat"] = {
                "max": float(np.nanmax(rhat)),
                "median": float(np.nanmedian(rhat)),
                "frac_above_1_01": float(np.nanmean(rhat > 1.01)),
            }
        except Exception as exc:          # noqa: BLE001
            out["r_hat"] = {"error": str(exc)}
        rng = np.random.default_rng(seed)
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
        if "diverging" in extra_fields:
            dv = np.asarray(extra_fields["diverging"])
            out["divergences"] = {
                "n": int(dv.sum()),
                "pct": float(100.0 * dv.mean()),
                "per_chain": [int(x) for x in dv.reshape(n_chains, -1).sum(axis=1)],
            }
        if "num_steps" in extra_fields:
            ns = np.asarray(extra_fields["num_steps"])
            out["leapfrog"] = {
                "total": int(ns.sum()),
                "mean_per_iter": float(ns.mean()),
                "pct_at_max_tree_depth": float(100.0 * (ns >= 1023).mean()),
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
