import numpyro
import numpyro.distributions as dist
from numpyro.diagnostics import split_gelman_rubin, effective_sample_size
import jax.numpy as jnp
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


def optimized_init(model, model_kwargs, rng_key, num_chains, steps=2000, lr=1e-2,
                   jitter=0.1):
    """Starting points for NUTS from a short optimisation of the log posterior.

    Each chain starts from its own init_to_median draw, where G is essentially
    empty and far from R_hat, so the likelihood gradient is steep. Descending
    the same potential in many small Adam steps first moves the chains into a
    plausible region before NUTS adapts its step size. The posterior is
    unchanged; only where the sampler begins is.

    Uses numpyro.optim.Adam, a wrapper over jax.example_libraries.optimizers,
    so it adds no dependency. Adam moves each coordinate by about lr per step
    whatever the size of its gradient, which is what keeps the descent from
    overshooting the way an untuned leapfrog step does. The last iterate is
    returned rather than the best one visited; see NOTES.md, 2026-10-02.

    The chains can end up close together, which would leave split R-hat little
    to compare, so `jitter` (in unconstrained units) is added to log(lam) and
    eps afterwards. The spike is left as optimised, since its prior scale is
    sigma0 = 1e-3. steps = 0 skips the optimisation and starts NUTS at the
    jittered init_to_median draws.

    Args:
        model: The numpyro model.
        model_kwargs (dict): Keyword arguments for the model.
        rng_key: PRNG key.
        num_chains (int): Number of starting points.
        steps (int): Adam steps; 0 for none.
        lr (float): Adam learning rate.
        jitter (float): Standard deviation of the noise added afterwards.

    Returns:
        z (dict): Unconstrained starting values with a leading chain axis,
            ready for MCMC.run(init_params=...).
        info (dict): Potential energy before and after, and rho(G) at the
            start of each chain.
    """
    import jax
    from numpyro import infer
    from numpyro.infer.util import initialize_model
    from numpyro.optim import Adam

    k_init, k_jit = jax.random.split(rng_key)
    mi = initialize_model(jax.random.split(k_init, num_chains), model,
                          model_kwargs=model_kwargs,
                          init_strategy=infer.init_to_median(num_samples=50))
    potential = mi.potential_fn
    grad = jax.grad(potential)
    opt = Adam(lr)

    def add_jitter(z, key):
        keys = jax.random.split(key, 2)
        for kk, name in zip(keys, ("lam", "eps")):
            if name in z:
                z[name] = z[name] + jitter * jax.random.normal(kk, z[name].shape)
        return z

    def optimise(z0):
        def step(state, _):
            return opt.update(grad(opt.get_params(state)), state), None
        state, _ = jax.lax.scan(step, opt.init(z0), None, length=steps)
        z = opt.get_params(state)
        return z, potential(z0), potential(z)

    z, pe_start, pe_end = jax.jit(jax.vmap(optimise))(mi.param_info.z)
    z = add_jitter(dict(z), k_jit)

    G = np.asarray(jax.vmap(mi.postprocess_fn)(z)["G"], dtype=np.float64)
    rho = np.abs(np.linalg.eigvals(G)).max(axis=1)
    info = {
        "potential_before": [float(x) for x in np.asarray(pe_start)],
        "potential_after": [float(x) for x in np.asarray(pe_end)],
        "rho_start_per_chain": [float(x) for x in rho],
    }
    return z, info


def matrix_model_spike_horseshoe(obs_data, pi0_ij, U_lower, V_lower, D, sigma0=0.001, tau=0.1,
                                 epsilon=1e-5, rho_penalty=None, rho_sigma=0.5,
                                 slab_width=None, mvn_factor=None, mvn_diag=None):
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
        epsilon (float): Regularization added to the diagonal of I - G before
            it is inverted.
        rho_penalty (str): None for no constraint on rho(G), or 'gaussian' for
            the Appendix H prior N(0, rho_sigma^2) on rho, estimated by
            `power_iteration_radius`.
        rho_sigma (float): Scale of the Gaussian penalty; 0.5 in Appendix H.
        slab_width (float): None for the horseshoe as published. Otherwise a
            regularised horseshoe (Piironen & Vehtari 2017) applied to each
            entry's whole slab scale s = (1 - pi0) tau lam, replaced by
            c s / sqrt(c^2 + s^2): an entry can still reach scale c through a
            large lam, but its tail beyond c is Gaussian rather than Cauchy,
            where the horseshoe's tail falls only as tau / t.
        mvn_factor (array): None for the matrix-normal likelihood. Otherwise
            the Table 7 ablation: vec(R_hat) ~ N(vec(R), F F^T + diag(mvn_diag))
            with F the D^2 x r factor of a truncated eigendecomposition of the
            full covariance S (equation 9); U_lower and V_lower are unused.
        mvn_diag (array): Diagonal added to F F^T, length D^2.
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
    # zeroed: the model has no self-loops and R_ii is fixed at 1 by definition.
    if slab_width is None:
        G = pi0_ij * spike_vals + (1. - pi0_ij) * slab_vals
    else:
        # Cap each entry's whole slab scale, (1 - pi0) tau lam, at the width.
        # Capping lam alone would cap the entry at (1 - pi0) c, which in rows
        # the SF prior shrinks hard leaves true edges needing a huge |eps|;
        # see NOTES.md, 2026-10-02. Capping the product keeps the horseshoe's
        # escape hatch (a large lam, at ~2 log lam) up to the width.
        s = (1. - pi0_ij) * tau * lam
        s = slab_width * s / jnp.sqrt(slab_width ** 2 + s ** 2)
        G = pi0_ij * spike_vals + s * eps
    G = G * (1. - jnp.eye(D))
    numpyro.deterministic("G", G)           # keep only G

    if rho_penalty:
        if rho_penalty != "gaussian":
            raise ValueError(f"unknown rho penalty '{rho_penalty}'")
        rho_est = power_iteration_radius(G)
        # recorded so each run can compare the estimate the penalty acted on
        # with the exact spectral radius computed afterwards
        numpyro.deterministic("rho_estimate", rho_est)
        numpyro.factor("spectral_radius_penalty",
                       dist.Normal(0.0, rho_sigma).log_prob(rho_est))

    # MatrixNormal mean R = (I - G)^-1, the sum over directed paths for a DAG.
    I = jnp.eye(D)
    I_minus_G = I - G + epsilon * I
    R_mean = jnp.linalg.solve(I_minus_G, I)

    if mvn_factor is not None:
        numpyro.sample(
            "R_hat_obs",
            dist.LowRankMultivariateNormal(loc=R_mean.reshape(-1),
                                           cov_factor=mvn_factor, cov_diag=mvn_diag),
            obs=obs_data[:D, :D].reshape(-1),
        )
        return

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

def draw_spectral_radius(G_draws):
    """Spectral radius rho(G) of each posterior draw.

    R = sum_d G^d converges only when rho(G) < 1; the model computes R as
    (I - G)^-1, which is defined for any G without an eigenvalue at 1, so rho
    is reported rather than used to reject draws.

    Parameters:
        G_draws (array): Posterior draws of G, shape (..., D, D).

    Returns:
        np.ndarray: rho of each draw, flattened over the leading axes.
    """
    g = np.asarray(G_draws)
    flat = g.reshape(-1, g.shape[-2], g.shape[-1])
    return np.abs(np.linalg.eigvals(flat)).max(axis=1)


def report_spectral_radius(rho, n_chains, pathological_rho=10.0):
    """Print a note on draws with rho >= 1, and warn when any is pathological.

    Every draw is kept; see NOTES.md, 2026-10-02. A note is printed when draws
    have rho >= 1 (R is then not a convergent path sum for them), and a
    RuntimeWarning raised when any exceeds `pathological_rho`, which flags a
    posterior that has run away to implausibly large cycles.

    Args:
        rho (np.ndarray): Spectral radius of each flattened draw.
        n_chains (int): Number of chains the draws come from, chain-major.
        pathological_rho (float): rho above which a warning is raised.
    """
    import warnings
    rho = np.asarray(rho)
    n = rho.size
    per_chain = n // n_chains
    ge1 = rho >= 1.0
    if ge1.any():
        counts = ge1.reshape(n_chains, -1).sum(axis=1)
        by_chain = ", ".join(f"chain {c}: {int(k)}/{per_chain}" for c, k in enumerate(counts))
        print(f"{int(ge1.sum())} of {n} draws ({100 * ge1.mean():.1f}%) have spectral "
              f"radius >= 1 ({by_chain}); median rho {np.median(rho):.3g}, max "
              f"{rho.max():.3g}. They are kept; for them R = (I - G)^-1 is not a "
              "convergent path sum.")
    if rho.max() > pathological_rho:
        warnings.warn(f"max spectral radius {rho.max():.3g} exceeds {pathological_rho}; "
                      "the posterior has likely run away from plausible graphs "
                      "(see spectral_radius in diagnostics.json).", RuntimeWarning)


def posterior_diagnostics(G_draws, rho, extra_fields=None,
                          max_ess_entries=5000, seed=0, max_tree_depth=12,
                          max_rhat_entries=20000):
    """Summarise sampler behaviour and the spectral radius of the draws for one run.

    Convergence statistics are computed over the entries of G. ESS is
    estimated on a random subsample of entries, since an autocorrelation
    estimate per entry is O(D^2) and D can be 500.

    Parameters:
        G_draws (array): Posterior draws, shape (chains, draws, D, D).
        rho (array): Spectral radius of each flattened draw.
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
    rho = np.asarray(rho)

    out = {
        "n_chains": int(n_chains),
        "n_draws_per_chain": int(n_draws),
        "n_draws_total": int(rho.size),
        "D": int(D),
    }

    # draws outside the region where R = sum_d G^d converges
    nonconv = rho >= 1.0
    out["nonconvergent"] = {
        "n": int(nonconv.sum()), "pct": float(100.0 * nonconv.mean()),
        "per_chain": [int(x) for x in nonconv.reshape(n_chains, n_draws).sum(axis=1)],
    }
    out["spectral_radius"] = {
        k: float(v) for k, v in zip(
            ["min", "median", "p95", "p99", "max"],
            np.percentile(rho, [0, 50, 95, 99, 100]),
        )
    }

    # convergence over the entries of G
    offdiag = ~np.eye(D, dtype=bool)
    if n_chains >= 2:
        # Index the entries before the chains: g[:, :, offdiag] would
        # materialise the whole off-diagonal block, which is 3 GB at D = 500 on
        # top of the draws themselves. Subsampling first keeps it to a few
        # hundred MB, and r_hat over a random subset of entries answers the
        # same question as r_hat over all of them.
        rng = np.random.default_rng(seed)
        off_idx = np.flatnonzero(offdiag)
        n_rhat = min(max_rhat_entries, off_idx.size)
        r_idx = np.sort(rng.choice(off_idx, size=n_rhat, replace=False))
        sub = g.reshape(n_chains, n_draws, D * D)[:, :, r_idx]
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
        out["r_hat"] = {"note": "split R-hat and ESS need at least two chains"}
        out["ess"] = {"note": "split R-hat and ESS need at least two chains"}

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
    Compute Local False Sign Rate for each edge from posterior samples.

    Parameters:
        flat_samples (array): Posterior samples of shape (samples, D, D)

    Returns:
        array: LFSR matrix of shape (D, D)
    """
    pos_probs = np.mean(flat_samples > 0, axis=0)
    neg_probs = np.mean(flat_samples < 0, axis=0)
    return np.minimum(pos_probs, neg_probs)
