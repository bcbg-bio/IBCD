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


def rho_log_penalty(rho, kind, sigma=0.5, start=1.0, width=0.2):
    """Log-density term that discourages a large spectral radius.

    'gaussian' is the Appendix H prior, log N(rho; 0, sigma^2). It acts at
    every rho, so it also shrinks cycles well inside rho < 1. 'barrier' is
    exactly zero below `start` and grows as ((rho - start) / width)^2 above
    it, so it leaves the interior alone. The hinge is squared so the gradient
    stays continuous, which NUTS needs.

    The defaults (start 1.0, width 0.2) are deliberately soft. At D = 150 a
    steeper barrier (0.9, 0.05) left chains that met it early in warmup with
    a collapsed step size, frozen wherever they were; the softer one lets a
    draw sit slightly past rho = 1 rather than pinning it there.
    """
    if kind == "gaussian":
        return dist.Normal(0.0, sigma).log_prob(rho)
    if kind == "barrier":
        return -jnp.square(jnp.maximum(rho - start, 0.0) / width)
    raise ValueError(f"unknown rho penalty '{kind}'")


def rhat_start(R_hat, se_hat, k=3.0):
    """Starting G from the total effects: R_hat soft-thresholded at k standard errors.

    sign(R_hat) * max(|R_hat| - k * S, 0), with a zero diagonal (R_hat_ii = 1 and
    S_ii = 0 exactly, so subtracting I and zeroing the diagonal agree). This is
    R - I, the total effects, not G: every node starts connected to all of its
    descendants, and the descent has to prune the indirect paths. What it does
    carry is orientation, from the asymmetry interventions create in R_hat (for
    a true edge i -> j, R_ij is nonzero and R_ji is not), and it uses only the
    model's own inputs, with no matrix inversion.

    k matters. A null entry clears one standard error about a third of the time,
    so at D = 500 k = 1 left ~77,000 spurious entries and ~14,000 pairs in both
    directions; k = 3 left ~600 and under ten, with 83% of the direct edges
    present and rho ~ 0.25. Higher k dropped too many true edges.

    Args:
        R_hat (np.ndarray): D x D total-effect estimates.
        se_hat (np.ndarray): D x D standard errors of R_hat.
        k (float): Threshold in standard errors.

    Returns:
        np.ndarray: D x D starting value for G.
    """
    R_hat = np.asarray(R_hat, dtype=float)
    G0 = np.sign(R_hat) * np.maximum(np.abs(R_hat) - k * np.asarray(se_hat, dtype=float), 0.0)
    np.fill_diagonal(G0, 0.0)
    return G0


def latents_from_G(G, pi0_ij, tau=0.1, slab_width=None):
    """Unconstrained values of the model's latents that reproduce a given G.

    The model does not sample G; it samples lam ~ HalfCauchy(1) (on the log
    scale), eps ~ N(0, 1) and spike ~ N(0, sigma0), and builds
    G = pi0 * spike + (1 - pi0) * tau * lam * eps. Only the product lam * eps is
    fixed by G, so each entry takes the most probable split under the prior,
    with spike at its mode of 0. "Most probable" is in the coordinates the
    sampler uses, log(lam), whose density carries the Jacobian term log(lam);
    in lam's own coordinates the HalfCauchy mode is at 0, which is degenerate.

    With m = (1 - pi0) tau and a = |G| / m, maximising
    -log(1 + lam^2) + log(lam) - eps^2 / 2 subject to lam * eps = a gives

        lam^2 = ((1 + a^2) + sqrt((1 + a^2)^2 + 4 a^2)) / 2,  eps = sign(G) a / lam,

    the only stationary point, a maximum. A zero entry gets lam = 1, eps = 0,
    the same values init_to_median gives; |eps| < 1 always. Entries the prior
    holds at the spike (pi0 = 1, so m = 0) get lam = 1, eps = 0 and stay at
    zero. Tied to the parameterisation of matrix_model_spike_horseshoe.

    With a slab width c the model caps each entry's scale s = (1 - pi0) tau lam
    at c via c s / sqrt(c^2 + s^2). lam is then taken from the same closed form
    and eps absorbs the cap, so G is still reproduced exactly; eps exceeds 1
    only for entries larger than about c (the edge scale), and the Adam
    descent rebalances the two.
    """
    G = np.asarray(G, dtype=float)
    m = (1.0 - np.asarray(pi0_ij, dtype=float)) * tau
    slab = m > 1e-12
    a = np.where(slab, np.abs(G) / np.where(slab, m, 1.0), 0.0)
    b = 1.0 + a * a
    lam = np.sqrt(0.5 * (b + np.sqrt(b * b + 4.0 * a * a)))
    if slab_width is None:
        eps = np.sign(G) * a / lam
    else:
        sc = m * lam
        sc = slab_width * sc / np.sqrt(slab_width ** 2 + sc ** 2)
        eps = np.where(slab, G / np.maximum(sc, 1e-300), 0.0)
    return {"lam": jnp.log(jnp.asarray(lam)), "eps": jnp.asarray(eps),
            "spike": jnp.zeros(G.shape)}


def optimized_init(model, model_kwargs, rng_key, num_chains, steps=2000, lr=1e-2,
                   jitter=0.1, start_G=None):
    """Starting points for NUTS from a short optimisation of the log posterior.

    From the default start, init_to_median, G is essentially empty and far from
    R_hat, so the likelihood gradient is steep. At D = 150 on SF graphs every
    chain started there left rho(G) < 1 within the first 20 warmup iterations
    and never came back, while chains started at the true G stayed at
    rho ~ 0.7. Descending the same potential from the same start in many small
    Adam steps reaches rho ~ 0.6 on every seed tested. The posterior is
    unchanged; only where the sampler begins is.

    At D = 500 the empty start is not enough: even a slow descent from it finds
    a wrong graph 3,000-5,500 nats below the optimum near the true G. Passing
    `start_G` (see rhat_start) begins the descent from the thresholded total
    effects instead, which lands within ~650-900 nats of it with valid rho.

    Uses numpyro.optim.Adam, a wrapper over jax.example_libraries.optimizers,
    so it adds no dependency. Adam moves each coordinate by about lr per step
    whatever the size of its gradient, which is what keeps the descent from
    overshooting the way an untuned leapfrog step does. At D >= 250 it does not
    fully settle at lr = 0.01 and the last iterate can sit a few thousand nats
    above the best one visited, but in the same basin (F1 within 0.005), which
    is negligible next to the ~3 D^2 / 2 nats NUTS climbs from the mode into
    the typical set, so the last iterate is returned.

    Without `start_G` each chain starts from its own init_to_median draw; with
    it, every chain starts from start_G. The chains end up close together,
    which would leave split R-hat little to compare, so `jitter` (in
    unconstrained units) is added to log(lam) and eps afterwards. The spike is left as optimised, since its
    prior scale is sigma0 = 1e-3. steps = 0 skips the optimisation and starts
    NUTS at the jittered start itself.

    Args:
        model: The numpyro model.
        model_kwargs (dict): Keyword arguments for the model.
        rng_key: PRNG key.
        num_chains (int): Number of starting points.
        steps (int): Adam steps; 0 for none.
        lr (float): Adam learning rate.
        jitter (float): Standard deviation of the noise added afterwards.
        start_G (np.ndarray): Optional D x D G to start the descent from.

    Returns:
        z (dict): Unconstrained starting values with a leading chain axis,
            ready for MCMC.run(init_params=...).
        info (dict): Potential energy before and after, and rho(G) at the
            start of each chain.
    """
    import inspect
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

    if start_G is None:
        z0 = mi.param_info.z
    else:
        tau = model_kwargs.get("tau", inspect.signature(model).parameters["tau"].default)
        one = latents_from_G(start_G, model_kwargs["pi0_ij"], tau,
                             slab_width=model_kwargs.get("slab_width"))
        # every chain descends from start_G itself; the jitter after the
        # descent separates them. Jittering before it as well put small
        # nonzero values on the null entries of the hub rows, which Adam at
        # lr 0.01 grew into cycles: from the true G at D = 500, rho went from
        # 0.38 to ~3 over 2000 steps, against 0.62 without it.
        z0 = {name: jnp.broadcast_to(one[name], (num_chains,) + one[name].shape)
              for name in mi.param_info.z}

    def optimise(z0):
        def step(state, _):
            return opt.update(grad(opt.get_params(state)), state), None
        state, _ = jax.lax.scan(step, opt.init(z0), None, length=steps)
        z = opt.get_params(state)
        return z, potential(z0), potential(z)

    z, pe_start, pe_end = jax.jit(jax.vmap(optimise))(z0)
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
                                 epsilon=1e-5, truncated_series=False, series_order=24,
                                 rho_penalty=None, rho_estimator="power", rho_sigma=0.5,
                                 rho_barrier_start=1.0, rho_barrier_width=0.2,
                                 slab_width=None):
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
        slab_width (float): None for the horseshoe as published. Otherwise a
            regularised horseshoe (Piironen & Vehtari 2017) applied to each
            entry's whole slab scale s = (1 - pi0) tau lam, replaced by
            c s / sqrt(c^2 + s^2): an entry can still reach scale c through a
            large lam, but its tail beyond c is Gaussian rather than Cauchy.
            The horseshoe's tail falls only as tau / t, so at D = 500 a prior
            draw has ~350 entries with |G| > 1.
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
    if slab_width is None:
        G = pi0_ij * spike_vals + (1. - pi0_ij) * slab_vals
    else:
        # Cap each entry's whole slab scale, (1 - pi0) tau lam, at the width.
        # Capping lam alone would cap the entry at (1 - pi0) c, which in rows
        # the SF prior shrinks hard (1 - pi0 down to ~1e-7) leaves true edges
        # needing |eps| in the thousands: at D = 500, 57-59% of true edges
        # would need |eps| > 3. Capping the product keeps the horseshoe's
        # escape hatch (a large lam, at ~2 log lam) up to the width.
        s = (1. - pi0_ij) * tau * lam
        s = slab_width * s / jnp.sqrt(slab_width ** 2 + s ** 2)
        G = pi0_ij * spike_vals + s * eps
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
