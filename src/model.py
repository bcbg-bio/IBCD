import numpyro
import numpyro.distributions as dist
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
