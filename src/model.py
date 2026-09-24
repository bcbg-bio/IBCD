import numpyro
import numpyro.distributions as dist
import jax.numpy as jnp
import arviz as az
import numpy as np

def matrix_model_spike_horseshoe(obs_data, pi0_ij, U_lower, V_lower, D, sigma0=0.001, tau=0.1, epsilon=1e-5):
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
        epsilon (float): Regularization for stability.
    """

    # Sample horseshoe local scales (HalfCauchy), shape (D, D)
    lam = numpyro.sample("lam", dist.HalfCauchy(1.).expand((D, D)), infer={"is_auxiliary": True})

    # Sample standard normal noise, for non-centered parameterization
    eps = numpyro.sample("eps", dist.Normal(0., 1.).expand((D, D)), infer={"is_auxiliary": True})

    # Slab values: non-centered horseshoe (continuous)
    slab_vals = tau * lam * eps

    # Spike values: zero-mean Normal with tiny variance
    spike_vals = numpyro.sample("spike", dist.Normal(0., sigma0).expand((D, D)), infer={"is_auxiliary": True})

    # Binary mixture via convex combination (no discrete Cat)
    G = pi0_ij * spike_vals + (1. - pi0_ij) * slab_vals
    numpyro.deterministic("G", G)           # keep only G

    # MatrixNormal mean
    I = jnp.eye(D)
    I_minus_G = I - G + epsilon * jnp.eye(D)
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

def spectral_radius(G, max_iter=50):
    v = jnp.ones(G.shape[0])
    for _ in range(max_iter):
        v = G @ v
        v = v / (jnp.linalg.norm(v) + 1e-8)
    return jnp.linalg.norm(G @ v)

    
def matrix_model_spike_horseshoe_spectral(obs_data, pi0_ij, U_lower, V_lower, D, sigma0=0.001, tau=0.1, epsilon=1e-5):
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
        epsilon (float): Regularization for stability.
    """

    # Sample horseshoe local scales (HalfCauchy), shape (D, D)
    lam = numpyro.sample("lam", dist.HalfCauchy(1.).expand((D, D)), infer={"is_auxiliary": True})

    # Sample standard normal noise, for non-centered parameterization
    eps = numpyro.sample("eps", dist.Normal(0., 1.).expand((D, D)), infer={"is_auxiliary": True})

    # Slab values: non-centered horseshoe (continuous)
    slab_vals = tau * lam * eps

    # Spike values: zero-mean Normal with tiny variance
    spike_vals = numpyro.sample("spike", dist.Normal(0., sigma0).expand((D, D)), infer={"is_auxiliary": True})

    # Binary mixture via convex combination (no discrete Cat)
    G = pi0_ij * spike_vals + (1. - pi0_ij) * slab_vals
    numpyro.deterministic("G", G)           # keep only G
    
    rho_G = spectral_radius(G)
    sigma_R = 0.5
    numpyro.factor(
        "spectral_radius_penalty",
        dist.Normal(0.0, sigma_R).log_prob(rho_G)
    )

    # beta = 0.1  
    # numpyro.factor(
    #     "spectral_radius_penalty",
    #     - jnp.log1p((rhofj_G / beta)**2)
    # )

    
    # MatrixNormal mean
    I = jnp.eye(D)
    I_minus_G = I - G + epsilon * jnp.eye(D)
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
