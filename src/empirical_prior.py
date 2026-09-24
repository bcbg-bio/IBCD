import pandas as pd
import numpy as np
import jax.numpy as jnp
from jax.scipy.stats import norm
import jax
from tqdm import tqdm


def load_R_and_SE_hat(r_hat_file, se_hat_file):
    """
    Load R_hat and SE_hat matrices, extract all off-diagonal entries.

    Returns:
        w  (jnp.ndarray): vector of selected R_hat entries.
        se (jnp.ndarray): vector of corresponding SE_hat entries.
    """
    R_hat = pd.read_csv(r_hat_file).values
    SE_hat = pd.read_csv(se_hat_file).values
    N = R_hat.shape[0]

    mask = ~np.eye(N, dtype=bool)

    w = R_hat[mask]
    se = SE_hat[mask]

    return w, se


def _project_to_budget(x, budget, max_iter=200, tol=1e-14):
    """Euclidean projection of x onto {y : sum(y) == budget, 0 <= y <= 1}.

    Both prior paths minimise a squared distance to xi subject to the box and
    a single constraint fixing the total slab mass. That is exactly this
    projection, whose solution is clip(x + c, 0, 1) for the unique offset c
    matching the budget. c is found by bisection; the bracket below is chosen
    so the clipped sum runs from 0 to len(x) across it.

    Args:
        x (np.ndarray): Point to project.
        budget (float): Required sum of the result.

    Returns:
        np.ndarray: The projection of x.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n == 0:
        return x.copy()

    budget = float(np.clip(budget, 0.0, n))
    lo = -float(x.max())       # sum(clip(x + lo, 0, 1)) == 0
    hi = 1.0 - float(x.min())  # sum(clip(x + hi, 0, 1)) == n
    for _ in range(max_iter):
        c = 0.5 * (lo + hi)
        if np.clip(x + c, 0.0, 1.0).sum() < budget:
            lo = c
        else:
            hi = c
        if hi - lo < tol:
            break
    return np.clip(x + 0.5 * (lo + hi), 0.0, 1.0)


def empirical_bayes_em(w, se_hat, K=50, sigma0=0.001, sigma_min=0.01, sigma_max=1.0, alpha_er=2, tol=1e-6, safety_limit=5_000_000):
    """
    Run the EM algorithm for Empirical Bayes estimation of spike-and-slab mixture.

    Parameters:
        w (jnp.ndarray): Off-diagonal entries of R_hat (observations).
        se_hat (jnp.ndarray): Standard error estimates for each entry.
        K (int): Number of slab components.
        sigma0 (float): Standard deviation of the spike component.
        sigma_min (float): Minimum variance for slab components.
        sigma_max (float): Maximum variance for slab components.
            sigma_k holds the corresponding standard deviations.
        alpha_er (float): Penalty hyperparameter (default: 2).
        tol (float): Convergence threshold for stopping criteria.
        safety_limit (int): Maximum number of EM iterations before stopping.

    Returns:
        pi_0 (float): Estimated proportion of the spike component.
        pi_k (jnp.ndarray): Estimated proportions of the slab components.
    """
    N = len(w)
    sigma_k = jnp.linspace(jnp.sqrt(sigma_min), jnp.sqrt(sigma_max), K)
    pi_0 = jnp.array(0.5)
    pi_k = jnp.full(K, (1 - pi_0) / K)

    @jax.jit
    def em_step(params):
        pi_0, pi_k = params
        f_0 = norm.pdf(w, loc=0, scale=jnp.sqrt(sigma0**2 + se_hat**2))
        f_k = norm.pdf(w[:, None], loc=0, scale=jnp.sqrt(sigma_k[None, :]**2 + se_hat[:, None]**2))
        numerator_0 = pi_0 * f_0
        numerator_k = pi_k * f_k
        denominator = numerator_0 + numerator_k.sum(axis=1)
        z_0 = numerator_0 / denominator
        z_k = numerator_k / denominator[:, None]
        n_0 = jnp.sum(z_0)
        n_k = jnp.sum(z_k, axis=0)
        pi_0_new = n_0 / (N + (alpha_er - 1) * K)
        pi_k_new = (n_k + (alpha_er - 1)) / (N + (alpha_er - 1) * K)
        return pi_0_new, pi_k_new
    
    pbar = tqdm(
        total=0,
        position=0,
        leave=True,
        dynamic_ncols=True,
        desc="EM",
        unit=" iterations"
    )
    
    iteration = 0
    while True:
        pi_0_new, pi_k_new = em_step((pi_0, pi_k))
        if float(jnp.max(jnp.abs(pi_k_new - pi_k))) < tol:
            pi_0, pi_k = pi_0_new, pi_k_new
            break
        pi_0, pi_k = pi_0_new, pi_k_new

        iteration += 1
        pbar.update(1) 

        if iteration >= safety_limit:
            print("EM did not reach tolerance. Stopped at safety limit.")
            break

    pbar.close()

    return pi_0, pi_k, sigma_k

def solve_spike_slab_diagonal_spike(
    xi,
    pi0,
    alpha_data=1.0,
    beta_global=1.0
):
    """
    Solve for edge-specific spike and slab weights under constraints.

    Args:
        xi (np.ndarray): Empirical Bayes statistic matrix (e.g., lfsr or xi scores).
        pi0 (float): Global spike proportion.
        alpha_data (float): Weight on data residual term.

    Returns:
        pi0_ij (np.ndarray): Edge-specific spike weight matrix.
        pi_k_ij (np.ndarray): Edge-specific slab weight matrix.
        val (float): Optimization objective value.
    """

    D = xi.shape[0]

    offdiag_mask = np.ones((D, D), dtype=bool)
    np.fill_diagonal(offdiag_mask, False)

    # Scale xi so that its total over the off-diagonal equals the slab budget
    # implied by the sparsity constraint, S = (1 - pi0) * (D^2 - D). The data
    # term then agrees with the constraint, so the projection below returns
    # xi_norm itself wherever the [0, 1] box does not bind.
    slab_budget = float(1.0 - pi0) * (D**2 - D)
    xi_offdiag_sum = float(xi[offdiag_mask].sum())
    if xi_offdiag_sum > 1e-12:
        xi_norm = xi * (slab_budget / xi_offdiag_sum)
    else:
        xi_norm = np.full((D, D), slab_budget / (D**2 - D))

    # pi0_ij + pi_k_ij == 1 off the diagonal, and the diagonal is all spike,
    # so only the off-diagonal slab weights are free.
    pik_sol = np.zeros((D, D))
    pi0_sol = np.ones((D, D))
    pik_sol[offdiag_mask] = _project_to_budget(xi_norm[offdiag_mask], slab_budget)
    pi0_sol[offdiag_mask] = 1.0 - pik_sol[offdiag_mask]

    # Objective value: with pi0 = 1 - pi_k the two residual terms are equal.
    resid = pik_sol[offdiag_mask] - xi_norm[offdiag_mask]
    val = float(alpha_data * 2.0 * np.sum(resid ** 2))

    return pi0_sol, pik_sol, val

def scale_free_degree(R):
    """Per-node spike proportions from the in- and out-strengths of R_hat.

    Let A = |R_hat|^2 with a zero diagonal, theta_i = sum_j A_ij the
    out-strength of node i and phi_j = sum_i A_ij the in-strength of node j.
    The scale-free prior matches an edge-probability matrix P in [0, 1] to
    those marginals after rescaling them so the largest becomes D - 1, the
    most edges a node can have. Only the row means of P are used downstream,
    and any P attaining the marginals has row sums theta_i * (D-1) / m with
    m = max(max theta, max phi), so

        pi0_i = 1 - (sum_j P_ij) / (D - 1) = 1 - theta_i / m

    which needs no solver. theta_i <= m by construction, so pi0_i is in
    [0, 1] without clipping.

    Args:
        R (np.ndarray): D x D matrix of estimated total causal effects.

    Returns:
        np.ndarray: Length-D vector of per-node spike proportions.
    """
    D = R.shape[0]

    A = np.abs(R)**2
    np.fill_diagonal(A, 0)
    theta = A.sum(axis=1)
    phi = A.sum(axis=0)

    m = max(theta.max(), phi.max())
    if m <= 0:
        return np.ones(D)
    return 1.0 - theta / m


def solve_edge_weights_rowwise(xi, pi0_i):
    """Solve for edge-specific spike and slab weights, one row at a time.

    Args:
        xi (np.ndarray): Symmetric D x D matrix of interaction strengths.
        pi0_i (np.ndarray): Length-D vector of per-node spike proportions.

    Returns:
        pi0_ij (np.ndarray): Edge-specific spike weight matrix.
        pi_k_ij (np.ndarray): Edge-specific slab weight matrix.
    """
    D = xi.shape[0]
    pi0_ij  = np.zeros((D, D))
    pi_k_ij = np.zeros((D, D))

    for i in range(D):
        # all off-diagonals in row i
        idx = np.r_[np.arange(0, i), np.arange(i+1, D)]
        n = len(idx)
        if n == 0:
            continue

        # Scale this row of xi so its total equals the row's slab budget,
        # d_i = (1 - pi0_i) * n, which is what the sparsity constraint fixes.
        row = xi[i, idx]
        row_budget = float(1.0 - pi0_i[i]) * n
        row_sum = float(row.sum())
        if row_sum > 1e-12:
            xnorm = row * (row_budget / row_sum)
        else:
            xnorm = np.full(n, row_budget / n)

        pk = _project_to_budget(xnorm, row_budget)
        p0 = 1.0 - pk

        pi0_ij[i, idx] = p0
        pi_k_ij[i, idx] = pk

    # diagonal always spike
    np.fill_diagonal(pi0_ij, 1.0)
    np.fill_diagonal(pi_k_ij, 0.0)

    return pi0_ij, pi_k_ij
