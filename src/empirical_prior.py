import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="cvxpy")

import pandas as pd
import numpy as np
import jax.numpy as jnp
from jax.scipy.stats import norm
import jax
import cvxpy as cp
from tqdm import tqdm


def load_R_and_SE_hat(r_hat_file, se_hat_file, causal_order: bool = False):
    """
    Load R_hat and SE_hat matrices, extract off-diagonal entries.

    If causal_order=False:
        use all off-diagonal entries (symmetric, undirected).
    If causal_order=True:
        use only upper-triangular off-diagonal entries (j < k).

    Returns:
        w  (jnp.ndarray): vector of selected R_hat entries.
        se (jnp.ndarray): vector of corresponding SE_hat entries.
    """
    R_hat = pd.read_csv(r_hat_file).values
    SE_hat = pd.read_csv(se_hat_file).values
    N = R_hat.shape[0]

    if causal_order:
        # upper triangle, excluding diagonal
        mask = np.triu(np.ones((N, N), dtype=bool), k=1)
    else:
        # all off-diagonal entries
        mask = ~np.eye(N, dtype=bool)

    w = R_hat[mask]
    se = SE_hat[mask]

    return w, se


def empirical_bayes_em(w, se_hat, K=50, sigma_min=0.01, sigma_max=1.0, alpha_er=2, tol=1e-6, safety_limit=5_000_000):
    """
    Run the EM algorithm for Empirical Bayes estimation of spike-and-slab mixture.

    Parameters:
        w (jnp.ndarray): Off-diagonal entries of R_hat (observations).
        se_hat (jnp.ndarray): Standard error estimates for each entry.
        K (int): Number of slab components.
        sigma_min (float): Minimum variance for slab components.
        sigma_max (float): Maximum variance for slab components.
        alpha (float): Penalty hyperparameter (default: 2).
        tol (float): Convergence threshold for stopping criteria.

    Returns:
        pi_0 (float): Estimated proportion of the spike component.
        pi_k (jnp.ndarray): Estimated proportions of the slab components.
    """
    N = len(w)
    sigma_k = jnp.linspace(jnp.sqrt(sigma_min), sigma_max, K)**2
    pi_0 = jnp.array(0.5)
    pi_k = jnp.full(K, (1 - pi_0) / K)

    @jax.jit
    def em_step(params):
        pi_0, pi_k = params
        f_0 = norm.pdf(w, loc=0, scale=jnp.sqrt(0.001) + se_hat)
        f_k = norm.pdf(w[:, None], loc=0, scale=jnp.sqrt(sigma_k)[None, :] + se_hat[:, None])
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

    pi0_var = cp.Variable((D,D), nonneg=True)
    pik_var = cp.Variable((D,D), nonneg=True)


    constraints = []

    offdiag_mask = np.ones((D, D), dtype=bool)
    np.fill_diagonal(offdiag_mask, False)

    constraints.append(pi0_var[offdiag_mask] + pik_var[offdiag_mask] == 1.0)
    constraints.append(cp.diag(pi0_var) == 1.0)
    constraints.append(cp.diag(pik_var) == 0.0)

    max_xi = xi.max()
    if max_xi < 1e-12:
        max_xi = 1e-12
    xi_norm = xi / max_xi


    data_resid = cp.sum_squares(
        cp.multiply(offdiag_mask, pik_var - xi_norm)
    ) + cp.sum_squares(
        cp.multiply(offdiag_mask, pi0_var - (1 - xi_norm))
    )

    # C) Global prior penalty: 
    #    We only want the sum over off-diagonal pi0_ij to be ~ pi0 * (D^2 - D).
    sum_pi0_offdiag = cp.sum(pi0_var) - cp.sum(cp.diag(pi0_var))
    #   i.e. pi0*(D^2 - D) is the target
    #global_penalty = cp.square(
    #    sum_pi0_offdiag - pi0*(D**2 - D)
    #)
    constraints.append(
        sum_pi0_offdiag == pi0*(D**2 - D)
    )
    # Weighted objective
    objective = alpha_data*data_resid # + beta_global*global_penalty

    # Solve
    prob = cp.Problem(cp.Minimize(objective), constraints)
    val = prob.solve()

    # Extract numeric solutions
    pi0_sol = pi0_var.value
    pik_sol = pik_var.value

    return pi0_sol, pik_sol, val

def scale_free_degree(R):
    """
    Given a matrix R (numpy array), solve for P using a degree-matching
    optimization, then return pi0 = 1 - P.

    Args:
        R (np.ndarray): Input square matrix of shape (D, D).

    Returns:
        pi0 (np.ndarray): Matrix of shape (D, D), with pi0 = 1 - P.
    """
    D = R.shape[0]

    # 1. Compute theta and phi
    A = np.abs(R)**2
    np.fill_diagonal(A, 0)
    theta_raw = A.sum(axis=1)
    phi_raw   = A.sum(axis=0)

    scale = (D-1) / max(theta_raw.max(), phi_raw.max())
    theta = theta_raw * scale
    phi   = phi_raw * scale

    # 2. Build constraint matrices A_out and A_in
    A_out = np.zeros((D, D*D))
    A_in = np.zeros((D, D*D))
    for i in range(D):
        for j in range(D):
            if i != j:
                idx = i * D + j
                A_out[i, idx] = 1
                A_in[j, idx] = 1

    # Remove diagonal entries
    keep = np.ones(D*D, dtype=bool)
    for i in range(D):
        keep[i * D + i] = False
    A_out = A_out[:, keep]
    A_in = A_in[:, keep]

    A_full = np.vstack([A_out, A_in])  # shape (2D, D^2 - D)
    b_full = np.concatenate([theta, phi])
    Nvar = A_full.shape[1]

    # 3. Solve quadratic program
    x = cp.Variable(Nvar)
    objective = cp.Minimize(cp.sum_squares(A_full @ x - b_full))
    constraints = [x >= 0, x <= 1]
    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.ECOS)

    P_flat = x.value
    if P_flat is None:
        raise RuntimeError("Solver failed: status %s" % prob.status)

    # 4. Reconstruct P with zero diagonal
    P = np.zeros((D, D))
    k = 0
    for i in range(D):
        for j in range(D):
            if i != j:
                P[i, j] = P_flat[k]
                k += 1

    # 5. Return pi0
    pi0 = 1 - P
    return pi0

def solve_edge_weights_rowwise(xi, pi0_i, alpha_sf=1.0, solver=cp.ECOS, causal_order=False, symmetric=False):
    D = xi.shape[0]
    pi0_ij  = np.zeros((D, D))
    pi_k_ij = np.zeros((D, D))

    for i in range(D):
        if causal_order:
            # Causal order: only j>i
            idx = np.arange(i+1, D)
        else:
            # non-causal order: all off-diagonals in row i
            idx = np.r_[np.arange(0, i), np.arange(i+1, D)]
        n = len(idx)
        if n == 0:
            continue

        row = xi[i, idx]
        m = row.max() if row.max() > 0 else 1.0
        xnorm = row / m

        p0 = cp.Variable(n, nonneg=True)
        pk = cp.Variable(n, nonneg=True)
        cons = [p0 + pk == 1, cp.sum(p0) == pi0_i[i] * n]
        obj = cp.sum_squares(pk - xnorm) + cp.sum_squares(p0 - (1 - xnorm))

        prob = cp.Problem(cp.Minimize(alpha_sf * obj), cons)
        prob.solve(solver=solver, verbose=False)

        pi0_ij[i, idx] = p0.value
        pi_k_ij[i, idx] = pk.value

        if not causal_order and symmetric:
            # optional: mirror to make undirected
            pi0_ij[idx, i] = p0.value
            pi_k_ij[idx, i] = pk.value

    # diagonal always spike
    np.fill_diagonal(pi0_ij, 1.0)
    np.fill_diagonal(pi_k_ij, 0.0)

    if causal_order:
        # only when enforcing DAG
        tril_i, tril_j = np.tril_indices(D)
        pi0_ij[tril_i, tril_j] = 1.0
        pi_k_ij[tril_i, tril_j] = 0.0

    return pi0_ij, pi_k_ij