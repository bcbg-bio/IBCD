"""Tests for IBCD's empirical prior construction and end-to-end pipeline.

Covers the invariants each prior path is supposed to satisfy (the constraints
in eq 18 and eq 19) and the well-formedness of the pipeline's outputs.

Run everything, including the slow end-to-end test:

    python tests/test_ibcd.py --slow

Run only the fast tests:

    python tests/test_ibcd.py

The functions are plain `test_*` functions, so `pytest tests/test_ibcd.py`
works too if pytest is available.
"""

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src"
sys.path.insert(0, str(SRC))

from model import (  # noqa: E402
    convergent_draws,
    gelfand_radius,
    latents_from_G,
    matrix_model_spike_horseshoe,
    optimized_init,
    posterior_diagnostics,
    power_iteration_radius,
    screen_draws,
    rhat_start,
    rho_log_penalty,
)
from empirical_prior import (  # noqa: E402
    _cap_and_redistribute,
    _project_to_budget,
    em_slab_scale,
    load_R_and_SE_hat,
    scale_free_degree,
    solve_edge_weights_rowwise,
    solve_spike_slab_diagonal_spike,
)


def _fixture(D=10, seed=0):
    """A small symmetric interaction matrix xi and a matching R_hat."""
    rng = np.random.default_rng(seed)
    A = rng.gamma(1.0, 1.0, size=(D, D))
    xi = (A + A.T) / 2.0  # xi is symmetric by construction in the pipeline
    np.fill_diagonal(xi, 0.0)
    xi /= np.sqrt((xi ** 2).sum())

    R = rng.normal(0, 0.3, size=(D, D))
    np.fill_diagonal(R, 1.0)
    return xi, R


# --------------------------------------------------------------------------
# load_R_and_SE_hat
# --------------------------------------------------------------------------

def test_load_R_and_SE_hat_returns_offdiagonal_entries():
    D = 7
    rng = np.random.default_rng(1)
    R = rng.normal(size=(D, D))
    SE = np.abs(rng.normal(size=(D, D))) + 0.1

    with tempfile.TemporaryDirectory() as tmp:
        r_path, se_path = Path(tmp) / "R.csv", Path(tmp) / "SE.csv"
        pd.DataFrame(R).to_csv(r_path, index=False)
        pd.DataFrame(SE).to_csv(se_path, index=False)
        w, se = load_R_and_SE_hat(str(r_path), str(se_path))

    offdiag = ~np.eye(D, dtype=bool)
    assert w.shape[0] == D * (D - 1)
    assert np.allclose(np.sort(w), np.sort(R[offdiag]))
    assert np.allclose(np.sort(se), np.sort(SE[offdiag]))


# --------------------------------------------------------------------------
# SF path: solve_edge_weights_rowwise
# --------------------------------------------------------------------------

def test_solve_edge_weights_rowwise_satisfies_its_constraints():
    D = 10
    xi, _ = _fixture(D)
    pi0_i = np.full(D, 0.8)
    pi0_ij, pi_k_ij = solve_edge_weights_rowwise(xi, pi0_i)

    offdiag = ~np.eye(D, dtype=bool)
    assert np.allclose(pi0_ij[offdiag] + pi_k_ij[offdiag], 1.0, atol=1e-6)
    assert np.allclose(np.diag(pi0_ij), 1.0)
    assert np.allclose(np.diag(pi_k_ij), 0.0)
    assert pi0_ij.min() >= -1e-6 and pi0_ij.max() <= 1 + 1e-6
    assert pi_k_ij.min() >= -1e-6 and pi_k_ij.max() <= 1 + 1e-6

    # eq 19: the per-row spike mass equals pi0_i * (number of candidates in the row)
    for i in range(D):
        row = pi0_ij[i, np.arange(D) != i]
        assert np.isclose(row.sum(), pi0_i[i] * (D - 1), atol=1e-5), (
            f"row {i}: spike mass {row.sum():.6f} != {pi0_i[i] * (D - 1):.6f}"
        )


def test_solve_edge_weights_rowwise_leaves_no_hard_spike():
    """No candidate edge may get pi0 = 1.

    A hard spike leaves G = spike ~ N(0, sigma0^2) with sigma0 = 1e-3, so the
    edge cannot be recovered however strong the evidence. The budget-matched
    normalisation of xi keeps every off-diagonal entry strictly interior.
    """
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0_i = np.full(D, 0.8)
    pi0_ij, _ = solve_edge_weights_rowwise(xi, pi0_i)

    offdiag = ~np.eye(D, dtype=bool)
    assert pi0_ij[offdiag].max() < 1.0 - 1e-6, (
        f"{(pi0_ij[offdiag] >= 1.0 - 1e-6).sum()} of {offdiag.sum()} candidates "
        "are hard-clamped to the spike"
    )


def test_solve_edge_weights_rowwise_matches_budget_scaled_xi():
    """With xi scaled to the row budget the QP solution is xi itself.

    The data term and the sparsity constraint then agree, so the optimal
    offset is zero and no mass is redistributed.
    """
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0_i = np.full(D, 0.8)
    _, pi_k_ij = solve_edge_weights_rowwise(xi, pi0_i)

    for i in range(D):
        idx = np.arange(D) != i
        row = xi[i, idx]
        target = row * ((1.0 - pi0_i[i]) * (D - 1) / row.sum())
        assert target.max() < 1.0, "fixture must not make the box bind"
        assert np.allclose(pi_k_ij[i, idx], target, atol=1e-5), (
            f"row {i}: solution departs from budget-scaled xi"
        )


def test_solve_edge_weights_rowwise_tracks_xi():
    """The QP fits pi_k to xi, so within a row the two should be co-monotone."""
    D = 12
    xi, _ = _fixture(D, seed=3)
    pi0_i = np.full(D, 0.7)
    _, pi_k_ij = solve_edge_weights_rowwise(xi, pi0_i)

    for i in range(D):
        j = np.arange(D) != i
        order_xi = np.argsort(xi[i, j])
        pk = pi_k_ij[i, j][order_xi]
        assert np.all(np.diff(pk) >= -1e-6), f"row {i}: pi_k is not increasing in xi"


# --------------------------------------------------------------------------
# ER path: solve_spike_slab_diagonal_spike
# --------------------------------------------------------------------------

def test_solve_spike_slab_diagonal_spike_satisfies_its_constraints():
    D = 10
    xi, _ = _fixture(D)
    pi0 = 0.8
    pi0_ij, pi_k_ij, _ = solve_spike_slab_diagonal_spike(xi, pi0=pi0)

    offdiag = ~np.eye(D, dtype=bool)
    assert np.allclose(pi0_ij[offdiag] + pi_k_ij[offdiag], 1.0, atol=1e-6)
    assert np.allclose(np.diag(pi0_ij), 1.0, atol=1e-6)
    assert np.allclose(np.diag(pi_k_ij), 0.0, atol=1e-6)
    # eq 18: total off-diagonal spike mass
    assert np.isclose(pi0_ij[offdiag].sum(), pi0 * (D ** 2 - D), atol=1e-4)


def test_solve_spike_slab_diagonal_spike_leaves_no_hard_spike():
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0_ij, _, _ = solve_spike_slab_diagonal_spike(xi, pi0=0.8)

    offdiag = ~np.eye(D, dtype=bool)
    assert pi0_ij[offdiag].max() < 1.0 - 1e-6


def test_solve_spike_slab_diagonal_spike_matches_budget_scaled_xi():
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0 = 0.8
    _, pi_k_ij, _ = solve_spike_slab_diagonal_spike(xi, pi0=pi0)

    offdiag = ~np.eye(D, dtype=bool)
    target = xi * ((1.0 - pi0) * (D ** 2 - D) / xi[offdiag].sum())
    assert target[offdiag].max() < 1.0, "fixture must not make the box bind"
    assert np.allclose(pi_k_ij[offdiag], target[offdiag], atol=1e-5)


def test_er_prior_is_symmetric():
    """xi is symmetric, so the ER prior carries no directional information.

    Documented behaviour, not a defect - it is why all orientation signal in
    the ER path comes from the likelihood.
    """
    D = 10
    xi, _ = _fixture(D)
    pi0_ij, _, _ = solve_spike_slab_diagonal_spike(xi, pi0=0.8)
    assert np.allclose(pi0_ij, pi0_ij.T, atol=1e-5)


# --------------------------------------------------------------------------
# Closed-form projection vs the quadratic program it replaced
# --------------------------------------------------------------------------

def test_project_to_budget_matches_known_solutions():
    # box slack: the offset is zero when x already sums to the budget
    x = np.array([0.1, 0.2, 0.3])
    assert np.allclose(_project_to_budget(x, x.sum()), x)
    # box binds above: [0, 0, 5] onto sum == 2 gives [0.5, 0.5, 1]
    assert np.allclose(_project_to_budget(np.array([0.0, 0.0, 5.0]), 2.0),
                       [0.5, 0.5, 1.0])
    # degenerate budgets
    assert np.allclose(_project_to_budget(np.array([0.3, 0.7]), 0.0), [0.0, 0.0])
    assert np.allclose(_project_to_budget(np.array([0.3, 0.7]), 2.0), [1.0, 1.0])


def test_rowwise_solver_matches_the_cvxpy_program():
    import cvxpy as cp
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0_i = np.full(D, 0.8)
    _, pi_k_ij = solve_edge_weights_rowwise(xi, pi0_i)

    for i in range(D):
        idx = np.arange(D) != i
        n = D - 1
        row = xi[i, idx]
        xnorm = row * ((1.0 - pi0_i[i]) * n / row.sum())
        p0 = cp.Variable(n, nonneg=True)
        pk = cp.Variable(n, nonneg=True)
        cons = [p0 + pk == 1, cp.sum(p0) == pi0_i[i] * n]
        obj = cp.sum_squares(pk - xnorm) + cp.sum_squares(p0 - (1 - xnorm))
        cp.Problem(cp.Minimize(obj), cons).solve(solver=cp.ECOS)
        # ECOS is iterative and solves only to its own tolerance, so compare
        # loosely on the iterates and strictly on the objective.
        assert np.allclose(pi_k_ij[i, idx], pk.value, atol=1e-4), f"row {i}"
        f = lambda v: np.sum((v - xnorm) ** 2) + np.sum(((1 - v) - (1 - xnorm)) ** 2)
        assert f(pi_k_ij[i, idx]) <= f(pk.value) + 1e-12, (
            f"row {i}: projection is worse than the QP solution"
        )


def test_global_solver_matches_the_cvxpy_program():
    import cvxpy as cp
    D = 12
    xi, _ = _fixture(D, seed=5)
    pi0 = 0.8
    pi0_ij, pi_k_ij, _ = solve_spike_slab_diagonal_spike(xi, pi0=pi0)

    offdiag = np.ones((D, D), dtype=bool)
    np.fill_diagonal(offdiag, False)
    xi_norm = xi * ((1.0 - pi0) * (D ** 2 - D) / xi[offdiag].sum())

    pi0_var = cp.Variable((D, D), nonneg=True)
    pik_var = cp.Variable((D, D), nonneg=True)
    cons = [pi0_var[offdiag] + pik_var[offdiag] == 1.0,
            cp.diag(pi0_var) == 1.0, cp.diag(pik_var) == 0.0,
            cp.sum(pi0_var) - cp.sum(cp.diag(pi0_var)) == pi0 * (D ** 2 - D)]
    obj = (cp.sum_squares(cp.multiply(offdiag, pik_var - xi_norm))
           + cp.sum_squares(cp.multiply(offdiag, pi0_var - (1 - xi_norm))))
    cp.Problem(cp.Minimize(obj), cons).solve()

    assert np.allclose(pi_k_ij[offdiag], pik_var.value[offdiag], atol=1e-4)
    assert np.allclose(pi0_ij[offdiag], pi0_var.value[offdiag], atol=1e-4)
    f = lambda v: 2.0 * np.sum((v - xi_norm[offdiag]) ** 2)
    assert f(pi_k_ij[offdiag]) <= f(pik_var.value[offdiag]) + 1e-12, (
        "projection is worse than the QP solution"
    )


# --------------------------------------------------------------------------
# scale_free_degree
# --------------------------------------------------------------------------

def _reference_scale_free_degree_qp(R):
    """The cvxpy program scale_free_degree replaced. Returns P and its targets.

    Fits an edge-probability matrix P in [0, 1] to the rescaled in- and
    out-strengths of |R_hat|^2 by least squares.
    """
    import cvxpy as cp

    D = R.shape[0]
    A = np.abs(R) ** 2
    np.fill_diagonal(A, 0)
    theta_raw, phi_raw = A.sum(axis=1), A.sum(axis=0)
    scale = (D - 1) / max(theta_raw.max(), phi_raw.max())
    theta, phi = theta_raw * scale, phi_raw * scale

    A_out = np.zeros((D, D * D))
    A_in = np.zeros((D, D * D))
    for i in range(D):
        for j in range(D):
            if i != j:
                A_out[i, i * D + j] = 1
                A_in[j, i * D + j] = 1
    keep = np.ones(D * D, dtype=bool)
    for i in range(D):
        keep[i * D + i] = False
    A_full = np.vstack([A_out[:, keep], A_in[:, keep]])

    x = cp.Variable(A_full.shape[1])
    cp.Problem(cp.Minimize(cp.sum_squares(A_full @ x - np.concatenate([theta, phi]))),
               [x >= 0, x <= 1]).solve(solver=cp.ECOS)

    P = np.zeros((D, D))
    k = 0
    for i in range(D):
        for j in range(D):
            if i != j:
                P[i, j] = x.value[k]
                k += 1
    return P, theta, phi


def test_scale_free_degree_returns_valid_probabilities():
    D = 10
    _, R = _fixture(D)
    pi0_i = scale_free_degree(R)
    assert pi0_i.shape == (D,)
    assert pi0_i.min() >= 0.0 and pi0_i.max() <= 1.0


def _exact_rho(G):
    return float(np.abs(np.linalg.eigvals(np.asarray(G, dtype=np.float64))).max())


def test_gelfand_radius_bounds_and_tracks_the_spectral_radius():
    rng = np.random.default_rng(0)
    for scale in (0.02, 0.1, 1.0, 50.0):                  # healthy through runaway
        G = rng.normal(0, scale, (30, 30))
        exact = _exact_rho(G)
        est = float(gelfand_radius(G))
        assert est >= exact * (1 - 1e-4), (scale, est, exact)
        assert est <= exact * 1.15, (scale, est, exact)


def test_gelfand_radius_vanishes_on_a_dag_with_a_finite_gradient():
    import jax
    rng = np.random.default_rng(1)
    G = np.triu(rng.normal(0, 0.3, (20, 20)), 1)          # nilpotent: rho = 0
    # the log-space floor leaves a small residue instead of exactly zero, but
    # the gradient must stay finite there
    assert float(gelfand_radius(G)) < 0.1
    g = np.asarray(jax.grad(lambda x: gelfand_radius(x))(G))
    assert np.isfinite(g).all()
    # sampled G are never exactly nilpotent: the spike noise alone lifts rho
    # well above zero through non-normal amplification, and the bound follows
    Gn = G + rng.normal(0, 1e-3, G.shape); np.fill_diagonal(Gn, 0.0)
    exact = _exact_rho(Gn)
    assert exact > 0.1
    assert exact * (1 - 1e-4) <= float(gelfand_radius(Gn)) <= exact * 1.15


def test_power_iteration_radius_matches_main():
    """Port of spectral_radius from 42660bd: right when the dominant
    eigenvalue is real and well separated."""
    rng = np.random.default_rng(2)
    Q = np.linalg.qr(rng.normal(size=(12, 12)))[0]
    G = Q @ np.diag(np.r_[0.8, np.linspace(0.3, 0.05, 11)]) @ Q.T
    assert abs(float(power_iteration_radius(G)) - 0.8) < 1e-3


def test_power_iteration_radius_is_unreliable_on_a_complex_pair():
    """Why it is tested against the Gelfand bound: when the leading
    eigenvalues are a complex-conjugate pair the iterate rotates, so the
    estimate need not approach rho."""
    th = 2.0
    blk = 5.0 * np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    G = np.zeros((6, 6)); G[:2, :2] = blk; G[:2, :2] += np.array([[0.0, -20.0], [0.0, 0.0]])   # non-normal, still complex
    G[2:, 2:] = np.diag([0.4, 0.3, 0.2, 0.1])
    ev = np.linalg.eigvals(G)
    assert abs(ev[np.argmax(np.abs(ev))].imag) > 1.0      # the premise of the test
    exact = _exact_rho(G)
    est = float(power_iteration_radius(G))
    assert abs(est / exact - 1.0) > 0.1, (est, exact)
    assert float(gelfand_radius(G)) >= exact * (1 - 1e-4)


def test_rho_barrier_is_zero_inside_and_grows_outside():
    kw = dict(start=0.9, width=0.05)
    assert float(rho_log_penalty(0.5, "barrier", **kw)) == 0.0
    assert float(rho_log_penalty(0.9, "barrier", **kw)) == 0.0
    a = float(rho_log_penalty(1.0, "barrier", **kw)); b = float(rho_log_penalty(1.5, "barrier", **kw))
    assert a < 0.0 and b < a
    assert np.isclose(a, -((1.0 - 0.9) / 0.05) ** 2)


def test_rho_barrier_defaults_are_the_soft_barrier():
    """start 1.0, width 0.2: free up to rho = 1, one nat at rho = 1.2."""
    assert float(rho_log_penalty(1.0, "barrier")) == 0.0
    assert np.isclose(float(rho_log_penalty(1.2, "barrier")), -1.0)


def test_rho_gaussian_is_the_appendix_h_prior():
    from scipy import stats
    for r in (0.0, 0.5, 2.0):
        assert np.isclose(float(rho_log_penalty(r, "gaussian", sigma=0.5)),
                          stats.norm(0, 0.5).logpdf(r))


def test_model_log_density_is_finite_under_each_rho_penalty():
    import jax
    import jax.numpy as jnp
    from numpyro.infer.util import initialize_model
    D = 6
    rng = np.random.default_rng(3)
    kw = dict(obs_data=np.eye(D) + rng.normal(0, 0.05, (D, D)),
              pi0_ij=np.full((D, D), 0.8), U_lower=jnp.eye(D) * 0.1,
              V_lower=jnp.eye(D), D=D)
    for pen, est in [(None, "power"), ("gaussian", "power"),
                     ("gaussian", "gelfand"), ("barrier", "gelfand")]:
        info = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe,
                                model_kwargs=dict(kw, rho_penalty=pen, rho_estimator=est))
        z = info.param_info.z
        pe = float(info.potential_fn(z))
        grads = jax.grad(info.potential_fn)(z)
        assert np.isfinite(pe), (pen, est)
        assert all(np.isfinite(np.asarray(v)).all() for v in grads.values()), (pen, est)


def _small_model_kwargs(D=6, seed=3):
    import jax.numpy as jnp
    rng = np.random.default_rng(seed)
    G = np.triu(rng.normal(0, 0.3, (D, D)) * (rng.random((D, D)) < 0.5), 1)
    R = np.linalg.inv(np.eye(D) - G)
    return dict(obs_data=R + rng.normal(0, 0.02, (D, D)), pi0_ij=np.full((D, D), 0.5),
                U_lower=jnp.eye(D) * 0.05, V_lower=jnp.eye(D), D=D)


def test_optimized_init_lowers_the_potential_and_gives_distinct_chains():
    import jax
    kw = _small_model_kwargs()
    z, info = optimized_init(matrix_model_spike_horseshoe, kw, jax.random.PRNGKey(0),
                             num_chains=3, steps=300, lr=1e-2, jitter=0.1)
    for name in ("lam", "eps", "spike"):
        assert z[name].shape == (3, 6, 6), name
        assert np.isfinite(np.asarray(z[name])).all(), name
    before, after = np.array(info["potential_before"]), np.array(info["potential_after"])
    assert (after < before).all(), (before, after)
    assert len(info["rho_start_per_chain"]) == 3
    # the jitter keeps the chains from starting at one point
    eps = np.asarray(z["eps"])
    assert not np.allclose(eps[0], eps[1])


def test_optimized_init_starts_are_accepted_by_nuts():
    import jax
    from numpyro.infer import MCMC, NUTS
    kw = _small_model_kwargs()
    z, _ = optimized_init(matrix_model_spike_horseshoe, kw, jax.random.PRNGKey(1),
                          num_chains=2, steps=100)
    m = MCMC(NUTS(matrix_model_spike_horseshoe, max_tree_depth=4), num_warmup=5,
             num_samples=5, num_chains=2, chain_method="vectorized", progress_bar=False)
    m.run(jax.random.PRNGKey(2), init_params=z, **kw)
    assert np.isfinite(np.asarray(m.get_samples()["G"])).all()


def test_rhat_start_soft_thresholds_at_k_standard_errors():
    R = np.array([[1.0, 0.5, -0.05], [0.02, 1.0, -0.4], [0.0, 0.1, 1.0]])
    S = np.array([[0.0, 0.1, 0.1], [0.1, 0.0, 0.1], [0.1, 0.1, 0.0]])
    G0 = rhat_start(R, S, k=3)
    assert np.allclose(np.diag(G0), 0.0)
    assert np.isclose(G0[0, 1], 0.2) and np.isclose(G0[1, 2], -0.1)     # shrunk, sign kept
    assert G0[0, 2] == 0.0 and G0[1, 0] == 0.0 and G0[2, 1] == 0.0      # below 3 SE


def test_latents_from_G_reproduces_G_through_the_model():
    import jax
    from numpyro.infer.util import initialize_model
    kw = _small_model_kwargs()
    rng = np.random.default_rng(4)
    G = rng.normal(0, 0.3, (6, 6)) * (rng.random((6, 6)) < 0.5); np.fill_diagonal(G, 0.0)
    mi = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe, model_kwargs=kw)
    z = latents_from_G(G, kw["pi0_ij"], tau=0.1)
    assert np.abs(np.asarray(z["eps"])).max() <= 1.0 + 1e-9
    assert np.allclose(np.asarray(mi.postprocess_fn(z)["G"]), G, atol=1e-6)


def test_latents_from_G_takes_the_most_probable_split():
    """For each entry the split of lam * eps maximises the prior density in
    the sampler's coordinates (log lam), and a zero entry gets the values
    init_to_median gives."""
    pi0 = np.full((3, 3), 0.5); tau = 0.1; m = 0.5 * tau
    G = np.array([[0.0, 0.03, -0.2], [1.5, 0.0, 0.0], [0.0, -0.001, 0.0]])
    z = latents_from_G(G, pi0, tau)
    lam = np.exp(np.asarray(z["lam"])); eps = np.asarray(z["eps"])
    assert np.allclose(lam[G == 0], 1.0) and np.allclose(eps[G == 0], 0.0)
    f = lambda l, e: -np.log1p(l * l) + np.log(l) - 0.5 * e * e
    for i, j in zip(*np.nonzero(G)):
        a = abs(G[i, j]) / m
        best = f(lam[i, j], eps[i, j])
        for scale in (0.9, 0.99, 1.01, 1.1):          # other splits with the same product
            l = lam[i, j] * scale
            assert f(l, np.sign(G[i, j]) * a / l) < best


def test_em_slab_scale_is_the_rms_of_the_slab():
    pi_k = np.array([0.3, 0.1, 0.0]); sigma_k = np.array([0.1, 0.5, 1.0])
    assert np.isclose(em_slab_scale(pi_k, sigma_k), np.sqrt((0.3 * 0.01 + 0.1 * 0.25) / 0.4))


def test_regularised_slab_saturates_at_its_width():
    """With a slab width c, an entry's whole slab scale saturates at c however
    large lam gets, which is what removes the Cauchy tail; every entry can
    still reach c, whatever its (1 - pi0)."""
    import jax
    import jax.numpy as jnp
    from numpyro.infer.util import initialize_model
    kw = dict(_small_model_kwargs(), slab_width=0.2)
    mi = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe, model_kwargs=kw)
    z = {"lam": jnp.full((6, 6), np.log(1e6)), "eps": jnp.ones((6, 6)), "spike": jnp.zeros((6, 6))}
    G = np.asarray(mi.postprocess_fn(z)["G"])
    off = ~np.eye(6, dtype=bool)
    assert np.allclose(G[off], 0.2, rtol=1e-3)
    # and without one, the same latents give an enormous G
    mi0 = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe, model_kwargs=_small_model_kwargs())
    assert np.asarray(mi0.postprocess_fn(z)["G"])[off].min() > 1e3


def test_latents_from_G_reproduces_G_with_a_slab_width():
    import jax
    from numpyro.infer.util import initialize_model
    kw = dict(_small_model_kwargs(), slab_width=0.2)
    rng = np.random.default_rng(6)
    G = rng.normal(0, 0.3, (6, 6)) * (rng.random((6, 6)) < 0.5); np.fill_diagonal(G, 0.0)
    mi = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe, model_kwargs=kw)
    z = latents_from_G(G, kw["pi0_ij"], tau=0.1, slab_width=0.2)
    assert np.allclose(np.asarray(mi.postprocess_fn(z)["G"]), G, atol=1e-6)


def test_optimized_init_from_a_start_G():
    import jax
    kw = _small_model_kwargs()
    G0 = rhat_start(kw["obs_data"], np.full((6, 6), 0.01), k=3)
    z, info = optimized_init(matrix_model_spike_horseshoe, kw, jax.random.PRNGKey(0),
                             num_chains=3, steps=200, start_G=G0)
    assert z["eps"].shape == (3, 6, 6)
    assert all(a <= b for a, b in zip(info["potential_after"], info["potential_before"]))


def test_optimized_init_with_no_steps_starts_at_start_G():
    import jax
    from numpyro.infer.util import initialize_model
    kw = _small_model_kwargs()
    G0 = rhat_start(kw["obs_data"], np.full((6, 6), 0.01), k=3)
    z, info = optimized_init(matrix_model_spike_horseshoe, kw, jax.random.PRNGKey(0),
                             num_chains=2, steps=0, jitter=0.0, start_G=G0)
    mi = initialize_model(jax.random.PRNGKey(0), matrix_model_spike_horseshoe, model_kwargs=kw)
    for c in range(2):
        Gc = np.asarray(mi.postprocess_fn(jax.tree_util.tree_map(lambda x: x[c], z))["G"])
        assert np.allclose(Gc, G0, atol=1e-6)
    assert np.allclose(info["potential_after"], info["potential_before"])


def test_cap_and_redistribute_preserves_the_total():
    b = np.array([10.0, 1.0, 1.0, 1.0])
    out = _cap_and_redistribute(b, cap=5.0)
    assert out.max() <= 5.0 + 1e-9
    assert np.isclose(out.sum(), b.sum())


def test_cap_and_redistribute_clips_when_the_total_cannot_fit():
    b = np.array([10.0, 10.0])
    out = _cap_and_redistribute(b, cap=3.0)
    assert np.allclose(out, 3.0)


def test_anchored_scale_free_degree_matches_the_eq_18_budget():
    """The EM level sets the total slab mass; theta only shares it out."""
    D = 12
    _, R = _fixture(D)
    pi0_global = 0.8
    pi0_i = scale_free_degree(R, pi0_global=pi0_global)
    budget = (1.0 - pi0_i) * (D - 1)
    assert np.isclose(budget.sum(), (1.0 - pi0_global) * (D * D - D))
    assert pi0_i.min() >= 0.0 and pi0_i.max() <= 1.0


def test_anchored_scale_free_degree_leaves_every_row_shrunk():
    """No row may be handed the whole row as budget, which is zero shrinkage."""
    D = 12
    _, R = _fixture(D)
    floor = 0.05
    pi0_i = scale_free_degree(R, pi0_global=0.1, pi0_floor=floor)
    assert pi0_i.min() >= floor - 1e-9


def test_anchored_scale_free_degree_is_not_hostage_to_the_largest_node():
    """The defect being fixed: under the max-normalisation pi0_i = 1 - theta_i/m,
    one dominant hub inflates m and drives every other node's budget toward
    zero while taking the whole of its own row. Anchoring the total to the EM
    level makes the other rows' budgets robust to that."""
    D = 12
    _, R = _fixture(D)
    hub = R.copy()
    hub[0, 1:] *= 50.0                      # one node with a huge out-strength

    legacy, legacy_hub = scale_free_degree(R), scale_free_degree(hub)
    anchored = scale_free_degree(R, pi0_global=0.9)
    anchored_hub = scale_free_degree(hub, pi0_global=0.9)

    def share_left_to_others(p):
        b = (1.0 - p) * (D - 1)
        return float(b[1:].sum() / b.sum())

    # legacy: the hub is handed its whole row, with no shrinkage at all, and
    # what is left for the other D-1 rows is a rounding error
    assert legacy_hub[0] < 1e-12
    assert share_left_to_others(legacy_hub) < 0.01
    # anchored: the hub keeps the pi0_floor of shrinkage, and the other rows
    # retain a usable share of a total that no longer depends on the hub
    assert anchored_hub[0] >= 0.05 - 1e-9
    assert share_left_to_others(anchored_hub) > 0.15
    assert share_left_to_others(anchored_hub) > 10 * share_left_to_others(legacy_hub)
    # the anchored total is the EM budget whether or not the hub is present
    assert np.isclose(((1.0 - anchored_hub) * (D - 1)).sum(),
                      ((1.0 - anchored) * (D - 1)).sum())


def test_anchored_scale_free_degree_responds_to_the_em_level():
    """Unlike the max-normalisation, which is invariant to R -> cR, the
    anchored version tracks the estimated sparsity."""
    D = 12
    _, R = _fixture(D)
    sparse = scale_free_degree(R, pi0_global=0.95)
    dense = scale_free_degree(R, pi0_global=0.5)
    assert sparse.mean() > dense.mean()


def test_posterior_diagnostics_bounds_the_rhat_subsample():
    rng = np.random.default_rng(0)
    D, n_chains, n_draws = 8, 3, 40
    g = rng.normal(0, 0.01, (n_chains, n_draws, D, D))
    keep, rho = convergent_draws(g.reshape(-1, D, D))
    out = posterior_diagnostics(g, rho, keep, seed=0, max_rhat_entries=17)
    assert out["r_hat"]["n_entries_used"] == 17
    assert out["r_hat"]["n_entries_total"] == D * D - D


def test_scale_free_degree_is_all_spike_when_R_carries_no_signal():
    R = np.eye(5)
    assert np.allclose(scale_free_degree(R), np.ones(5))


def test_scale_free_degree_matches_the_cvxpy_program():
    D = 10
    _, R = _fixture(D)
    pi0_i = scale_free_degree(R)
    P, theta, phi = _reference_scale_free_degree_qp(R)

    # only the row means of P feed the SF prior
    assert np.allclose(pi0_i, 1.0 - P.sum(axis=1) / (D - 1), atol=1e-4)

    # the closed form attains the out-strength targets exactly; the QP, being
    # iterative, does not. ECOS misses them badly on some real inputs.
    cf_rows = (1.0 - pi0_i) * (D - 1)
    assert np.allclose(cf_rows, theta, atol=1e-12)
    assert (np.sum((cf_rows - theta) ** 2)
            <= np.sum((P.sum(axis=1) - theta) ** 2) + 1e-12)


# --------------------------------------------------------------------------
# Truncated path sum and draw convergence
# --------------------------------------------------------------------------

def keep_mask_count(diag):
    return diag["n_draws_total"] - diag["excluded"]["n"]


def test_posterior_diagnostics_reports_what_the_run_did():
    rng = np.random.default_rng(0)
    C, N, D = 2, 30, 6
    draws = np.stack([np.stack([_dag(D, seed=c * N + n) for n in range(N)])
                      for c in range(C)])
    draws[1, :5] *= 0.0                      # make a few draws trivially distinct
    keep = np.ones(C * N, dtype=bool)
    keep[[3, 7]] = False                     # pretend two draws were rejected
    rho = rng.random(C * N)
    extra = {"diverging": np.zeros((C, N), dtype=bool),
             "num_steps": np.full((C, N), 7),
             "accept_prob": np.full((C, N), 0.83),
             "adapt_state.step_size": np.full((C, N), 3e-3)}
    extra["diverging"][0, :4] = True

    d = posterior_diagnostics(draws, rho, keep, extra_fields=extra)
    assert d["n_chains"] == C and d["n_draws_total"] == C * N and d["D"] == D
    assert d["excluded"]["n"] == 2                       # what the filter removed
    assert d["nonconvergent"]["n"] == int((np.asarray(rho) >= 1).sum())   # rho >= 1, kept or not
    assert d["divergences"]["n"] == 4 and d["divergences"]["per_chain"] == [4, 0]
    assert d["leapfrog"]["total"] == C * N * 7
    assert np.isclose(d["step_size"]["median"], 3e-3)
    assert d["step_size"]["per_chain"] == [3e-3, 3e-3]
    assert np.isclose(d["accept_prob"]["mean"], 0.83)

    # the saturation cap follows max_tree_depth rather than being hardcoded
    assert d["leapfrog"]["max_leapfrog_steps"] == 2 ** 12 - 1   # the default
    extra_sat = dict(extra, num_steps=np.full((C, N), 2 ** 12 - 1))
    d10 = posterior_diagnostics(draws, rho, keep, extra_fields=extra_sat,
                                max_tree_depth=10)
    d12 = posterior_diagnostics(draws, rho, keep, extra_fields=extra_sat,
                                max_tree_depth=12)
    assert d10["leapfrog"]["max_leapfrog_steps"] == 1023
    assert d12["leapfrog"]["max_leapfrog_steps"] == 4095
    # 4095 steps saturates at depth 12 and also exceeds the depth-10 cap
    assert d10["leapfrog"]["pct_at_max_tree_depth"] == 100.0
    assert d12["leapfrog"]["pct_at_max_tree_depth"] == 100.0
    # ... but 1023 steps only saturates at depth 10
    extra_mid = dict(extra, num_steps=np.full((C, N), 1023))
    assert posterior_diagnostics(draws, rho, keep, extra_fields=extra_mid,
                                 max_tree_depth=12)["leapfrog"]["pct_at_max_tree_depth"] == 0.0
    assert d["spectral_radius"]["max"] <= 1.0
    import json as _json
    _json.dumps(d)                            # must be serialisable


def _dag(D, seed=0, density=0.15, scale=0.25):
    """A strictly upper-triangular G, i.e. a DAG in the given variable order."""
    rng = np.random.default_rng(seed)
    G = np.triu(rng.normal(0, scale, size=(D, D)), 1)
    G *= rng.random((D, D)) < density
    return G


def _series(G, order):
    D = G.shape[0]
    R = np.eye(D)
    for _ in range(order):
        R = np.eye(D) + G @ R
    return R


def test_series_equals_inverse_at_sufficient_order():
    """For a DAG the two agree exactly once the order reaches the longest path."""
    D = 20
    G = _dag(D, seed=2)
    exact = np.linalg.inv(np.eye(D) - G)
    # G is nilpotent: G^D is identically zero, so order D-1 suffices
    assert np.allclose(_series(G, D - 1), exact, atol=1e-12)
    assert np.allclose(np.linalg.matrix_power(G, D), 0.0)


def test_series_stays_bounded_where_the_inverse_blows_up():
    """Near the singularity of (I - G) the inverse explodes; the sum does not.

    This is the numerical reason to prefer the series: the inverse has a pole
    the sampler can approach, and its gradients scale as ||(I - G)^-1||^2.
    """
    D = 20
    G = _dag(D, seed=2)
    G = G + 3.0 * G.T                                     # make it cyclic
    G = G / (np.abs(np.linalg.eigvals(G)).max() * 1.001)  # push rho just under 1

    inv_norm = np.linalg.norm(np.linalg.inv(np.eye(D) - G), 2)
    series_norm = np.linalg.norm(_series(G, 24), 2)
    assert np.isfinite(series_norm)
    assert inv_norm > 20 * series_norm, (
        f"inverse {inv_norm:.1f} vs series {series_norm:.1f}"
    )


def test_convergent_draws_separates_runaway_draws():
    D = 8
    good = np.stack([_dag(D, seed=s) for s in range(5)])
    bad = good * 500.0  # far outside the region where sum_d G^d converges
    draws = np.concatenate([good, bad])

    keep, rho = convergent_draws(draws)
    assert keep.shape == (10,) and rho.shape == (10,)
    # a nilpotent G stays nilpotent under scaling, so rho is 0 for all of these
    assert keep.all() and np.allclose(rho, 0.0)

    # a genuinely cyclic draw with rho >= 1 must be dropped
    cyc = np.stack([g + 3.0 * g.T for g in good]) * 5.0
    keep2, rho2 = convergent_draws(np.concatenate([good, cyc]))
    assert keep2[:5].all()
    assert not keep2[5:].any(), f"cyclic draws must be dropped, rho={rho2[5:]}"


def test_convergent_draws_threshold_is_the_spectral_radius():
    D = 6
    G = np.zeros((D, D))
    G[0, 1] = G[1, 0] = 0.4          # 2-cycle, rho = 0.4
    keep, rho = convergent_draws(G[None])
    assert np.isclose(rho[0], 0.4) and keep[0]
    G[0, 1] = G[1, 0] = 1.2          # rho = 1.2
    keep, rho = convergent_draws(G[None])
    assert np.isclose(rho[0], 1.2) and not keep[0]


# --------------------------------------------------------------------------
# combine_chains
# --------------------------------------------------------------------------

def _write_chain(path, draws, colnames=None):
    """Lay out a directory the way a single-chain ibcd.py run does."""
    path.mkdir(parents=True, exist_ok=True)
    np.save(path / "G_draws.npy", draws.astype(np.float32))
    if colnames is not None:
        mean = draws.reshape(-1, draws.shape[-2], draws.shape[-1]).mean(axis=0)
        pd.DataFrame(mean, columns=colnames, index=colnames).to_csv(path / "G.csv")
    return path


def _chain_of_dags(n_draws, D, seed):
    return np.stack([_dag(D, seed=seed * 1000 + n) for n in range(n_draws)])[None]


def _combine(chain_dirs, out, epsilon=0.05, max_spectral_radius=None):
    import argparse
    import combine_chains
    combine_chains.main(argparse.Namespace(
        chain_dirs=[str(c) for c in chain_dirs], output_dir=str(out),
        epsilon=epsilon, seed=42, max_spectral_radius=max_spectral_radius,
    ))
    return {n: pd.read_csv(out / f"{n}.csv", index_col=0)
            for n in ("G", "pip", "lfsr")}, json.load((out / "diagnostics.json").open())


def test_load_chains_stacks_and_recovers_names():
    import combine_chains
    D, N = 6, 4
    names = [f"V{i+1}" for i in range(D)]
    with tempfile.TemporaryDirectory() as tmp:
        dirs = [_write_chain(Path(tmp) / f"chain{c}", _chain_of_dags(N, D, c), names)
                for c in range(3)]
        draws, colnames = combine_chains.load_chains([str(d) for d in dirs])
    assert draws.shape == (3, N, D, D)
    assert colnames == names


def test_load_chains_rejects_inconsistent_or_missing_input():
    import combine_chains
    with tempfile.TemporaryDirectory() as tmp:
        a = _write_chain(Path(tmp) / "a", _chain_of_dags(4, 6, 0))
        b = _write_chain(Path(tmp) / "b", _chain_of_dags(4, 8, 1))   # different D
        try:
            combine_chains.load_chains([str(a), str(b)])
        except ValueError as exc:
            assert "disagree" in str(exc)
        else:
            raise AssertionError("mismatched shapes must raise")

        try:
            combine_chains.load_chains([str(a), str(Path(tmp) / "nope")])
        except FileNotFoundError:
            pass
        else:
            raise AssertionError("a missing G_draws.npy must raise")


def test_combine_chains_matches_concatenated_draws():
    D, N = 6, 5
    names = [f"V{i+1}" for i in range(D)]
    with tempfile.TemporaryDirectory() as tmp:
        chains = [_chain_of_dags(N, D, c) for c in range(3)]
        dirs = [_write_chain(Path(tmp) / f"chain{c}", chains[c], names)
                for c in range(3)]
        out, diag = _combine(dirs, Path(tmp) / "combined")

        flat = np.concatenate(chains, axis=0).reshape(-1, D, D).astype(np.float32)
        assert np.allclose(out["G"].values, flat.mean(axis=0), atol=1e-6)
        assert np.allclose(out["pip"].values, (np.abs(flat) > 0.05).mean(axis=0),
                           atol=1e-6)
        for df in out.values():
            assert list(df.columns) == names and list(df.index) == names
        assert diag["n_chains"] == 3 and diag["n_draws_per_chain"] == N
        assert diag["nonconvergent"]["n"] == 0


def _runaway_chains(tmp, D, N):
    good = [_chain_of_dags(N, D, c) for c in range(2)]
    bad = _chain_of_dags(N, D, 2).copy()
    bad[:, :, 0, 1] = 2.0          # a 2-cycle of weight 2 gives rho = 2
    bad[:, :, 1, 0] = 2.0
    dirs = [_write_chain(Path(tmp) / f"chain{c}", ch) for c, ch in enumerate(good + [bad])]
    return good, bad, dirs


def test_combine_chains_keeps_every_draw_by_default():
    """Without --max_spectral_radius nothing is excluded; rho >= 1 is reported."""
    D, N = 6, 5
    with tempfile.TemporaryDirectory() as tmp:
        good, bad, dirs = _runaway_chains(tmp, D, N)
        out, diag = _combine(dirs, Path(tmp) / "combined")
        assert diag["excluded"]["n"] == 0
        assert diag["nonconvergent"]["per_chain"] == [0, 0, N]
        allv = np.concatenate(good + [bad], axis=0).reshape(-1, D, D).astype(np.float32)
        assert np.allclose(out["G"].values, allv.mean(axis=0), atol=1e-6)


def test_screen_draws_keeps_everything_by_default_and_cuts_when_asked():
    import warnings
    rho = np.array([0.5, 0.8, 1.5, 3.0, 0.2, 0.9])
    assert screen_draws(rho, 2).all()
    assert (screen_draws(rho, 2, max_spectral_radius=1.0) == (rho < 1.0)).all()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        screen_draws(np.array([0.5, 12.0]), 1)
        assert any("exceeds 10" in str(x.message) for x in w)


def test_combine_chains_excludes_nonconvergent_draws():
    """With --max_spectral_radius 1, a runaway chain is dropped from the summaries."""
    D, N = 6, 5
    with tempfile.TemporaryDirectory() as tmp:
        good = [_chain_of_dags(N, D, c) for c in range(2)]
        bad = _chain_of_dags(N, D, 2).copy()
        bad[:, :, 0, 1] = 2.0          # a 2-cycle of weight 2 gives rho = 2
        bad[:, :, 1, 0] = 2.0
        assert (np.abs(np.linalg.eigvals(bad[0])).max(axis=1) >= 1).all()
        dirs = [_write_chain(Path(tmp) / f"chain{c}", ch)
                for c, ch in enumerate(good + [bad])]
        out, diag = _combine(dirs, Path(tmp) / "combined", max_spectral_radius=1.0)

        assert diag["excluded"]["n"] == N
        assert diag["excluded"]["per_chain"] == [0, 0, N]
        assert diag["nonconvergent"]["per_chain"] == [0, 0, N]
        # the combined mean must use only the two healthy chains
        kept = np.concatenate(good, axis=0).reshape(-1, D, D).astype(np.float32)
        assert np.allclose(out["G"].values, kept.mean(axis=0), atol=1e-6)


# --------------------------------------------------------------------------
# End to end (slow)
# --------------------------------------------------------------------------

def _run_pipeline(truncated_series, save_diagnostics=True, rho_penalty="none",
                  rho_estimator="power", init_strategy="median", sf_anchor="em",
                  slab_width="none", init_G_matrix=None):
    """Run the real pipeline on a small subset of the shipped example data."""
    import argparse
    import ibcd

    data = pd.read_csv(REPO / "data" / "input" / "data.csv")
    keep = [f"V{i}" for i in range(1, 11)]
    sub = data[data["target"].isin(["control"] + keep)][keep + ["target"]]

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "data.csv"
        sub.to_csv(path, index=False)
        out = Path(tmp) / "out"
        init_G = None
        if init_G_matrix is not None:
            init_G = str(Path(tmp) / "G_start.csv")
            pd.DataFrame(init_G_matrix, columns=keep).to_csv(init_G, index=False)
        ibcd.main(argparse.Namespace(
            data=str(path), prior="sf", output_dir=str(out),
            alpha_er=2.0, pi0_floor=0.05, sf_anchor=sf_anchor,
            num_warmup=20, num_samples=40, num_chains=1, epsilon=0.05,
            chain_method="vectorized",
            target_accept_prob=0.7, max_tree_depth=10,
            truncated_series=truncated_series, series_order=24, seed=42,
            save_diagnostics=save_diagnostics,
            rho_penalty=rho_penalty, rho_estimator=rho_estimator, rho_sigma=0.5,
            rho_barrier_start=1.0, rho_barrier_width=0.2,
            init_strategy=init_strategy, init_opt_steps=200, init_opt_lr=0.01,
            init_jitter=0.1, init_rhat_k=3.0, init_G=init_G, slab_width=slab_width,
            max_spectral_radius=None,
        ))

        pip = pd.read_csv(out / "pip.csv", index_col=0)
        G = pd.read_csv(out / "G.csv", index_col=0)
        lfsr = pd.read_csv(out / "lfsr.csv", index_col=0)
        diag_path = out / "diagnostics.json"
        assert diag_path.exists() == save_diagnostics, (
            "diagnostics.json presence must follow --save_diagnostics"
        )
        diag = json.load(diag_path.open()) if save_diagnostics else None
    return keep, pip, G, lfsr, diag


def _check_outputs(keep, pip, G, lfsr, diag):
    D = len(keep)
    for name, df in [("pip", pip), ("G", G), ("lfsr", lfsr)]:
        assert df.shape == (D, D), f"{name}.csv has shape {df.shape}"
        assert list(df.columns) == keep, f"{name}.csv columns not preserved"
        assert list(df.index) == keep, f"{name}.csv index not preserved"
        assert np.isfinite(df.values).all(), f"{name}.csv contains non-finite values"

    assert pip.values.min() >= 0.0 and pip.values.max() <= 1.0
    assert lfsr.values.min() >= 0.0 and lfsr.values.max() <= 0.5
    # no self-loops: the diagonal of G is zeroed in the model
    assert np.allclose(np.diag(G.values), 0.0)
    assert np.allclose(np.diag(pip.values), 0.0)
    # the posterior should not be degenerate: some edges get real support
    assert pip.values.max() > 0.5, "no edge reached PIP > 0.5"

    if diag is None:
        return
    # when requested, the run must leave a usable diagnostics record
    for key in ("n_chains", "n_draws_total", "nonconvergent", "spectral_radius",
                "divergences", "leapfrog", "step_size", "accept_prob",
                "runtime_seconds", "config"):
        assert key in diag, f"diagnostics.json missing {key}"
    assert diag["config"]["seed"] == 42
    assert diag["config"]["target_accept_prob"] == 0.7
    assert diag["config"]["max_tree_depth"] == 10
    assert diag["excluded"]["n"] + int(keep_mask_count(diag)) == diag["n_draws_total"]
    if diag["config"]["max_spectral_radius"] is None:      # the default keeps every draw
        assert diag["excluded"]["n"] == 0


def test_end_to_end_outputs_are_well_formed():
    _check_outputs(*_run_pipeline(truncated_series=False, save_diagnostics=True))


def test_end_to_end_outputs_are_well_formed_with_truncated_series():
    # also covers the default, where no diagnostics file is written
    _check_outputs(*_run_pipeline(truncated_series=True, save_diagnostics=False))


def test_end_to_end_outputs_are_well_formed_with_a_rho_penalty():
    keep, pip, G, lfsr, diag = _run_pipeline(
        truncated_series=False, save_diagnostics=True,
        rho_penalty="barrier", rho_estimator="gelfand")
    _check_outputs(keep, pip, G, lfsr, diag)
    assert diag["config"]["rho_penalty"] == "barrier"
    assert diag["config"]["rho_estimator"] == "gelfand"
    assert diag["rho_estimate"]["median"] >= 0.0
    # gelfand is an upper bound, so it should not sit far below the exact rho
    assert diag["rho_estimate"]["median_ratio_to_exact"] > 0.9


def test_end_to_end_outputs_are_well_formed_with_an_optimized_start():
    keep, pip, G, lfsr, diag = _run_pipeline(
        truncated_series=False, save_diagnostics=True,
        rho_penalty="barrier", rho_estimator="gelfand", init_strategy="optimized")
    _check_outputs(keep, pip, G, lfsr, diag)
    assert diag["config"]["init_strategy"] == "optimized"
    assert diag["config"]["rho_barrier_start"] == 1.0
    assert len(diag["init"]["rho_start_per_chain"]) == diag["n_chains"]
    assert all(a < b for a, b in zip(diag["init"]["potential_after"],
                                     diag["init"]["potential_before"]))


def test_end_to_end_outputs_are_well_formed_without_the_sf_anchor():
    keep, pip, G, lfsr, diag = _run_pipeline(
        truncated_series=False, save_diagnostics=True, sf_anchor="none")
    _check_outputs(keep, pip, G, lfsr, diag)
    assert diag["config"]["sf_anchor"] == "none"


def test_end_to_end_outputs_are_well_formed_with_the_rhat_start():
    keep, pip, G, lfsr, diag = _run_pipeline(
        truncated_series=False, save_diagnostics=True, init_strategy="rhat")
    _check_outputs(keep, pip, G, lfsr, diag)
    assert diag["config"]["init_strategy"] == "rhat"
    assert diag["config"]["init_rhat_k"] == 3.0
    assert len(diag["init"]["rho_start_per_chain"]) == diag["n_chains"]


def test_end_to_end_outputs_are_well_formed_with_a_slab_width_and_a_file_start():
    rng = np.random.default_rng(7)
    G0 = np.triu(rng.normal(0, 0.2, (10, 10)) * (rng.random((10, 10)) < 0.3), 1)
    keep, pip, G, lfsr, diag = _run_pipeline(
        truncated_series=False, save_diagnostics=True, init_strategy="file",
        slab_width="em", init_G_matrix=G0)
    _check_outputs(keep, pip, G, lfsr, diag)
    assert diag["config"]["init_strategy"] == "file"
    assert diag["config"]["slab_width"] == "em"
    assert 0.1 <= diag["config"]["slab_width_value"] <= 1.0      # within the EM's grid


SLOW = {"test_end_to_end_outputs_are_well_formed_with_a_slab_width_and_a_file_start",
        "test_end_to_end_outputs_are_well_formed_with_the_rhat_start",
        "test_end_to_end_outputs_are_well_formed_without_the_sf_anchor",
        "test_end_to_end_outputs_are_well_formed_with_an_optimized_start",
        "test_end_to_end_outputs_are_well_formed",
        "test_end_to_end_outputs_are_well_formed_with_truncated_series",
        "test_end_to_end_outputs_are_well_formed_with_a_rho_penalty"}


def _main():
    run_slow = "--slow" in sys.argv or os.environ.get("IBCD_SLOW") == "1"
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    failed = skipped = 0
    for name, fn in tests:
        if name in SLOW and not run_slow:
            print(f"SKIP  {name}  (pass --slow to run)")
            skipped += 1
            continue
        try:
            fn()
            print(f"PASS  {name}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"FAIL  {name}: {type(exc).__name__}: {exc}")
    total = len(tests) - skipped
    print(f"\n{total - failed}/{total} passed"
          + (f", {skipped} skipped" if skipped else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_main())
