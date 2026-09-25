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

from model import convergent_draws, posterior_diagnostics  # noqa: E402
from empirical_prior import (  # noqa: E402
    _project_to_budget,
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
    return diag["n_draws_total"] - diag["nonconvergent"]["n"]


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
             "num_steps": np.full((C, N), 7)}
    extra["diverging"][0, :4] = True

    d = posterior_diagnostics(draws, rho, keep, extra_fields=extra)
    assert d["n_chains"] == C and d["n_draws_total"] == C * N and d["D"] == D
    assert d["nonconvergent"]["n"] == 2
    assert d["divergences"]["n"] == 4 and d["divergences"]["per_chain"] == [4, 0]
    assert d["leapfrog"]["total"] == C * N * 7
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
# End to end (slow)
# --------------------------------------------------------------------------

def _run_pipeline(truncated_series):
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
        ibcd.main(argparse.Namespace(
            data=str(path), prior="sf", output_dir=str(out),
            alpha_er=2.0,
            num_warmup=20, num_samples=40, num_chains=1, epsilon=0.05,
            truncated_series=truncated_series, series_order=24, seed=42,
        ))

        pip = pd.read_csv(out / "pip.csv", index_col=0)
        G = pd.read_csv(out / "G.csv", index_col=0)
        lfsr = pd.read_csv(out / "lfsr.csv", index_col=0)
        with open(out / "diagnostics.json") as fh:
            diag = json.load(fh)
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

    # the run must leave a usable diagnostics record
    for key in ("n_chains", "n_draws_total", "nonconvergent", "spectral_radius",
                "divergences", "leapfrog", "runtime_seconds", "config"):
        assert key in diag, f"diagnostics.json missing {key}"
    assert diag["config"]["seed"] == 42
    assert diag["nonconvergent"]["n"] + int(keep_mask_count(diag)) == diag["n_draws_total"]


def test_end_to_end_outputs_are_well_formed():
    _check_outputs(*_run_pipeline(truncated_series=False))


def test_end_to_end_outputs_are_well_formed_with_truncated_series():
    _check_outputs(*_run_pipeline(truncated_series=True))


SLOW = {"test_end_to_end_outputs_are_well_formed",
        "test_end_to_end_outputs_are_well_formed_with_truncated_series"}


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
