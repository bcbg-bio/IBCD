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

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src"
sys.path.insert(0, str(SRC))

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

def test_scale_free_degree_returns_valid_probabilities():
    D = 10
    _, R = _fixture(D)
    pi0 = scale_free_degree(R)
    assert pi0.shape == (D, D)
    assert pi0.min() >= -1e-6 and pi0.max() <= 1 + 1e-6


# --------------------------------------------------------------------------
# End to end (slow)
# --------------------------------------------------------------------------

def test_end_to_end_outputs_are_well_formed():
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
            alpha_sf=1.0, alpha_er=2.0,
            num_warmup=20, num_samples=40, num_chains=1, epsilon=0.05,
        ))

        pip = pd.read_csv(out / "pip.csv", index_col=0)
        G = pd.read_csv(out / "G.csv", index_col=0)
        lfsr = pd.read_csv(out / "lfsr.csv", index_col=0)

    D = len(keep)
    for name, df in [("pip", pip), ("G", G), ("lfsr", lfsr)]:
        assert df.shape == (D, D), f"{name}.csv has shape {df.shape}"
        assert list(df.columns) == keep, f"{name}.csv columns not preserved"
        assert list(df.index) == keep, f"{name}.csv index not preserved"
        assert np.isfinite(df.values).all(), f"{name}.csv contains non-finite values"

    assert pip.values.min() >= 0.0 and pip.values.max() <= 1.0
    assert lfsr.values.min() >= 0.0 and lfsr.values.max() <= 0.5
    # the posterior should not be degenerate: some edges get real support
    assert pip.values.max() > 0.5, "no edge reached PIP > 0.5"


SLOW = {"test_end_to_end_outputs_are_well_formed"}


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
