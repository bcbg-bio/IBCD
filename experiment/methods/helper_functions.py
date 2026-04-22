import conditional_independence
from causaldag import igsp
import numpy as np
import pandas as pd
import scipy.linalg as slin
import scipy.optimize as sopt
import sys

def loss_notears(W, X, loss_type='l2'):
    """Evaluate value and gradient of loss."""
    M = X @ W
    if loss_type == 'l2':
        R = X - M
        loss = 0.5 / X.shape[0] * (R ** 2).sum()
        G_loss = - 1.0 / X.shape[0] * X.T @ R
    elif loss_type == 'logistic':
        loss = 1.0 / X.shape[0] * (np.logaddexp(0, M) - X * M).sum()
        G_loss = 1.0 / X.shape[0] * X.T @ (sigmoid(M) - X)
    elif loss_type == 'poisson':
        S = np.exp(M)
        loss = 1.0 / X.shape[0] * (S - X * M).sum()
        G_loss = 1.0 / X.shape[0] * X.T @ (S - X)
    else:
        raise ValueError('unknown loss type')
#         print(loss, G_loss)
    return loss, G_loss	


def notears_linear(X, lambda1, loss_type='l2', max_iter=100, h_tol=1e-8, rho_max=1e+16, w_threshold=0.3):
    """Solve min_W L(W; X) + lambda1 ‖W‖_1 s.t. h(W) = 0 using augmented Lagrangian.
    Args:
        X (np.ndarray): [n, d] sample matrix
        lambda1 (float): l1 penalty parameter
        loss_type (str): l2, logistic, poisson
        max_iter (int): max num of dual ascent steps
        h_tol (float): exit if |h(w_est)| <= htol
        rho_max (float): exit if rho >= rho_max
        w_threshold (float): drop edge if |weight| < threshold
    Returns:
        W_est (np.ndarray): [d, d] estimated DAG
    """
    def _loss(W):
        """Evaluate value and gradient of loss."""
        M = X @ W
        if loss_type == 'l2':
            R = X - M
            loss = 0.5 / X.shape[0] * (R ** 2).sum()
            G_loss = - 1.0 / X.shape[0] * X.T @ R
        elif loss_type == 'logistic':
            loss = 1.0 / X.shape[0] * (np.logaddexp(0, M) - X * M).sum()
            G_loss = 1.0 / X.shape[0] * X.T @ (sigmoid(M) - X)
        elif loss_type == 'poisson':
            S = np.exp(M)
            loss = 1.0 / X.shape[0] * (S - X * M).sum()
            G_loss = 1.0 / X.shape[0] * X.T @ (S - X)
        else:
            raise ValueError('unknown loss type')
            
#         print(loss, G_loss)
        return loss, G_loss

    def _h(W):
        """Evaluate value and gradient of acyclicity constraint."""
        E = slin.expm(W * W)  # (Zheng et al. 2018)
        h = np.trace(E) - d
        #     # A different formulation, slightly faster at the cost of numerical stability
        #     M = np.eye(d) + W * W / d  # (Yu et al. 2019)
        #     E = np.linalg.matrix_power(M, d - 1)
        #     h = (E.T * M).sum() - d
        G_h = E.T * W * 2
        return h, G_h

    def _adj(w):
        """Convert doubled variables ([2 d^2] array) back to original variables ([d, d] matrix)."""
        return (w[:d * d] - w[d * d:]).reshape([d, d])

    def _func(w):
        """Evaluate value and gradient of augmented Lagrangian for doubled variables ([2 d^2] array)."""
        W = _adj(w)
        loss, G_loss = _loss(W)
        h, G_h = _h(W)
        obj = loss + 0.5 * rho * h * h + alpha * h + lambda1 * w.sum()
        G_smooth = G_loss + (rho * h + alpha) * G_h
        g_obj = np.concatenate((G_smooth + lambda1, - G_smooth + lambda1), axis=None)
        return obj, g_obj

    n, d = X.shape
    w_est, rho, alpha, h = np.zeros(2 * d * d), 1.0, 0.0, np.inf  # double w_est into (w_pos, w_neg)
    bnds = [(0, 0) if i == j else (0, None) for _ in range(2) for i in range(d) for j in range(d)]
    if loss_type == 'l2':
        X = X - np.mean(X, axis=0, keepdims=True)
    for i in range(max_iter):
        w_new, h_new = None, None
        while rho < rho_max:
            sol = sopt.minimize(_func, w_est, method='L-BFGS-B', jac=True, bounds=bnds)
            w_new = sol.x
            h_new, _ = _h(_adj(w_new))
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        # print(i, sol.fun, h)
        if h <= h_tol or rho >= rho_max:
            break
    W_est = _adj(w_est)
    W_est[np.abs(W_est) < w_threshold] = 0
    return W_est


def make_dotears_data2(X, targets):
  data_dotears = {}
  for key in np.unique(targets):
    if key != "control":
      data_dotears[int(key.replace("V", ""))-1] = X[targets==key,]
    else:
      data_dotears["obs"] = X[targets==key,]
  return data_dotears


def make_dotears_data(X, targets, colnames):
    data_dotears = {}
    colnames = np.array(colnames)

    for key in np.unique(targets):
        if key == "control":
            data_dotears["obs"] = X[targets == key, :]
        else:
            col_idx = np.where(colnames == key)[0][0]
            data_dotears[col_idx] = X[targets == key, :]

    return data_dotears

def run_igsp(X_dotears, alpha=0.001, alpha_inv=0.001):
    obs_data = X_dotears['obs']
    p = obs_data.shape[1]
    nodes = list(range(p))

    # iv_samples_list is a list of n x p ndarrays
    inv_samples_from_data = [v for k, v in X_dotears.items() if k != 'obs']

    # setting_list is list of dicts
    # each dict is key 'intervention' to a list of nodes
    settings_from_data = [dict(interventions=[k]) for k, v in X_dotears.items() if k != 'obs']

    obs_suffstat = conditional_independence.partial_correlation_suffstat(obs_data)
    invariance_suffstat = conditional_independence.gauss_invariance_suffstat(
        obs_data, inv_samples_from_data)

    ci_tester = conditional_independence.MemoizedCI_Tester(
        conditional_independence.partial_correlation_test, obs_suffstat, alpha=alpha)
    invariance_tester = conditional_independence.MemoizedInvarianceTester(
        conditional_independence.gauss_invariance_test, invariance_suffstat, alpha=alpha_inv)

    est_dag = igsp(settings_from_data, nodes, ci_tester, invariance_tester)
    W_igsp = est_dag.to_amat()[0]
    return(W_igsp)

