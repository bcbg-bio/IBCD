import numpy as np

def calc_metrics(X: np.ndarray, X_true: np.ndarray, eps: float = 1e-8) -> dict:
    """
    Compute various evaluation metrics between a predicted matrix X and ground truth X_true.

    Parameters:
    - X (np.ndarray): Predicted square matrix of shape (D, D).
    - X_true (np.ndarray): Ground truth square matrix of same shape.
    - eps (float): Threshold below which values are treated as zero to avoid noise.

    Returns:
    dict with keys:
      - precision, recall, F1, rmse, mae, shd, weight_acc, TP, FP, TN, FN
    """
    # Ensure inputs are numpy arrays
    X = np.array(X, dtype=float)
    X_true = np.array(X_true, dtype=float)

    # Dimension
    D = X.shape[1]

    # Remove diagonal entries
    mask = ~np.eye(D, dtype=bool)
    x = X[mask]
    xt = X_true[mask]

    # Zero out near-zero values
    x[np.abs(x) < eps] = 0
    xt[np.abs(xt) < eps] = 0

    # Error metrics
    rmse = np.sqrt(np.mean((x - xt) ** 2))
    mae = np.mean(np.abs(x - xt))

    # Signs of entries
    sign_x = np.sign(x)
    sign_xt = np.sign(xt)

    # True Negatives: both are zero
    TN = int(np.sum((sign_x == 0) & (sign_xt == 0)))
    # Total matches minus TN = True Positives
    total_matches = int(np.sum(sign_x == sign_xt))
    TP = total_matches - TN

    # False Negatives: predicted zero, true non-zero
    FN = int(np.sum((sign_x == 0) & (sign_xt != 0)))
    # False Positives: total mismatches minus FN
    total_mismatches = int(np.sum(sign_x != sign_xt))
    FP = total_mismatches - FN

    # Structural Hamming Distance
    shd = total_mismatches

    # Class counts in ground truth
    N_pos = int(np.sum(xt > 0))
    N_neg = int(np.sum(xt < 0))
    N_zero = int(np.sum(xt == 0))

    # Weighted accuracy
    weights = np.zeros_like(xt, dtype=float)
    if N_pos > 0:
        weights[sign_xt > 0] = 1 / N_pos
    if N_neg > 0:
        weights[sign_xt < 0] = 1 / N_neg
    if N_zero > 0:
        weights[sign_xt == 0] = 1 / N_zero
    weight_acc = np.sum((sign_x == sign_xt) * weights) / np.sum(weights)

    # Precision and recall
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

    # F1 score
    F1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        'precision': precision,
        'recall': recall,
        'F1': F1,
        'rmse': rmse,
        'mae': mae,
        'shd': shd,
        'weight_acc': weight_acc,
        'TP': TP,
        'FP': FP,
        'TN': TN,
        'FN': FN
    }
