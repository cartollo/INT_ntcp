from ntcp_common_imports import *
from ntcp_optimization_function import *

def predict_logistic(data, params):
    """
    Computes logistic model predictions by applying a sigmoid to the computed scores.

    Parameters
    ----------
    data : dict
        Container holding the dose-related inputs and risk factors required by the model.

    params : array-like
        Model parameters in the following order:
        n, a_b, b_0, b_eud, followed by optional risk-factor coefficients.

    Returns
    -------
    numpy.ndarray
        Predicted probabilities from the logistic model.
    """
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()

    n = params[0]
    a_b = params[1]
    b_0 = params[2]
    b_eud = params[3]
    b_risk = params[4:] if risk_matrix.shape[1] > 0 else []

    corrected_doses = EQD2(original_bins, n_fractions, a_b)       # EQD2-corrected dose bins
    EUDs = compute_EUDs(surfaces_matrix, corrected_doses, n)      # Compute EUD values
    scores = b_0 + b_eud * EUDs                                   # Base logistic score
    if risk_matrix is not None and len(b_risk) > 0:
        scores += np.dot(risk_matrix, np.array(b_risk))           # Add risk-factor contribution
    return sigmoid(scores)                                         # Convert scores to probabilities

def sigmoid(x): 
    # Implements the sigmoid function
    x = np.clip(x, -500, 500)
    return 1 / (1 + np.exp(-x))

def logistic_to_lyman(params):
    # Implements the tranformation from logistic parameters (input: params) to lyman parameters
    n = params[0]
    a_b = params[1]
    b0 = params[2]
    b_eud = params[3]
    b_risks = params[4:]

    td50 = - b0/b_eud
    m = -1/(b0 + sum(b_risks))
    dmfs = [1 + b/b0 for b in b_risks]

    return np.array([n, a_b, td50, m] + dmfs)