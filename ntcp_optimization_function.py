from ntcp_common_imports import *
from ntcp_global_variables import *

# region Optimization Functions
def fit_model(data, N_jk=1, p_jk=1, verbose=1):
    """
    Fits the logistic NTCP model, optionally using Jack-Knife resampling. The resampling is stratified based on data['outcome'], 
    keeping the same proportion of events and non-events.

    Parameters
    ----------
    data : dict
        Input data used for model fitting.

    N_jk : int, optional
        Number of Jack-Knife iterations. If >1, Jack-Knife resampling is applied.

    p_jk : float, optional
        Fraction of the dataset sampled in each Jack-Knife iteration.

    verbose : int, optional
        Controls the verbosity of printed output.

    Returns
    -------
    dict
        If N_jk > 1:
            DataFrame containing the N_jk final solutions after refinement (N_jk rows)
        Else:
            {
                "best_refined_results": DataFrame,
                "all_refined_results": DataFrame,
                "best_bf_results": DataFrame
            }  
            - best_refined_results: the final solution after refinement (1 row)  
            - all_refined_results: all initial guesses after refinement  
            - best_bf_results: selected initial guesses
    """
    if verbose > 0:
        if N_jk > 1:
            print(f"\nStart Jack-Knife Fitting ({N_jk:.0f} iterations, {p_jk*100:.0f}% sampling)...")
        else:
            print(f"\nStart Fitting (No Jack-Knife)...")

    risk_factors_names = list(data["risk_factors"].columns)
    param_names = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in risk_factors_names] + ["-log(L)"]

    if N_jk > 1: # Jack-Knife fitting: each core is assigned a batch of iterations, and each optimization is performed without parallelization
        best_refined_results = Parallel(n_jobs=-1, verbose=verbose)(
            delayed(_jackknife_single)(
                seed, data, p_jk, parallelize=False, return_only_best=True
            ) for seed in range(N_jk)
        ) # calls _jackknife_single, that performs a single Jack-Knife optimization

        best_refined_results = np.array(best_refined_results)
        df_results = pd.DataFrame(best_refined_results, columns=param_names).reset_index(drop=True)
        return df_results
    else: 
        # No Jack-Knife: N_jk = 1 and p_jk = 1 (no resampling), and the optimization is performed parallelizing both grid search and refining steps
        best_refined_results, all_refined_results, best_bf_results = _jackknife_single(
            1, data, 1, parallelize=True, return_only_best=False
        )
        best_refined_results = best_refined_results[np.newaxis, :]

        df_best_refined = pd.DataFrame(best_refined_results, columns=param_names).reset_index(drop=True)
        df_all_refined  = pd.DataFrame(all_refined_results, columns=param_names).reset_index(drop=True)
        df_best_bf      = pd.DataFrame(best_bf_results, columns=param_names).reset_index(drop=True)

        return {
            "best_refined_results": df_best_refined,
            "all_refined_results": df_all_refined,
            "best_bf_results": df_best_bf
        }

def _jackknife_single(seed, data, p_jk, parallelize, return_only_best=False):
    """
    Performs a single Jack-Knife iteration of the NTCP model fitting. The optimization procedure consists in a 
    broad grid search (using brute-force approach) with refinement of the most promising solutions to ensure global convergence.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility of the resampling.

    data : dict
        Input dataset containing at least the keys "outcome", "n_fractions", "DSH", "original_bins", and "risk_factors".

    p_jk : float
        Fraction of the dataset to sample for the Jack-Knife iteration, stratified by outcome.

    parallelize : bool
        Whether to run the optimization steps in parallel when possible.

    return_only_best : bool, optional
        If True, only the best refined solution is returned; 
        if False, returns also all refined solutions and selected brute-force solutions.

    Returns
    -------
    If return_only_best is True:
        numpy.ndarray
            The single best refined solution from this Jack-Knife iteration.
    Else:
        tuple of numpy.ndarray
            - best_refined_results: The best refined solution.
            - all_refined_results: All refined solutions from the initial guesses.
            - best_bf_results: The best solutions selected from the brute-force grid search.
    """

    # Extract data and perform Jack-Knife Sampling
    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()

    np.random.seed(seed)
    pos_idx = np.where(labels == 1)[0]
    neg_idx = np.where(labels == 0)[0]
    n_pos_samples = int(p_jk * len(pos_idx))
    n_neg_samples = int(p_jk * len(neg_idx))
    sampled_pos_idx = np.random.choice(pos_idx, size=n_pos_samples, replace=False)
    sampled_neg_idx = np.random.choice(neg_idx, size=n_neg_samples, replace=False)
    idx = np.concatenate([sampled_pos_idx, sampled_neg_idx])
    np.random.shuffle(idx)
    labels_jk = labels[idx]
    surfaces_matrix_jk = surfaces_matrix[idx]
    n_fractions_jk = n_fractions[idx]
    if risk_matrix is not None and len(risk_matrix) > 0 and len(idx) > 0:
        risk_matrix_jk = risk_matrix[idx]
    else:
        risk_matrix_jk = np.empty((0, risk_matrix.shape[1] if risk_matrix is not None else 0))

    data_jk = {
      "outcome": labels_jk,
      "n_fractions": n_fractions_jk,
      "DSH": surfaces_matrix_jk,
      "original_bins": original_bins,
      "risk_factors": pd.DataFrame(risk_matrix_jk)
      }

    # First step of the optimization: find best n_best_brute_force points from grid search
    all_results_bf = find_opt_brute_force(
        data_jk, parallelize
    )

    # Second step of the optimization: refine the most promising solutions
    # Initial guesses are selected as the top `n_initial_guesses` solutions with the lowest objective values, 
    # retaining only those whose first two parameters (n and alpha/beta ratio) are unique to ensure diversity.

    idxs = np.argpartition(all_results_bf[:, -1], INITIAL_GUESSES_BRUTE_FORCE-1)[:INITIAL_GUESSES_BRUTE_FORCE]
    subset = all_results_bf[idxs]
    initial_guesses = subset[:, :-1]

    initial_guesses_kept = []
    best_idx = np.argmin(subset[:, -1])
    best_guess = initial_guesses[best_idx]
    initial_guesses_kept.append(best_guess)

    for i, g in enumerate(initial_guesses):
        if i == best_idx:
            continue

        keep = True
        for kept in initial_guesses_kept:
            if not np.all(g[:2] != kept[:2]):
                keep = False
                break

        if keep:
            initial_guesses_kept.append(g)
    
    if not return_only_best:
        initial_guesses_kept_with_NLL = []

        for g in initial_guesses_kept:
            idx_in_subset = np.where((initial_guesses == g).all(axis=1))[0][0]
            nll_val = subset[idx_in_subset, -1]
            initial_guesses_kept_with_NLL.append(np.append(g, nll_val))

        initial_guesses_kept_with_NLL = np.array(initial_guesses_kept_with_NLL)

    all_refined_results = refine_solutions(
        data_jk, initial_guesses_kept, parallelize
    )
    best_refined_results = all_refined_results[np.argmin(all_refined_results[:, -1]), :]

    if return_only_best:
        return best_refined_results
    else:
        return best_refined_results, all_refined_results, initial_guesses_kept_with_NLL

def find_opt_brute_force(data, parallelize, huge_grid=False, verbose=0):
    """
    Performs a brute-force grid search to evaluate initial candidate solutions for the NTCP model.

    Parameters
    ----------
    data : dict
        Input dataset containing at least the keys "outcome" and "risk_factors".

    parallelize : bool
        Whether to evaluate the grid points in parallel.

    huge_grid : bool, optional
        If True, uses a significantly broader and coarser parameter grid. This extended grid is
        intended only for convergence checking, to test the robustness of the optimization
        process, and is not used during the standard optimization workflow.

    verbose : int, optional
        Verbosity level for diagnostic output.

    Returns
    -------
    numpy.ndarray
        Array containing the evaluated solutions for all points in the grid. Each row corresponds
        to a parameter set and its associated log-likelihood value.
    """
    # Define search grid for initial points
    p = np.mean(data['outcome'])
    n_risk = data['risk_factors'].shape[1]

    if huge_grid:
        grid_dict = GRID_CHECK_CONVERGENCE(p, n_risk) # Huge grid, used just for checking convergence
    else:
        grid_dict = GRID_OPTIMIZATION_BF(p, n_risk) # Usedd for fitting

    # Evaluate the points efficiently using a custom function
    results = evaluate_loglikelihood_grid(data, grid_dict, parallelize, verbose = verbose)
    return results

def evaluate_loglikelihood_grid(data, grid_dict, parallelize=True, threshold=np.inf, verbose = 0):
    """
    Evaluates the log-likelihood of the logistic NTCP model over a grid of parameter values.

    Parameters
    ----------
    data : dict
        Input dataset containing at least the keys "outcome", "n_fractions", "DSH", "original_bins", and "risk_factors".

    grid_dict : dict
        Dictionary defining the parameter grids. Expected keys include:
        "n", "a_b", "beta_0", "beta_EUD", and "beta_risk".

    parallelize : bool, optional
        If True, evaluations over the a_b grid are done in parallel.

    threshold : float, optional
        Maximum allowed NLL value; points exceeding this threshold are not returned. 
        Useful when computing simultaneous confidence intervals to avoid extremely large outputs.

    Returns
    -------
    numpy.ndarray
        Array of evaluated points. Each row corresponds to a parameter set along with its 
        associated NLL value. The number of columns depends on the number of parameters.
    """
    if parallelize:
        results = Parallel(n_jobs=-1, verbose = verbose)(
        delayed(evaluate_grid_ab)(             # Parallelize over alpha/beta values
            a_b, grid_dict, data, threshold)
            for a_b in grid_dict.get("a_b")
            )

        results = [r for r in results if len(r) > 0]
        if len(results) > 0:
            results = np.vstack(results)
        else:
            results = np.empty((0, len(grid_dict.get("n")) + len(grid_dict.get("beta_risk")) + 2))
    else:
        labels = data['outcome']
        n_fractions = data['n_fractions']
        surfaces_matrix = data['DSH']
        original_bins = data['original_bins']
        risk_matrix = data['risk_factors']

        n_values = grid_dict.get("n")
        a_b_values = grid_dict.get("a_b")
        b_0_values = grid_dict.get("beta_0")
        b_eud_values = grid_dict.get("beta_EUD")
        b_risk_values = grid_dict.get("beta_risk")

        # Evaluate the grid hierarchically, since operations involving alpha/beta and n are more expensive 
        # than those involving the beta coefficients. First, compute the EQD2 correction (which depends only on alpha/beta), 
        # then compute the EUDs (which depend only on n). After that, evaluate all combinations of the beta coefficients 
        # using the fixed EUD values.

        results = []
        for a_b in a_b_values:
            corrected_bins = EQD2(original_bins, n_fractions, a_b)           
            for n in n_values:
                EUDs = compute_EUDs(surfaces_matrix, corrected_bins, n)
                for b_0 in b_0_values:
                    for b_eud in b_eud_values:
                        if risk_matrix.shape[1] == 0:
                            scores = b_0 + b_eud * EUDs
                            fun_value = np.mean(np.logaddexp(0, scores) - labels * scores)
                            if fun_value <= threshold:
                                results.append([n, a_b, b_0, b_eud, fun_value])
                        else:
                            for b_risk_comb in product(*b_risk_values):
                                scores = b_0 + b_eud * EUDs + np.dot(risk_matrix, np.array(b_risk_comb))
                                fun_value = np.mean(np.logaddexp(0, scores) - labels * scores)
                                if fun_value <= threshold:
                                    results.append([n, a_b, b_0, b_eud] + list(b_risk_comb) + [fun_value])

    return np.array(results)

def evaluate_grid_ab(a_b, grid_dict, data, threshold):
    """
    Evaluates the log-likelihood of the logistic NTCP model for a fixed alpha/beta value (a_b) 
    over a grid of other parameters. Useful only when the grid search is parallelized.

    Parameters
    ----------
    a_b : float
        Fixed alpha/beta value for which the evaluation is performed.

    grid_dict : dict
        Dictionary defining the parameter grids. Expected keys include "n", "beta_0", "beta_EUD", and "beta_risk".

    data : dict
        Input dataset containing at least the keys "outcome", "n_fractions", "DSH", "original_bins", and "risk_factors".

    threshold : float
        Maximum allowed log-likelihood value; points exceeding this threshold are not returned.

    Returns
    -------
    list of list
        Each inner list contains a combination of parameter values along with its mean log-likelihood value.
        The structure of each list is [n, a_b, beta_0, beta_EUD, (...beta_risk...), fun_value].
    """
    n_values = grid_dict.get("n")
    b_0_values = grid_dict.get("beta_0")
    b_eud_values = grid_dict.get("beta_EUD")
    b_risk_values = grid_dict.get("beta_risk")

    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()

    corrected_bins = EQD2(original_bins, n_fractions, a_b)
    temp_results = []
    for n in n_values:
        EUDs = compute_EUDs(surfaces_matrix, corrected_bins, n)
        for b_0 in b_0_values:
            for b_eud in b_eud_values:
                if risk_matrix.shape[1] == 0:
                    scores = b_0 + b_eud * EUDs
                    loss = np.logaddexp(0, scores) - labels * scores
                    fun_value = np.mean(loss)
                    if fun_value <= threshold:
                        temp_results.append([n, a_b, b_0, b_eud, fun_value])
                else:
                    for b_risk_comb in product(*b_risk_values):
                        scores = b_0 + b_eud * EUDs + np.dot(risk_matrix, np.array(b_risk_comb))
                        loss = np.logaddexp(0, scores) - labels * scores
                        fun_value = np.mean(loss)
                        if fun_value <= threshold:
                            temp_results.append([n, a_b, b_0, b_eud] + list(b_risk_comb) + [fun_value])
    return temp_results

def refine_solutions(data, initial_guesses, parallelize, verbose = 0):
    """
    Refines a set of initial candidate solutions for the NTCP model using an independent optimization 
    for each starting point.

    Parameters
    ----------
    data : dict
        Input dataset used for the optimization.

    initial_guesses : array-like
        Array of initial parameter sets to refine.

    parallelize : bool
        If True, optimizations for each initial guess are run in parallel.

    Returns
    -------
    numpy.ndarray
        Array containing the refined solutions for all initial guesses. Each row corresponds 
        to a single optimized parameter set.
    """
    # Run and independent optimization for each starting point
    if parallelize:
        results = Parallel(n_jobs=-1, verbose = verbose)(
            delayed(refine_single_point)(params, data)
            for params in initial_guesses
        )
    else:
        results = []
        for params in initial_guesses:
            res = refine_single_point(params, data)
            results.append(res)

    results = np.array(results)
    return results

def refine_single_point(params, data):
    """
    Performs a local optimization of the NTCP model parameters starting from a single set of initial parameters.

    Parameters
    ----------
    params : array-like
        Initial parameter values for the optimization. Expected order: [n, a_b, beta_0, beta_EUD, ...beta_risk].

    data : dict
        Input dataset used for the optimization.

    Returns
    -------
    tuple
        Optimized parameter values followed by the corresponding negative log-likelihood value.
        The structure is (n, a_b, beta_0, beta_EUD, ...beta_risk, fun_value).
    """

    res = minimize(
        neg_loglikelihood,
        x0=list(params),
        args=data,
        method=OPTIMIZATION_ALGORITHM,
        bounds= OPTIMIZATION_BOUNDS(len(params) - 4),
        options = REFINEMENT_OPTIMIZATION_OPTIONS
    )

    return (*res.x, res.fun)

def neg_loglikelihood(params, data):
    """
    Computes the normalized negative log-likelihood (NLL) of the logistic NTCP model for a given set of parameters.

    Parameters
    ----------
    params : array-like
        Model parameters in the order [n, a_b, beta_0, beta_EUD, ...beta_risk].

    data : dict
        Input dataset containing at least the keys "outcome", "n_fractions", "DSH", "original_bins", and "risk_factors".

    Returns
    -------
    float
        normalized negative log-likelihood
    """
    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()

    n = params[0]
    a_b = params[1]
    b_0 = params[2]
    b_eud = params[3]
    b_risk = params[4:] if risk_matrix.shape[1] > 0 else []

    corrected_doses = EQD2(original_bins, n_fractions, a_b)       # Correct the dose bins using EQD2 correction
    EUDs = compute_EUDs(surfaces_matrix, corrected_doses, n)      # Compute the EUDs
    scores = b_0 + b_eud * EUDs
    if risk_matrix is not None and len(b_risk) > 0:
        scores += np.dot(risk_matrix, np.array(b_risk))           # Compute the scores

    neg_ll = np.logaddexp(0, scores) - labels * scores
    return np.mean(neg_ll)                                        # Return the normalized negative log-likelihood (divided by the number of patients)

def compute_EUDs(DSH, corrected_doses, n):
    """
    Computes the Equivalent Uniform Dose (EUD) for each patient based on the dose-surface histogram.

    Parameters
    ----------
    DSH : array-like
        Matrix where each row represents a patient's dose-surface histogram (fractions of organ surface).

    corrected_doses : array-like
        Corrected dose values (EQD2) corresponding to the bins of the dose-surface histogram.

    n : float
        Tissue-specific parameter used in the EUD calculation.

    Returns
    -------
    numpy.ndarray
        Array of EUD values, one per patient.
    """
    doses_1_n = np.power(corrected_doses, 1 / n)
    S = np.sum(DSH * doses_1_n, axis=1)
    return np.power(S, n)

def EQD2(original_bins, n_fractions, a_b):
    """
    Computes the EQD2 (Equivalent Dose in 2 Gy fractions) correction for a set of dose bins.

    Parameters
    ----------
    original_bins : array-like
        Original dose bins for each patient.

    n_fractions : array-like
        Number of fractions in which the dose was delivered to each patient.

    a_b : float
        Tissue-specific alpha/beta ratio.

    Returns
    -------
    numpy.ndarray
        Corrected doses for each patient and each bin, after applying the EQD2 formula.
    """
    corrected_doses = np.zeros((len(n_fractions), len(original_bins)))
    for i in range(len(n_fractions)):
        corrected_doses[i, :] = original_bins * (a_b + original_bins / n_fractions[i]) / (a_b + 2)
    return corrected_doses