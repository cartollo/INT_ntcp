from ntcp_gian_common_imports import *
from ntcp_gian_optimization_function import *
from ntcp_gian_functions import *

def summarize_fit(results_jk, coverage=0.68, plot_histograms=False):
    """
    Summarizes the fit of a logistic NTCP model, optionally using Jack-Knife results.

    This function computes parameter estimates for both the logistic and Lyman models,
    including derived quantities TD50, m, and DMFs for binary risk factors. If multiple
    Jack-Knife iterations are present, it calculates medians and confidence intervals
    at the specified coverage. Optionally, it can plot histograms of the parameters
    across Jack-Knife iterations, highlighting median and confidence intervals.

    Parameters
    ----------
    results_jk : DataFrame
        DataFrame containing model parameter estimates. Columns should include:
        - "n", "a_b", "beta_0", "beta_EUD"
        - "beta_<risk factor>" for any binary risk factors
        - optionally "-log(L)" for negative log-likelihoods

    coverage : float, optional
        Confidence level for Jack-Knife estimates (default is 0.68).

    plot_histograms : bool, optional
        Whether to plot histograms of the parameters across Jack-Knife iterations.

    Returns
    -------
    df_out : DataFrame
        Summary table with parameter medians and confidence intervals (if Jack-Knife).

    parameters_logistic_model : DataFrame
        Estimated parameters for the logistic NTCP model: [n, a_b, beta_0, beta_EUD, beta_<risk factors>].

    parameters_lyman_model : DataFrame
        Estimated parameters for the Lyman model: [n, a_b, TD50, m, DMFs].
    """
    # Define the name expected for each parameter
    n_name = "n"
    a_b_name = "a_b"
    beta0_name = "beta_0"
    betaEUD_name = "beta_EUD"
    beta_prefix = "beta_"
    TD50_name = "TD50"
    m_name = "m"
    dmf_prefix = "dmf_"
    loglik_name = "-log(L)"

    # Infer risk factor names from columns in results_jk
    risk_factors_names = [col.replace(beta_prefix, "") for col in results_jk.columns 
                          if col.startswith(beta_prefix) and col not in [beta0_name, betaEUD_name]]
    
    param_base = [n_name, a_b_name, beta0_name, betaEUD_name]
    param_risk = [beta_prefix + rf for rf in risk_factors_names]
    all_params = param_base + param_risk
    logistic_cols = param_base + [beta_prefix + rf for rf in risk_factors_names]
    lyman_cols = [n_name, a_b_name, "TD50", "m"] + [dmf_prefix + rf for rf in risk_factors_names]

    single_row = len(results_jk) == 1
    df_params = results_jk[all_params].copy()

    if single_row: # No Jack-Knife: no uncertainty estimates are computed
        rows_df = []
        for name in all_params:
            rows_df.append({"Parameter": name, "Estimate": df_params[name].iloc[0]})

        # TD50
        TD50_val = -df_params[beta0_name].iloc[0] / df_params[betaEUD_name].iloc[0]

        # m
        if param_risk:
            m_val = -1 / (df_params[beta0_name].iloc[0] + df_params[param_risk].sum(axis=1).iloc[0])
        else:
            m_val = -1 / df_params[beta0_name].iloc[0]

        rows_df.append({"Parameter": TD50_name, "Estimate": TD50_val})
        rows_df.append({"Parameter": m_name, "Estimate": m_val})

        # dmfs
        for rf in risk_factors_names:
            val = 1 + df_params[beta_prefix + rf].iloc[0] / df_params[beta0_name].iloc[0]
            rows_df.append({"Parameter": dmf_prefix + rf, "Estimate": val})

        df_out = pd.DataFrame(rows_df)
        print(f"{'Parameter':<30} {'Estimate':>20}")
        for r in rows_df:
            print(f"{r['Parameter']:<30} {r['Estimate']:>20.6f}")

        parameters_logistic_model = pd.DataFrame({name: df_out.loc[df_out["Parameter"] == name, "Estimate"].values[0] 
                                                for name in logistic_cols}, index=[0])

        parameters_lyman_model = pd.DataFrame({name: df_out.loc[df_out["Parameter"] == name, "Estimate"].values[0] 
                                            for name in lyman_cols}, index=[0])

        return df_out, parameters_logistic_model, parameters_lyman_model

    # Jack-Knife: medians, confidence intervals and histograms are computed
    alpha = 1 - coverage

    # Compute the population of TD50 among the Jack-Knife iterations
    TD50_vals = -df_params[beta0_name] / df_params[betaEUD_name]

    # Same for m and dmfs
    if param_risk:
        m_vals = -1 / (df_params[beta0_name] + df_params[param_risk].sum(axis=1))
    else:
        m_vals = -1 / df_params[beta0_name]

    dmf_vals = {rf: 1 + df_params[beta_prefix + rf] / df_params[beta0_name] 
                for rf in risk_factors_names}

    # Compute confidence intervals using quantiles [q_{alpha/2}, median, q_{1-alpha/2}]
    def get_ci(vals):
        med = np.median(vals)
        lo = np.quantile(vals, alpha/2)
        hi = np.quantile(vals, 1 - alpha/2)
        return med, lo, hi

    # Organize results and print
    rows_print = []
    rows_df = []

    def add_param(name, vals):
        med, lo, hi = get_ci(vals)
        rows_print.append({"Parameter": name, "Estimate": f"{med:.3f}", "CI": f"[{lo:.3f}, {hi:.3f}]"})
        rows_df.append({"Parameter": name, "Estimate": med,
                        f"Low(CI JK {int(coverage*100)}%)": lo,
                        f"High(CI JK {int(coverage*100)}%)": hi})

    for name in all_params:
        add_param(name, df_params[name])

    add_param("TD50", TD50_vals)
    add_param("m", m_vals)

    for rf in risk_factors_names:
        add_param(dmf_prefix + rf, dmf_vals[rf])

    if loglik_name in results_jk.columns:
        add_param(loglik_name, results_jk[loglik_name].to_numpy())

    print(f"{'Parameter':<30} {'Estimate (Median)':>20} {f'CI (JK {int(coverage*100)}%)':>20}")
    for r in rows_print:
        print(f"{r['Parameter']:<30} {r['Estimate']:>20} {r['CI']:>20}")

    df_out = pd.DataFrame(rows_df)

    parameters_logistic_model = pd.DataFrame({name: df_out.loc[df_out["Parameter"] == name, "Estimate"].values[0] 
                                            for name in logistic_cols}, index=[0])

    parameters_lyman_model = pd.DataFrame({name: df_out.loc[df_out["Parameter"] == name, "Estimate"].values[0] 
                                        for name in lyman_cols}, index=[0])

    # Plot histograms for all quantities
    if plot_histograms:
        cols_to_plot = all_params + ["TD50", "m"] + [dmf_prefix + rf for rf in risk_factors_names]
        if loglik_name in results_jk.columns:
            cols_to_plot.append(loglik_name)

        n_cols = 4
        n_rows = (len(cols_to_plot) + n_cols - 1) // n_cols
        plt.figure(figsize=(5 * n_cols, 4 * n_rows))

        for i, col in enumerate(cols_to_plot):
            if col in df_params.columns:
                vals = df_params[col]
            elif col == "TD50":
                vals = TD50_vals
            elif col == "m":
                vals = m_vals
            elif col.startswith(dmf_prefix):
                key = col.replace(dmf_prefix, "")
                vals = dmf_vals.get(key, None)
            elif col == loglik_name:
                vals = results_jk[loglik_name].to_numpy()
            else:
                vals = None

            if vals is None:
                continue

            n_bins = min(max(len(vals)//10, 10), 50)
            v_min, v_max = vals.min(), vals.max()
            bin_edges = np.linspace(v_min, v_max, n_bins + 1)

            plt.subplot(n_rows, n_cols, i+1)
            plt.hist(vals, bins=bin_edges, color="steelblue", alpha=0.7)
            med, lo, hi = get_ci(vals)
            plt.axvline(med, color="red", linestyle="--")
            plt.axvline(lo, color="green", linestyle="--")
            plt.axvline(hi, color="green", linestyle="--")
            plt.title(col)
            if i == 0:
                plt.legend(["Median", "CI low", "CI high"])

        plt.tight_layout()
        plt.show()

    return df_out, parameters_logistic_model, parameters_lyman_model

def likelihood_ratio_test(model_H1, data, model_H0=None):
    """
    Performs a likelihood ratio test (LRT) comparing two nested logistic NTCP models.

    This function tests whether adding additional factors (model_H1) to a null model
    (model_H0) significantly improves the fit. If a null model `model_H0` 
    is provided, the LRT compares it to `model_H1`. If `model_H0` is not provided,
    the test compares the full model against an intercept-only model.

    The function prints the following:
    - Parameters tested (automatically inferred from the structures of model_H1 and model_H0)
    - Specification of H0 and H1
    - Log-likelihoods of H0 and H1
    - Likelihood ratio statistic, p-value, and significance stars

    Parameters
    ----------
    model_H1 : DataFrame
        Fitted model results for the alternative hypothesis. Must contain columns:
        - "Parameter": name of each parameter ("n", "a_b", "beta_0", "beta_EUD", "beta_<risk factor>")
        - "Estimate": fitted value for the parameter

    data : dict
        Usual data object

    model_H0 : DataFrame, optional
        Fitted model results for the null hypothesis. If not provided, the null model is
        intercept-only (no EUD or risk factors).

    Returns
    -------
    dict
        Dictionary containing:
        - "LR_stat": likelihood ratio statistic
        - "p_value": p-value of the test

    Notes
    -----
    The LR statistic is computed as 2 * (NLL_H0 - NLL_H1) * N, where NNL denotes the normalized
    negative log-likelihood of the model. Significance is indicated with stars:
    '***' p<0.001, '**' p<0.01, '*' p<0.05, '.' p<0.1, 'ns' otherwise.
    """
    labels = data["outcome"]
    risk_matrix = data["risk_factors"]

    def extract_params(model, cols):
        p = [
            model["n"].values[0],
            model["a_b"].values[0],
            model["beta_0"].values[0],
            model["beta_EUD"].values[0]
        ]
        for c in cols:
            name = "beta_" + c
            if name in model.columns:
                p.append(model[name].values[0])
        return np.array(p)

    if model_H0 is not None:
        beta_risk_H1 = [c.replace("beta_", "") for c in model_H1.columns if c.startswith("beta_") and c not in ["beta_0", "beta_EUD"]]
        beta_risk_H0 = [c.replace("beta_", "") for c in model_H0.columns if c.startswith("beta_") and c not in ["beta_0", "beta_EUD"]]

        tested_factors = [p for p in beta_risk_H1 if p not in beta_risk_H0]
        kept_factors = [p for p in beta_risk_H0]

        data_H0 = deepcopy(data)
        data_H0["risk_factors"] = data_H0["risk_factors"][kept_factors]

        params_H1 = extract_params(model_H1, list(risk_matrix.columns))
        params_H0 = extract_params(model_H0, kept_factors)

        nll0 = neg_loglikelihood(params_H0, data_H0)
        nll1 = neg_loglikelihood(params_H1, data)

        LR_stat = 2 * (nll0 - nll1) * len(labels)
        df_test = len(tested_factors)
        p = 1 - chi2.cdf(LR_stat, df=df_test)
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '.' if p < 0.1 else 'ns'

        print(f"Testing {df_test} parameter" + ("s" if df_test>1 else "") + f": {', '.join(tested_factors)}")
        print(f"H0: EUD" + "".join([f" + {rf}" for rf in kept_factors]) if len(kept_factors)>0 else "H0: base_model(EUD)")
        print(f"H1: EUD" + "".join([f" + {p}" for p in beta_risk_H1]))
        print(f"NLL(H0)={nll0:.6f}, NLL(H1)={nll1:.6f}, LR={LR_stat:.4f}, p={p:.4e} {stars}")

        return {"LR_stat": LR_stat, "p_value": p}
    else:
        beta_risk_H1 = [c.replace("beta_", "") for c in model_H1.columns if c.startswith("beta_") and c not in ["beta_0", "beta_EUD"]]

        params_H1 = extract_params(model_H1, beta_risk_H1)

        p_mean = np.mean(labels)
        nll0 = -(p_mean*np.log(p_mean) + (1-p_mean)*np.log(1-p_mean))
        nll1 = neg_loglikelihood(params_H1, data)

        LR_stat = 2 * (nll0 - nll1) * len(labels)
        df_test = len(params_H1) - 3
        p = 1 - chi2.cdf(LR_stat, df=df_test)
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '.' if p < 0.1 else 'ns'

        print(f"Testing {df_test} parameter" + ("s" if df_test>1 else "") + ": EUD" + (" + " + " + ".join(beta_risk_H1) if len(beta_risk_H1)>0 else ""))
        print("H0: intercept-only")
        print("H1: EUD" + "".join([f" + {rf}" for rf in beta_risk_H1]))
        print(f"NLL(H0)={nll0:.6f}, NLL(H1)={nll1:.6f}, LR={LR_stat:.4f}, p={p:.4f} {stars}")

        return {"LR_stat": LR_stat, "p_value": p}

def simultaneous_CI(data, parameters_logistic_model, fit_parameters=None, coverage=0.68, n_points=50):
    '''
    Calculates simultaneous confidence intervals for logistic NTCP model parameters using a grid-based likelihood evaluation.

    Parameters
    ----------
    data : dict
        Usual data object.
        
    parameters_logistic_model : array-like
        DataFrame of fitted logistic model parameters
    fit_parameters : list of str or None, optional
        List of parameter names for which confidence intervals should be calculated. Parameters not listed
        are held fixed at their fitted values. Default is None, meaning all parameters are considered.
        Example: fit_parameters = ["n", "a_b"] computes the 2-dimensional simultaneous CI for n and alpha/beta only.

    coverage : float, optional
        Desired confidence level for the intervals (default 0.68, corresponding to ~1 standard deviation).

    n_points : int, optional
        Number of points per parameter to use in the evaluated grid (default 20).

    Returns
    -------
    pandas.DataFrame
        A DataFrame with the following columns:
        - 'Parameter': name of the parameter.
        - 'Estimate': fitted value of the parameter.
        - 'Low(CI Sim coverage%)': lower bound of the simultaneous confidence interval at the specified coverage.
        - 'High(CI Sim coverage%)': upper bound of the simultaneous confidence interval at the specified coverage.

    Notes
    -----
    - Confidence intervals are computed using a likelihood threshold approach:
    a grid of parameter values is evaluated, and points with negative log-likelihood
    within a threshold from the fitted value are considered valid.
   - The threshold (`target_value`) is computed as the negative log-likelihood at the fitted
    parameters plus 0.5 * chi2.ppf(coverage, df=number_of_parameters) / N, where N is the number
    of patients.
    - The function automatically expands the grid if the valid region touches the initial edges,
    ensuring proper coverage of the confidence region.

    '''
    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in list(data["risk_factors"].columns)]
    params = parameters_logistic_model[param_order].to_numpy()[0]
    N = len(data['outcome'])
    df = len(param_order) if fit_parameters is None else len(fit_parameters)
    target_value = neg_loglikelihood(params, data) + 0.5 * chi2.ppf(coverage, df=df) / N

    param_config = {
        "n": {"min": 0.01, "max": 1.0, "step": 0.01, "window": 0.1},
        "a_b": {"min": 0.5, "max": 10.0, "step": 0.1, "window": 3.0},
        "beta_0": {"min": -np.inf, "max": np.inf, "step": 0.1, "window": 4.0},
        "beta_EUD": {"min": 0.0, "max": np.inf, "step": 0.01, "window": 0.1},
    }
    for rf in data["risk_factors"].columns:
        param_config[f"beta_{rf}"] = {"min": -np.inf, "max": np.inf, "step": 0.1, "window": 1.0}

    risk_idx = {rf: i for i, rf in enumerate(data["risk_factors"].columns)}

    grid_dict = {}
    beta_risk = []

    for i, name in enumerate(param_order):
        median = params[i]
        if fit_parameters is not None and name not in fit_parameters:
            if name.startswith("beta_") and name not in ["beta_0", "beta_EUD"]:
                beta_risk.append(np.array([median]))
            else:
                grid_dict[name] = np.array([median])
        else:
            if name.startswith("beta_") and name not in ["beta_0", "beta_EUD"]:
                lower = median - param_config[name]["window"]
                upper = median + param_config[name]["window"]
                beta_risk.append(np.linspace(lower, upper, n_points))
            else:
                lower = max(param_config[name]["min"], median - param_config[name]["window"])
                upper = min(param_config[name]["max"], median + param_config[name]["window"])
                grid_dict[name] = np.linspace(lower, upper, n_points)

    grid_dict["beta_risk"] = beta_risk

    valid_points = evaluate_loglikelihood_grid(data, grid_dict, parallelize=True, threshold=target_value)
    valid_points = np.array(valid_points)

    expanded = True
    while expanded:
        expanded = False
        for i, name in enumerate(param_order):
            if fit_parameters is not None and name not in fit_parameters:
                continue

            step = param_config[name]["step"]
            new_points = []
            n_points_expansion = int(step/param_config[name]["window"] * n_points)

            if name.startswith("beta_") and name not in ["beta_0", "beta_EUD"]:
                rf = name.split("_", 1)[1]
                idx = risk_idx[rf]
                grid_vals = grid_dict["beta_risk"][idx]
                vals = valid_points[:, i]
                q_low = np.min(vals)
                q_high = np.max(vals)

                if q_low != param_config[name]["min"] and q_low == grid_vals[0]:
                    bound_min = max(param_config[name]["min"], grid_vals[0] - step)
                    new_points = np.linspace(bound_min, grid_vals[0], n_points_expansion, endpoint=False)
                    grid_dict["beta_risk"][idx] = np.concatenate((new_points, grid_vals))
                    expanded = True

                elif q_high != param_config[name]["max"] and q_high == grid_vals[-1]:
                    bound_max = min(param_config[name]["max"], grid_vals[-1] + step)
                    new_points = np.linspace(grid_vals[-1], bound_max, n_points_expansion + 1, endpoint=False)[1:]
                    grid_dict["beta_risk"][idx] = np.concatenate((grid_vals, new_points))
                    expanded = True

                if len(new_points) > 0:
                    new_grid_dict = grid_dict.copy()
                    new_grid_dict["beta_risk"] = grid_dict["beta_risk"][:]
                    new_grid_dict["beta_risk"][idx] = new_points
                    new_valid = evaluate_loglikelihood_grid(data, new_grid_dict, parallelize=True, threshold=target_value)
                    if len(new_valid) > 0:
                        valid_points = np.vstack([valid_points, np.array(new_valid)])

            else:
                grid_vals = grid_dict[name]
                vals = valid_points[:, i]
                q_low = float(np.min(vals))
                q_high = float(np.max(vals))

                if q_low != param_config[name]["min"] and q_low == grid_vals[0]:
                    bound_min = max(param_config[name]["min"], grid_vals[0] - step)
                    new_points = np.linspace(bound_min, grid_vals[0], n_points_expansion, endpoint=False)
                    grid_dict[name] = np.concatenate((new_points, grid_vals))
                    expanded = True

                elif q_high != param_config[name]["max"] and q_high == grid_vals[-1]:
                    bound_max = min(param_config[name]["max"], grid_vals[-1] + step)
                    new_points = np.linspace(grid_vals[-1], bound_max, n_points_expansion + 1, endpoint=False)[1:]
                    grid_dict[name] = np.concatenate((grid_vals, new_points))
                    expanded = True

                if len(new_points) > 0:
                    new_grid_dict = grid_dict.copy()
                    new_grid_dict[name] = new_points
                    new_valid = evaluate_loglikelihood_grid(data, new_grid_dict, parallelize=True, threshold=target_value)
                    if len(new_valid) > 0:
                        valid_points = np.vstack([valid_points, np.array(new_valid)])

    results_list = []
    for i, name in enumerate(param_order):
        median = params[i]
        if fit_parameters is not None and name not in fit_parameters:
            results_list.append({
                "Parameter": name,
                "Estimate": float(median),
                f"Low(CI Sim {int(coverage*100)}%)": np.nan,
                f"High(CI Sim {int(coverage*100)}%)": np.nan
            })
        else:
            vals = valid_points[:, i]
            q_low = float(np.min(vals))
            q_high = float(np.max(vals))
            results_list.append({
                "Parameter": name,
                "Estimate": float(median),
                f"Low(CI Sim {int(coverage*100)}%)": q_low,
                f"High(CI Sim {int(coverage*100)}%)": q_high
            })

    return pd.DataFrame(results_list)

def predict(data, parameters_logistic_model, threshold=0.5):
    """
    Computes and prints standard classification metrics for a logistic model.

    Parameters
    ----------
    data : dict
        Usual data object

    parameters_logistic_model : DataFrame
        Parameters of the logistic model used to generate predictions

    threshold : float, optional
        Threshold to convert predicted probabilities into binary labels. Default is 0.5.

    Returns
    -------
    DataFrame
        Predicted probabilities for each observation.
    
    Notes
    -----
    Prints metrics independent of threshold (AUC, Average Precision) and metrics dependent on the 
    chosen threshold (Accuracy, Precision, Recall, F1 Score), along with the confusion matrix.
    """
    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in list(data["risk_factors"].columns)]
    parameters_logistic_model = parameters_logistic_model[param_order].to_numpy()[0]
    predictions = predict_logistic(data, parameters_logistic_model)
    labels = data["outcome"]
    pred_labels = (predictions >= threshold).astype(int)

    auc = roc_auc_score(labels, predictions)
    avg_prec = average_precision_score(labels, predictions)

    precision = precision_score(labels, pred_labels, zero_division=0)
    recall = recall_score(labels, pred_labels, zero_division=0)
    f1 = f1_score(labels, pred_labels, zero_division=0)
    accuracy = accuracy_score(labels, pred_labels)

    cm = confusion_matrix(labels, pred_labels)
    cm_width = max(len(str(cm.max())), 2)

    print("=== Classification Metrics ===\n")
    print("Metrics independent of threshold:")
    print(f"AUC (ROC):         {auc:.4f}")
    print(f"Average Precision: {avg_prec:.4f}\n")

    print(f"Metrics using threshold = {threshold}:")
    print(f"Accuracy:  {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1 Score:  {f1:.4f}\n")

    print(f"Confusion Matrix using threshold = {threshold}:")
    header = f"Predicted → | {'0'.rjust(cm_width)} | {'1'.rjust(cm_width)} |"
    print(header)

    for i, row in enumerate(cm):
        print(f"Actual {i}    | {str(row[0]).rjust(cm_width)} | {str(row[1]).rjust(cm_width)} |")
    
    return pd.DataFrame({
        'id': data['id'],
        'predicted_prob': predictions
    })
 
def hosmer_lemeshow_test(labels, predictions, g=10, bins_as_quantiles=False):
    """
    Performs the Hosmer-Lemeshow goodness-of-fit test for a logistic regression model.

    Parameters
    ----------
    labels : array-like
        True binary outcomes (0 or 1) for each observation (data['outcome'])

    predictions : array-like
        Predicted probabilities from the logistic regression model for each observation.

    g : int, optional
        Number of groups (bins) to use for the test. Default is 10.

    bins_as_quantiles : bool, optional
        If True, bins are defined based on quantiles of the predicted probabilities.
        If False, bins are uniformly spaced between the minimum and maximum predicted probabilities.

    Returns
    -------
    dict
        Dictionary with the following keys:
        - "hl_stat": Hosmer-Lemeshow test statistic.
        - "p_value": p-value of the test.
    """
    labels = np.array(labels)
    predictions = np.array(predictions)

    if len(labels) != len(predictions):
        raise ValueError("labels e predictions devono avere la stessa lunghezza")

    if bins_as_quantiles:
        quantiles = np.percentile(predictions, np.linspace(0, 100, g + 1))
        bin_method = "quantiles"
    else:
        quantiles = np.linspace(np.min(predictions), np.max(predictions), g + 1)
        bin_method = "uniform"

    groups = np.digitize(predictions, bins=quantiles, right=True)

    observed = np.zeros(g)
    expected = np.zeros(g)
    n_group = np.zeros(g)

    for i in range(1, g + 1):
        mask = groups == i
        n_group[i-1] = np.sum(mask)
        observed[i-1] = np.sum(labels[mask])
        expected[i-1] = np.sum(predictions[mask])

    valid = n_group > 0
    hl_stat = np.sum(((observed[valid] - expected[valid])**2) /
                     (expected[valid] * (1 - expected[valid]/n_group[valid])))

    p = 1 - chi2.cdf(hl_stat, g - 2)

    print(f"Hosmer-Lemeshow test (g={g}, bins selected as {bin_method}):")
    print(f"HL statistic: {hl_stat:.4f}")
    print(f"p-value: {p:.4f}")

    return {"hl_stat": hl_stat, "p_value": p}

def plot_ntcp_stratified(data, parameters_logistic_model, coverage=0.68, min_patients=5, n_quantiles = 2):
    """
    Plots NTCP curves stratified by binary risk factors, including Wilson confidence intervals.

    ## Parameters

    data:  dict
        Usual data object

    parameters_logistic_model : DataFrame
        Logistic model parameters

    coverage : float, optional
        Confidence level for Wilson intervals (default 0.68).

    min_pts : int, optional
        Minimum number of patients per subgroup/bin to display a point and confidence interval (default 5).

    n_quantiles : int, optional
        Number of quantiles used to group patients for confidence interval computation (default 2).

    ## Returns

    None
    Displays NTCP curves with stratified confidence intervals; does not return values.
    """
    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in list(data["risk_factors"].columns)]
    params = parameters_logistic_model[param_order].to_numpy()[0]
    n, a_b = params[[0,1]]
    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()
    names = list(data["risk_factors"].columns)

    EUDs = compute_EUDs(surfaces_matrix, EQD2(original_bins, n_fractions, a_b), n)
    labels_array = np.array(labels).reshape(-1)
    dose_values = np.linspace(np.min(EUDs), np.max(EUDs), 100)

    fig, ax = plt.subplots(1, 1, figsize=(7,5))

    if risk_matrix.shape[1] == 0:  # No risk factors, simple calibration plot
        score = params[2] + dose_values * params[3]
        curve = sigmoid(score)
        ax.plot(dose_values, curve, color='blue', linewidth=2, label='NTCP')
        if len(EUDs) >= min_patients:
            plot_ci(EUDs, labels_array, 'blue', min_patients, coverage, n_quantiles, ax=ax)
        ax.set_title("NTCP calibration plot")
        ax.set_xlabel("EUD [Gy]")
        ax.set_ylabel("NTCP")
        ax.legend(fontsize=8)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.show()
        return

    # Stratified plotting when risk factors exist
    K = risk_matrix.shape[1]
    combos = list(itertools.product([0,1], repeat=K))
    labels_comb = []
    for comb in combos:
        parts = []
        for val, name in zip(comb, names):
            parts.append(name if val == 1 else f"No {name}")
        labels_comb.append("/".join(parts))

    cmap = colormaps['RdBu_r']
    colors = [cmap(i/(len(combos)-1)) for i in range(len(combos))]
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']

    for j, combo in enumerate(combos):
        mask = np.ones(len(EUDs), dtype=bool)
        rf_term = 0
        for k, val in enumerate(combo):
            mask &= (risk_matrix[:, k] == val)
            rf_term += val * params[4 + k]

        score = params[2] + dose_values * params[3] + rf_term
        curve = sigmoid(score)

        ax.plot(dose_values, curve,
                color=colors[j],
                linestyle=linestyles[j % len(linestyles)],
                linewidth=2,
                label=labels_comb[j])

        if np.sum(mask) >= min_patients:
            plot_ci(EUDs[mask], labels_array[mask],
                    colors[j], min_patients, coverage, n_quantiles, ax=ax)

    ax.set_title("NTCP curves stratified")
    ax.set_xlabel("EUD [Gy]")
    ax.set_ylabel("NTCP")
    ax.legend(fontsize=8)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.show()
   
def calibration_plot(data, parameters_logistic_model, bin_width=0.1, coverage=0.68, min_pts=5):
    """
    Plots a calibration plot comparing predicted probabilities from a logistic model to observed proportions in outcomes.

    Parameters
    ----------
    data : dict
        Usual data object

    parameters_logistic_model : DataFrame
        Logistic model parameters

    bin_width : float, optional
        Width of probability bins for calibration (default is 0.1).

    coverage : float, optional
        Confidence level for the Wilson intervals around observed proportions (default is 0.68).

    min_pts : int, optional
        Minimum number of patients required in a probability bin to display a confidence interval (default is 5).

    Notes
    -----
    - Predicted probabilities are divided into bins of width `bin_width`.
    - For each bin with at least `min_pts` observations, the observed proportion of events and the Wilson
      confidence interval are computed.
    - The plot includes:
        * Dashed diagonal line representing perfect calibration.
        * Observed proportions with error bars for the confidence interval.
        * Histogram bars indicating the number of events (red) and non-events (black) per bin.
    - Axes are scaled to show all bins and confidence intervals.
    """ 
    hist_max_height = 0.35 # Maximum height of the histograms, only for visualization

    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in list(data["risk_factors"].columns)]
    params = parameters_logistic_model[param_order].to_numpy()[0]

    labels = data["outcome"]
    pred_probs = predict_logistic(data, params)
    labels = np.array(labels).flatten()

    x_max = np.max(pred_probs)
    bins_edges = np.arange(0, x_max + bin_width, bin_width)
    bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2

    obs_points = []
    obs_centers = []
    ci_lowers = []
    ci_uppers = []
    n_events_bin = []
    n_total_bin = []

    all_bins_centers = []
    all_n_events_bin = []
    all_n_total_bin = []

    for i in range(len(bins_centers)):
        in_bin = (pred_probs >= bins_edges[i]) & (pred_probs < bins_edges[i+1])
        n_tot = np.sum(in_bin)
        n_events = np.sum(labels[in_bin])

        all_bins_centers.append(bins_centers[i])
        all_n_events_bin.append(n_events)
        all_n_total_bin.append(n_tot)

        if n_tot >= min_pts:
            lower, center, upper = wilson_ci(n_events, n_tot, coverage) # Compute wilson confidence interval for each bin
            obs_points.append(center)
            obs_centers.append(bins_centers[i])
            ci_lowers.append(lower)
            ci_uppers.append(upper)
            n_events_bin.append(n_events)
            n_total_bin.append(n_tot)

    scale = hist_max_height / max(all_n_total_bin) if all_n_total_bin else 1

    plt.figure(figsize=(8,8))
    plt.plot([0,1],[0,1], 'k--', label='Ideal calibration')

    if len(obs_centers) > 0:
        plt.errorbar(obs_centers, obs_points,
                     yerr=[np.array(obs_points)-np.array(ci_lowers), np.array(ci_uppers)-np.array(obs_points)],
                     fmt='o', color='blue', markersize=8, markerfacecolor='blue', elinewidth=2, capsize=4)

    for x, ne, nt in zip(all_bins_centers, all_n_events_bin, all_n_total_bin):
        width = bin_width*0.9
        plt.bar(x, ne*scale, width=width, color='red', align='center')
        plt.bar(x, -(nt - ne)*scale, width=width, color='black', align='center')

    plt.axhline(0, color='grey', linestyle='--', linewidth=1)
    plt.xlabel('Predicted probability', fontsize=14)
    plt.ylabel('Observed frequency', fontsize=14)
    plt.title('', fontsize=16)

    plt.legend(handles=[
        plt.Line2D([0],[0], color='k', linestyle='--', label='Ideal calibration'),
        plt.Line2D([0],[0], color='blue', marker='o', linewidth=0, label='Observed proportion (68% CI)'),
        plt.Rectangle((0,0),1,1, color='red', label='Events'),
        plt.Rectangle((0,0),1,1, color='black', label='Non-events')
    ], fontsize=12)

    plt.grid(True)

    xlim_right = min(1.0, np.max(pred_probs) * 1.05)
    plt.xlim(0, xlim_right)

    if len(ci_lowers) > 0 and len(ci_uppers) > 0:
        ylim_lower = -hist_max_height * 1.2
        ylim_upper = max(ci_uppers) * 1.05
    else:
        ylim_lower = -hist_max_height * 1.2
        ylim_upper = 1.0

    plt.ylim(ylim_lower, ylim_upper)
    plt.show()

def likelihood_surface_plots(data, lyman_params, coverage=0.68, n_points=100, bounds=None):
    """
    Generates contour plots of the negative log-likelihood surfaces for Lyman model parameters
    transformed into the logistic NTCP framework. White contours represent the 2D confidence
    region for the specified coverage. Returns also confidence intervals and 4 extreme vertices.

    Parameters
    ----------
    data : dict
        Standard data object containing outcomes and risk factors.
    lyman_params : list or array
        Fitted Lyman model parameters in order:
        [n, alpha/beta, TD50, m, DMF1, DMF2, ..., DMFk]
    coverage : float, optional
        Confidence level for the plotted region (default 0.68)
    n_points : int, optional
        Number of points per axis in the grid (default 50)
    bounds : dict, optional
        Dictionary specifying min and max ranges for parameters.

    Returns
    -------
    dict
        For each pair of parameters:
        - 'ci_vertices': list of 4 tuples (x, y, NLL) corresponding to extreme points
        - 'ci_1d': dictionary with 1D confidence intervals for each parameter
    """       
    param_order = ["n", "a_b", "TD50", "m"] + [f"dmf_{rf}" for rf in list(data["risk_factors"].columns)]
    lyman_params = lyman_params[param_order].to_numpy()[0]
    n_fit, a_b_fit, TD50_fit, m_fit, *dmfs_fit = lyman_params
    dmfs_fit = np.array(dmfs_fit)
    N = len(data["outcome"])
    threshold = compute_Lyman_negll(lyman_params, data) + 0.5 * chi2.ppf(coverage, df=2) / N

    def grid(fit, ws, n_points, vmin=None, vmax=None):
        lower = max(fit - ws, vmin) if vmin is not None else fit - ws
        upper = min(fit + ws, vmax) if vmax is not None else fit + ws
        return np.sort(np.concatenate((np.linspace(lower, upper, n_points), [fit])))

    if bounds is None:
        n_values = grid(n_fit, 0.2, n_points, vmin=0.01, vmax=1)
        a_b_values = grid(a_b_fit, 3, n_points, vmin=0.5, vmax=10)
        T_values = grid(TD50_fit, 15, n_points)
    else:
        n_lower, n_upper = bounds.get("n", (max(n_fit-0.2,0.01), min(n_fit+0.2,1)))
        n_values = np.sort(np.concatenate((np.linspace(n_lower, n_upper, n_points), [n_fit]))) if n_lower <= n_fit <= n_upper else np.linspace(n_lower, n_upper, n_points)

        a_b_lower, a_b_upper = bounds.get("a_b", (max(a_b_fit-3,0.5), min(a_b_fit+3,10)))
        a_b_values = np.sort(np.concatenate((np.linspace(a_b_lower, a_b_upper, n_points), [a_b_fit]))) if a_b_lower <= a_b_fit <= a_b_upper else np.linspace(a_b_lower, a_b_upper, n_points)

        T_lower, T_upper = bounds.get("TD50", (TD50_fit-15, TD50_fit+15))
        T_values = np.sort(np.concatenate((np.linspace(T_lower, T_upper, n_points), [TD50_fit]))) if T_lower <= TD50_fit <= T_upper else np.linspace(T_lower, T_upper, n_points)

    output = {}

    def plot_surface(eval_fun, x_vals, y_vals, xlabel, ylabel, param_x, param_y, fixed_params, fitted_point):
        results = np.array(
            Parallel(n_jobs=-1)(
                delayed(eval_fun)([y_val, x_val]) if param_x=='alpha/beta' else delayed(eval_fun)([x_val, y_val])
                for x_val in x_vals for y_val in y_vals
            )
        ).reshape(len(x_vals), len(y_vals))

        mask = results <= threshold

        X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
        plt.figure(figsize=(8,6))
        cf = plt.contourf(X, Y, results, levels=int(n_points**2 / 5), cmap="jet")
        plt.contour(X, Y, mask, levels=[0.5], colors='white', linewidths=2)
        plt.plot(fitted_point[0], fitted_point[1], 'g*', markersize=12, label='Fitted')
        plt.plot(X[results.argmin()//len(y_vals), results.argmin()%len(y_vals)],
                 Y[results.argmin()//len(y_vals), results.argmin()%len(y_vals)], 'r*', markersize=15, label='Minimum')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(f"$-\\frac{{1}}{{N}} \\log(L)$ depends on {param_x}, {param_y}; {fixed_params} fixed\nmin: {fitted_point[0]:.2f}, {fitted_point[1]:.2f}")
        cbar = plt.colorbar(cf)
        cbar.set_label('Negative log-likelihood')
        plt.legend()
        plt.show()

        x_masked = X[mask]
        y_masked = Y[mask]
        results_masked = results[mask]

        min_x, max_x = round(float(x_masked.min()),3), round(float(x_masked.max()),3)
        min_y, max_y = round(float(y_masked.min()),3), round(float(y_masked.max()),3)
        ci_1d = {param_x: (min_x, max_x), param_y: (min_y, max_y)}

        idx_y_min_x = (x_masked == x_masked.min()).argmax()
        idx_y_max_x = (x_masked == x_masked.max()).argmax()
        idx_x_min_y = (y_masked == y_masked.min()).argmax()
        idx_x_max_y = (y_masked == y_masked.max()).argmax()

        ci_vertices = [
            (min_x, round(float(y_masked[idx_y_min_x]),3), round(float(results_masked[idx_y_min_x]),3)),
            (max_x, round(float(y_masked[idx_y_max_x]),3), round(float(results_masked[idx_y_max_x]),3)),
            (round(float(x_masked[idx_x_min_y]),3), min_y, round(float(results_masked[idx_x_min_y]),3)),
            (round(float(x_masked[idx_x_max_y]),3), max_y, round(float(results_masked[idx_x_max_y]),3))
        ]

        output[f"{param_x}_vs_{param_y}"] = {
            "ci_vertices": ci_vertices,
            "ci_1d": ci_1d
        }

    plot_surface(
        lambda params_xy: compute_Lyman_negll([params_xy[0], params_xy[1], TD50_fit, m_fit] + list(dmfs_fit), data),
        a_b_values, n_values,
        r'$\alpha/\beta$ [Gy]', r'$n$',
        'alpha/beta', 'n',
        'TD50, m, DMFs',
        (a_b_fit, n_fit)
    )

    plot_surface(
        lambda params_xy: compute_Lyman_negll([n_fit, params_xy[1], params_xy[0], m_fit] + list(dmfs_fit), data),
        a_b_values, T_values,
        r'$\alpha/\beta$ [Gy]', r'$EUD_{50}$ [Gy]',
        'alpha/beta', 'TD50',
        'n, m, DMFs',
        (a_b_fit, TD50_fit)
    )

    plot_surface(
        lambda params_xy: compute_Lyman_negll([params_xy[0], a_b_fit, params_xy[1], m_fit] + list(dmfs_fit), data),
        n_values, T_values,
        r'$n$', r'$EUD_{50}$ [Gy]',
        'n', 'TD50',
        'alpha/beta, m, DMFs',
        (n_fit, TD50_fit)
    )

    output["Likelihood_threshold"] = round(float(threshold), 5)
    return output

def wilson_ci(n_events, n_tot, coverage):
    """
    Computes the Wilson score confidence interval for a binomial proportion.

    Parameters
    ----------
    n_events : int
        Number of observed events (successes).

    n_tot : int
        Total number of trials.

    coverage : float
        Desired confidence level (e.g., 0.95 for a 95% confidence interval).

    Returns
    -------
    tuple of floats
        - lower : float, lower bound of the confidence interval
        - p : float, observed proportion n_events / n_tot
        - upper : float, upper bound of the confidence interval

    Notes
    -----
    If n_tot is 0, returns (0.0, 0.0, 0.0).
    """
    if n_tot == 0:
        return 0.0, 0.0
    p = n_events / n_tot
    alpha = 1 - coverage
    z = norm.ppf(1 - alpha/2)
    denom = 1 + z**2 / n_tot
    centre = p + z**2 / (2 * n_tot)
    half_width = z * np.sqrt(p*(1 - p) / n_tot + z**2 / (4 * n_tot**2))
    lower = (centre - half_width) / denom
    upper = (centre + half_width) / denom
    return lower, p, upper

def plot_ci(EUD_subset, labels_subset, color, min_pts, coverage, n_quantiles, ax=None):
    """
    Support function for plot_ntcp_stratified function
    
    Divides data in bins (computed with quantiles). 
    Then, for each bin, plots the observed proportion of positive events with error bars representing 
    the Wilson confidence interval, and annotates the bin with the number of events over total.

    Note: Variables have the suffix "subset" because they represent subgroups of patients, stratified by clinical factors.

    Parameters
    ----------
    EUD_subset : array-like
        EUD values used to define bins.

    labels_subset : array-like
        Binary outcomes corresponding to the EUD values.

    color : str
        Color used for markers and error bars.

    min_pts : int
        Minimum number of observations required in a bin to plot the confidence interval.

    coverage : float
        Confidence level for the Wilson interval, in [0, 1] range (e.g., 0.95 for 95% CI).

    n_quantiles : int
        Number of quantile-based bins to divide the EUD values into.

    ax : matplotlib.axes.Axes, optional
        Axes object on which to plot. If None, uses the current axes.
    """
    if ax is None:
        ax = plt.gca()
    if len(EUD_subset) < 3:
        return
    bins_edges = np.quantile(EUD_subset, np.linspace(0, 1, n_quantiles + 1))
    bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
    for i in range(len(bins_centers)):
        if i == len(bins_centers)-1:
            in_bin = (EUD_subset >= bins_edges[i]) & (EUD_subset <= bins_edges[i+1])
        else:
            in_bin = (EUD_subset >= bins_edges[i]) & (EUD_subset < bins_edges[i+1])
        n_tot = np.sum(in_bin)
        if n_tot >= min_pts:
            N_events = np.sum(labels_subset[in_bin])
            lower, center, upper = wilson_ci(N_events, n_tot, coverage)
            ax.errorbar(
                bins_centers[i], center,
                yerr=[[center-lower],[upper-center]],
                fmt='o', color=color, capsize=5,
                markersize=10, markerfacecolor=color, elinewidth=2
            )
            ax.text(
                bins_centers[i] + 0.6, center,
                f'{int(N_events)}/{int(n_tot)}',
                color=color, fontsize=13, verticalalignment='center'
            )

def compute_Lyman_negll(params, data):
    """
    Computes the negative log-likelihood for the logistic NTCP model using Lyman model parameters.

    Parameters
    ----------
    params : array-like
        Vector of Lyman model parameters in order:
        [n, alpha/beta, TD50, m, DMF1, DMF2, ..., DMFk]

    data : dict
        Standard data object with outcomes, dose distributions, fractions, and risk factors.

    Returns
    -------
    negll : float
        Mean negative log-likelihood of the model predictions.
    """
    n, a_b, TD50, m, *dmfs = params
    dmfs = np.array(dmfs)

    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()

    if len(dmfs) == 0:
        b0 = -1 / m
        b_eud = -b0 / TD50
        b_risk = np.array([])
    else:
        b0 = 1 / (m * (len(dmfs) - 1 - np.sum(dmfs)))
        b_risk = b0 * (dmfs - 1)
        b_eud = -b0 / TD50

    beta = np.concatenate(([b0, b_eud], b_risk))
    EUDs = compute_EUDs(surfaces_matrix, EQD2(original_bins, n_fractions, a_b), n)

    if risk_matrix.size == 0:
        scores_logistic = beta[0] + beta[1] * EUDs
    else:
        scores_logistic = beta[0] + beta[1] * EUDs + risk_matrix @ beta[2:]

    preds = sigmoid(scores_logistic)
    negll = -np.mean(labels * np.log(preds) + (1 - labels) * np.log(1 - preds))
    return negll