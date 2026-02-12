from ntcp_gian_common_imports import *
from ntcp_gian_optimization_function import *

def check_convergence(data, coverage=0.68, do_plot=False):
    '''
    The purpose of this procedure is to verify the robustness of the optimization algorithm. By exploring a substantially broader parameter 
    grid, retaining only the initial guesses that satisfy the likelihood-ratio criterion (alpha = 1 - coverage) and with unique values of 
    n and alpha/beta ratio, the refinement step begins from widely separated and potentially suboptimal starting points. 
    If the algorithm consistently converges to the same optimum despite these distant initial conditions, this indicates that the optimization
    landscape is well-behaved and that the identified solution is stable and reliable.

    Parameters
    ----------
    data : dict
        Usual data object
    coverage : float, optional
        Statistical coverage used in the chi-square confidence threshold that defines which
        brute-force solutions are kept. Default = 0.68.
    do_plot : bool, optional
        If True, displays a scatter plot of the brute-force negative log-likelihood values with
        the threshold line (time-consuming).

    Returns
    -------
    best_refined_results : ndarray
        Best solution found.
    df_all_refined_results : pandas.DataFrame
        DataFrame containing all refined solutions sorted by NLL.
    df_initial_kept_with_NLL : pandas.DataFrame
        DataFrame containing the initial brute-force solutions that passed the threshold,
        ordered in the same NLL-sorted sequence as df_all_refined_results.
    all_results_bf : ndarray
        Array of all brute-force results before filtering.
    '''

    print("Starting grid search...")
    # Call find_opt_brute_force using huge_grid = True
    all_results_bf = find_opt_brute_force(data, parallelize=True, huge_grid=True, verbose=10) 
    
    print("Grid search completed:")
    # Compute the likelihood threshold
    threshold = np.min(all_results_bf[:, -1]) + 0.5 * chi2.ppf(coverage, df = 4 + data['risk_factors'].shape[1]) / len(data['outcome'])

    if do_plot:
        # Plot the NLL column, sorted
        last_col_sorted = np.sort(all_results_bf[:, -1])
        plt.figure(figsize=(10, 6))
        plt.scatter(range(len(last_col_sorted)), last_col_sorted, 
                    c=['red' if v <= threshold else 'blue' for v in last_col_sorted])
        plt.axhline(y=threshold, color='r', linestyle='--', linewidth=2)
        plt.ylabel('Negative log-likelihood')
        plt.title('Brute force grid search results')
        plt.show()

    # The NLL has to respect the threshold criterion
    initial_guesses_with_NLL = all_results_bf[all_results_bf[:, -1] <= threshold, :]
    min_idx = np.argmin(initial_guesses_with_NLL[:, -1])
    initial_guesses_kept_with_NLL = [initial_guesses_with_NLL[min_idx]]

    # The values of (n, alpha/beta) must be unique
    selected_pairs = {tuple(initial_guesses_with_NLL[min_idx, :2])}

    for i, row in enumerate(initial_guesses_with_NLL):
        if i == min_idx:
            continue
        pair = tuple(row[:2])
        if pair not in selected_pairs:
            initial_guesses_kept_with_NLL.append(row)
            selected_pairs.add(pair)

    initial_guesses_kept_with_NLL = np.array(initial_guesses_kept_with_NLL)
    initial_guesses_kept = initial_guesses_kept_with_NLL[:, :-1]

    print(f"Selected {initial_guesses_kept.shape[0]} out of {all_results_bf.shape[0]} " +
          f"solutions ({initial_guesses_kept.shape[0]/all_results_bf.shape[0] * 100:.2f} %), using likelihood threshold = {threshold:.5f}")

    # Refine solutions
    all_refined_results = refine_solutions(data, initial_guesses_kept, parallelize=True, verbose=5)
    best_refined_results = all_refined_results[np.argmin(all_refined_results[:, -1]), :]

    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in data["risk_factors"].columns]

    df_all_refined_results = pd.DataFrame(all_refined_results, columns=param_order + ["NLL"])
    df_initial_kept_with_NLL = pd.DataFrame(initial_guesses_kept_with_NLL, columns=param_order + ["NLL"])

    order = np.argsort(df_all_refined_results["NLL"].values)
    df_all_refined_results = df_all_refined_results.iloc[order].reset_index(drop=True)
    df_initial_kept_with_NLL = df_initial_kept_with_NLL.iloc[order].reset_index(drop=True)

    return best_refined_results, df_all_refined_results, df_initial_kept_with_NLL, all_results_bf

def plot_refined_vs_initial(df_init, df_refined, N, x_col=None, y_col=None,
                            s_base=40, fontsize=14, conf=0.68, verbose=1):
    """
    Plots and compares clustered initial solutions versus clustered refined solutions for two
    selected parameters, and prints a summary table of the refined representative clusters.

    The function groups nearby solutions in both the initial and refined sets using parameter-
    specific tolerances, identifies a representative solution for each refined cluster (the one
    with minimal NLL), and visualizes the clustering in two side-by-side scatter plots:
    - left panel: clustered initial solutions
    - right panel: clustered refined solutions with Likelihood-Ratio p-values

    The goal is to visually assess whether distant initial guesses converge toward the same
    regions in parameter space after refinement, providing insight into the stability and
    robustness of the optimization procedure.

    Parameters
    ----------
    df_init : pandas.DataFrame
        DataFrame containing the initial solutions (before refinement), including parameter
        values and NLL.

    df_refined : pandas.DataFrame
        DataFrame containing the refined solutions, including parameter values and NLL.

    N : int
        Number of samples used in the likelihood (needed for LR test computation).
        Typically equal to len(data['outcome']).

    x_col : str or None, optional
        Name of the parameter to plot on the x-axis. If None, the first column of df_refined
        is used.

    y_col : str or None, optional
        Name of the parameter to plot on the y-axis. If None, the second column of df_refined
        is used.

    s_base : float, optional
        Base size of the scatter points. Cluster size scales with the number of points.

    fontsize : int, optional
        Font size for annotations and axis labels.

    conf : float, optional
        Confidence level for the Likelihood-Ratio test when computing the p-values of refined
        solutions. Default is 0.68.

    verbose : int, optional
        Controls the amount of textual information shown on the plots.
        - verbose = 0: no annotations
        - verbose = 1: show only the index of the representative solution
        - verbose = 2: show index, cluster size (if >1), NLL value and LR-test p-value

    Returns
    -------
    None

        The clusters are defined based on the following distance condition:

        xi and xj belong to the same cluster if
        sqrt( ((xi - xj) / tol_x)^2 + ((yi - yj) / tol_y)^2 ) <= 1

        Where tol_x and tol_y are defined using the dictionary:

        tols = {"n": 0.01,
                "a_b": 0.2,
                "beta_0": 0.5,
                "beta_EUD": 0.005,
                "beta_risk": 0.1}
    """
    if x_col is None:
        x_col = df_refined.columns[0]
    if y_col is None:
        y_col = df_refined.columns[1]

    tols = {"n":0.01, "a_b": 0.2, "beta_0": 0.5, "beta_EUD": 0.005, "beta_risk": 0.1}

    if x_col in tols:
        tol_x = tols[x_col]
    elif x_col.startswith("beta_") and x_col not in ["beta_0", "beta_EUD"]:
        tol_x = tols["beta_risk"]
    else:
        tol_x = None

    if y_col in tols:
        tol_y = tols[y_col]
    elif y_col.startswith("beta_") and y_col not in ["beta_0", "beta_EUD"]:
        tol_y = tols["beta_risk"]
    else:
        tol_y = None

    x_min = min(df_init[x_col].min(), df_refined[x_col].min())
    x_max = max(df_init[x_col].max(), df_refined[x_col].max())
    y_min = min(df_init[y_col].min(), df_refined[y_col].min())
    y_max = max(df_init[y_col].max(), df_refined[y_col].max())
    dx = (x_max - x_min) * 0.05
    dy = (y_max - y_min) * 0.05

    param_cols = [c for c in df_refined.columns if c not in ["NLL", "LR_prob"]]
    param_count = len(param_cols)
    NLL_min = df_refined["NLL"].min()

    probs = []
    for idx, row in df_refined.iterrows():
        if row["NLL"] == NLL_min:
            probs.append(1.0)
        else:
            LR = 2 * N * (row["NLL"] - NLL_min)
            p = 1 - chi2.cdf(LR, param_count)
            probs.append(p)
    df_refined["LR_prob"] = probs

    points_to_plot = []
    coords = df_refined[[x_col, y_col]].values
    assigned = np.zeros(len(df_refined), dtype=bool)
    for i, (xi, yi) in enumerate(coords):
        if assigned[i]:
            continue
        group_idx = [i]
        for j in range(i+1, len(coords)):
            if assigned[j]:
                continue
            xj, yj = coords[j]
            d = np.sqrt(((xi - xj)/tol_x)**2 + ((yi - yj)/tol_y)**2)
            if d <= 1:
                group_idx.append(j)
        assigned[group_idx] = True
        group = df_refined.iloc[group_idx]
        idx_min = group["NLL"].idxmin()
        points_to_plot.append({
            "x": group[x_col].mean(),
            "y": group[y_col].mean(),
            "idx_min": idx_min,
            "count": len(group),
            "indices": list(group.index)
        })

    indices_to_show_right = [p["idx_min"] for p in points_to_plot]

    points_init_to_plot = []
    coords_init = df_init[[x_col, y_col]].values
    assigned_init = np.zeros(len(df_init), dtype=bool)
    for i, (xi, yi) in enumerate(coords_init):
        if assigned_init[i]:
            continue
        group_idx = [i]
        for j in range(i+1, len(coords_init)):
            if assigned_init[j]:
                continue
            xj, yj = coords_init[j]
            d = np.sqrt(((xi - xj)/tol_x)**2 + ((yi - yj)/tol_y)**2)
            if d <= 1:
                group_idx.append(j)
        assigned_init[group_idx] = True
        group = df_init.iloc[group_idx]
        points_init_to_plot.append({
            "x": group[x_col].mean(),
            "y": group[y_col].mean(),
            "count": len(group),
            "indices": list(group.index)
        })

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    for p in points_to_plot:
        size = s_base * min(p["count"], 10)
        axes[1].scatter(p["x"], p["y"], s=size, c='gray')
        nll_val = df_refined.loc[p["idx_min"], "NLL"]
        pval = df_refined.loc[p["idx_min"], "LR_prob"]
        color = 'green' if pval >= (1-conf) else 'red'

        if verbose == 1:
            axes[1].annotate(f"{p['idx_min']}",
                             (p["x"] + dx*0.03, p["y"] + dy*0.03),
                             fontsize=fontsize, color=color)

        if verbose == 2:
            if p["count"] == 1:
                axes[1].annotate(f"{p['idx_min']} (NLL={nll_val:.4f}, p={pval:.2f})",
                                 (p["x"] + dx*0.03, p["y"] + dy*0.03),
                                 fontsize=fontsize, color=color)
            else:
                axes[1].annotate(f"{p['idx_min']} ({p['count']} pts, NLL={nll_val:.4f}, p={pval:.2f})",
                                 (p["x"] + dx*0.03, p["y"] + dy*0.03),
                                 fontsize=fontsize, color=color)

    axes[1].set_xlabel(x_col, fontsize=fontsize+2)
    axes[1].set_ylabel(y_col, fontsize=fontsize+2)
    axes[1].set_title(f"Refined solutions: {x_col} vs {y_col}", fontsize=fontsize+4)
    axes[1].set_xlim(x_min - dx, x_max + dx)
    axes[1].set_ylim(y_min - dy, y_max + dy)

    for p in points_init_to_plot:
        size = s_base * min(p["count"], 10)
        axes[0].scatter(p["x"], p["y"], s=size, c='gray')
        idx_in_right = [i for i in p["indices"] if i in indices_to_show_right]

        if verbose == 1 and idx_in_right:
            idx = idx_in_right[0]
            nll_val = df_refined.loc[idx, "NLL"]
            pval = df_refined.loc[idx, "LR_prob"]
            color = 'green' if pval >= (1-conf) else 'red'
            axes[0].annotate(f"{idx}",
                             (p["x"] + dx*0.03, p["y"] + dy*0.03),
                             fontsize=fontsize, color=color)

        if verbose == 2 and idx_in_right:
            label = ", ".join(map(str, idx_in_right))
            if p["count"] > 1:
                axes[0].annotate(f"{label}, {p['count']} pts",
                                 (p["x"] + dx*0.03, p["y"] + dy*0.03),
                                 fontsize=fontsize, color='black')
            else:
                axes[0].annotate(f"{label}",
                                 (p["x"] + dx*0.03, p["y"] + dy*0.03),
                                 fontsize=fontsize, color='black')

    axes[0].set_xlabel(x_col, fontsize=fontsize+2)
    axes[0].set_ylabel(y_col, fontsize=fontsize+2)
    axes[0].set_title(f"Initial solutions in the low-NLL region: {x_col} vs {y_col}", fontsize=fontsize+4)
    axes[0].set_xlim(x_min - dx, x_max + dx)
    axes[0].set_ylim(y_min - dy, y_max + dy)

    plt.tight_layout()
    plt.show()

    table = []
    for p in points_to_plot:
        idx = p["idx_min"]
        nll_val = df_refined.loc[idx, "NLL"]
        param_vals = df_refined.loc[idx, param_cols]
        if p["count"] > 1:
            idx_str = f"{idx}, {p['count']} pts"
        else:
            idx_str = f"{idx}"
        row = [idx_str] + [f"{v:.2f}" for v in param_vals] + [f"{nll_val:.4f}"]
        table.append(row)

    headers = ["indice"] + param_cols + ["NLL"]
    col_widths = [max(len(str(row[i])) for row in ([headers] + table))
                  for i in range(len(headers))]
    header_row = "  ".join(f"{h:<{col_widths[i]}}" for i, h in enumerate(headers))
    print(header_row)
    print("-" * len(header_row))
    for row in table:
        print("  ".join(f"{row[i]:<{col_widths[i]}}" for i in range(len(row))))
