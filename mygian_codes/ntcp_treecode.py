from ntcp_common_imports import *
from ntcp_optimization_function import *
from ntcp_functions import *

def reorder_coefficients(coefficients, risk_factors, risk_factors_sorted):
    #  Reorders logistic model coefficients to match the sorted risk factor order.
    beta_0 = coefficients[0]
    beta_EUD = coefficients[1]
    betas_dict = dict(zip(risk_factors, coefficients[2:]))
    reordered = [beta_0, beta_EUD] + [betas_dict[rf] for rf in risk_factors_sorted]
    return np.array(reordered)

def sort_risk_factors(df_bin, risk_factors, order=True):
    # Sorts risk factors by their frequency deviation from 0.5 to prioritize splits in the tree.
    if not order:
        return risk_factors
    freq = df_bin[risk_factors].mean()
    return (freq - 0.5).abs().sort_values().index.tolist()

def compute_predictions(EUDs, risk_matrix, coefficients):
    # Computes predicted probabilities using the logistic NTCP model for all patients.
    scores = coefficients[0] + EUDs * coefficients[1] + risk_matrix @ coefficients[2:]
    return np.array([sigmoid(s) for s in scores])

def build_tree(G, df, risk_factors_sorted, predictions, labels, depth=0, parent='root', path=""):
    # Recursively builds the hierarchical tree by splitting patients at each binary risk factor.
    if not risk_factors_sorted:
        return

    current_factor = risk_factors_sorted[0]

    for val in [0, 1]:
        mask = df[current_factor] == val
        if parent != 'root':
            mask &= df['mask_parent']
        n_patients = int(mask.sum())
        if n_patients == 0:
            continue
        n_events = int(labels[mask].sum())
        predicted_ratio = np.median(predictions[mask]) if n_patients > 0 else np.nan

        node_label = "NO" if val == 0 else "YES"
        node_name = f"{path}{node_label}"

        G.add_node(node_name,
                   n_patients=n_patients,
                   n_events=n_events,
                   factor=current_factor,
                   value=val,
                   depth=depth + 1,
                   predicted_ratio=predicted_ratio)
        G.add_edge(parent, node_name)

        df_next = df.copy()
        df_next['mask_parent'] = mask
        build_tree(G, df_next, risk_factors_sorted[1:], predictions, labels,
                   depth + 1, node_name, node_name + "/")

def layout_hierarchy_compact(G, node, pos, x=0, dx=1.0):
    # Assigns x/y positions to nodes for a compact visualization of the tree.
    children = list(G.successors(node))
    depth = G.nodes[node]['depth']
    pos[node] = (x, -depth * 2.0)
    if not children:
        return
    width = dx / len(children)
    for i, child in enumerate(children):
        layout_hierarchy_compact(G, child, pos, x - dx / 2 + width / 2 + i * width, dx=width)

def draw_risk_tree(G, pos, risk_factors_sorted, color_pred=True, tolerance=0.05):
    # Plots the tree using networkx/Matplotlib, adding patient counts, event rates, predicted probabilities, and node labels.
    total_width = max(x for x, _ in pos.values())
    plt.figure(figsize=(total_width / 2, len(risk_factors_sorted) * 4))
    nx.draw_networkx_nodes(G, pos, node_shape='s', node_size=22000, node_color='lightblue')
    nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-|>', arrowsize=18, edge_color='gray')

    for n in G.nodes:
        node_data = G.nodes[n]
        n_pat = node_data['n_patients']
        n_evt = node_data['n_events']
        rate = 100 * n_evt / n_pat if n_pat > 0 else 0
        pred = 100 * node_data.get('predicted_ratio', np.nan)
        x, y = pos[n]

        if color_pred and not np.isnan(pred):
            diff = pred - rate
            if diff > tolerance * 100:
                col = "red"
            elif diff < -tolerance * 100:
                col = "blue"
            else:
                col = "darkgreen"
        else:
            col = "black"

        val_text = ""
        if n != "root":
            val_text = "NO" if node_data.get("value", 0) == 0 else "YES"

        label = f"{val_text}\n{n_pat} patients\n{n_evt} events\nObs = {rate:.1f}%"
        plt.text(x, y + 0.125, label, fontsize=18, ha='center', va='center', color='black', linespacing=1.5)
        if not np.isnan(pred):
            plt.text(x, y - 0.62, f"Pred = {pred:.1f}%", color=col,
                     fontsize=18, ha='center', va='center', linespacing=1.5)

    for depth, rf in enumerate(risk_factors_sorted, start=1):
        y_pos = -depth * 2.0
        plt.text(-5, y_pos, rf, fontsize=22, fontweight='bold',
                 ha='right', va='center', color='black')
        plt.hlines(y_pos - 1, -2, total_width + 2, colors='gray',
                   linestyles='dashed', linewidth=1.2)

    plt.axis('off')
    plt.show()

def risk_tree(data, parameters_logistic_model, order=True, color_pred=True, tolerance=0.05):
    """
    Constructs and visualizes a decision tree of binary risk factors for logistic NTCP predictions.

    The function builds a hierarchical tree representing all risk factors, splitting patients 
    based on risk factor presence/absence. 
    Each node displays the number of patients, events, observed event rate, and optionally predicted event 
    probability from the logistic NTCP model. Predicted probabilities are colored based on the difference between 
    observed and predicted rates:
    (green: predicted and observed rates are within tolerance, blue: predicted rates are lower than observed, red: predicted rates are higher than observed)
    
    Parameters
    ----------
    data : dict
        Dataset containing:
        - "outcome": binary labels (0/1)
        - "n_fractions": number of fractions delivered
        - "DSH": dose-surface histogram matrix
        - "original_bins": original dose bins
        - "risk_factors": DataFrame with binary risk factors

    parameters_logistic_model : DataFrame
        DataFrame containing fitted logistic NTCP parameters: [n, a_b, beta_0, beta_EUD, beta_<risk factors>].

    order : bool, optional
        Whether to sort risk factors based on the balance between events and non-events in each branch before building the tree (default is True).

    color_pred : bool, optional
        Whether to color predicted rates text based on the difference between observed and predicted event rates (default is True).

    tolerance : float, optional
        Threshold (%) for coloring nodes if predicted rates differ from observed rates (default is 0.05).

    Returns
    -------
    None
        Displays a hierarchical tree plot of risk factors with patient counts, events, observed rates, and predicted rates.
    """
    labels = data["outcome"]
    n_fractions = data["n_fractions"]
    surfaces_matrix = data["DSH"]
    original_bins = data["original_bins"]
    risk_matrix = data["risk_factors"].to_numpy()
    risk_factors = list(data["risk_factors"].columns)

    param_order = ["n", "a_b", "beta_0", "beta_EUD"] + [f"beta_{rf}" for rf in list(data["risk_factors"].columns)]
    parameters_logistic_model = parameters_logistic_model[param_order].to_numpy()[0]
    
    n, a_b = parameters_logistic_model[[0, 1]]
    coefficients = parameters_logistic_model[2:]
    
    if len(risk_factors) == 0 or risk_matrix.size == 0:
        return None

    df = pd.concat([pd.DataFrame({"label": labels}), pd.DataFrame(risk_matrix, columns=risk_factors)], axis=1)
    risk_factors_sorted = sort_risk_factors(df, risk_factors, order)  # assumes this sorts based on event balance
    coefficients = reorder_coefficients(coefficients, risk_factors, risk_factors_sorted)
    EUDs = compute_EUDs(surfaces_matrix, EQD2(original_bins, n_fractions, a_b), n)
    predictions = sigmoid(coefficients[0] + coefficients[1] * EUDs + risk_matrix @ coefficients[2:])

    G = nx.DiGraph()
    root = "root"
    G.add_node(root, n_patients=len(df), n_events=int(df["label"].sum()), depth=0)
    G.nodes[root]['predicted_ratio'] = np.median(predictions)

    df['mask_parent'] = np.ones(len(df), dtype=bool)
    build_tree(G, df, risk_factors_sorted, predictions, labels)

    pos = {}
    total_width = max(30, (2 ** len(risk_factors) + len(risk_factors)**2) * 3)
    layout_hierarchy_compact(G, root, pos, x=total_width/2, dx=total_width)

    draw_risk_tree(G, pos, risk_factors_sorted, color_pred=color_pred, tolerance=tolerance)