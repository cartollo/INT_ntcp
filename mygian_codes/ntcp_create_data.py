import pandas as pd
import numpy as np

def create_data(outcome, cumulative_dsh, n_fractions, risk_factors=None, impute_risk_factors=True):
    """
    Prepares and aligns clinical and dosimetric data for NTCP modeling.

    This function aligns outcome labels, cumulative dose-volume histograms,
    fractionation data, and optional clinical risk factors using patient IDs
    as the unique key. The cumulative DVH is converted into a Differential
    Dose-Volume Histogram (DSH).

    Missing values in risk factors can either be imputed with zeros or cause
    the corresponding patients to be removed from all datasets.

    Parameters
    ----------
    outcome : pandas.DataFrame
        DataFrame containing patient outcomes. The first column must contain
        patient IDs, and the second column the binary outcome (0/1).

    cumulative_dsh : pandas.DataFrame
        DataFrame containing cumulative dose-volume histograms. The first column
        must contain patient IDs, and the remaining columns cumulative DVH values.

    n_fractions : int, float, or pandas.DataFrame
        Number of fractions per patient. If a scalar is provided, the same value
        is assigned to all patients. If a DataFrame is provided, the first column
        must be patient IDs and the second the number of fractions.

    risk_factors : pandas.DataFrame, optional
        DataFrame containing additional clinical risk factors. The first column
        must contain patient IDs. Remaining columns are treated as covariates.

    impute_risk_factors : bool, default=True
        If True, missing values in risk factors are imputed with 0.
        If False, patients with missing risk factors are removed and all other
        data are realigned accordingly.

    Returns
    -------
    dict
        Dictionary containing:
            - "outcome": numpy array of binary outcomes
            - "DSH": numpy array of Differential Dose-Volume Histograms
            - "n_fractions": numpy array of fraction numbers
            - "original_bins": numpy array of DSH bin starting points
            - "risk_factors": pandas DataFrame of aligned risk factors
            - "id": numpy array of patient IDs
    """
    # --- Standardize column names and copy inputs ---
    outcome = outcome.copy().rename(columns={outcome.columns[0]: "id"})
    cumulative_dsh = cumulative_dsh.copy().rename(columns={cumulative_dsh.columns[0]: "id"})

    # --- Standardize n_fractions (rename columns and create data frame if a single value is given) ---
    if isinstance(n_fractions, (int, float)):
        n_fractions = pd.DataFrame({
            "id": outcome["id"],
            "n_fractions": n_fractions
        })
    else:
        n_fractions = n_fractions.copy()
        n_fractions = n_fractions.rename(columns={
            n_fractions.columns[0]: "id",
            n_fractions.columns[1]: "n_fractions"
        })

    # --- Standardize risk factors (rename columns) ---
    if risk_factors is not None and risk_factors.shape[1] > 0:
        risk_factors = risk_factors.copy().rename(columns={risk_factors.columns[0]: "id"})
        if any(c is None or c == "" for c in risk_factors.columns[1:]):
            risk_factors.columns = ["id"] + [
                f"risk_{i+1}" for i in range(risk_factors.shape[1] - 1)
            ]

    # --- Build master table defining valid IDs and their order ---
    data = outcome.merge(n_fractions, on="id", how="inner")

    if risk_factors is not None and risk_factors.shape[1] > 1:
        data = data.merge(risk_factors, on="id", how="left")
        rf_cols = risk_factors.columns[1:]

        if impute_risk_factors:
            data[rf_cols] = data[rf_cols].fillna(0)
        else:
            data = data.dropna(subset=rf_cols)

    # --- Align cumulative DSH to the final set and order of IDs ---
    valid_ids = data["id"]
    cumulative_dsh = (
        cumulative_dsh[cumulative_dsh["id"].isin(valid_ids)]
        .set_index("id")
        .loc[valid_ids]
        .reset_index()
    )

    # --- Compute Differential DSH ---
    DSH, original_bins = compute_DSH(cumulative_dsh)
    DSH = DSH.iloc[:, 1:]

    # --- Extract arrays for modeling ---
    outcome_array = data.iloc[:, 1].to_numpy()
    n_frac = data["n_fractions"].to_numpy()
    id_array = data["id"].to_numpy()

    if risk_factors is not None and risk_factors.shape[1] > 1:
        risk_factors_df = data[rf_cols].reset_index(drop=True)
    else:
        risk_factors_df = pd.DataFrame(index=data.index)

    return {
        "outcome": outcome_array,
        "DSH": DSH.to_numpy(),
        "n_fractions": n_frac,
        "original_bins": original_bins,
        "risk_factors": risk_factors_df,
        "id": id_array
    }

def compute_DSH(cumulative_dsh):
    """
    Computes the Differential Dose-Volume Histogram (DSH) from a cumulative dose-volume histogram.

    This function converts a cumulative dose distribution into a differential histogram by taking 
    the difference between consecutive dose bins. It also normalizes values if they are in percentage 
    format (0-100) and returns the corresponding original bin edges.

    Parameters
    ----------
    cumulative_dsh : pandas.DataFrame
        DataFrame containing cumulative dose distributions. The first column should be patient IDs,
        subsequent columns contain cumulative dose-volume values for each dose bin.
        IMPORTANT: Each column name must represent the corresponding dose level.
    Returns
    -------
    DSH : pandas.DataFrame
        DataFrame containing the differential dose-volume histogram. The first column is patient IDs,
        and the remaining columns correspond to differential dose bins in the format "[bin_start, bin_end)".

    original_bins : numpy.ndarray
        Array of the starting values of each dose bin, extracted from the column names of the DSH.
    """
    cumulative_dsh = cumulative_dsh.copy()
    cumulative_dsh.fillna(0, inplace=True)
    if cumulative_dsh.iloc[:, 1:].max().max() > 1:
        cols = cumulative_dsh.columns[1:]
        cumulative_dsh[cols] = cumulative_dsh[cols].to_numpy(dtype=float) / 100.0
    DSH = cumulative_dsh.iloc[:, [0]]
    diff_cols = {}
    dsh_cols = []
    for i in range(1, len(cumulative_dsh.columns) - 1):
        col1 = cumulative_dsh.columns[i]
        col2 = cumulative_dsh.columns[i + 1]
        new_col_name = f"[{col1}, {col2})"
        dsh_cols.append(new_col_name)
        diff_cols[new_col_name] = cumulative_dsh[col1] - cumulative_dsh[col2]
    DSH = pd.concat([DSH, pd.DataFrame(diff_cols)], axis=1)
    last_col = cumulative_dsh.columns[-1]
    DSH[last_col] = cumulative_dsh[last_col]
    dsh_cols.append(last_col)
    original_bins = [float(str(col).strip('[]()').split(',')[0]) for col in dsh_cols]
    return DSH, np.array(original_bins)

def data_to_dataframe(data):
    """
    Converts the output of create_data into a single pandas DataFrame.

    The resulting DataFrame contains one row per patient and includes:
    - patient identifier
    - outcome label
    - number of fractions
    - clinical risk factors (if present)
    - differential DSH values

    DSH columns are named using the corresponding values in
    data["original_bins"], ensuring a direct mapping between dose bins
    and histogram values.

    Parameters
    ----------
    data : dict
        Dictionary returned by create_data, with keys:
        "id", "outcome", "n_fractions", "DSH",
        "original_bins", and "risk_factors".

    Returns
    -------
    pandas.DataFrame
        DataFrame with all patient-level information and DSH bins
        as separate columns.
    """
    df = pd.DataFrame({
        "id": data["id"],
        "outcome": data["outcome"],
        "n_fractions": data["n_fractions"]
    })

    if data["risk_factors"] is not None and data["risk_factors"].shape[1] > 0:
        rf = data["risk_factors"].reset_index(drop=True)
        df = pd.concat([df, rf], axis=1)

    dsh_cols = [str(b) for b in data["original_bins"]]
    dsh_df = pd.DataFrame(data["DSH"], columns=dsh_cols)

    df = pd.concat([df, dsh_df], axis=1)

    return df