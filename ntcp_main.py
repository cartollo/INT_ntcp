from ntcp_common_imports import *
from ntcp_check_convergence import *
from ntcp_create_data import *
from ntcp_functions import *
from ntcp_global_variables import *
from ntcp_optimization_function import *
from ntcp_treecode import *
from ntcp_visualization_cl import *


# # Import Data
path = r'Synthetic_Data'
outcome = pd.read_excel(os.path.join(path, "outcome.xlsx")) # ID + outcome variable
cumulative_dsh = pd.read_excel(os.path.join(path, "cumulative_dsh.xlsx")) # ID + cumulative DSH
n_fractions = pd.read_excel(os.path.join(path, "n_fractions.xlsx")) # ID + number of fractions
clinical_factors = pd.read_excel(os.path.join(path, "clinical_factors.xlsx")) # ID + columns of risk factors

data = create_data(outcome, cumulative_dsh, n_fractions, clinical_factors) # Aligns data and impute missing values in risk factors

data_to_dataframe(data)

# Fit the model and visualize the results
results_jk = fit_model(data, N_jk = 2, p_jk=0.9) # Fits the model, using 1000 Jack-Knife iterations with 90% resampling 
# takes ~ 18 minutes on 8 processors including one risk factor
results_summarized, parameters_logistic_model, parameters_lyman_model = summarize_fit(results_jk, coverage = 0.68, plot_histograms = True)

# Summarize jackknife results and visualize final parameter estimates with histograms of the jackknife distributions
print(parameters_logistic_model)
print("")
print(parameters_lyman_model)
res_LRT = likelihood_ratio_test(model_H1=parameters_logistic_model, data=data, model_H0=None)

# Compute the likelihood ratio test against the intercept-only model

# Suppose you have the EUD-only model; then you can perform a likelihood ratio test for the TURP factor only
EUD_only_model_parameters = pd.DataFrame([{
    'n': 0.01,
    'a_b': 0.75,
    'beta_0': -6.5,
    'beta_EUD': 0.06
}])
res_LRT_turp = likelihood_ratio_test(model_H1=parameters_logistic_model, data=data, model_H0=EUD_only_model_parameters)

# Compute the simultaneous confidence intervals for parameters n, a_b and beta_EUD, with global coverage of 68%
sim_CI = simultaneous_CI(data, parameters_logistic_model = parameters_logistic_model, fit_parameters = ["n", "a_b", "beta_EUD"], coverage = 0.68)
sim_CI

# # Predictions and model calibration
predictions = predict(data, parameters_logistic_model)
# Compute the predicted probabilities
predictions

res_HL = hosmer_lemeshow_test(data['outcome'], predictions['predicted_prob'], g = 10, bins_as_quantiles=False)

# Goodness of fit Hosmer-Lemeshow test
risk_tree(data, parameters_logistic_model)

# Visualize the risk tree
plot_ntcp_stratified(data, parameters_logistic_model, coverage=0.68, min_patients=5, n_quantiles=2)

# Plot the NTCP curves stratifies dy risk factors
calibration_plot(data, parameters_logistic_model, coverage=0.68, min_pts = 5)

# Plot predicted vs observed probabilities

# Convergence Diagnostic
bounds_plot = {"TD50": (80, 120), "n": (0.01, 0.11), "a_b": (0.5, 3.5)}
results_LSP = likelihood_surface_plots(data, parameters_lyman_model, coverage=0.68, bounds=bounds_plot)

# Likelihood plots for TD50, n and a_b parameters (~ 40s on 8 processors)
best, df_refined, df_init, all_results_bf = check_convergence(data)

# Checks optimization robustness by exploring a broad grid and diverse starting points to confirm solution stability
print("\nn vs alpha/beta:")
plot_refined_vs_initial(df_init, df_refined, N=len(data['outcome']), x_col="n", y_col="a_b")

print("\nn vs beta_0:")
plot_refined_vs_initial(df_init, df_refined, N=len(data['outcome']), x_col="n", y_col="beta_EUD")

print("\nbeta_0 vs beta_EUD:")
plot_refined_vs_initial(df_init, df_refined, N=len(data['outcome']), x_col="beta_0", y_col="beta_EUD")

print(f"\nbeta_{data['risk_factors'].columns[0]} vs beta_EUD:")
plot_refined_vs_initial(df_init, df_refined, N=len(data['outcome']), x_col=f"beta_{data['risk_factors'].columns[0]}", y_col="beta_EUD")



