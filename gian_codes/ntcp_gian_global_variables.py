import numpy as np

# region Global Variables
def OPTIMIZATION_BOUNDS(n_risk):
    return [(0.01, 1), (0.5, 10), (-np.inf, np.inf), (0, np.inf)] + [(-np.inf, np.inf)] * n_risk
# n ∈ [0.01, 1], α/β ∈ [0.5, 10], beta_0 ∈ ℝ, beta_eud > 0, beta_risk ∈ ℝ ∀ risk factor

def GRID_OPTIMIZATION_BF(p, n_risk):
    bounds = OPTIMIZATION_BOUNDS(n_risk)

    n_min, n_max = bounds[0]
    ab_min, ab_max = bounds[1]
    beta0_min, beta0_max = bounds[2]

    step_n = 0.05
    step_ab = 0.5

    n_steps = int(round(n_max / step_n))
    ab_steps = int(round(ab_max / step_ab))

    n_vals = np.array([n_min] + [step_n * i for i in range(1, n_steps + 1)])
    a_b_vals = np.array([step_ab * i for i in range(1, ab_steps + 1)])

    logit_p = np.log(p / (1 - p))
    beta_0_vals = np.linspace(
        max(0.15 * logit_p, beta0_min),
        min(4.0 * logit_p, beta0_max),
        15
    )

    beta_EUD_vals = [0.01 * i for i in range(8)]
    beta_risk_vals = [[0] for _ in range(n_risk)]

    return {
        "n": n_vals,
        "a_b": a_b_vals,
        "beta_0": beta_0_vals,
        "beta_EUD": beta_EUD_vals,
        "beta_risk": beta_risk_vals
    }

def GRID_CHECK_CONVERGENCE(p, n_risk):
    bounds = OPTIMIZATION_BOUNDS(n_risk)

    n_min, n_max = bounds[0]
    ab_min, ab_max = bounds[1]

    step_n = 0.05
    step_ab = 0.5

    n_steps = int(round(n_max / step_n))
    ab_steps = int(round(ab_max / step_ab))

    n_vals = np.array([n_min] + [step_n * i for i in range(1, n_steps + 1)])
    a_b_vals = np.array([ab_min] + [step_ab * i for i in range(2, ab_steps + 1)])

    beta_0_vals = np.linspace(0, 10 * np.log(p / (1 - p)), 20)

    alpha_eud = np.linspace(0, 1, 20) ** 2
    beta_EUD_vals = np.log(1.0 + (1.5 - 1.0) * alpha_eud)

    beta_risk_vals = []
    for _ in range(n_risk):
        alpha_r = np.linspace(0, 1, 15) ** 2
        beta_risk_vals.append(np.log(0.3 + (3.0 - 0.3) * alpha_r))

    return {
        "n": n_vals,
        "a_b": a_b_vals,
        "beta_0": beta_0_vals,
        "beta_EUD": beta_EUD_vals,
        "beta_risk": beta_risk_vals
    }
