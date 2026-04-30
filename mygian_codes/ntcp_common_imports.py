OPTIMIZATION_ALGORITHM = 'L-BFGS-B'
INITIAL_GUESSES_BRUTE_FORCE = 5
# REFINEMENT_OPTIMIZATION_OPTIONS = {'maxiter': 1000, 'ftol': 1e-6, 'gtol': 1e-8}
REFINEMENT_OPTIMIZATION_OPTIONS = None 
# Default options: {'maxiter': 15000, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05}

# maxiter: Maximum number of algorithm iterations.
# ftol: The iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol
# gtol: The iteration will stop when max{|proj g_i | i = 1, ..., n} <= gtol where proj g_i is the i-th component of the projected gradient.

# region Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
from joblib import Parallel, delayed    # Necessary for parallelization
import networkx as nx                   # Only for classification tree
from copy import deepcopy
from itertools import product
from matplotlib import colormaps
import os
from scipy.optimize import minimize
from scipy.stats import (
    chi2,
    norm)
from sklearn.metrics import (
    roc_auc_score,
    precision_score,
    recall_score,
    f1_score,
    average_precision_score,
    accuracy_score,
    confusion_matrix
)
