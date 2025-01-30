#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to perform a multi-variable linear regression based on the effective 
thermal retardation values (R_eff) determined for all scenarios considered in 
this study. The regression is performed for each normalized travel time step, 
where R_eff is treated as the response variable, and log-conductivity variance 
and thermal Peclet number serve as predictor variables. 

This script also determines the normalized time step, where the intercept of the 
regression curve equals the apparent thermal retardation factor (R_app), in order 
to derive a formula for estimating R_eff as a function of the log-conductivity 
variance and the Peclet number.
 
Author: H. Gebhardt
"""
import pandas as pd
import numpy as np
import statsmodels.api as sm

### --- Define functions --- ###
################################

# find nearest index
def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

# load Reff data
def load_Reff_data(path, scenarios, n_var, n_t):
    Reff_scenarios = {scenario: np.zeros((n_var, n_t)) for scenario in scenarios}
    for scenario in scenarios:
        for i in range(1, 4):
            Reff = pd.read_csv(f"{path}{scenario}/var{i}/Reff.csv")
            Reff_scenarios[scenario][i-1, :] = Reff['R_eff']
    return Reff_scenarios

# calculate normalized time
def calculate_normalized_time(params, ratio_c):
    return {
        scenario: t_q * q * ratio_c / (L * l)
        for scenario, (t_q, q, L, l) in params.items()
    }

### --- Set path and parameters --- ###
#######################################

# initialize output directory
path_to_file = "../results/"

# number of variances and time steps
n_var, n_t = 3, 450

# scenarios of this study
scenarios = ["L1q1", "L1q2", "L1q3", "L2q2", "L3q2"]

# time steps for each considered mean Darcy velocity q_0
t_q1, t_q2, t_q3 = np.arange(3e8, 1.3501e11, 3e8), np.arange(1e8, 4.51e10, 1e8), np.arange(6e7, 2.706e10, 6e7)

c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
n = 0.25                                    # porosity (-)
c_b = c_s*(1-n) + c_f*n                     # volumetric heat capacity of the aquifer (J/m^3*K)
ratio_c = c_f / c_b
R_app = c_b/(n*c_f)                         # apparent thermal retardation factor (-)

L = [100,125,150]                           # list containing the correlation lengths L_x (m) considered in this study
q = [1.95e-8,5.85e-8,9.8e-8]                # list containing the mean Darcy velocities q (m/s) considered in this study 
e = [0.5,0.4,1/3]                           # list containing the anisotropy ratios e = L_z/L_x considered in this study
Pe = [4.98, 14.94, 25.03, 18.67, 22.41]     # Peclet numbers considered in this study
var = [1, 2, 3]                             # variances of the ln K field considered in this study

# scenario parameter
scenario_params = {
    "L1q1": (t_q1, q[0], L[0], e[0]),
    "L1q2": (t_q2, q[1], L[0], e[0]),
    "L1q3": (t_q3, q[2], L[0], e[0]),
    "L2q2": (t_q2, q[1], L[1], e[1]),
    "L3q2": (t_q2, q[1], L[2], e[2]),
}

### --- Load data and calculate normalized time --- ###
#######################################################

Reff_scenarios = load_Reff_data(path_to_file, scenarios, n_var, n_t)
t_dl_scenarios = calculate_normalized_time(scenario_params, ratio_c)

### --- Perform multi-variable linear regression --- ###
########################################################
t_dl = np.arange(0.5, 10, 0.1)
results = []
for t in t_dl:
    indices = {scenario: find_nearest(t_dl_scenarios[scenario], t) for scenario in scenarios}
    Reff_values = [
        Reff_scenarios[scenario][var_idx, indices[scenario]]
        for var_idx in range(3)
        for scenario in scenarios
    ]
    df = pd.DataFrame({'var': np.repeat(var, len(scenarios)), 'Pe': np.tile(Pe, 3), 'Reff': Reff_values})
    X, y = sm.add_constant(df[['var', 'Pe']]), df['Reff']
    model = sm.OLS(y, X).fit()
    results.append([t, model.params[0], model.params[1], model.params[2], model.rsquared, model.rsquared_adj, np.sqrt(np.mean(model.resid**2))])

### --- Save regression results --- ###
#######################################

columns = ['normalized_time', 'intercept', 'var_coeff', 'Pe_coeff', 'R_squared', 'R_squared_adjusted', 'RMSE']
lin_reg = pd.DataFrame(results, columns=columns)
lin_reg.to_csv(path_to_file + 'linear_regression.csv', index=False)

### --- Find index for which the intercept is equal to R_app and save results --- ###
#####################################################################################

i = find_nearest(lin_reg['intercept'], round(R_app,2)) 
indices = {scenario: find_nearest(t_dl_scenarios[scenario], t_dl[i]) for scenario in scenarios}
print(i)

Reff_values = [
    Reff_scenarios[scenario][var_idx, indices[scenario]]
    for var_idx in range(3)
    for scenario in scenarios
]
lin_reg_Rapp = pd.DataFrame({
    'var': np.repeat(var, len(scenarios)), 
    'Pe': np.tile(Pe, 3), 
    'Reff': Reff_values, 
    'normalized_time': np.repeat(round(t_dl[i],1), len(scenarios)*len(var)),
    'norm_time_idx': np.repeat(i, len(scenarios)*len(var))}) 

# save R_eff values from direct simulations for which the intercept equals R_app 
lin_reg_Rapp.to_csv(path_to_file + 'lin_reg_identify_Rapp.csv', index=False)






