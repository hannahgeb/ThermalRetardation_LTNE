#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S7 of the supporting information containing the effective 
thermal retardation (R_eff) as function of normalized travel time for different 
values of the log-conductivity variance and the thermal Peclet number (all scenarios
considered in this study). Figure S7 also shows the apparent thermal retardation
factor (R_app).

Author: H. Gebhardt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

# set path to load results and save figures
path_to_file = "../results/"

# number of variances and time steps
n_var, n_t = 3, 450
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

scenario_params = {
    "L1q1": (t_q1, q[0], L[0], e[0]),
    "L1q2": (t_q2, q[1], L[0], e[0]),
    "L1q3": (t_q3, q[2], L[0], e[0]),
    "L2q2": (t_q2, q[1], L[1], e[1]),
    "L3q2": (t_q2, q[1], L[2], e[2]),
}

### --- Load and prepare data --- ###
##################################### 

# load Reff data and calculate normalized time
Reff_scenarios = load_Reff_data(path_to_file, scenarios, n_var, n_t)
t_dl_scenarios = calculate_normalized_time(scenario_params, ratio_c)

### --- Plot Figure S7 --- ###
##############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19, 10))

ax1 = plt.subplot2grid((2, 3), (0, 0))
ax2 = plt.subplot2grid((2, 3), (0, 1))
ax3 = plt.subplot2grid((2, 3), (0, 2))
ax4 = plt.subplot2grid((2, 3), (1, 0))
ax5 = plt.subplot2grid((2, 3), (1, 1))

plt.subplots_adjust(hspace=0.05, wspace=0.05)

for ax, scenario, Pe_value in zip(
    [ax1, ax2, ax3, ax4, ax5],
    ["L1q1", "L1q2", "L2q2", "L3q2", "L1q3"],
    [4.98, 14.94, 18.67, 22.41, 25.03]
):
    t_dl = t_dl_scenarios[scenario]
    Reff = Reff_scenarios[scenario]

    for i, color in enumerate(['red', 'green', 'blue']):
        ax.scatter(t_dl, Reff[i, :], label=f'$\\sigma_{{\\ln K}}^2$ = {i+1}', facecolors='None', edgecolors=color)

    ax.hlines(R_app, xmin=0, xmax=10, linestyles='--', colors='black', label=f'$R_{{\\mathrm{{app}}}}$ = {round(R_app, 2)}')
    ax.text(0.25, 2.71, f'$Pe$ = {Pe_value}', fontsize=18)

    if ax in [ax1, ax4]:
        ax.set_ylabel('$R_{\\mathrm{eff}}$', fontsize=18)
        ax.yaxis.set_ticks(np.arange(2.0, 2.85, 0.1))

    if ax in [ax3, ax4, ax5]:
        ax.set_xlabel(r'$\frac{\rho_f c_f}{\rho_b c_b}$ $\frac{q_0 t}{L_x e}$', fontsize=22)

for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
for ax in [ax2, ax3, ax5]:
    plt.setp(ax.get_yticklabels(), visible=False)
for ax in fig.get_axes():
    ax.set_xlim(0, 10.0)
    ax.set_ylim(2, 2.8)
    
# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S7_Reff_normalized_time.pdf', bbox_inches='tight')






