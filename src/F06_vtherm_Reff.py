#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 6 of the manuscript containing the temporal evolution 
of the thermal velocity (Fig. 6a) and effective thermal retardation (R_eff) 
(Fig. 6b) for different values of the log-conductivity variance (scenario L1q2) 
and the apparent thermal retardation factor (R_app).

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# path to load results and save figures
path_to_file = '../results/'

# define time steps for scenario L1q2
t_q2 = np.arange(1E8, 4.51E10, 1E8) / 1e10

# number of time steps
n_t = 450  

# Define parameter values
lambda_b = 1.6369               # Bulk thermal conductivity (W/m*K)
c_f = 4.18e6                    # Volumetric heat capacity of water (J/m^3*K)
c_s = 880 * 2650                # Volumetric heat capacity of the solid (J/m^3*K)
n = 0.25                        # Porosity (-)
c_b = c_s * (1 - n) + c_f * n   # Volumetric heat capacity of the aquifer (J/m^3*K)
R_app = c_b / (n * c_f)         # Apparent thermal retardation factor (-)

### --- Load and prepare data --- ###
#####################################

# initialize arrays for results
thermal_velocity_avg = np.zeros((n_t, 3))   # Average thermal velocity for each variance
thermal_velocity_moe = np.zeros((n_t, 3))   # Margin of error of thermal velocity for each variance
Reff = np.zeros((n_t, 3, 4))                # Effective thermal retardation factors for scenario L1q2

# load results for all variances
for i in range(1, 4):
    # load thermal velocity for scenario L1q2
    vt = pd.read_csv(f"{path_to_file}L1q2/var{i}/v_therm.csv")
    thermal_velocity_avg[:, i - 1] = vt['vt_average'] / 1e-7
    thermal_velocity_moe[:, i - 1] = vt['vt_margin_of_error'] / 1e-7
    
    # load effective thermal retardation factors
    for j in range(4):
        R_eff = pd.read_csv(f"{path_to_file}L1q{i}/var{j}/Reff.csv")
        Reff[:, i - 1, j] = R_eff['R_eff']
        
# thermal velocity of homogeneous case for scenario L1q2
vt_var0 = pd.read_csv(f"{path_to_file}L1q2/var0/v_therm.csv")['vt'] / 1e-7

### --- Plot Figure 6 --- ###
############################# 

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(14, 5))
gs = fig.add_gridspec(1, 2, hspace=0.1, wspace=0.25)
(ax1, ax2) = gs.subplots()

# plot thermal velocity (Figure 6a)
ax1.scatter(t_q2[::2], vt_var0[::2], marker='x', color='black', label='hom.', zorder=6)

# loop through the variances and plot thermal velocity with margin of error
colors = ['red', 'green', 'blue']
for i, color in enumerate(colors):
    ax1.fill_between(t_q2[1::2], thermal_velocity_avg[1::2, i] + thermal_velocity_moe[1::2, i],
                     thermal_velocity_avg[1::2, i] - thermal_velocity_moe[1::2, i], alpha=0.08, color=color, zorder=2)
    ax1.scatter(t_q2[1::2], thermal_velocity_avg[1::2, i], facecolor='None', edgecolor=color)
    ax1.plot(t_q2[1::2], thermal_velocity_avg[1::2, i] + thermal_velocity_moe[1::2, i], color=color, linewidth=0.1, alpha=0.2, zorder=1)
    ax1.plot(t_q2[1::2], thermal_velocity_avg[1::2, i] - thermal_velocity_moe[1::2, i], color=color, linewidth=0.1, alpha=0.2, zorder=1)

ax1.set_ylabel('$v_{t}$ (10$^{-7}$ m/s)', fontsize=18)
ax1.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax1.set_ylim(0.7, 1.22)
ax1.set_xlim(0, 1.4)
ax1.text(0.1, 0.75, '$Pe$ = 14.94', fontsize=18)
ax1.text(0.02, 1.235, 'a)', weight='bold', fontsize=20)

# plot effective thermal retardation (Figure 6b)
ax2.scatter(t_q2[::2], Reff[::2, 1, 0], marker='x', color='black', label='hom.')
colors = ['red', 'green', 'blue']
labels = ['$\sigma_{\ln K}^2$ = 1', '$\sigma_{\ln K}^2$ = 2', '$\sigma_{\ln K}^2$ = 3']
for i, (color, label) in enumerate(zip(colors, labels)):
    ax2.scatter(t_q2[::2], Reff[::2, 1, i + 1], facecolor='None', edgecolor=color, label=label)

ax2.hlines(round(R_app, 2), xmin=0, xmax=1.1e10, linestyles='--', colors='black', label='$R_{\mathrm{app}}$ = 2.67')
ax2.vlines(x=0.1, ymin=2, ymax=3, alpha=0.4, color='grey', linewidth=7)
ax2.vlines(x=0.3, ymin=2, ymax=3, alpha=0.4, color='grey', linewidth=7)
ax2.vlines(x=0.5, ymin=2, ymax=3, alpha=0.4, color='grey', linewidth=7)

# annotations for time markers
ax2.text(0.03, 2.03, 't$_1$', fontsize=15)
ax2.text(0.23, 2.03, 't$_2$', fontsize=15)
ax2.text(0.43, 2.03, 't$_3$', fontsize=15)

ax2.set_ylabel('$R_{\mathrm{eff}}$', fontsize=18)
ax2.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax2.set_ylim(2, 3)
ax2.set_xlim(0, 1.4)
ax2.legend(loc='lower right', frameon=False)
ax2.text(0.02, 3.03, 'b)', weight='bold', fontsize=20)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_paper/6_vt_Reff.pdf',  bbox_inches='tight', dpi=200)

