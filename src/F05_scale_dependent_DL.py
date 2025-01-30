#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 5 of the manuscript containing the scale-dependent 
longitudinal dispersion coefficient normalized by the bulk thermal diffusivity 
for different values of  the log-conductivity variance for scenario L1q2. 
Figure 5 also contains the asymptotic values for each variance based on the 
closed-form expression for asymptotic longitudinal macrodispersivity provided 
by Chang & Yeh (2012).

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import pandas as pd

### --- Define functions --- ###
################################

# closed-form expression for asymptotic longitudinal thermal dispersivity provided by Chang & Yeh (2012)
def a_asymptotic(L, var, q):
    P = c_f * q * L / lambda_b
    beta = var * L * (np.sqrt(np.pi)/2 + 16*np.sqrt(np.pi)/P**4 - 16/P**3 - 8/3*1/P 
                       - 16*np.sqrt(np.pi)/P**4* np.exp(P**2/4)* special.erfc(P/2))
    return beta

### --- Set path and parameters --- ###
#######################################

# path to load and plot results
path_to_file = ('../results/')

# first node for which BTC fitting is conducted
start_node = 21

# number of nodes in x-direction
n_nodes_x = 95

# x-coordinate 
x_coord = np.concatenate([
    np.arange(-200, -40, 20, dtype=np.float32),
    np.arange(-50, -20, 10, dtype=np.float32),
    np.arange(-25, 15, 5, dtype=np.float32),
    np.arange(20, 50, 10, dtype=np.float32),
    np.arange(60, 1520, 20, dtype=np.float32)
])

n = 0.25                                    # porosity (-)
lambda_b = 1.6369                           # bulk thermal conductivity (W/m*K)
lambda_f = 0.6                              # thermal conductivity of water (W/m*K)
lambda_s = (lambda_b - lambda_f*n)/(1-n)    # thermal conductivity of solid (W/m*K)
c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
c_b = (c_s*(1-n) + n*c_f)                   # volumetric heat capacity of aquifer (J/m^3*K)
D_diff = lambda_b / c_b                     # bulk thermal diffusivity (m^2/s)
q2 = 5.85e-8                                # mean Darcy velocity (m/s)

### --- Load and prepare data --- ###
#####################################

# initialize arrays for loading D_L and alpha_L
DL_arr = np.zeros((4, n_nodes_x - start_node))
aL_arr = np.zeros((4, n_nodes_x - start_node))
for i in range(4):
    DL = pd.read_csv(path_to_file + 'L1q2/var{}/D_L.csv'.format(i))
    DL_arr[i, :] = DL['D_L']
    aL_arr[i, :] = DL['a_L']
    
# normalize D_L
DL_norm = DL_arr/D_diff    
    
# calculate asymptotic D_L/D_diff based on asymptotic macrodispersivity
DL_norm_inf_var1 = (a_asymptotic(100, 1, q2)*q2 + D_diff)/D_diff
DL_norm_inf_var2 = (a_asymptotic(100, 2, q2)*q2 + D_diff)/D_diff
DL_norm_inf_var3 = (a_asymptotic(100, 3, q2)*q2 + D_diff)/D_diff
    
### --- Plot Figure 5 --- ###
#############################    
    
fig = plt.figure(figsize=(7,5))
scatter1 = plt.scatter(x_coord[start_node::2], DL_norm[0,::2], label='hom.', marker = 'x', color='black', zorder=2)
scatter2 = plt.scatter(x_coord[start_node::2], DL_norm[1,::2], label='$\sigma_{\ln K}^2$ = 1', facecolor='None', edgecolor='red', zorder=2)
scatter3 = plt.scatter(x_coord[start_node::2], DL_norm[2,::2], label='$\sigma_{\ln K}^2$ = 2', facecolor='None', edgecolor='green', zorder=2)
scatter4 = plt.scatter(x_coord[start_node::2], DL_norm[3,::2], label='$\sigma_{\ln K}^2$ = 3', facecolor='None', edgecolor='blue', zorder=2)

plt.hlines(DL_norm_inf_var1, 0, 1400, color='red', linestyles = '--', alpha=0.3)
plt.hlines(DL_norm_inf_var2, 0, 1400, color='green', linestyles = '--', alpha=0.3)
plt.hlines(DL_norm_inf_var3, 0, 1400, color='blue', linestyles = '--', alpha=0.3)
hlines4 = plt.hlines(0, 0, 0, color='black', linestyles = '--', label =r'$D_L^{\infty}/D_{\mathrm{diff}}$')

plt.hlines(1.0,0,1500, linestyles=(0,(1,1)), color='grey', linewidth=2, zorder=1)

plt.ylabel('$D_L / D_{\mathrm{diff}}$', fontsize=18)
plt.xlabel('Distance from heat source (m)', fontsize=18)
plt.text(980,3.2,'$Pe$ = 14.94', fontsize = 18) 
plt.text(770,22.8,'based on Chang&Yeh (2012)', fontsize = 12, zorder=6)  

scatter_handles = [scatter1, scatter2, scatter3, scatter4]
line_handles = [hlines4]
legend1 = plt.legend(handles=scatter_handles, loc='upper left', bbox_to_anchor=(-0.03, 1.005), ncol=1, frameon=False, fontsize=15)
legend2 = plt.legend(handles=line_handles, loc='upper left', bbox_to_anchor=(0.25, 1.005), ncol=1, frameon=False, fontsize=15)
plt.gca().add_artist(legend1)

plt.xlim(0,1400)
plt.ylim(0,25)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_paper/5_scale_dependent_DL.pdf', bbox_inches='tight')



