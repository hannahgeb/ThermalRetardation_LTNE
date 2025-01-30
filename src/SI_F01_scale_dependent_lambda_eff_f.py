#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S1 of the supporting information containing the
scale-dependent effective thermal conductivity of the fluid for different values
of the log-conductivity variance and a Peclet number of 14.94 (scenario L1q2).

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = ('../results/')

# node where heat source is located
source_node = 16

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

### --- Load and prepare data --- ###
#####################################  
    
# initialize arrays for results
lambda_f_arr = np.zeros((3,n_nodes_x-source_node)) 

# load effective thermal conductivity of the fluid   
for i in range(1,4,1):
    lambda_eff_f = pd.read_csv(path_to_file + 'L1q2/var{}/lambda_eff_f.csv'.format(i))
    lambda_f_arr[i-1, :] = lambda_eff_f['lambda_eff_f']
    
### --- Plot Figure S1 --- ###
##############################   

fig = plt.figure(figsize=(7,5))
plt.hlines(0.6,0,1400, linestyles='--', color='darkgrey', linewidth=2, zorder=1, label='hom.: $\lambda_{\mathrm{eff}, f} = 0.6$ Wm$^{-1}$K$^{-1}$')
plt.scatter(x_coord[source_node:], lambda_f_arr[0,:], label='$\sigma_{\ln K}^2$ = 1', facecolor='None', edgecolor='red', zorder=2)
plt.scatter(x_coord[source_node:], lambda_f_arr[1,:], label='$\sigma_{\ln K}^2$ = 2', facecolor='None', edgecolor='green', zorder=2)
plt.scatter(x_coord[source_node:], lambda_f_arr[2,:], label='$\sigma_{\ln K}^2$ = 3', facecolor='None', edgecolor='blue', zorder=2)
plt.ylabel(r'$\lambda_{\mathrm{eff}, f}$ (Wm$^{-1}$K$^{-1}$)', fontsize=18)
plt.xlabel('Distance from heat source (m)', fontsize=18)
plt.xlim(0,1400)
plt.ylim(-5,150)
legend = plt.legend(frameon=False, ncol=1, loc='upper left', fontsize=15)  #bbox_to_anchor=(0.31, 0.47)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S1_scale_dependent_lambda_eff_f.pdf',  bbox_inches='tight')