#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S3 of the supporting information containing the normalized
transverse dispersion coefficient as a function of the distance from the heat source
for different values of the log-conductivity variance and the thermal Peclet number.

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = '../results/' 

# x-coordinate 
x_coord = np.concatenate([
    np.arange(-200, -40, 20, dtype=np.float32),
    np.arange(-50, -20, 10, dtype=np.float32),
    np.arange(-25, 15, 5, dtype=np.float32),
    np.arange(20, 50, 10, dtype=np.float32),
    np.arange(60, 1520, 20, dtype=np.float32)])

# number of nodes
n_nodes_x = 95

# parameter values
lambda_b = 1.6369                           # bulk thermal conductivity (W/m*K)
c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
D = lambda_b / c_f

### --- Load and prepare data --- ###
#####################################  

# initialize arrays for results
D_T = np.zeros((n_nodes_x,3,4))

for i in range(1, 4, 1):
    for j in range(0, 4, 1):
        # transverse dispersion coefficient 
        DT = pd.read_csv(path_to_file + 'L1q{}/var{}/D_T.csv'.format(i,j))
        D_T[:,i-1,j] = DT['D_T']

# normalize D_T
DT_norm = D_T/D

### --- Plot Figure S3 --- ###
##############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
(ax1, ax2, ax3) = gs.subplots(sharex='col', sharey='row')

ax1.scatter(x_coord[20::2],DT_norm[20::2,0,0], label='hom.', marker='x', color='black', zorder=2)
ax1.scatter(x_coord[20::2],DT_norm[20::2,0,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red', zorder=2)
ax1.scatter(x_coord[20::2],DT_norm[20::2,0,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green', zorder=2)
ax1.scatter(x_coord[20::2],DT_norm[20::2,0,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue', zorder=2)
ax1.hlines(1.0,0,1400,linestyles=(0,(1,1)), color='grey', linewidth=1, zorder=1)
ax1.text(100,0.25,'$Pe$ = 4.98', fontsize = 18)
ax1.set_ylabel('$D_T/D^{\mathrm{SS}}_{\mathrm{diff}}$', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 1300, 200))
ax1.legend(ncol=1,loc='upper right', frameon=False, fontsize=16) 

ax2.scatter(x_coord[20::2],DT_norm[20::2,1,0], label='hom.', marker='x', color='black', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue', zorder=2)
ax2.hlines(1.0,0,1400, linestyles=(0,(1,1)), color='grey', linewidth=1, zorder=1)
ax2.text(100,0.25,'$Pe$ = 14.94', fontsize = 18) 
ax2.xaxis.set_ticks(np.arange(0, 1300, 200))

ax3.scatter(x_coord[20::2],DT_norm[20::2,2,0], label='hom.', marker='x', color='black', zorder=2)
ax3.scatter(x_coord[20::2],DT_norm[20::2,2,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red', zorder=2)
ax3.scatter(x_coord[20::2],DT_norm[20::2,2,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green', zorder=2)
ax3.scatter(x_coord[20::2],DT_norm[20::2,2,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue', zorder=2)
ax3.hlines(1.0,0,1400, linestyles=(0,(1,1)), color='grey', linewidth=1, zorder=1)
ax3.text(100,0.25,'$Pe$ = 25.03', fontsize = 18)
ax3.xaxis.set_ticks(np.arange(0, 1500, 200))

for ax in fig.get_axes():
    ax.label_outer()
    ax.set_xlim(0,1400)
    ax.set_ylim(0,3)
    ax.set_xlabel('Distance from heat source (m)', fontsize=18)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S3_DT.pdf',  bbox_inches='tight')








