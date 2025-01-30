#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 4 of the manuscript containing the steady-state 
transverse plume extents (Fig. 4a) and the corresponding transverse dispersion 
coefficients (Fig. 4b) for different log-conductivity values and scenario L1q2.


Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# path to load and plot results
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
lambda_b = 1.6369       # bulk thermal conductivity (W/m*K)
c_f = 4.18e6            # volumetric heat capacity of water (J/m^3*K)
D = lambda_b / c_f 
q2 = 5.85e-8            # mean Darcy velocity (m/s) for scenario L1q2, L2q2 and L3q2

# theoretical value of the lateral plume extent (kyy) for the homogeneous case 
# (Hidalgo et al. 2009)
kyy_theoret_q2 = ((2* D * x_coord) / q2) / 1e3

### --- Load and prepare data --- ###
#####################################

# initialize arrays for loading kyy and D_T
kyy = np.zeros((n_nodes_x,3,4))
D_T = np.zeros((n_nodes_x,3,4))

for i in range(1, 4, 1):
    for j in range(0, 4, 1):
        # steady-state lateral plume extent 
        k_yy = pd.read_csv(path_to_file + 'L1q{}/var{}/kyy_x_sst.csv'.format(i,j))
        kyy[:,i-1,j] = k_yy['kyy'] / 1e3
        
        # transverse dispersion coefficient 
        DT = pd.read_csv(path_to_file + 'L1q{}/var{}/D_T.csv'.format(i,j))
        D_T[:,i-1,j] = DT['D_T']
        
# normalize D_T
DT_norm = D_T/D

### --- Plot Figure 4 --- ###
#############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(14,5))
gs = fig.add_gridspec(1, 2, hspace=0.1, wspace=0.25)
(ax1, ax2) = gs.subplots()

ax1.plot(x_coord,kyy_theoret_q2, label=r'$2D_\mathrm{diff}^{\mathrm{SS}}x/q_0$', color='black', linewidth=1)
ax1.scatter(x_coord[::2],kyy[::2,1,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(x_coord[::2],kyy[::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.scatter(x_coord[::2],kyy[::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax1.scatter(x_coord[::2],kyy[::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax1.text(900,2.5,'$Pe$ = 14.94', fontsize = 18)
ax1.set_xlim(0,1400)
ax1.set_ylim(0,40)
ax1.set_ylabel('$\kappa^{\mathrm{eff}}_{yy}$ (10$^3$ m$^2$)', fontsize=18)
ax1.set_xlabel('Distance from heat source (m)', fontsize=18)
ax1.legend(ncol=1,loc='upper left', frameon=False, fontsize=16)
ax1.text(10,41.5,'a)', weight='bold', fontsize=20)

ax2.scatter(x_coord[20::2],DT_norm[20::2,1,0], label='hom.', marker='x', color='black', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green', zorder=2)
ax2.scatter(x_coord[20::2],DT_norm[20::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue', zorder=2)
ax2.hlines(1.0,0,1400, linestyles=(0,(1,1)), color='grey', linewidth=1, zorder=1)
ax2.text(10,3.1,'b)', weight='bold', fontsize=20)

ax2.set_ylabel('$D_T/D^{\mathrm{SS}}_{\mathrm{diff}}$', fontsize=18)
ax2.set_xlim(0,1400)
ax2.set_ylim(0,3)
ax2.set_xlabel('Distance from heat source (m)', fontsize=18)

# save figure as PDF
#fig.savefig(path_to_file + 'Figures_paper/4_kyy_DT.pdf', bbox_inches='tight')








