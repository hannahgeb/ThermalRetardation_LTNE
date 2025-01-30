#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S2 of the supporting information containing the 
steady-state transverse plume extent as a function of the distance from the heat 
source for different values of the log-conductivity variance and the thermal 
Peclet number (scenario L1q1 - L1q3), and the theoretical values for the 
homogeneous case.

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
q1 = 1.95e-8                                # mean Darcy velocity (m/s) for scenario L1q1
q2 = 5.85e-8                                # mean Darcy velocity (m/s) for scenario L1q2, L2q2 and L3q2
q3 = 9.80e-8                                # mean Darcy velocity (m/s) for scenario L1q3

# theoretical value of the lateral plume extent (kyy) for the homogeneous case
kyy_theoret_q1 = (2* D * x_coord) / q1
kyy_theoret_q2 = (2* D * x_coord) / q2
kyy_theoret_q3 = (2* D * x_coord) / q3

### --- Load and prepare data --- ###
#####################################  

# initialize arrays for results
kyy = np.zeros((n_nodes_x,3,4))

for i in range(1, 4, 1):    
    for j in range(0, 4, 1):
        # steady-state lateral plume extent
        k_yy = pd.read_csv(path_to_file + 'L1q{}/var{}/kyy_x_sst.csv'.format(i,j))
        kyy[:,i-1,j] = k_yy['kyy'] / 1e3

### --- Plot Figure S2 --- ###
##############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
(ax1, ax2, ax3) = gs.subplots(sharex='col', sharey='row')

ax1.plot(x_coord,kyy_theoret_q1/1e3, label=r'$2D_\mathrm{diff}^\mathrm{SS}x/q_0$', color='black', linewidth=1)
ax1.scatter(x_coord[::2],kyy[::2,0,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(x_coord[::2],kyy[::2,0,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.scatter(x_coord[::2],kyy[::2,0,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax1.scatter(x_coord[::2],kyy[::2,0,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax1.text(100,70,'$Pe$ = 4.98', fontsize = 18)
ax1.xaxis.set_ticks(np.arange(0, 1300, 200))
ax1.set_ylabel('$\kappa^{\mathrm{eff}}_{yy}$ (10$^3$ m$^2$)', fontsize=18)

ax2.plot(x_coord,kyy_theoret_q2/1e3, label=r'$2D_\mathrm{diff}^\mathrm{SS}x/q_0$', color='black', linewidth=1)
ax2.scatter(x_coord[::2],kyy[::2,1,0], label='hom.',  marker = 'x', color='black')
ax2.scatter(x_coord[::2],kyy[::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax2.scatter(x_coord[::2],kyy[::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax2.scatter(x_coord[::2],kyy[::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax2.text(100,70,'$Pe$ = 14.94', fontsize = 18) 
ax2.xaxis.set_ticks(np.arange(0, 1300, 200))

ax3.plot(x_coord,kyy_theoret_q3/1e3, label=r'$2D_\mathrm{diff}^\mathrm{SS}x/q_0$', color='black', linewidth=1)
ax3.scatter(x_coord[::2],kyy[::2,2,0], label='hom.',  marker = 'x', color='black')
ax3.scatter(x_coord[::2],kyy[::2,2,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax3.scatter(x_coord[::2],kyy[::2,2,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax3.scatter(x_coord[::2],kyy[::2,2,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax3.text(100,70,'$Pe$ = 25.03', fontsize = 18)
ax3.xaxis.set_ticks(np.arange(0, 1500, 200))
ax3.legend(ncol=1,loc='upper right', frameon=False, fontsize=16)

for ax in fig.get_axes():
    ax.label_outer()
    ax.set_xlim(0,1400)
    ax.set_ylim(0,80)
    ax.set_xlabel('Distance from heat source (m)', fontsize=18)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S2_kyy.pdf', bbox_inches='tight')







