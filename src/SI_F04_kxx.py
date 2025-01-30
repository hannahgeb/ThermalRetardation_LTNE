#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S4 of the supporting information containing the temporal
evolution of the longitudinal plume extent for different values of the log-conductivity 
variance and the thermal Peclet number.

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = '../results/' 

# time steps for each considered mean Darcy velocity q_0
t_q1 = np.arange(3E8, 1.3501E11, 3E8) /1e10
t_q2 = np.arange(1E8, 4.51E10, 1E8) /1e10
t_q3 = np.arange(6e7, 2.706e10, 6e7) /1e10

# number of time steps
n_t = 450

### --- Load and prepare data --- ###
#####################################  

# initialize arrays for results
kxx = np.zeros((n_t,3,4))

for i in range(1, 4, 1):    
    for j in range(0, 4, 1):
        # longitudinal plume extent averaged in transverse direction for heterogeneous cases
        k_xx = pd.read_csv(path_to_file + 'L1q{}/var{}/kxx_t.csv'.format(i,j))
        kxx[:,i-1,j] = k_xx['kxx'] /1e3

### --- Plot Figure S4 --- ###
##############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
(ax1, ax2, ax3) = gs.subplots(sharex='col', sharey='row')

ax1.scatter(t_q1[::4],kxx[::4,0,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(t_q1[::4],kxx[::4,0,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.scatter(t_q1[::4],kxx[::4,0,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax1.scatter(t_q1[::4],kxx[::4,0,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax1.text(0.5,110,'$Pe$ = 4.98', fontsize = 18)
ax1.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax1.set_xlim(0,13.7)
ax1.set_ylim(0,120)
ax1.set_ylabel('$\kappa^{\mathrm{eff}}_{xx}$ (10$^3$ m$^2$)', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 13, 2.5))
ax1.legend(ncol=1,loc='lower right', frameon=False, fontsize=16)

ax2.scatter(t_q2[::4],kxx[::4,1,0], label='hom.',  marker = 'x', color='black')
ax2.scatter(t_q2[::4],kxx[::4,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax2.scatter(t_q2[::4],kxx[::4,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax2.scatter(t_q2[::4],kxx[::4,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax2.text(0.15,110,'$Pe$ = 14.94', fontsize = 18) 
ax2.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax2.set_xlim(0,4.57)
ax2.xaxis.set_ticks(np.arange(0, 5, 1))

ax3.scatter(t_q3[::4],kxx[::4,2,0], label='hom.',  marker = 'x', color='black')
ax3.scatter(t_q3[::4],kxx[::4,2,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax3.scatter(t_q3[::4],kxx[::4,2,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax3.scatter(t_q3[::4],kxx[::4,2,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax3.text(0.1,110,'$Pe$ = 25.03', fontsize = 18)
ax3.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax3.set_xlim(0,2.74)

for ax in fig.get_axes():
    ax.label_outer()

# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S4_kxx.pdf', bbox_inches='tight')


