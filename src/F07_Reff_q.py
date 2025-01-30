#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 7 of the manuscript containing the temporal evolution
of the effective thermal retardation (R_eff) for three Darcy velocities (scenario 
L1q1 - L1q3), constant correlation length of L=100 m for all cases and three 
values of the log-conductivity variance. Figure 7 also contains the apparent 
thermal retardation factor.

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

# time steps for each considered mean Darcy velocity q_0
t_q1 = np.arange(3E8, 1.3501E11, 3E8) /1e10
t_q2 = np.arange(1E8, 4.51E10, 1E8) /1e10
t_q3 = np.arange(6e7, 2.706e10, 6e7) /1e10

# number of time steps
n_t = 450

lambda_b = 1.6369                           # bulk thermal conductivity (W/m*K)
c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
n = 0.25                                    # porosity (-)
c_b = c_s*(1-n) + c_f*n                     # volumetric heat capacity of the aquifer (J/m^3*K)
R_app = c_b/(n*c_f)                         # apparent thermal retardation factor (-)

### --- Load and prepare data --- ###
#####################################

# initialize arrays for results
Reff = np.zeros((n_t,3,4))

for i in range(1, 4, 1):    
    for j in range(0, 4, 1):
        # effective thermal retardation factors, scenarios L1q1, L1q2 and L1q3
        R_eff = pd.read_csv(path_to_file + 'L1q{}/var{}/Reff.csv'.format(i,j))
        Reff[:,i-1,j] = R_eff['R_eff']

### --- Plot Figure 8 --- ###
#############################   

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
((ax1, ax2, ax3)) = gs.subplots(sharex='col', sharey='row')

ax1.scatter(t_q1[::],Reff[::,0,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(t_q1[::],Reff[::,0,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.scatter(t_q1[::],Reff[::,0,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax1.scatter(t_q1[::],Reff[::,0,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax1.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black', label='$R_{\mathrm{app}}$ =' + ' {}'.format(round(R_app,2)))
ax1.text(1.4,2.1,'$Pe$ = 4.98', fontsize = 18)
ax1.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax1.legend(ncol=1,loc='upper left', frameon=False, fontsize=16, bbox_to_anchor=(0.01, 1.03))
ax1.set_ylabel('$R_{\mathrm{eff}}$', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 1.7, 0.5))

ax2.scatter(t_q2[::2],Reff[::2,1,0], label='hom.',  marker = 'x', color='black')
ax2.scatter(t_q2[::2],Reff[::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax2.scatter(t_q2[::2],Reff[::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax2.scatter(t_q2[::2],Reff[::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax2.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax2.text(1.35,2.1,'$Pe$ = 14.94', fontsize = 18) 
ax2.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax2.xaxis.set_ticks(np.arange(0, 1.7, 0.5))

ax3.scatter(t_q3[::2],Reff[::2,2,0], label='hom.',  marker = 'x', color='black')
ax3.scatter(t_q3[::2],Reff[::2,2,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax3.scatter(t_q3[::2],Reff[::2,2,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax3.scatter(t_q3[::2],Reff[::2,2,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax3.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax3.text(1.35,2.1,'$Pe$ = 25.03', fontsize = 18)
ax3.set_xlabel('Time (10$^{10}$ s)', fontsize=18)

for ax in fig.get_axes():
    ax.label_outer()
    ax.set_xlim(0,2)
    ax.set_ylim(2,3.6)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_paper/7_Reff_q.pdf',  bbox_inches='tight')

    




