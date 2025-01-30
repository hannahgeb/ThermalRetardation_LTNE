#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 8 of the manuscript containing temporal evolution of 
the effective thermal retardation (R_eff) for the three heterogeneous cases compared 
to the homogeneous case with and without thermal macrodispersion (scenario L1q2) 
and the apparent thermal retardation factor (R_app).

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = '../results/' 

# time steps for scenario L1q2
t_q1 = np.arange(3E8, 1.3501E11, 3E8) /1e10
t_q2 = np.arange(1E8, 4.51E10, 1E8) /1e10

# number of time steps
n_t = 450

# parameter values
lambda_b = 1.6369                           # bulk thermal conductivity (W/m*K)
c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
n = 0.25                                    # porosity (-)
c_b = c_s*(1-n) + c_f*n                     # volumetric heat capacity of the aquifer (J/m^3*K)
R_app = c_b/(n*c_f)                         # apparent thermal retardation factor (-)

### --- Load and prepare data --- ###
#####################################

# initialize arrays for results
Reff_disp = np.zeros((n_t,3))
Reff = np.zeros((n_t,3,4))

for i in range(1, 4, 1):
    # effective thermal retardation for homogeneous cases with added dispersion 
    R_eff = pd.read_csv(path_to_file + 'L1q2/disp_var{}/Reff.csv'.format(i))
    Reff_disp[:,i-1] = R_eff['R_eff']
    
    for j in range(0, 4, 1):
        # effective thermal retardation factors, scenarios L1q1, L1q2 and L1q3
        R_eff = pd.read_csv(path_to_file + 'L1q{}/var{}/Reff.csv'.format(i,j))
        Reff[:,i-1,j] = R_eff['R_eff']

### --- Plot Figure 7 --- ###
############################# 

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
((ax1, ax2, ax3)) = gs.subplots(sharey='row', sharex='col')

ax1.scatter(t_q2[::2],Reff[::2,1,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(t_q2[::2],Reff_disp[::2,0], label=r'hom. + macrodisp.', marker = 'x', color='red')
ax1.scatter(t_q2[::2],Reff[::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black', label='$R_{\mathrm{app}}$ ='+' {}'.format(round(R_app,2)))
ax1.text(0.95,2.1,'$Pe$ = 14.94', fontsize = 18)
ax1.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax1.legend(ncol=1,loc='upper left', frameon=False, fontsize=16)
ax1.xaxis.set_ticks(np.arange(0, 1.3, 0.2))
ax1.set_ylabel('$R_{\mathrm{eff}}$', fontsize=18)

ax2.scatter(t_q2[::2],Reff[::2,1,0], label='hom.',  marker = 'x', color='black')
ax2.scatter(t_q2[::2],Reff_disp[::2,1], label=r'hom. + macrodisp.', marker = 'x', color='green')
ax2.scatter(t_q2[::2],Reff[::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax2.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax2.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax2.legend(ncol=1,loc='upper left', frameon=False, fontsize=16)
ax2.xaxis.set_ticks(np.arange(0, 1.3, 0.2))

ax3.scatter(t_q2[::2],Reff[::2,1,0], label='hom.',  marker = 'x', color='black')
ax3.scatter(t_q2[::2],Reff_disp[::2,2], label=r'hom. + macrodisp.', marker = 'x', color='blue')
ax3.scatter(t_q2[::2],Reff[::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax3.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax3.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax3.legend(ncol=1,loc='upper left', frameon=False, fontsize=16)
ax3.xaxis.set_ticks(np.arange(0, 1.5, 0.2))

for ax in fig.get_axes():
    ax.label_outer()
    ax.set_xlim(0,1.4)
    ax.set_ylim(2,3.6)

# save figure as PDF
#fig.savefig(path_to_file + 'Figures_paper/8_Reff_add_disp.pdf', bbox_inches='tight')





