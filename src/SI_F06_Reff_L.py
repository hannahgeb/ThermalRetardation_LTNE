#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S6 of the supporting information containing the temporal
evolution of the effective thermal retardation (R_eff) for varying correlation
lengths, a constant Darcy velocity (scenarios L1q2 - L3q2) and different values 
of the log-conductivity variance. Figure S6 also shows the apparent thermal 
retardation factor (R_app).

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = '../results/' 

# time steps for scenario L1q2, L2q2 and L3q2
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
Reff = np.zeros((n_t,3,4)) 

for i in range(1, 4, 1):
    for j in range(0, 4, 1):
        # effective thermal retardation factors (Figure S5), scenarios L1q2, L2q2 and L3q2 
        R_eff = pd.read_csv(path_to_file + 'L{}q2/var{}/Reff.csv'.format(i,j))
        Reff[:,i-1,j] = R_eff['R_eff']
        
### --- Plot Figure S6 --- ###
##############################

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(19,5))
gs = fig.add_gridspec(1, 3, hspace=0.1, wspace=0.03)
((ax4, ax5, ax6)) = gs.subplots(sharex='col', sharey='row')

ax4.scatter(t_q2[::2],Reff[::2,0,0], label='hom.',  marker = 'x', color='black')
ax4.scatter(t_q2[::2],Reff[::2,0,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax4.scatter(t_q2[::2],Reff[::2,0,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax4.scatter(t_q2[::2],Reff[::2,0,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax4.text(0.95,2.1,'$Pe$ = 14.94', fontsize = 18)
ax4.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black', label='$R_{\mathrm{app}}$ ='+' {}'.format(round(R_app,2))) 
ax4.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax4.set_ylabel('$R_{\mathrm{eff}}$', fontsize=18)
ax4.legend(ncol=1,loc='upper left', frameon=False, fontsize=16, bbox_to_anchor=(-0.01, 1.02))
ax4.xaxis.set_ticks(np.arange(0, 1.3, 0.2))

ax5.scatter(t_q2[::2],Reff[::2,1,0], label='hom.',  marker = 'x', color='black')
ax5.scatter(t_q2[::2],Reff[::2,1,1], label='$\sigma_{\ln K}^2$ = 1',  facecolors='None', edgecolors='red')
ax5.scatter(t_q2[::2],Reff[::2,1,2], label='$\sigma_{\ln K}^2$ = 2',  facecolors='None', edgecolors='green')
ax5.scatter(t_q2[::2],Reff[::2,1,3], label='$\sigma_{\ln K}^2$ = 3',  facecolors='None', edgecolors='blue')
ax5.text(0.95,2.1,'$Pe$ = 18.67', fontsize = 18)
ax5.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax5.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax5.xaxis.set_ticks(np.arange(0, 1.3, 0.2))

ax6.scatter(t_q2[::2],Reff[::2,2,0], label='hom.',  marker = 'x', color='black')
ax6.scatter(t_q2[::2],Reff[::2,2,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax6.scatter(t_q2[::2],Reff[::2,2,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax6.scatter(t_q2[::2],Reff[::2,2,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax6.text(0.95,2.1,'$Pe$ = 22.41', fontsize = 18)
ax6.hlines(round(R_app,2),xmin=0,xmax=2,linestyles='--',colors='black')
ax6.set_xlabel('Time (10$^{10}$ s)', fontsize=18)

for ax in fig.get_axes():
    ax.label_outer()
    ax.set_xlim(0,1.4)
    ax.set_ylim(2,3.6)

# save figure as PDF    
fig.savefig(path_to_file + 'Figures_SI/S6_Reff_L.pdf',  bbox_inches='tight')







