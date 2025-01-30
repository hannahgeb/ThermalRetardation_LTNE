#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure S5 of the supporting information containing the steady-state
transverse plume extent (Fig. S5a) and the temporal evolution of the longitudinal
plume extent (Fig. S5b) for three heterogeneous cases compared to the homogeneous 
case with and without thermal macrodispersion (scenario L1q2).

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

# time steps for scenario L1q2
t_q2 = np.arange(1E8, 4.51E10, 1E8) /1e10

# number of time steps
n_t = 450

# number of nodes
n_nodes_x = 95

### --- Load and prepare data --- ###
#####################################  

# initialize arrays for results
kyy_disp = np.zeros((n_nodes_x,3))
kxx_disp = np.zeros((n_t,3))
kyy = np.zeros((n_nodes_x,3,4))
kxx = np.zeros((n_t,3,4))

for i in range(1, 4, 1):

    # steady-state lateral plume extent for homogeneous cases with added dispersion (Figure S4a)
    k_yy_disp = pd.read_csv(path_to_file + 'L1q2/disp_var{}/kyy_x_sst.csv'.format(i))
    kyy_disp[:,i-1] = k_yy_disp['kyy'] / 1e3
    
    # longitudinal plume extent averaged in transverse direction for homogeneous cases with added dispersion (Figure S4b)
    k_xx_disp = pd.read_csv(path_to_file + 'L1q2/disp_var{}/kxx_t.csv'.format(i))
    kxx_disp[:,i-1] = k_xx_disp['kxx'] / 1e3
    
    for j in range(0, 4, 1):
        # steady-state lateral plume extent (Figure S1)
        k_yy = pd.read_csv(path_to_file + 'L1q{}/var{}/kyy_x_sst.csv'.format(i,j))
        kyy[:,i-1,j] = k_yy['kyy'] / 1e3
        
        # longitudinal plume extent averaged in transverse direction for heterogeneous cases of L = 100 m (Figure S3)
        k_xx = pd.read_csv(path_to_file + 'L1q{}/var{}/kxx_t.csv'.format(i,j))
        kxx[:,i-1,j] = k_xx['kxx'] / 1e3
      
### --- Plot Figure S5 --- ###
##############################              

plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(14,5))
gs = fig.add_gridspec(1, 2, hspace=0.1, wspace=0.25)
(ax1, ax2) = gs.subplots()

ax1.scatter(x_coord[::2],kyy[::2,1,0], label='hom.',  marker = 'x', color='black')
ax1.scatter(x_coord[::2],kyy_disp[::2,0], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 1$', marker = 'x', color='red', zorder=2)
ax1.scatter(x_coord[::2],kyy_disp[::2,1], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 2$', marker = 'x', color='green', zorder=2)
ax1.scatter(x_coord[::2],kyy_disp[::2,2], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 3$', marker = 'x', color='blue', zorder=2)

ax1.scatter(x_coord[::2],kyy[::2,1,1], label='$\sigma_{\ln K}^2$ = 1', facecolors='None', edgecolors='red')
ax1.scatter(x_coord[::2],kyy[::2,1,2], label='$\sigma_{\ln K}^2$ = 2', facecolors='None', edgecolors='green')
ax1.scatter(x_coord[::2],kyy[::2,1,3], label='$\sigma_{\ln K}^2$ = 3', facecolors='None', edgecolors='blue')
ax1.text(900,2.5,'$Pe$ = 14.94', fontsize = 18)
ax1.set_xlim(0,1400)
ax1.set_ylim(0,40)
ax1.set_ylabel('$\kappa^{\mathrm{eff}}_{yy}$ (10$^3$ m$^2$)', fontsize=18)
ax1.set_xlabel('Distance from heat source (m)', fontsize=18)
ax1.legend(ncol=1,loc='upper left', frameon=False, fontsize=14, bbox_to_anchor=(-0.03, 1.03))
ax1.text(10,41.5,'a)', weight='bold', fontsize=20)

ax2.scatter(t_q2[::4],kxx[::4,1,0], label='hom.',  marker = 'x', color='black')
ax2.scatter(t_q2[::4],kxx_disp[::4,0], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 1$', marker = 'x', color='red', zorder=2)
ax2.scatter(t_q2[::4],kxx_disp[::4,1], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 2$', marker = 'x', color='green', zorder=2)
ax2.scatter(t_q2[::4],kxx_disp[::4,2], label=r'hom. + macrodisp., $\sigma_{\ln K}^2 = 3$', marker = 'x', color='blue', zorder=2)
ax2.scatter(t_q2[::4],kxx[::4,1,1],  facecolors='None', edgecolors='red', zorder=1)
ax2.scatter(t_q2[::4],kxx[::4,1,2],  facecolors='None', edgecolors='green', zorder=1)
ax2.scatter(t_q2[::4],kxx[::4,1,3],  facecolors='None', edgecolors='blue', zorder=1)
ax2.set_xlabel('Time (10$^{10}$ s)', fontsize=18)
ax2.set_xlim(0,3)
ax2.set_ylabel('$\kappa^{\mathrm{eff}}_{xx}$ (10$^3$ m$^2$)', fontsize=18)
ax2.set_ylim(0,90)
ax2.text(0.05,93,'b)', weight='bold', fontsize=20)

# save figure as PDF
fig.savefig(path_to_file + 'Figures_SI/S5_compare_moments.pdf', bbox_inches='tight')

