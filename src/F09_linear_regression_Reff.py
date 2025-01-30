#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script reproducing Figure 9 of the manuscript containing the results of multi-variable
linear regression to derive an approximate formula for the effective thermal
retardation (R_eff) as a function of the log-conductivity variance and the thermal
Peclet number. Fig. 9a contains the identification of coefficients for a normalized
travel time where the intercept of the regression curve equals the apparent thermal
retardation (R_app), Fig. 9b shows the corresponding regression curves for each 
log-conductivity variance and R_eff from the direct simulations.

Author: H. Gebhardt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

### --- Set path and parameters --- ###
#######################################

# set path to load results and save figures
path_to_file = "../results/"

c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
n = 0.25                                    # porosity (-)
c_b = c_s*(1-n) + c_f*n                     # volumetric heat capacity of the aquifer (J/m^3*K)
ratio_c = c_f / c_b
R_app = c_b/(n*c_f)                         # apparent thermal retardation factor (-)

# Peclet numbers considered in this study
Pe = np.array([4.98, 14.94, 25.03, 18.67, 22.41])  

# variances considered in this study    
var = np.array([1, 2, 3])    
    
### --- Load and prepare data --- ###
#####################################

# load regression results
lin_reg = pd.read_csv(path_to_file + 'linear_regression.csv')

t_dl = lin_reg['normalized_time']
interc_arr = lin_reg['intercept']
Pe_coeff_arr = lin_reg['Pe_coeff']
var_coeff_arr = lin_reg['var_coeff']
R_squared_arr = lin_reg['R_squared']

lin_reg_Rapp = pd.read_csv(path_to_file + 'lin_reg_identify_Rapp.csv')
t_dl_Rapp = lin_reg_Rapp['normalized_time'][0]
t_dl_Rapp_idx = lin_reg_Rapp['norm_time_idx'][0]
y = lin_reg_Rapp['Reff']

# calculate regression curves
reg_curve_var1 = interc_arr[t_dl_Rapp_idx]+Pe_coeff_arr[t_dl_Rapp_idx]*Pe+var_coeff_arr[t_dl_Rapp_idx]*var[0]
reg_curve_var2 = interc_arr[t_dl_Rapp_idx]+Pe_coeff_arr[t_dl_Rapp_idx]*Pe+var_coeff_arr[t_dl_Rapp_idx]*var[1]
reg_curve_var3 = interc_arr[t_dl_Rapp_idx]+Pe_coeff_arr[t_dl_Rapp_idx]*Pe+var_coeff_arr[t_dl_Rapp_idx]*var[2]

### --- Plot Figure 9 --- ###
############################# 

plt.rcParams.update({'font.size': 14})

fig = plt.figure()
fig.set_figheight(9)
fig.set_figwidth(17)

ax1 = plt.subplot2grid(shape=(6, 6), loc=(0, 0),rowspan=2, colspan=2)
ax2 = plt.subplot2grid(shape=(6, 6), loc=(2, 0),rowspan=2, colspan=2, sharex=ax1)
ax3 = plt.subplot2grid(shape=(6, 6), loc=(4, 0),rowspan=2, colspan=2, sharex=ax1)
ax4 = plt.subplot2grid((6, 6), (1, 2), rowspan=4,colspan=3)

plt.subplots_adjust(hspace=0.1,wspace=0.9)
 
ax1.scatter(t_dl[:-2:4],interc_arr[:-2:4], color = 'black', zorder=2)
ax1.set_ylabel(r'$\beta_0$', fontsize=18, labelpad=34)
ax1.hlines(y=R_app, xmin=0, xmax=t_dl_Rapp, color='grey', linewidth=1, linestyles ='dashed', zorder=1)
ax1.vlines(x=t_dl_Rapp, ymin=0, ymax=R_app, color='grey', linewidth=1, linestyles ='dashed', zorder=1)
ax1.set_ylim(2.5,2.77)
ax1.yaxis.set_ticks(np.arange(2.5, 2.77, 0.05))
ax1.set_xlim(0,10)
ax1.text(1, 2.79, r'$R_{\mathrm{eff}}$ = $\beta_0$ - $\beta_1 \cdot \sigma_{\ln K}^2$ - $\beta_2 \cdot Pe$', weight='bold', fontsize=16)
ax1.text(0.2, 2.679, r'$\beta_0$(6.6) = 2.67 = $R_{\mathrm{app}}$ ', color='black',  fontsize=14)
ax1.text(-0.1, 2.8, 'a)', weight='bold', fontsize=20)

ax2.scatter(t_dl[:-2:4],var_coeff_arr[:-2:4], color = 'black', zorder=2)
ax2.set_ylabel(r'$\beta_1$', fontsize=18, labelpad=12)
ax2.hlines(y=-0.071, xmin=0, xmax=t_dl_Rapp, color='grey', linestyles='dashed', linewidth=1, zorder=1)
ax2.vlines(x=t_dl_Rapp, ymin=-0.15, ymax=var_coeff_arr[t_dl_Rapp_idx], color='grey', linestyles='dashed', linewidth=1, zorder=1)
ax2.set_ylim(-0.15,0.025)
ax2.yaxis.set_ticks(np.arange(-0.15, 0.03, 0.025))
ax2.set_xlim(0,10)
ax2.text(0.2, -0.066, r'$\beta_1(6.6)$ = -0.071', color='black',  fontsize=14)

ax3.scatter(t_dl[:-2:4],Pe_coeff_arr[:-2:4], color = 'black', zorder=2)
ax3.set_ylabel(r'$\beta_2$', fontsize=18, labelpad=3)
ax3.set_xlabel(r'$\frac{\rho_f c_f}{\rho_b c_b}$ $\frac{q_0 t}{L_x e}$', fontsize=18)
ax3.hlines(y=-0.0056, xmin=0, xmax=t_dl_Rapp, color='grey', linewidth=1, linestyles ='dashed',zorder=1)
ax3.vlines(x=t_dl_Rapp, ymin=-0.0125, ymax=Pe_coeff_arr[t_dl_Rapp_idx], color='grey', linewidth=1, linestyles ='dashed',zorder=1)
ax3.set_ylim(-0.0125,0)
ax3.set_xlim(0,10)
ax3.text(0.2, -0.0052, r'$\beta_2(6.6)$ = -0.0056', color='black',  fontsize=14)

ax4.scatter(Pe, y[:5], color = 'red', label='Direct ($\sigma_{\ln K}^2$ = 1)')
ax4.scatter(Pe, y[5:10], color = 'green', label='Direct ($\sigma_{\ln K}^2$ = 2)')
ax4.scatter(Pe, y[10:15], color = 'blue', label='Direct ($\sigma_{\ln K}^2$ = 3)')
ax4.plot(Pe, reg_curve_var1, color='red')
ax4.plot(Pe, reg_curve_var2, color='green')
ax4.plot(Pe, reg_curve_var3, color='blue')   
ax4.legend(frameon=False, fontsize=16, loc='lower left')
ax4.set_ylabel('$R_{\mathrm{eff}}$', fontsize=18)
ax4.set_xlabel('$Pe$', fontsize=18)
ax4.text(20, 2.222, r'$\frac{\rho_f c_f}{\rho_b c_b}$ $\frac{q_0 t}{L_x e}$ = 6.6', color='black', fontsize=16)
ax4.text(11.7,2.59,'R$^2$ = {}'.format(round(R_squared_arr[t_dl_Rapp_idx],2)), fontsize=16)
ax4.text(11.7,2.62,r'$R_{\mathrm{eff}}$ = '+'{}{}'.format(round(interc_arr[t_dl_Rapp_idx],2),round(var_coeff_arr[t_dl_Rapp_idx],4)) 
         + '$\sigma_{\ln K}^2$'+'{}$Pe$'.format(round(Pe_coeff_arr[t_dl_Rapp_idx],4)),fontsize=16)
ax4.set_ylim(2.2,2.65)
ax4.xaxis.set_ticks(np.arange(5, 30, 5))
ax4.yaxis.set_ticks(np.arange(2.2, 2.66, 0.05))
ax4.set_xticklabels(ax4.get_xticks(), fontsize=16)
ax4.set_yticklabels(ax4.get_yticks(), fontsize=16)
ax4.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax4.text(4, 2.67, 'b)', weight='bold', fontsize=20)

for ax in [ax1, ax2]:
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_yticks(ax.get_yticks()[1:])

# save figure as PDF
fig.savefig(path_to_file + 'Figures_paper/9_linear_regression_Reff.pdf',  bbox_inches='tight')




