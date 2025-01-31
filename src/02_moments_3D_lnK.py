#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script calculating and saving spatial moments based on the temperature fields 
(numpy arrays) resulting from MOOSE simulations ('01_run_sim.py' has to be 
executed before).

NOTE:
This script will not run for general users as the temperature fields (numpy arrays) 
are not provided in the repository.

Author: H. Gebhardt
"""

import pandas as pd
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from scipy.stats import norm

### --- Set path and parameters --- ###
#######################################

# set scenario and variance
var = 3     # variance of ln K field
L_idx = 1   # index of considered correlation length L_x = [100, 125, 150]
q_idx = 2   # index of mean Darcy velocity q_0 = [1.95e-8, 5.85e-8, 9.8e-8]

# set path to load temperature fields and to save results from post-processing
dT_arrays_path = '../L1q2/var{}/dT_arrays/'.format(var) 

save_path = ('../results/')
path_to_save_results = (save_path + 'L{}q{}/var{}/'.format(L_idx, q_idx, var))

# number of realizations
n_real = 1  #100 when heterogeneous realizations are loaded

# number of nodes
n_nodes_x, n_nodes_y, n_nodes_z = 95, 79, 26

# end node for integration in x direction, exclude last 100 m from integration
end = -5     

# number of time steps
time_steps = 450

# time steps for each scenario
time_arr = [
    np.arange(3E8, 1.3501E11, 3E8),
    np.arange(1E8, 4.51E10, 1E8),
    np.arange(6e7, 2.706e10, 6e7)
]

# set time step based on considered scenario
t = time_arr[q_idx-1]

# parameters
q_arr = [1.95e-8, 5.85e-8, 9.8e-8]              # mean Darcy velocity q_0 (m/s)
q = q_arr[q_idx-1]                             
n = 0.25                                        # porosity (-)
v_a = q / n                                     # seepage velocity (m/s)
c_f = 4.18e6                                    # volumetric heat capacity of water (J/m^3*K)
lambda_b = 1.6369                               # bulk thermal conductivity (W/m*K)
D = lambda_b / c_f

### --- Discretization and weights for averaging --- ###
########################################################

# x-coordinates
x_coord = np.concatenate([
    np.arange(-200, -40, 20, dtype=np.float32),
    np.arange(-50, -20, 10, dtype=np.float32),
    np.arange(-25, 15, 5, dtype=np.float32),
    np.arange(20, 50, 10, dtype=np.float32),
    np.arange(60, 1520, 20, dtype=np.float32)])
x_coord = np.tile(x_coord, (n_nodes_z,1)) 

# y-coordinates 
y_coord = np.concatenate([
    np.arange(-650,-90,20,dtype=np.float32), 
    np.arange(-100,0,10,dtype=np.float32), 
    np.arange(-5,15,5,dtype=np.float32), 
    np.arange(20,120,10,dtype=np.float32), 
    np.arange(130,670,20,dtype=np.float32)])
y_coord = np.tile(y_coord, (n_nodes_z,1))

# array of weights for averaging in y-direction
width_y = 1300
w_y = np.concatenate([
    np.full(1, 10 / width_y),
    np.full(26, 20 / width_y),
    np.full(1, 15 / width_y),
    np.full(9, 10 / width_y),
    np.full(1, 7.5 / width_y),
    np.full(3, 5 / width_y),
    np.full(1, 7.5 / width_y),
    np.full(9, 10 / width_y),
    np.full(1, 15 / width_y),
    np.full(26, 20 / width_y),
    np.full(1, 10 / width_y)
])

### --- Main loop for spatial moment analysis --- ###
#####################################################

# initialize storage for k_yy, k_xx, and v_t
k_yy_n = []         # transverse plume extent
k_xx_n = []         # longitudinal plume extent, averaged in transverse direction
v_therm = []        # thermal velocity

for k in range(n_real): 
    # load temperature fields as arrays
    dT_file = (dT_arrays_path + 'dT_arr_{}.npy'.format(k))
    dT_arr = np.load(dT_file)
    dT_arr = np.where(dT_arr < 1e-4,0,dT_arr) # eliminate numerical artifacts

    # create arrays to store transverse and longitudinal spatial moments
    # transverse direction
    k_yy = np.zeros((time_steps,n_nodes_z,n_nodes_x))
    # longitudinal direction 
    k_xx = np.zeros((time_steps,n_nodes_z,n_nodes_y))
    zeroth_x, first_x = np.zeros_like(k_xx), np.zeros_like(k_xx)

    for i in range(time_steps):
        for j in range(n_nodes_x):
            
            dT = dT_arr[:,:,j,i]
            zeroth = spi.trapezoid(dT, x=y_coord, axis=1)
            dTy = dT * y_coord
            first = spi.trapezoid(dTy, x=y_coord, axis=1)
            dTyy = dT * y_coord * y_coord
            second = spi.trapezoid(dTyy, x=y_coord, axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                kyy = np.where(zeroth > 0, (second / zeroth) - (first / zeroth)**2, 0)
            k_yy[i,:,j] = kyy
                
        for j in range(n_nodes_y):
            
            dT = dT_arr[:,j,:end,i]
            zeroth = spi.trapezoid(dT, x=x_coord[:,:end], axis=1)
            zeroth_x[i,:,j]= zeroth
            dTx = dT * x_coord[:,:end]
            first = spi.trapezoid(dTx, x=x_coord[:,:end], axis=1)
            first_x[i,:,j]= first
            dTxx = dT * x_coord[:,:end] * x_coord[:,:end]
            second = spi.trapezoid(dTxx ,x=x_coord[:,:end], axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                kxx = np.where(zeroth > 0, (second / zeroth) - (first / zeroth)**2, 0)
            k_xx[i,:,j] = kxx
    
    # average vertically and store transverse plume extent for each realization       
    k_yy_n.append(np.average(k_yy, axis=1))
    
    # calculate time derivative of first longitudinal plume extent and average in z- and y-direction 
    dt_first_x = np.gradient(first_x,t, axis=0)
    dt_first_x = np.average(dt_first_x, axis=1)
    dt_first_x = np.average(dt_first_x, axis=1, weights=w_y)
    
    # average zeroth moment vertically and in transversely
    zeroth_x = np.average(zeroth_x, axis=1)
    zeroth_x = np.average(zeroth_x, axis=1, weights=w_y)
    
    # calculate thermal velocity and store for each realization
    v_t = dt_first_x/zeroth_x
    with np.errstate(divide='ignore', invalid='ignore'):
        v_t = np.where(zeroth_x > 0, dt_first_x / zeroth_x, np.nan)   
    v_therm.append(v_t)
    
    # average longitudinal second centered moment vertically and transversely and store for each realization
    k_xx = np.average(k_xx, axis=1)
    k_xx = np.average(k_xx, axis=1, weights=w_y)
    k_xx_n.append(k_xx)

### --- Avergaing moments over all realizations --- ###
#######################################################

kyy_x = np.average(k_yy_n, axis=0)
kxx_t = np.average(k_xx_n, axis=0)
vt = np.average(v_therm, axis=0)

### --- Save results as CSV --- ###
###################################
def save_csv(data, filename, headers):
    df = pd.DataFrame(data, columns=headers)
    df.to_csv(filename, index=False)
    
# save transverse plume as CSV
save_csv(np.column_stack([x_coord[0, :], kyy_x[-1, :]]), '{}kyy_x_sst.csv'.format(path_to_save_results), ['x_coord', 'kyy'])

# calculate transverse dispersion coefficient for steady-state (Hidalgo et al. 2009)
D_T = 0.5 * np.gradient(kyy_x[-1,:],x_coord[0,:]) * q

# save transverse dispersion coefficient as csv
save_csv(np.column_stack([x_coord[0, :], D_T]), '{}D_T.csv'.format(path_to_save_results), ['x_coord', 'D_T'])

# save longitudinal plume extent averaged in transverse direction as CSV
save_csv(np.column_stack([t, kxx_t]), '{}kxx_t.csv'.format(path_to_save_results), ['time', 'kxx'])

### --- Calculate 95% confidence interval of thermal velocity and save results --- ###
######################################################################################
# (only for heterogeneous cases)
"""
vt_std_dev = np.std(v_therm, ddof=1, axis=0)
standard_error = vt_std_dev / np.sqrt(n_real)
confidence_level = 0.95
z_score = norm.ppf(1 - (1 - confidence_level) / 2)
margin_of_error= z_score * standard_error
"""
# save average, standard deviation and margin of error of the thermal velocity as csv

combined_array = np.column_stack([t, vt])
#combined_array = np.column_stack([t, vt_av, vt_std_dev, margin_of_error])

headers = ['time', 'vt']
#headers = ['time', 'vt_average', 'vt_std_dev', 'vt_margin_of_error']

save_csv(combined_array, '{}v_therm.csv'.format(path_to_save_results), headers)

# calculate effective thermal retardation factor and save as CSV
with np.errstate(divide='ignore', invalid='ignore'):
    R_eff = np.where(vt > 0, v_a / vt, np.nan)
save_csv(np.column_stack([t, R_eff]), '{}Reff.csv'.format(path_to_save_results), ['time', 'R_eff'])

### --- Plotting to check results --- ###
#########################################

plt.rcParams.update({'font.size': 16})

# plot steady-state transverse plume extent
y_theoret = (2* D * x_coord[0,:]) / q   # theoretical value of kyy for the homogeneous case
fig = plt.figure(figsize=(8,5.5))
plt.scatter(x_coord[0,11::2],kyy_x[-1,11::2]/1e3)
plt.plot(x_coord[0,11::2],y_theoret[11::2]/1e3, label=r'$2Dx/q_0$', color='black', linewidth=1)
plt.ylabel(r'$\kappa^{\mathrm{eff}}_{yy}$ (10$^3$ m$^2$)')
plt.xlabel('Distance from heat source (m)')
plt.xlim(0,1.5e3)
plt.legend(ncol=1,loc='lower right', frameon=False)
plt.show()

# plot effective thermal retardation factor
fig = plt.figure(figsize=(8,5.5))
plt.scatter(t,R_eff)
plt.hlines(2.67,xmin=0,xmax=1.5e10,linestyles='--',colors='black', label='$R_{app}$=2.67')
plt.ylabel('$R_{eff}$')
plt.xlabel('time (s)')
plt.ylim(2,3.6)
plt.xlim(-10000,1.5e10)
plt.show()
