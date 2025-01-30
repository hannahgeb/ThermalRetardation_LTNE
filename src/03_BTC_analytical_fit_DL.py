#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to determine the scale-dependent longitudinal dispersion coefficient (D_L(x)) 
for scenario L1q2 by fitting a 1D analytical solution (van Genuchten & Alves, 1982) 
to the numerical BTCs.

NOTE:
This script will not run for general users as the temperature fields (numpy arrays) 
are not provided in the repository. A single temperature field of an exemplary 
heterogeneous realization (scenario L1q2, log-conductivity variance of 3) is 
provided in the repository so that the script can be executed.

Author: H. Gebhardt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import scipy.integrate as spi
from scipy.optimize import curve_fit
import pandas as pd

### --- Define functions --- ###
################################

# 1D analytical solution (van Genuchten & Alves, 1982) used to fit the numerical BTCs
def analytical_T(t, v_t, D_L, x):
    term1 = 0.5 * special.erfc((x - v_t * t) / np.sqrt(4 * D_L * t))
    term2 = np.sqrt((t * v_t**2) / (np.pi * D_L)) * np.exp(-(x - v_t * t)**2 / (4 * D_L * t))
    term3 = 0.5 * (1 + (v_t * x / D_L) + (v_t**2 * t / D_L)) * np.exp(v_t * x / D_L) * special.erfc((x + v_t * t) / (2 * np.sqrt(D_L * t)))
    return term1 + term2 - term3

# function to fit analytical solution
def fit_function(t, v_t, D_L):
    return analytical_T(t, v_t, D_L, x)  

# find nearest index
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

### --- Set path and parameters --- ###
#######################################

# define variance for which the effective longitudinal dispersion coefficient D_L is calculated
var = 3

# set path to load temperature fields and to save results from post-processing
#dT_arrays_path = '../L1q2/var{}/dT_arrays/'.format(var)  
dT_arrays_path = ('../data/Example_L1q2_var3/') # path for example realization, scenario L1q2, variance of 3

save_path = ('../results/')
#path_to_save_results = (save_path + 'L1q2/var{}/'.format(var))
path_to_save_results = (save_path + 'Example_L1q2_var3/') # path for example realization, scenario L1q2, variance of 3

# number of realizations
n_real = 1

# number of nodes
n_nodes_x, n_nodes_y, n_nodes_z = 95, 79, 26

# start and end node in x direction considered for BTC fitting
source_node = 16
start_node =  21
end_node = 95

# parameter values
n = 0.25                                    # porosity (-)
lambda_b = 1.6369                           # bulk thermal conductivity (W/m*K)
lambda_f = 0.6                              # thermal conductivity of water (W/m*K)
lambda_s = (lambda_b - lambda_f*n)/(1-n)    # thermal conductivity of solid (W/m*K)
c_f = 4.18e6                                # volumetric heat capacity of water (J/m^3*K)
c_s = 880*2650                              # volumetric heat capacity of the solid (J/m^3*K)
c_b = (c_s*(1-n) + n*c_f)
D_diff = lambda_b / c_b                     # bulk thermal diffusivity (m^2/s)
q2 = 5.85e-8                                # Darcy flux (m/s)

### --- Discretization --- ###
##############################

# time steps for scenario L1q2
dt = 1E8
t_q2 = np.arange(1E8, 4.51E10, dt)

# x-coordinate 
x_coord = np.concatenate([
    np.arange(-200, -40, 20, dtype=np.float32),
    np.arange(-50, -20, 10, dtype=np.float32),
    np.arange(-25, 15, 5, dtype=np.float32),
    np.arange(20, 50, 10, dtype=np.float32),
    np.arange(60, 1520, 20, dtype=np.float32)])

# y-coordinate
y_coord = np.concatenate([
    np.arange(-650,-90,20,dtype=np.float32), 
    np.arange(-100,0,10,dtype=np.float32), 
    np.arange(-5,15,5,dtype=np.float32), 
    np.arange(20,120,10,dtype=np.float32), 
    np.arange(130,670,20,dtype=np.float32)])

### --- Main loop for fitting D_L --- ###
#####################################################

# create array to store values of fitted longitudinal dispersion coefficient, 
# average RMSE and effective thermal conductivty of the fluid
DL_arr = np.zeros((n_real,n_nodes_z,end_node-start_node))
RMSE_arr = np.zeros((n_real,n_nodes_z,end_node-start_node))
lambda_eff_f = np.zeros((end_node-source_node))

for i in range(n_real):  
    
    dT_file = (dT_arrays_path + 'dT_arr_{}.npy'.format(i))
    dT_arr = np.load(dT_file)
    dT_arr = np.where(dT_arr<1e-4,0,dT_arr) # eliminate numerical artifacts

    for j in range(n_nodes_z): 
        for k in range(start_node,end_node,1):
            
            dT = dT_arr[j,:,k,-1]
    
            # find transverse horizontal plume center of mass        
            zeroth_y = spi.trapezoid(dT, x=y_coord)    
            dTy = dT * y_coord
            first_y = spi.trapezoid(dTy, x=y_coord)
            com = first_y/zeroth_y
            com_pos = find_nearest(y_coord,com)
       
            # normalized BTC at transverse position of plume center of mass 
            dT_com = dT_arr[j, com_pos, k, :]  
            dT_max = np.max(dT_com)
            dT_norm = dT_com / dT_max
            
            # calculate thermal velocity as initial guess for curve fitting
            grad_dT = np.gradient(dT_norm, 1e8)
            t_arg = np.argmax(grad_dT)
            t_char = t_q2[t_arg] 
            x = x_coord[k] 
            vt = x / t_char
            
            # initial guesses for thermal velocity and D_L
            initial_guess_vt = vt
            initial_guess_D_L = 1e-6
            
            lower_bound_v_th = vt * 1e-2
            upper_bound_v_th = vt * 1e2 
            lower_bound_D_L = 1e-8
            upper_bound_D_L = 5e-5
            bounds = ([lower_bound_v_th, lower_bound_D_L], [upper_bound_v_th,upper_bound_D_L])
            
            # perform curve fitting
            params, covariance = curve_fit(
                fit_function, t_q2, dT_norm, p0=[initial_guess_vt, initial_guess_D_L], maxfev=500, #bounds=bounds 
            )
            fit_v_t, fit_D_L = params

            # calculate predicted values and residuals
            predicted_values = fit_function(t_q2, fit_v_t, fit_D_L)
            residuals = dT_norm - predicted_values

            # calculate RSS (Residual Sum of Squares)
            RSS = np.sum(residuals**2)

            # calculate RMSE (Root Mean Squared Error)
            rmse = np.sqrt(RSS / len(dT_norm))
            RMSE_arr[i, j, k-start_node] = rmse

            # store fitted dispersion coefficient
            DL_arr[i, j, k-start_node] = fit_D_L
            
            # plot curve fitting (comparison of analytical solution and numerical BTC)
            """
            T_analytical = analytical_T(t=t_q2, v_t=fit_v_t, D_L=fit_D_L, x=x) 
                        
            fig = plt.figure(figsize=(7.5,5.5))
            plt.scatter(t_q2[::3]/1e10,dT_norm[::3],label ='numerical simulation', s=50)
            plt.scatter(t_q2[::3]/1e10,T_analytical[::3], marker='x', label='1D analytical solution')
            plt.ylabel('Normalized temperature (-)', fontsize=16)
            plt.xlabel('Time (10$^{10}$ s)', fontsize=16)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.title('realization {}, layer {}'.format(i,j), fontsize=14)
            plt.legend(fontsize=14,frameon=False)
            plt.show()
            """
### --- Average vertically and across realizations --- ###
##########################################################

DL_av = np.nanmean(DL_arr, axis=0)
DL_av = np.average(DL_av, axis=0)

RMSE_av = np.nanmean(RMSE_arr,axis=0)
RMSE_av = np.average(RMSE_av,axis=0)


### --- Calculate longitudinal macrodispersivity from D_L  --- ###
##################################################################
           
alphaL = (DL_av - lambda_b / c_b) * 1 / (((c_f / c_b) * q2))   


### --- Save scale-dependent dispersion coefficient, macrodispersivity and RMSE  --- ###
########################################################################################

data = {
    "x_coord": x_coord[start_node:end_node],
    "D_L": DL_av,
    "RMSE": RMSE_av,
    "a_L": alphaL,
}
df = pd.DataFrame(data)
df = df.to_csv(path_to_save_results + 'D_L.csv', index=False)

RMSE_av = np.average(RMSE_av, axis=0) 
print(f'average RMSE: {RMSE_av}')

### --- Plotting to check results --- ###
#########################################

# plot scale-dependent dispersion coefficient
fig = plt.figure(figsize=(7.5,5.5))
plt.scatter(x_coord[start_node:end_node],DL_av/D_diff, s=50)
plt.ylabel('$D_L/D_{diff}$', fontsize=16)
plt.xlabel('Distance from heat source', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.show()

### --- Calculation of effective thermal conductivity of the fluid that incorporates macrodispersion --- ###
############################################################################################################

lambda_eff_f[5:] = ((lambda_b + alphaL * c_f * q2)  - (1 - n) * lambda_s) / n 
    
# at the location of the heat source, the effective thermal conductivity of the fluid is 0.6 W/m*K
lambda_eff_f[0] = lambda_f
    
# for the distances < 40 m from the heat source, the effective thermal conductivity is determined by linear interpolation
slope = (lambda_eff_f[5] - lambda_f) / (x_coord[start_node] - x_coord[source_node])
lambda_eff_f[1:5] = slope * x_coord[source_node + 1 : start_node] + lambda_f 

### --- Save effective thermal conductivty of the fluid  --- ###
################################################################

data = {
    "x_coord": x_coord[source_node:end_node],
    "lambda_eff_f": lambda_eff_f
}
df = pd.DataFrame(data)
df = df.to_csv(path_to_save_results + 'lambda_eff_f.csv', index=False)

# save effective thermal conductivity of the fluid as text file for MOOSE input
# only for heterogeneous cases
"""
output_file = '../data/MOOSE_input_files/lambda_eff_f_var{}.txt'.format(var)
with open(output_file, 'w') as f:
    for h, v in zip(np.array(x_coord[source_node:], np.int16), lambda_eff_f):
        formatted_value = f"{v:.2f}"  # Ensure exactly 2 decimal places 
        f.write(f"d_{h}m = {formatted_value}\n")  
"""