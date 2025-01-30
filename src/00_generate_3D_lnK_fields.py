#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script generating random 3D heterogeneous hydraulic conductivity (K) fields 
using the Python package GSTools (Müller et al. 2022). The hydraulic conductivity
field follows a log-normal distribution (ln K) and a Gaussian covariance model is used.
One exemplary realization for a log-conductivity variance of 1 is generated.
The generated ln K field is transformed to permeability, which is saved as a 
data input file (K{}.data) for the MOOSE simulation. 

Author: H. Gebhardt
"""

import gstools as gs
import numpy as np
from scipy.stats import gmean
import pandas as pd

### --- Set output path and parameters --- ###
##############################################

# output directory
output_dir = '../data/MOOSE_input_files/'

# number of realizations
n_real = 1

# settings for field generation
L = 100             # correlation length L_x (m)
e = 0.5             # anisotropy ratio e = L_z / L_x (-)
var = 1             # log-conductivity variance
mean = -10.1        # mean

# parameters
nu = 1e-3           # dynamic viscosity of water (Ns/m^2)
g = 9.81            # gravitational acceleration (m/s^2)
rho = 1000          # density of water (kg/m^3)

### --- Discretization and weights for averaging --- ###
########################################################

# define dimensions and discretization of the field
x = np.concatenate([
    np.arange(0, 160, 20, dtype=np.float32),
    np.arange(150, 180, 10, dtype=np.float32),
    np.arange(175, 215, 5, dtype=np.float32),
    np.arange(220, 250, 10, dtype=np.float32),
    np.arange(260, 1720, 20, dtype=np.float32)
])

y = np.concatenate([
    np.arange(0, 560, 20, dtype=np.float32),
    np.arange(550, 650, 10, dtype=np.float32),
    np.arange(645, 665, 5, dtype=np.float32),
    np.arange(670, 770, 10, dtype=np.float32),
    np.arange(780, 1320, 20, dtype=np.float32)
])

z = np.arange(0, 520, 20, dtype=np.float32)

# weights
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
    np.full(1, 10 / width_y),
])

width_x = 1700
w_x = np.concatenate([
    np.full(1, 10 / width_x),
    np.full(6, 20 / width_x),
    np.full(1, 15 / width_x),
    np.full(2, 10 / width_x),
    np.full(1, 7.5 / width_x),
    np.full(7, 5 / width_x),
    np.full(1, 7.5 / width_x),
    np.full(2, 10 / width_x),
    np.full(1, 15 / width_x),
    np.full(72, 20 / width_x),
    np.full(1, 10 / width_x),
])

### --- Generate multi-gaussian hydraulic conductivity fields --- ###
#####################################################################

# storage for geometric means and seeds
seed_arr = np.zeros(n_real, dtype=np.int64)
K_g = np.zeros(n_real)

# realization loop
for i in range(n_real):
    # generate seed and store it
    seed = np.random.randint(1, 100_000_000)
    seed_arr[i] = seed

    # generate random Gaussian field
    model = gs.Gaussian(dim=3, var=var, len_scale=L, anis=e)
    srf = gs.SRF(model, mean=mean, seed=seed)
    K = srf((x, y, z), mesh_type='structured')
    K = gs.transform.array_force_moments(K, mean=mean, var=var)
    K = gs.transform.array_to_lognormal(K)

    # compute geometric mean of hydraulic conductivity 
    Kg = gmean(gmean(gmean(K, axis=2), axis=1, weights=w_y), axis=0, weights=w_x)
    K_g[i] = Kg
    print('Geometric mean for realization {}: {}'.format(i, Kg))

    # transform hydraulic conductivity to permeability
    k = K * nu / (rho * g)

    # write data to file
    data_path = '{}/K{}.data'.format(output_dir, i)
    with open(data_path, 'w') as f:
        f.write('AXIS X\n')
        np.savetxt(f, x, fmt='%.2f', newline=' ')
        f.write('\n\nAXIS Y\n')
        np.savetxt(f, y, fmt='%.2f', newline=' ')
        f.write('\n\nAXIS Z\n')
        np.savetxt(f, z, fmt='%.2f', newline=' ')
        f.write('\n\nDATA\n')
        np.savetxt(f, k.flatten(order='F'))

    # adjust line endings (if Windows system is used for generating K data)
    with open(data_path, 'rb') as f:
        content = f.read().replace(b'\r\n', b'\n')
    with open(data_path, 'wb') as f:
        f.write(content)

# save geometric mean to CSV
# pd.DataFrame(K_g, columns=['K_g']).to_csv('{}/Kg.csv'.format(output_dir), index=False)

# save seeds in case they will be loaded and used again for other variances
# np.save('{}/seed_arr.npy'.format(output_dir), seed_arr)