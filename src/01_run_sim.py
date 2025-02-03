#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script running the MOOSE (Gaston et al. 2009, Permann et al. 2020) simulation 
for 100 heterogeneous realizations, which have been generated previously generated 
using script '00_generate_3D_lnK_fields.py'. This script also processes the temperature 
output (CSV files) resulting from each MOOSE simulation to numpy arrays.

NOTE: 
This script can only be executed if a MOOSE application has been created before, 
with the PorousFlow module library (Wilkins et al. 2020, 2021) enabled in the Makefile. 
The name of the application used in this study is 'hmgf' (Heat transport in 
Multi-Gaussian Fields). This Python script must be executed in the directory of the
created application (in this case in '../hmgf/'). The directory '../hmgf/problems/' 
should contain the MOOSE input file ('heat_trans_mgf.i') and the K{}.data files.

For system requirements and MOOSE installation visit 
https://mooseframework.inl.gov/getting_started/installation/index.html 
For creating an application visit 
https://mooseframework.inl.gov/getting_started/new_users.html

Author: H. Gebhardt
"""

import subprocess
import pandas as pd
import numpy as np

### --- Define functions --- ###
################################

# function to modify the K field in the input file for each realization
# the K{}.data files need to be saved in the same directory as the MOOSE input files
def update_input_file(k):
    input_file = path_input_file + 'heat_trans_mgf.i'
    if k > 0:
        with open(input_file, 'r') as file:
            filedata = file.read()
        filedata = filedata.replace('K{}.data'.format(k-1), 'K{}.data'.format(k))
        with open(input_file, 'w') as file:
            file.write(filedata)

# function to read and process CSV containing temperature values for a time step
def process_csv(file_index, time_step):
    filename = path_input_file + 'heat_trans_mgf_output_T_{:04}.csv'.format(file_index)
    data_frame = pd.read_csv(filename).drop(columns=['id'])
    data_frame = data_frame.sort_values(['z', 'y', 'x'])
    temperature_array = data_frame['temperature'].to_numpy()
    reshaped_array = temperature_array.reshape(n_nodes) - 293
    dT_arr[..., time_step] = reshaped_array
    
### --- Set path and parameters --- ###
#######################################

# application name
app_name = 'hmgf'

# path to folder of the application called 'hmgf' (Heat transport in Multi-Gaussian Fields)
path_to_app = ('../{}/'.format(app_name)) 			
# '../hmgf/problems/' should contain the MOOSE input file and the K{}.data files
path_input_file = ('../{}/problems/'.format(app_name))  	

# number of CPUs used for parallel simulation
n_cpu = 48

# number of realizations
n_real = 100

# number of time steps
time_steps = 450

# number of nodes
n_nodes = (26, 79, 95)  # (z, y, x)

# time multiplier
time_mult = 4

# create array to store temperature field
dT_arr = np.zeros((*n_nodes, time_steps))

### --- Main loop to run MOOSE application 'hmgf' for each realization --- ###
##############################################################################
for k in range(n_real):
    update_input_file(k)

    # run the MOOSE application: activate moose and run application 'hmgf'
    subprocess.run('. $CONDA_PREFIX/etc/profile.d/conda.sh && conda activate moose; mpiexec -n {} {}/{}-opt -i problems/heat_trans_mgf.i'.format(n_cpu, path_to_app, app_name), shell=True)
    
    # process temperature outputs for each time step
    for step in range(time_steps):
        file_index = (step + 1) * time_mult
        process_csv(file_index, step)

    # save the temperature field as a numpy array
    #np.save(path_to_save_dT + 'dT_arr_{}.npy'.format(k), dT_arr)



