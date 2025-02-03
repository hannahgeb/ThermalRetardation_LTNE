

# Effective Thermal Retardation and Dispersion in 3D Heterogeneous Aquifers 

This repository accompanies the manuscript:

Gebhardt, H., Zech, A., Rau, G.C., Bayer, P. (202X): Effective thermal retardation in aquifers of heterogeneous hydraulic conductivity

The repository provides scripts to perform transient heat transport simulations in 3D heterogeneous hydraulic conductivity fields using the Multiphysics Object-Oriented Simulation Environment ([MOOSE](https://mooseframework.inl.gov/)), an open-source, parallel finite element framework for numerical modeling (Gaston et al., 2009; Permann et al., 2020).
The heat transport model is created using the MOOSE library [PorousFlow](https://mooseframework.inl.gov/modules/porous_flow/index.html) (Wilkins et al., 2020, 2021), which enables the simulation of transport and flow in porous media.
This repository further provides routines to calculate thermal dispersion and effective thermal retardation based on the temperature fields resulting from heat transport simulations.

To run the heat transport model of this study, [MOOSE](https://mooseframework.inl.gov/) must be [installed](https://mooseframework.inl.gov/getting_started/installation/index.html) and an [application created](https://mooseframework.inl.gov/getting_started/new_users.html), with the [PorousFlow](https://mooseframework.inl.gov/modules/porous_flow/index.html) library enabled in the Makefile. 

## Structure

- `README.md` - description of the repository
- `LICENSE` - the default license is MIT
- `requirements.txt` - requirements for [pip](https://pip.pypa.io/en/stable/user_guide/#requirements-files) to install all needed packages (see below)
- `data/` - contains MOOSE input files and data required to run heat transport simulations for heterogeneous and homogeneous cases
	+ `heat_trans.i` - MOOSE input file for the homogeneous cases (all scenarios)
	+ `heat_trans_add_disp.i` - MOOSE input file for the homogeneous case with added macrodispersion (scenario 'L1q2')
	+ `heat_trans_mgf.i` - MOOSE input file for the heterogeneous cases (all scenarios)
	+ `K0.data` - three-dimensional heterogeneous permeability field generated with the script './src/00_generate_3D_lnK_fields.py'
	+ `lambda_eff_f_var1.txt` - scale-dependent effective longitudinal thermal conductivity of fluid for a log-conductivity variance of 1 
	+ `lambda_eff_f_var2.txt` - scale-dependent effective longitudinal thermal conductivity of fluid for a log-conductivity variance of 2 
	+ `lambda_eff_f_var3.txt` - scale-dependent effective longitudinal thermal conductivity of fluid for a log-conductivity variance of 3
	+ `mesh.e` - mesh for the homogeneous cases with added macrodispersion
- `results/` - contains results of the processed data (dispersion and retardation values) and the figures of the manuscript
	- `Figures_paper/` - contains figures of results as displayed in the main manuscript 
		+ `4_kyy_DT.pdf`
		+ `5_scale_dependent_DL.pdf`
		+ `6_vt_Reff.pdf`
		+ `7_Reff_q.pdf`
		+ `8_Reff_add_disp.pdf`
		+ `9_linear_regression_Reff.pdf`
	- `Figures_SI/` - contains figures of results as displayed in the supporting information
		+ `S1_scale_dependent_lambda_eff_f.pdf`
		+ `S2_kyy.pdf`
		+ `S3_DT.pdf`
		+ `S4_kxx.pdf`
		+ `S5_compare_moments.pdf`
		+ `S6_Reff_L.pdf`
		+ `S7_Reff_normalized_time.pdf`
	- `L1q1/` - contains results for scenario 'L1q1'
	- `L1q2/` - contains results for scenario 'L1q2'
		- `disp_var1/` - contains results for the homogeneous case with macrodispersion for a log-conductivity variance of 1
		- `disp_var2/` - contains results for the homogeneous case with macrodispersion for a log-conductivity variance of 2
		- `disp_var3/` - contains results for the homogeneous case with macrodispersion for a log-conductivity variance of 3
		- `var0/` - contains results for the homogeneous case
		- `var1/` - contains results for the ensemble average of the heterogeneous case for a log-conductivity variance of 1
			+ `D_L.csv` - longitudinal dispersion coefficient and macrodispersivity for each distance from the heat source, average RMSE for fitting procedure
			+ `D_T.csv` - transverse dispersion coefficient for each distance from the heat source 
			+ `kxx_t.csv` - longitudinal second centered moments as a function of time  
			+ `kyy_x_sst.csv` - steady-state transverse second centered moments for each distance from the heat source
			+ `lambda_eff_f.csv` - effective longitudinal thermal conductivity of the fluid for each distance from the heat source
			+ `Reff.csv` - effective thermal retardation factor as a function of time
			+ `v_therm.csv` - average thermal velocity, standard deviation and margin of error as function of time 
		- `var2/` - contains results for the ensemble average of the heterogeneous case for a log-conductivity variance of 2
		- `var3/` - contains results for the ensemble average of the heterogeneous case for a log-conductivity variance of 3
	- `L1q3/` - contains results for scenario 'L1q3'
	- `L2q2/` - contains results for scenario 'L2q2'
	- `L3q2/` - contains results for scenario 'L3q2'
	+ `lin_reg_identify_Rapp.csv` - identification of normalized time step and corresponding effective thermal retardation factors
   	+ `linear_regression.csv` - multi-variable linear regression results (intercept, coefficients, R^2, RMSE) based on effective thermal retardation values for each normalized time step
- `src/` - contains all scripts used for running simulations, data analyses and plotting of results
	+ `00_generate_3D_lnK_fields.py` - script generating random 3D heterogeneous hydraulic conductivity fields
	+ `01_run_sim.py` - script running MOOSE simulations for heterogeneous realizations generated with '00_generate_3D_lnK_fields.py' and storing temperature fields
	+ `02_moments_3D_lnK.py` - script calculating and saving spatial moments (dispersion and retardation) based on stored temperature fields
	+ `03_BTC_analytical_fit_DL.py` - script to determine the scale-dependent longitudinal dispersion coefficient for scenario 'L1q2' based on stored temperature fields
	+ `04_linear_regression_Reff.py` - script performing a multi-variable linear regression based on the effective thermal retardation values 
	+ `F04_kyy_DT.py` - reproducing Figure 4 of the manuscript
	+ `F05_scale_dependent_DL.py` - reproducing Figure 5 of the manuscript
	+ `F06_vtherm_Reff.py` - reproducing Figure 6 of the manuscript
	+ `F07_Reff_q.py` - reproducing Figure 7 of the manuscript
	+ `F08_Reff_add_disp.py`- reproducing Figure 8 of the manuscript
	+ `F09_linear_regression_Reff.py` - reproducing Figure 9 of the manuscript
	+ `SI_F01_scale_dependent_lambda_eff_f.py` - reproducing Figure S1 of the supporting information
	+ `SI_F02_kyy.py` - reproducing Figure S2 of the supporting information
	+ `SI_F03_DT.py` - reproducing Figure S3 of the supporting information
	+ `SI_F04_kxx.py` - reproducing Figure S4 of the supporting information
	+ `SI_F05_compare_kxx_kyy.py` - reproducing Figure S5 of the supporting information
	+ `SI_F06_Reff_L.py` - reproducing Figure S6 of the supporting information
	+ `SI_F07_Reff_normalized_time.py` - reproducing Figure S7 of the supporting information

## References

Gaston, D., Newman, C., Hansen, G., & Lebrun-Grandié, D. (2009). Moose: A parallel computational framework for coupled systems of nonlinear equations. Nuclear Engineering and Design, 239(10), 1768–1778. doi: 10.1016/j.nucengdes.2009.05.021

Permann, C. J., Gaston, D. R., Andrs, D., Carlsen, R. W., Kong, F., Lindsay, A. D., Miller, J. M., Peterson, J. W., Slaughter, A. E., Stogner, R. H., & Martineau, R. C. (2020). Moose: Enabling massively parallel multiphysics simulation. SoftwareX, 11(10), 100430. doi: 10.1016/j.softx.2020.100430

Wilkins, A., Green, C., & Ennis-King, J. (2020). Porousflow: A multiphysics simulation code for coupled problems in porous media. Journal of Open Source Software, 5(55), 2176. doi: 10.21105/joss.02176

Wilkins, A., Green, C. P., & Ennis-King, J. (2021). An open-source multiphysics simulation code for coupled problems in porous media. Computers & Geosciences, 154(2a), 104820. doi: 10.1016/j.cageo.2021.104820

## Contact

You can contact us via <hannah.gebhardt@geo.uni-halle.de>.


## License

MIT © 2020
