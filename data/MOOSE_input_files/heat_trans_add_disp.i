### --- MOOSE input file for the homogeneous case with added dispersion --- ####
################################################################################

### Heat transport model with homogeneous permeability values for the MOOSE application 
### using the PorousFlow library. Thermal macrodispersion is added based on the
### effective dispersion coefficients resulting from the corresponding heterogeneous 
### cases. The transverse and longitudinal effective dispersion coefficients are 
### implemented through adjustment of the thermal conductivity of the fluid.

### NOTE: To run this input file, a MOOSE application must be created with the
### PorousFlow module enabled in the Makefile.


## Longitudinal thermal conductivity 'lambda_xx' of the fluid (W/mK) for each 
## distance from the heat source is given by the following values, see files 
## 'lambda_f_var{i}.txt' with i = 1, 2, 3

# for log-conductivity variance of 1:
d_0m = 0.60
d_5m = 0.92
d_10m = 1.23
d_20m = 1.86
d_30m = 2.50
d_40m = 3.13
d_60m = 4.69
d_80m = 5.98
d_100m = 7.22
d_120m = 8.42
d_140m = 9.64
d_160m = 10.85
d_180m = 12.05
d_200m = 13.22
d_220m = 14.37
d_240m = 15.53
d_260m = 16.66
d_280m = 17.77
d_300m = 18.85
d_320m = 19.90
d_340m = 20.85
d_360m = 21.81
d_380m = 22.76
d_400m = 23.69
d_420m = 24.64
d_440m = 25.48
d_460m = 26.32
d_480m = 27.06
d_500m = 27.75
d_520m = 28.42
d_540m = 29.03
d_560m = 29.51
d_580m = 30.00
d_600m = 30.42
d_620m = 30.83
d_640m = 31.26
d_660m = 31.73
d_680m = 32.23
d_700m = 32.77
d_720m = 33.27
d_740m = 33.77
d_760m = 34.24
d_780m = 34.71
d_800m = 35.18
d_820m = 35.64
d_840m = 36.07
d_860m = 36.50
d_880m = 36.88
d_900m = 37.26
d_920m = 37.60
d_940m = 37.89
d_960m = 38.16
d_980m = 38.41
d_1000m = 38.66
d_1020m = 38.92
d_1040m = 39.17
d_1060m = 39.43
d_1080m = 39.71
d_1100m = 40.06
d_1120m = 40.43
d_1140m = 40.82
d_1160m = 41.22
d_1180m = 41.58
d_1200m = 41.90
d_1220m = 42.24
d_1240m = 42.52
d_1260m = 42.72
d_1280m = 42.93
d_1300m = 43.12
d_1320m = 43.31
d_1340m = 43.47
d_1360m = 43.61
d_1380m = 43.74
d_1400m = 43.83
d_1420m = 43.96
d_1440m = 44.05
d_1460m = 44.26
d_1480m = 44.21
d_1500m = 44.93

###############

## transverse thermal conductivity 'lambda_yy' of the fluid (W/mK)

lambda_yy = 1.9 	# for a variance of 1
#lambda_yy = 3.7 	# for a variance of 2
#lambda_yy = 5.2 	# for a variance of 3


## import mesh from file
## the mesh is divided into subdomains for each distance from the heat source
[Mesh]
   type = FileMesh
   file = 'mesh.e'
[]


[GlobalParams]
    PorousFlowDictator = dictator 
[]


[Variables]
    [porepressure]
        initial_condition = 9810
    []
    [temperature]
        initial_condition = 293.0
        scaling = 1E-8
    []
[]



[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure temperature'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[Kernels]
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = porepressure
        gravity = '0 0 0'
    []
    [energy_dot]
        type = PorousFlowEnergyTimeDerivative
        variable = temperature
    []
    [convection]
        type = PorousFlowFullySaturatedHeatAdvection
        variable = temperature
        gravity = '0 0 0'
    []
    [conduction]
        type = PorousFlowHeatConduction
        variable = temperature
    []
[]



[BCs]
    [flow_intlet]
        type = PorousFlowSink 
        variable = porepressure
        flux_function = -5.85e-5  	# for q2
        boundary = 'left'
    []
    [temp_line_source]
        type = DirichletBC
        variable = temperature
        value = 313
	#function = 313
        boundary = line_source
    []
    [const_temperature_inlet]
        type = DirichletBC
        variable = temperature
        value = 293.0
        boundary = 'left'
    []
    [outflow_heat]
        type = PorousFlowOutflowBC
        boundary = 'right'
        flux_type = heat
        variable = temperature
        gravity = '0 0 0'
        include_relperm = false
    []
    [outflow_water]
        type = PorousFlowOutflowBC
        boundary = 'right'
        flux_type = fluid
        variable = porepressure
        gravity = '0 0 0'
        include_relperm = false 
    []
[]




[FluidProperties]
    [simple_fluid]
        type = SimpleFluidProperties
        bulk_modulus = 1E15
        viscosity = 1.0E-3
        density0 = 1000.0
        thermal_expansion = 0 
        cp = 4194
        cv = 4186
        porepressure_coefficient = 0
    []
[]


[Materials]
    [temperature]
        type = PorousFlowTemperature
        temperature = temperature
    []
    [porosity]
        type = PorousFlowPorosityConst
        porosity = 0.25
    []
    [simple_fluid]
        type = PorousFlowSingleComponentFluid
        fp = simple_fluid
        phase = 0
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '4.18e-12 0 0   0 4.18e-12 0   0 0 4.18e-12'
    []
    [rock_internal_energy]
        type = PorousFlowMatrixInternalEnergy
        density = 2650.0
        specific_heat_capacity = 880.0
    []
    [lambda_upstream]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_0m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = '1 2 3'
    []
    [lambda_5m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_5m} 0 0  0 ${lambda_yy} 0  0 0 0.6'

        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 4
    []
    [lambda_10m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_10m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 5
    []
    [lambda_20m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_20m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 6
    []
    [lambda_30m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_30m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 7
    []
    [lambda_40m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_40m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 8
    []
    [lambda_60m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_60m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 9
    []
    [lambda_80m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_80m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 10
    []
    [lambda_100m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_100m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 11
    []
    [lambda_120m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_120m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 12
    []
    [lambda_140m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_140m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 13
    []
    [lambda_160m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_160m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 14
    []
    [lambda_180m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_180m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 15
    []
    [lambda_200m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_200m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 16
    []
    [lambda_220m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_220m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 17
    []
    [lambda_240m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_240m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 18
    []
    [lambda_260m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_260m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 19
    []
    [lambda_280m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_280m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 20
    []
    [lambda_300m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_300m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 21
    []
    [lambda_320m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_320m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 22
    []
    [lambda_340m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_340m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 23
    []
    [lambda_360m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_360m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 24
    []
    [lambda_380m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_380m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 25
    []
    [lambda_400m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_400m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 26
    []
    [lambda_420m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_420m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 27
    []
    [lambda_440m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_440m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 28
    []
    [lambda_460m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_460m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 29
    []
    [lambda_480m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_480m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 30
    []
    [lambda_500m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_500m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 31
    []
    [lambda_520m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_520m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 32
    []
    [lambda_540m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_540m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 33
    []
    [lambda_560m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_560m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 34
    []
    [lambda_580m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_580m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 35
    []
    [lambda_600m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_600m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 36
    []
    [lambda_620m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_620m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 37
    []
    [lambda_640m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_640m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 38
    []
    [lambda_660m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_660m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 39
    []
    [lambda_680m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_680m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 40
    []
    [lambda_700m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_700m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 41
    []
    [lambda_720m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_720m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 42
    []
    [lambda_740m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_740m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 43
    []
    [lambda_760m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_760m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 44
    []
    [lambda_780m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_780m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 45
    []
    [lambda_800m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_800m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 46
    []
    [lambda_820m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_820m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 47
    []
    [lambda_840m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_840m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 48
    []
    [lambda_860m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_860m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 49
    []
    [lambda_880m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_880m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 50
    []
    [lambda_900m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_900m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 51
    []
    [lambda_920m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_920m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 52
    []
    [lambda_940m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_940m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 53
    []
    [lambda_960m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_960m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 54
    []
    [lambda_980m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_980m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 55
    []
    [lambda_1000m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1000m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 56
    []
    [lambda_1020m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1020m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 57
    []
    [lambda_1040m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1040m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 58
    []
    [lambda_1060m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1060m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 59
    []
    [lambda_1080m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1080m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 60
    []
    [lambda_1100m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1100m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 61
    []
    [lambda_1120m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1120m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 62
    []
    [lambda_1140m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1140m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 63
    []
    [lambda_1160m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1160m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 64
    []
    [lambda_1180m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1180m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 65
    []
    [lambda_1200m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1200m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 66
    []
    [lambda_1220m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1220m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 67
    []
    [lambda_1240m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1240m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 68
    []
    [lambda_1260m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1260m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 69
    []
    [lambda_1280m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1280m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 70
    []
    [lambda_1300m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1300m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 71
    []
    [lambda_1320m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1320m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 72
    []
    [lambda_1340m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1340m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 73
    []
    [lambda_1360m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1360m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 74
    []
    [lambda_1380m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1380m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 75
    []
    [lambda_1400m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1400m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 76
    []
    [lambda_1420m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1420m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 77
    []
    [lambda_1440m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1440m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 78
    []
    [lambda_1460m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1460m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 79
    []
    [lambda_1480m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1480m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = 80
    []
    [lambda_1500m]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '${d_1500m} 0 0  0 ${lambda_yy} 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
	block = '81'
    []
    [PS]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
    []
    [massfrac]
        type = PorousFlowMassFraction
    []
    [relperm]
        type = PorousFlowRelativePermeabilityConst
        phase = 0
    []
[]


[Preconditioning]
  active = basic
  [basic]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'bcgs jacobi'  
  []
[]



[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    end_time = 4.5e10  		# for q2
    dt = 5e6   
    nl_abs_tol = 1E-8
    l_abs_tol = 1E-9
[]


[VectorPostprocessors]
    [T]
        type = NodalValueSampler
        sort_by = x
        variable = temperature
        outputs = output
    []
[]



# define output interval
[Times]
  [times]
    type = TimeIntervalTimes
    time_interval = 1E8		# for q2
  []
[]



[Outputs]
    [exodus]
        type = Exodus
	sync_only = true
	sync_times_object = times
    []
    [output] # temperature output
        type = CSV
	sync_only = true
	sync_times_object = times
    []
[]





