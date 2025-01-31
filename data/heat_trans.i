### --- MOOSE input file for the homogeneous cases --- ####
###########################################################

### Heat transport model with homogeneous permeability values for the MOOSE application 
### using the PorousFlow library. 

### NOTE: To run this input file, a MOOSE application must be created with the
### PorousFlow module enabled in the Makefile.


## Define time and inflow settings for respective scenario

# simulation time (s) 
#end_time = 1.35E11	# for q1
end_time = 4.5E10 	# for q2
#end_time = 2.7E10 	# for q3 

# output time interval (s)
#dt_out = 3E8		# for q1
dt_out = 1E8		# for q2
#dt_out = 6E7		# for q3

# incoming flux at the inflow boundary (kg/m^2s)
#flux = -1.95E-5	# for q1
flux = -5.85E-5		# for q2
#flux = -9.8E-5		# for q3



[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 3
        dx = '140.0 30.0 40.0 30.0 1460.0'   
        dy = '540.0 100.0 20.0 100.0 540.0'  
        dz = '500.0'
        ix = '7 3 8 3 73'     
        iy ='27 10 4 10 27'
        iz = '25' 
    []
    [line_source]
        input = cmg
        type = ExtraNodesetGenerator
        new_boundary = 'line_source'
        coord = '200 650 0; 200 650 20; 200 650 40; 200 650 60; 200 650 80; 
                 200 650 100; 200 650 120; 200 650 140; 200 650 160; 200 650 180; 
                 200 650 200; 200 650 220; 200 650 240; 200 650 260; 200 650 280;
                 200 650 300; 200 650 320; 200 650 340; 200 650 360; 200 650 380;
                 200 650 400; 200 650 420; 200 650 440; 200 650 460; 200 650 480;
                 200 650 500'
    []
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
        flux_function = ${flux}  
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
	#save_in = nodal_outflow 
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
    [thermal_conductivity]
        type = PorousFlowThermalConductivityFromPorosity
        lambda_f = '0.6 0 0  0 0.6 0  0 0 0.6'
        lambda_s = '1.98 0 0  0 1.98 0  0 0 1.98'
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
    end_time = ${end_time}  
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
    time_interval = ${dt_out}
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







