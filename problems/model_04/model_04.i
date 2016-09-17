[Mesh]
  file = /home/yechuda/projects/metamoq/problems/model_04/model_04.e
[]

[Variables]
  active = 'voltage density'

  [./voltage]
    order = FIRST
    family = LAGRANGE
  [../]

  [./density]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./E_field_x]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./E_field_y]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./E_field_magnitude]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./body_force_x]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./body_force_y]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./E_field_magnitude_Laplace]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff_voltage source_voltage diff_density drift_density time_voltage time_density'

  [./time_voltage]
    type = TimeDerivative
    variable = voltage
  [../]

  [./time_density]
    type = TimeDerivative
    variable = density
  [../]

  [./diff_voltage]
    type = Diffusion
    variable = voltage
  [../]

  [./source_voltage]
    type = CoupledSpaceChargeDensity
    variable = voltage
    space_charge_density = density
    permittivity_reciprocal = -1.12940907e11
  [../]

  [./diff_density]
    type = DensityDiffusion
    variable = density
    charge_diffusion_coefficient = 2.66e-5
  [../]

  [./drift_density]
    type = CoupledPotentialGradient
    variable = density
    potential = voltage
    mobility = 2.1e-4
  [../]
[]

[AuxKernels]
  [./E_field_x_aux]
    type = VariableGradientComponent
    variable = E_field_x
    component = x
    gradient_variable = voltage
  [../]

  [./E_field_y_aux]
    type = VariableGradientComponent
    variable = E_field_y
    component = y
    gradient_variable = voltage
  [../]

  [./E_field_magnitude_aux]
    type = VectorMagnitudeAux
    variable = E_field_magnitude
    x = E_field_x
    y = E_field_y
  [../]

  [./body_force_x_aux]
    type = BodyForceComponent
    variable = body_force_x
    component = x
    potential = voltage
    space_charge_density = density
  [../]

  [./body_force_y_aux]
    type = BodyForceComponent
    variable = body_force_y
    component = y
    potential = voltage
    space_charge_density = density
  [../]

  [./E_field_magnitude_Laplace_aux]
    type = SolutionAux
    solution = Laplace_solution
    variable = E_field_magnitude_Laplace
    from_variable = E_field_magnitude
    execute_on = initial
  [../]
[]

[UserObjects]
  [./Laplace_solution]
    type = SolutionUserObject
    mesh = model_04_Laplace_Exodus.e
    system_variables = 'E_field_magnitude'
    timestep = 'LATEST'
  [../]
[]

[BCs]
  active = 'cathode_voltage anode_voltage anode_tip_density anode_passive_density drift_flux_density'

  [./cathode_voltage]
    type = DirichletBC
    variable = voltage
    boundary = 'cathode'
    value = 0
  [../]

  [./anode_voltage]
    type = DirichletBC
    variable = voltage
    boundary = 'anode_tip anode_passive'
    value = 20000
  [../]

  [./anode_tip_density]
    type = InjectionDirichletBC
    variable = density
    boundary = 'anode_tip'
    coef = 2.24e-9
    E0 = 2.2196e7
    E_magnitude_Laplace = E_field_magnitude_Laplace
  [../]

  [./anode_passive_density]
    type = DirichletBC
    variable = density
    boundary = 'anode_passive'
    value = 0
  [../]

  [./drift_flux_density]
    type = DriftFluxBC
    variable = density
    boundary = 'inlet outlet cathode'
    mobility = 2.1e-4
    potential = voltage
  [../]
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  kernel_coverage_check = false
[]

[Postprocessors]
  active = 'E0_tip E0_corona'

  [./E0_tip]
    type = SideAverageValue
    boundary = 'anode_tip'
    variable = E_field_magnitude
  [../]

  [./E0_corona]
    type = ConditionalSideAverageValue
    boundary = 'anode_tip'
    variable = E_field_magnitude
    E0 = 2.2196e7
    E_magnitude_Laplace = E_field_magnitude_Laplace
  [../]
[]

[Preconditioning]
  active = 'smp'

  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  solve_type = 'PJFNK'
  type = Transient
  end_time = 5e-2
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dtmin = 1e-12
  trans_ss_check = true
  ss_check_tol = 5e-5
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    dt = 1e-7
    growth_factor = 1.2
    optimal_iterations = 15
  [../]
[]

[Outputs]
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
  show_var_residual = 'voltage density'
[]
