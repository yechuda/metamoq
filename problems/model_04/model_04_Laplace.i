[Mesh]
  file = /home/yechuda/projects/metamoq/problems/model_04/model_04.e
[]

[Variables]

  [./voltage]
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
[]

[Kernels]
  active = 'diff_voltage'

  [./diff_voltage]
    type = Diffusion
    variable = voltage
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

  [./E_field_magnitude]
    type = VectorMagnitudeAux
    variable = E_field_magnitude
    x = E_field_x
    y = E_field_y
  [../]
[]

[BCs]
  active = 'cathode_voltage anode_voltage'

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
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  kernel_coverage_check = false
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
  type = Steady
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  [./Exodus]
    type = Exodus
    elemental_as_nodal = true
    execute_elemental_on = none
  [../]

[]

[Debug]
  show_var_residual_norms = true
  show_var_residual = 'voltage'
[]
