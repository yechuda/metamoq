/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElectricFieldBodyForceExplicitTimeRamp.h"

template<>
InputParameters validParams<ElectricFieldBodyForceExplicitTimeRamp>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D

  // Required parameters
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addParam<Real>("initial_scaling", 1.0e-06, "Initial value of scaling factor for body force");
  params.addParam<Real>("ramp_factor", 0.04, "Multiplication factor for ramp pattern");

  return params;
}



ElectricFieldBodyForceExplicitTimeRamp::ElectricFieldBodyForceExplicitTimeRamp(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _body_force_x(coupledValue("body_force_x")),
  _body_force_y(coupledValue("body_force_y")),
  _body_force_z(coupledValue("body_force_z")),

  // Variable numberings
  _body_force_x_var_number(coupled("body_force_x")),
  _body_force_y_var_number(coupled("body_force_y")),
  _body_force_z_var_number(coupled("body_force_z")),

  // Required parameters
  _component(getParam<unsigned>("component")),
  _initial_scaling(getParam<Real>("initial_scaling")),
  _ramp_factor(getParam<Real>("ramp_factor"))

{
}

Real ElectricFieldBodyForceExplicitTimeRamp::computeQpResidual()
{
  // Real _scaling = _initial_scaling * std::pow(_ramp_factor, _t - 1);
  Real _scaling = _initial_scaling + (_t - 1) * _ramp_factor;

  if (_scaling > 1.0)
    _scaling = 1.0;

  RealVectorValue body_force(_body_force_x[_qp], _body_force_y[_qp], _body_force_z[_qp]);
  return -_scaling * body_force(_component) * _test[_i][_qp];
}

Real ElectricFieldBodyForceExplicitTimeRamp::computeQpJacobian()
{
  return 0;
}

Real ElectricFieldBodyForceExplicitTimeRamp::computeQpOffDiagJacobian(unsigned jvar)
{
  if ((jvar == _body_force_x_var_number) || (jvar == _body_force_y_var_number) || (jvar == _body_force_z_var_number))
  {
    // Real _scaling = _initial_scaling * std::pow(_ramp_factor, _t - 1);
    Real _scaling = _initial_scaling + (_t - 1) * _ramp_factor;

    if (_scaling > 1.0)
      _scaling = 1.0;

    return -_scaling * _phi[_j][_qp] * _test[_i][_qp];
  }

  else
    return 0;
}
