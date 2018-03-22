/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "InjectionOnsetDampedBC.h"

template<>
InputParameters validParams<InjectionOnsetDampedBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<Real>("initial_value", "Initial value of space charge density on the corona electrode");
  params.addRequiredCoupledVar("E_x", "Electric field strength x component");
  params.addCoupledVar("E_y", 0, "Electric field strength y component"); // only required in 2D and 3D
  params.addCoupledVar("E_z", 0, "Electric field strength z component"); // only required in 3D
  params.addRequiredCoupledVar("E_onset_x", "Corona onset electric field strength x component");
  params.addCoupledVar("E_onset_y", 0, "Corona onset  electric field strength y component"); // only required in 2D and 3D
  params.addCoupledVar("E_onset_z", 0, "Corona onset  electric field strength z component"); // only required in 3D
  params.addRequiredParam<Real>("damping_exp", "Damping exponent");

  return params;
}

InjectionOnsetDampedBC::InjectionOnsetDampedBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _init_val(getParam<Real>("initial_value")),
    _E_x_var(coupled("E_x")),
    _E_x(coupledValue("E_x")),
    _E_y_var(coupled("E_y")),
    _E_y(coupledValue("E_y")),
    _E_z_var(coupled("E_z")),
    _E_z(coupledValue("E_z")),
    _E_onset_x(coupledValue("E_onset_x")),
    _E_onset_y(coupledValue("E_onset_y")),
    _E_onset_z(coupledValue("E_onset_z")),
    _u_old(valueOld()),
    _b(getParam<Real>("damping_exp"))
{}

Real InjectionOnsetDampedBC::computeQpResidual()
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];

  RealVectorValue _E_onset_vector(_E_onset_x[_qp], _E_onset_y[_qp], _E_onset_z[_qp]);
  Real _E0 = _E_onset_vector * _normals[_qp];

  if (_t_step == 1)
    return _p * _test[_i][_qp] * (_u[_qp] - _init_val);
  else
    return _p * _test[_i][_qp] * (_u[_qp] - (std::pow(_E_n / _E0, _b) * _u_old[_qp]));
}

Real InjectionOnsetDampedBC::computeQpJacobian()
{
  return _p*_phi[_j][_qp]*_test[_i][_qp];
}

Real InjectionOnsetDampedBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];

  RealVectorValue _E_onset_vector(_E_onset_x[_qp], _E_onset_y[_qp], _E_onset_z[_qp]);
  Real _E0 = _E_onset_vector * _normals[_qp];

  if (jvar == _E_x_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](0) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else if (jvar == _E_y_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](1) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else if (jvar == _E_z_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](2) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else
    return 0.0;
}
