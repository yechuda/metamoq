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
#include "InjectionPeekVariableDampedBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionPeekVariableDampedBC>()
{
  MooseEnum geometry("sphere=0 cylinder=1");
  MooseEnum r_axis("x=0 y=1 z=2");
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<MooseEnum>("r_axis", r_axis, "The axis along which the corona electrode radius is measured");
  params.addRequiredParam<Real>("initial_value", "Initial value of space charge density on the corona electrode");
  params.addRequiredParam<MooseEnum>("geometry", geometry, "Geometry type of the corona electrode");
  params.addRequiredCoupledVar("E_x", "Electric field strength x component");
  params.addCoupledVar("E_y", 0, "Electric field strength y component"); // only required in 2D and 3D
  params.addCoupledVar("E_z", 0, "Electric field strength z component"); // only required in 3D
  params.addRequiredParam<Real>("damping_exp", "Damping exponent");

  return params;
}

InjectionPeekVariableDampedBC::InjectionPeekVariableDampedBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _r_axis(getParam<MooseEnum>("r_axis")),
    _init_val(getParam<Real>("initial_value")),
    _E_x_var(coupled("E_x")),
    _E_x(coupledValue("E_x")),
    _E_y_var(coupled("E_y")),
    _E_y(coupledValue("E_y")),
    _E_z_var(coupled("E_z")),
    _E_z(coupledValue("E_z")),
    _u_old(valueOld()),
    _geometry(getParam<MooseEnum>("geometry")),
    _b(getParam<Real>("damping_exp"))
{}

Real
InjectionPeekVariableDampedBC::computeQpResidual()
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];
  Real _E0;
  Real _r0 = _q_point[_qp](_r_axis);

  if (_geometry == 0)
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(0.5 * 100 * _r0)));
  else
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(100 * _r0)));

  if (_t_step == 1)
    return _p * _test[_i][_qp] * (_u[_qp] - _init_val);
  else
    return _p * _test[_i][_qp] * (_u[_qp] - (std::pow(_E_n / _E0, _b) * _u_old[_qp]));
}

Real
InjectionPeekVariableDampedBC::computeQpJacobian()
{
  return _p*_phi[_j][_qp]*_test[_i][_qp];
}

Real
InjectionPeekVariableDampedBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];
  Real _E0;
  Real _r0 = _q_point[_qp](_r_axis);

  if (_geometry == 0)
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(0.5 * 100 * _r0)));
  else
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(100 * _r0)));

  if (jvar == _E_x_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](0) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else if (jvar == _E_y_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](1) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else if (jvar == _E_z_var)
    return -_p * _test[_i][_qp] * (_phi[_j][_qp] * _normals[_qp](2) / _E0 * _u_old[_qp]) * _b * std::pow(_E_n / _E0, _b - 1);
  else
    return 0.0;
}
