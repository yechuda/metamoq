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
#include "InjectionPeekConstantOffsetBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionPeekConstantOffsetBC>()
{
  MooseEnum geometry("sphere=0 cylinder=1");
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<Real>("r0", "Radius of the corona electrode");
  params.addRequiredParam<Real>("initial_value", "Initial value of space charge density on the corona electrode");
  params.addRequiredParam<MooseEnum>("geometry", geometry, "Geometry type of the corona electrode");
  params.addRequiredCoupledVar("E_x", "Electric field strength x component");
  params.addCoupledVar("E_y", 0, "Electric field strength y component"); // only required in 2D and 3D
  params.addCoupledVar("E_z", 0, "Electric field strength z component"); // only required in 3D
  params.addRequiredParam<Real>("offset", "The offset of the space charge density value");

  return params;
}

InjectionPeekConstantOffsetBC::InjectionPeekConstantOffsetBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _r0(getParam<Real>("r0")),
    _init_val(getParam<Real>("initial_value")),
    _E_x_var(coupled("E_x")),
    _E_x(coupledValue("E_x")),
    _E_y_var(coupled("E_y")),
    _E_y(coupledValue("E_y")),
    _E_z_var(coupled("E_z")),
    _E_z(coupledValue("E_z")),
    _u_old(valueOld()),
    _geometry(getParam<MooseEnum>("geometry")),
    _offset(getParam<Real>("offset"))
{}

Real
InjectionPeekConstantOffsetBC::computeQpResidual()
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];
  Real _E0;

  if (_geometry == 0)
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(0.5 * 100 * _r0)));
  else
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(100 * _r0)));

  if (_t_step == 1)
    return _p * _test[_i][_qp] * (_u[_qp] - _init_val + _offset);
  else
    return _p * _test[_i][_qp] * (_u[_qp] - (_E_n / _E0 * (_u_old[_qp] + _offset)) + _offset);
}

Real
InjectionPeekConstantOffsetBC::computeQpJacobian()
{
  return _p*_phi[_j][_qp]*_test[_i][_qp];
}

Real
InjectionPeekConstantOffsetBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real _E0;

  if (_geometry == 0)
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(0.5 * 100 * _r0)));
  else
    _E0 = 3.1e06 * (1 + (0.308/std::sqrt(100 * _r0)));

  if (jvar == _E_x_var)
    return -_p*_test[_i][_qp]*(_phi[_j][_qp]*_normals[_qp](0) / _E0 * _u_old[_qp]);
  else if (jvar == _E_y_var)
    return -_p*_test[_i][_qp]*(_phi[_j][_qp]*_normals[_qp](1) / _E0 * _u_old[_qp]);
  else if (jvar == _E_z_var)
    return -_p*_test[_i][_qp]*(_phi[_j][_qp]*_normals[_qp](2) / _E0 * _u_old[_qp]);
  else
    return 0.0;
}
