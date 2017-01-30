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
#include "InjectionPenaltyNormalBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionPenaltyNormalBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<Real>("coef", "Space charge density scaling coefficient");
  params.addRequiredParam<Real>("E0", "Electric field strength at corona onset");
  params.addRequiredCoupledVar("E_x", "Electric field strength x component");
  params.addCoupledVar("E_y", 0, "Electric field strength y component"); // only required in 2D and 3D
  params.addCoupledVar("E_z", 0, "Electric field strength z component"); // only required in 3D

  return params;
}

InjectionPenaltyNormalBC::InjectionPenaltyNormalBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _coef(getParam<Real>("coef")),
    _E0(getParam<Real>("E0")),
    _E_x_var(coupled("E_x")),
    _E_x(coupledValue("E_x")),
    _E_y_var(coupled("E_y")),
    _E_y(coupledValue("E_y")),
    _E_z_var(coupled("E_z")),
    _E_z(coupledValue("E_z"))
{}

Real
InjectionPenaltyNormalBC::computeQpResidual()
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];
  Real _difference = _E_n - _E0;

  if (_difference > 0)
    return _p*_test[_i][_qp]*(_coef*(_E0-_E_n) + _u[_qp]);
  else
    return _p*_test[_i][_qp]*_u[_qp];
}

Real
InjectionPenaltyNormalBC::computeQpJacobian()
{
  return _p*_phi[_j][_qp]*_test[_i][_qp];
}

Real
InjectionPenaltyNormalBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue _E_vector(_E_x[_qp], _E_y[_qp], _E_z[_qp]);
  Real _E_n = _E_vector * _normals[_qp];
  Real _difference = _E_n - _E0;
  if (jvar == _E_x_var)
    if (_difference > 0)
      return -_p*_test[_i][_qp]*_coef*_phi[_j][_qp]*_normals[_qp](0);
    else
      return 0.0;
  else if (jvar == _E_y_var)
    if (_difference > 0)
      return -_p*_test[_i][_qp]*_coef*_phi[_j][_qp]*_normals[_qp](1);
    else
      return 0.0;
  else if (jvar == _E_z_var)
    if (_difference > 0)
      return -_p*_test[_i][_qp]*_coef*_phi[_j][_qp]*_normals[_qp](3);
    else
      return 0.0;
  else
    return 0.0;
}
