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

#include "DriftFluxOffsetBC.h"

template<>
InputParameters validParams<DriftFluxOffsetBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  params.addRequiredParam<Real>("offset", "The offset of the space charge density value");
  return params;
}

DriftFluxOffsetBC::DriftFluxOffsetBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _coef(getParam<Real>("mobility")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential")),
    _offset(getParam<Real>("offset"))
{}

Real DriftFluxOffsetBC::computeQpResidual()
{
  Real coefficient = _coef;
  return (-coefficient * _u[_qp] * _grad_potential[_qp] * _normals[_qp] - coefficient * _offset * _grad_potential[_qp] * _normals[_qp]) * _test[_i][_qp];
}

Real DriftFluxOffsetBC::computeQpJacobian()
{
  Real coefficient = _coef;
  return (-coefficient*_phi[_j][_qp]*_grad_potential[_qp]*_normals[_qp])*_test[_i][_qp];
}

Real DriftFluxOffsetBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real coefficient = _coef;
  if (jvar == _potential_var)
    return (-coefficient * _u[_qp] * _grad_phi[_j][_qp] * _normals[_qp] - coefficient * _offset * _grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];
  return 0.0;
}
