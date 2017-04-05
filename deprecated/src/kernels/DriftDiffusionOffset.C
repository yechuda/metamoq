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

#include "DriftDiffusionOffset.h"

template<>
InputParameters validParams<DriftDiffusionOffset>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  params.addRequiredParam<Real>("offset", "The offset of the space charge density value");
  return params;
}

DriftDiffusionOffset::DriftDiffusionOffset(const InputParameters & parameters) :
    Kernel(parameters),
    _coef(getParam<Real>("mobility")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential")),
    _offset(getParam<Real>("offset"))
{}

Real DriftDiffusionOffset::computeQpResidual()
{
  Real coefficient = _coef;
  return coefficient * _offset * _grad_potential[_qp] * _grad_test[_i][_qp];
}

Real DriftDiffusionOffset::computeQpJacobian()
{
  return 0.0;
}

Real DriftDiffusionOffset::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real coefficient = _coef;
  if (jvar == _potential_var)
    return coefficient * _offset * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  return 0.0;
}
