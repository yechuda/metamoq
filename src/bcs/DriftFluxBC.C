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

#include "DriftFluxBC.h"

template<>
InputParameters validParams<DriftFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addParam<Real>("mobility", 1.0, "Ion mobility coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  return params;
}

DriftFluxBC::DriftFluxBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _coef(getParam<Real>("mobility")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

Real DriftFluxBC::computeQpResidual()
{
  Real coefficient = _coef;
  return (-coefficient*_u[_qp]*_grad_potential[_qp]*_normals[_qp])*_test[_i][_qp];
}

Real DriftFluxBC::computeQpJacobian()
{
  Real coefficient = _coef;
  return (-coefficient*_phi[_j][_qp]*_grad_potential[_qp]*_normals[_qp])*_test[_i][_qp];
}

Real DriftFluxBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real coefficient = _coef;
  if (jvar == _potential_var)
    return (-coefficient*_u[_qp]*_grad_phi[_j][_qp]*_normals[_qp])*_test[_i][_qp];
  return 0.0;
}
