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

#include "DriftFluxLogBC.h"

template<>
InputParameters validParams<DriftFluxLogBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addParam<Real>("mobility", 1.0, "Ion mobility coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  return params;
}

DriftFluxLogBC::DriftFluxLogBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _mu(getParam<Real>("mobility")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

Real DriftFluxLogBC::computeQpResidual()
{
  return -_mu * std::exp(_u[_qp]) * _grad_potential[_qp] * _normals[_qp] * _test[_i][_qp];
}

Real DriftFluxLogBC::computeQpJacobian()
{
  return -_mu * std::exp(_u[_qp]) * _phi[_j][_qp] * _grad_potential[_qp] * _normals[_qp] * _test[_i][_qp];
}

Real DriftFluxLogBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
  {
    return -_mu * std::exp(_u[_qp]) * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
