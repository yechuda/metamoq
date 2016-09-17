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

#include "CoupledPotentialGradient.h"

template<>
InputParameters validParams<CoupledPotentialGradient>()
{
  InputParameters params = validParams<Kernel>();
  
  params.addParam<Real>("mobility", 0.0, "Ion mobility coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  return params;
}

CoupledPotentialGradient::CoupledPotentialGradient(const InputParameters & parameters) :
    Kernel(parameters),
    _coef(getParam<Real>("mobility")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

Real CoupledPotentialGradient::computeQpResidual()
{
  Real coefficient = _coef;
  return coefficient*_u[_qp]*_grad_potential[_qp]*_grad_test[_i][_qp];
}

Real CoupledPotentialGradient::computeQpJacobian()
{
  Real coefficient = _coef;
  return coefficient*_phi[_j][_qp]*_grad_potential[_qp]*_grad_test[_i][_qp];
}

Real CoupledPotentialGradient::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real coefficient = _coef;
  if (jvar == _potential_var)
    return coefficient*_u[_qp]*_grad_phi[_j][_qp]*_grad_test[_i][_qp];
  return 0.0;
}
