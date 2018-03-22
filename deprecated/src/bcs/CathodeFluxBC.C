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

#include "CathodeFluxBC.h"

template<>
InputParameters validParams<CathodeFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  return params;
}

CathodeFluxBC::CathodeFluxBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _mobility(getParam<Real>("mobility")),
    _diffusivity(getParam<Real>("charge_diffusion_coefficient")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

Real CathodeFluxBC::computeQpResidual()
{
  return -_mobility * _u[_qp] * _grad_potential[_qp] * _normals[_qp] * _test[_i][_qp] - _diffusivity * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
}

Real CathodeFluxBC::computeQpJacobian()
{
  return -_mobility * _phi[_j][_qp] * _grad_potential[_qp] * _normals[_qp] * _test[_i][_qp] - _diffusivity * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
}

Real CathodeFluxBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
    return -_mobility * _u[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
  else
    return 0.0;
}
