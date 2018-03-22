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
#include "FaresSchroderNoBCBC.h"

template<>
InputParameters validParams<FaresSchroderNoBCBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosiyty");
  params.addParam<Real>("sigma", 1.2, "sigma parameter");

  return params;
}

FaresSchroderNoBCBC::FaresSchroderNoBCBC(const InputParameters & parameters) :
    IntegratedBC(parameters),

    // Parameters
    _rho(getParam<Real>("rho")),
    _mu_mol(getParam<Real>("mu_mol")),
    _sigma(getParam<Real>("sigma"))

{}

Real
FaresSchroderNoBCBC::computeQpResidual()
{
  Real nu_mol = _mu_mol / _rho;
  return -nu_mol * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp] -
          _sigma * _u[_qp] * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
}

Real
FaresSchroderNoBCBC::computeQpJacobian()
{
  Real nu_mol = _mu_mol / _rho;
  return -nu_mol * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp] -
          _sigma * _phi[_j][_qp] * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp] -
          _sigma * _u[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
}
