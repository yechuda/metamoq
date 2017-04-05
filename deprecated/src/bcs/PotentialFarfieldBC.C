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

#include "PotentialFarfieldBC.h"

template<>
InputParameters validParams<PotentialFarfieldBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  return params;
}

PotentialFarfieldBC::PotentialFarfieldBC(const InputParameters & parameters) :
    IntegratedBC(parameters)

{}

Real PotentialFarfieldBC::computeQpResidual()
{
  return _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
}

Real PotentialFarfieldBC::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
}
