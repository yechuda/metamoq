/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LaplacianRZ.h"

template<>
InputParameters validParams<LaplacianRZ>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

LaplacianRZ::LaplacianRZ(const InputParameters & parameters) :
    Kernel(parameters)
{
}

Real
LaplacianRZ::computeQpResidual()
{
  Real r = _q_point[_qp](0);
  return -_grad_u[_qp](0) / r * _test[_i][_qp];
}

Real
LaplacianRZ::computeQpJacobian()
{
  Real r = _q_point[_qp](0);
  return -_grad_phi[_j][_qp](0) / r * _test[_i][_qp];
}
