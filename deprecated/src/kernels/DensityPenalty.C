/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DensityPenalty.h"

template<>
InputParameters validParams<DensityPenalty>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("penalty", "Penalty coefficient");
  return params;
}

DensityPenalty::DensityPenalty(const InputParameters & parameters) :
    Kernel(parameters),
    _penalty(getParam<Real>("penalty"))
{
}

Real DensityPenalty::computeQpResidual()
{
  if (_u[_qp] < 0)
  {
    return -_penalty * _u[_qp] * _test[_i][_qp];
  }
  else
    return 0.0;
}

Real DensityPenalty::computeQpJacobian()
{
  if (_u[_qp] < 0)
  {
    return -_penalty * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
    return 0.0;
}
