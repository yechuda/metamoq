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

#include "PoissonOffset.h"

template<>
InputParameters validParams<PoissonOffset>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<Real>("permittivity_reciprocal", "The reciprocal of the product of free space permittivity and relative permittivity");
  params.addRequiredParam<Real>("offset", "The offset of the space charge density value");

  return params;
}

PoissonOffset::PoissonOffset(const InputParameters & parameters) :
    Kernel(parameters),
    _coef(getParam<Real>("permittivity_reciprocal")),
    _offset(getParam<Real>("offset"))
{
}

Real
PoissonOffset::computeQpResidual()
{
  Real coefficient = _coef;
  return -coefficient * _offset * _test[_i][_qp];
}

Real
PoissonOffset::computeQpJacobian()
{
  return 0.0;
}

Real
PoissonOffset::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
