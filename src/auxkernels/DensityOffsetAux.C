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

#include "DensityOffsetAux.h"

template<>
InputParameters validParams<DensityOffsetAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("density","Space charge density solution.");
  params.addRequiredParam<Real>("offset", "The offset of the space charge density value");
  return params;
}

DensityOffsetAux::DensityOffsetAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _density(coupledValue("density")),
    _offset(getParam<Real>("offset"))
{
}

Real DensityOffsetAux::computeValue()
{
  return _density[_qp] + _offset;
}
