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

#include "AbsDensityAux.h"

template<>
InputParameters validParams<AbsDensityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("density","Space charge density solution.");
  return params;
}

AbsDensityAux::AbsDensityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _density(coupledValue("density"))
{
}

Real AbsDensityAux::computeValue()
{
  if (_density[_qp] < 0.0)
    return -_density[_qp];
  else
    return _density[_qp];
}
