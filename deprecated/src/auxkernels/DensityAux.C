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

#include "DensityAux.h"

template<>
InputParameters validParams<DensityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("log_density","The naural logarithm of the space charge density.");
  return params;
}

DensityAux::DensityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _log_density(coupledValue("log_density"))
{
}

Real DensityAux::computeValue()
{
  return std::exp(_log_density[_qp]);
}
