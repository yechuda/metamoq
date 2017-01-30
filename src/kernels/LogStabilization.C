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

#include "LogStabilization.h"

template<>
InputParameters validParams<LogStabilization>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("offset","The offset parameter that goes into the exponential function");
  return params;
}


LogStabilization::LogStabilization(const InputParameters & parameters) :
    Kernel(parameters),
    _offset(getParam<Real>("offset"))
{
}

LogStabilization::~LogStabilization()
{
}

Real
LogStabilization::computeQpResidual()
{
  return -_test[_i][_qp] * std::exp(-(_offset + _u[_qp]));
}

Real
LogStabilization::computeQpJacobian()
{
  return -_test[_i][_qp] * std::exp(-(_offset + _u[_qp])) * -_phi[_j][_qp];
}
