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

#include "ConditionalSideAverageValue.h"

template<>
InputParameters validParams<ConditionalSideAverageValue>()
{
  InputParameters params = validParams<ConditionalSideIntegralVariablePostprocessor>();
  return params;
}

ConditionalSideAverageValue::ConditionalSideAverageValue(const InputParameters & parameters) :
    ConditionalSideIntegralVariablePostprocessor(parameters),
    _volume(0)
{}

void
ConditionalSideAverageValue::initialize()
{
  ConditionalSideIntegralVariablePostprocessor::initialize();
  _volume = 0;
}

void
ConditionalSideAverageValue::execute()
{
  ConditionalSideIntegralVariablePostprocessor::execute();
  if (_E_mag[_qp] > _E0)
    _volume += volume();
}

Real
ConditionalSideAverageValue::getValue()
{
  Real integral = ConditionalSideIntegralVariablePostprocessor::getValue();
  gatherSum(_volume);
  return integral / _volume;
}

Real
ConditionalSideAverageValue::volume()
{
  return _current_side_volume;
}

void
ConditionalSideAverageValue::threadJoin(const UserObject & y)
{
  ConditionalSideIntegralVariablePostprocessor::threadJoin(y);
  const ConditionalSideAverageValue & pps = static_cast<const ConditionalSideAverageValue &>(y);
  _volume += pps._volume;
}
