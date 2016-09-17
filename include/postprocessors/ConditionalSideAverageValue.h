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

#ifndef CONDITIONALSIDEAVERAGEVALUE_H
#define CONDITIONALSIDEAVERAGEVALUE_H

#include "ConditionalSideIntegralVariablePostprocessor.h"

//Forward Declarations
class ConditionalSideAverageValue;

template<>
InputParameters validParams<ConditionalSideAverageValue>();

class ConditionalSideAverageValue : public ConditionalSideIntegralVariablePostprocessor
{
public:
  ConditionalSideAverageValue(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  virtual Real volume();
  Real _volume;
};

#endif
