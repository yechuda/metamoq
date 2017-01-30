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
#ifndef BODYFORCEFILTERAUX_H
#define BODYFORCEFILTERAUX_H

#include "AuxKernel.h"

class BodyForceFilterAux;

template<>
InputParameters validParams<BodyForceFilterAux>();

class BodyForceFilterAux : public AuxKernel
{
public:
  BodyForceFilterAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const Real _treshold;
  const VariableValue & _body_force_x;
  const VariableValue & _body_force_y;
  const VariableValue & _body_force_z;
  unsigned _component;
};

#endif /* BODYFORCEFILTERAUX_H */
