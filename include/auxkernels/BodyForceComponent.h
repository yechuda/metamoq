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

#ifndef BODYFORCECOMPONENT_H
#define BODYFORCECOMPONENT_H

#include "AuxKernel.h"

class BodyForceComponent;

template<>
InputParameters validParams<BodyForceComponent>();

class BodyForceComponent : public AuxKernel
{
public:

  BodyForceComponent(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:

  const VariableGradient & _grad_potential;
  const VariableValue & _space_charge_density;
  int _component;
};

#endif // BODYFORCECOMPONENT_H
