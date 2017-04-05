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

#ifndef CURRENTDENSITYMAGNITUDEAUX_H
#define CURRENTDENSITYMAGNITUDEAUX_H

#include "AuxKernel.h"

class CurrentDensityMagnitudeAux;

template<>
InputParameters validParams<CurrentDensityMagnitudeAux>();

class CurrentDensityMagnitudeAux : public AuxKernel
{
public:
  CurrentDensityMagnitudeAux(const InputParameters & parameters);

  virtual ~CurrentDensityMagnitudeAux() {}

protected:

  virtual Real computeValue();

  const VariableGradient & _grad_potential;
  const VariableValue & _density;
  const VariableGradient & _grad_density;
  const Real _mu;
  const Real _D;
};

#endif //CURRENTDENSITYMAGNITUDEAUX_H
