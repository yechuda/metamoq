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

#ifndef PECLETAUX_H
#define PECLETAUX_H

#include "AuxKernel.h"

class PecletAux;

template<>
InputParameters validParams<PecletAux>();

class PecletAux : public AuxKernel
{
public:
  PecletAux(const InputParameters & parameters);

  virtual ~PecletAux() {}

protected:

  virtual Real computeValue();

  const VariableGradient & _grad_potential;
  Real _D;
  Real _mu;
};

#endif //PECLETAUX_H
