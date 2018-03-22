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

#ifndef OMEGAAUX_H
#define OMEGAAUX_H

#include "AuxKernel.h"

class omegaAux;

template<>
InputParameters validParams<omegaAux>();

class omegaAux : public AuxKernel
{
public:
  omegaAux(const InputParameters & parameters);

  virtual ~omegaAux() {}

protected:

  virtual Real computeValue();

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Parameters
  Real _beta_c_star;
};

#endif //OMEGAAUX_H
