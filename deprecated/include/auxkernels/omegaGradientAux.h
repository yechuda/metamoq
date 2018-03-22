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

#ifndef OMEGAGRADIENTAUX_H
#define OMEGAGRADIENTAUX_H

#include "AuxKernel.h"

class omegaGradientAux;

template<>
InputParameters validParams<omegaGradientAux>();

class omegaGradientAux : public AuxKernel
{
public:
  omegaGradientAux(const InputParameters & parameters);

  virtual ~omegaGradientAux() {}

protected:

  virtual Real computeValue();

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Second derivative tensors
  const VariableSecond & _second_u_vel;
  const VariableSecond & _second_v_vel;
  const VariableSecond & _second_w_vel;

  // Parameters
  unsigned _component;
  Real _beta_c_star;
};

#endif //OMEGAGRADIENTAUX_H
