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

#ifndef CURRENTPOSTPROCESSOR_H
#define CURRENTPOSTPROCESSOR_H

// MOOSE includes
#include "CurrentPostprocessor.h"

// Forward Declarations
class CurrentPostprocessor;

template<>
InputParameters validParams<CurrentPostprocessor>();

class CurrentPostprocessor : public SideIntegralVariablePostprocessor
{
public:
  CurrentPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableGradient & _grad_potential;
  const Real _mu;
  const Real _diffusivity;
  const Real _diam;
};

#endif // CURRENTPOSTPROCESSOR_H
