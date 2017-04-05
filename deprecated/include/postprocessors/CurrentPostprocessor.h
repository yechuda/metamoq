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

#include "SideIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

//Forward Declarations
class CurrentPostprocessor;

template<>
InputParameters validParams<CurrentPostprocessor>();

class CurrentPostprocessor :
  public SideIntegralPostprocessor,
  public MooseVariableInterface
{
public:
  CurrentPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableGradient & _grad_potential;
  const Real _mu;
  const Real _diffusivity;
  const Real _diam;
};

#endif
