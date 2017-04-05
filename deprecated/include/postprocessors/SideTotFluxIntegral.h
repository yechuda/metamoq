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

#ifndef SIDETOTFLUXINTEGRAL_H
#define SIDETOTFLUXINTEGRAL_H

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class SideTotFluxIntegral;

template<>
InputParameters validParams<SideTotFluxIntegral>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class SideTotFluxIntegral : public SideIntegralVariablePostprocessor
{
public:
  SideTotFluxIntegral(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const VariableGradient & _grad_potential;
  const Real _mu;
  const Real _diffusivity;
  const Real _d;
};

#endif // SIDETOTFLUXINTEGRAL_H
