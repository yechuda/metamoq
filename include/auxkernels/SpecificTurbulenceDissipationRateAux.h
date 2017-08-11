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

#ifndef SPECIFICTURBULENCEDISSIPATIONRATEAUX_H
#define SPECIFICTURBULENCEDISSIPATIONRATEAUX_H

#include "AuxKernel.h"

class SpecificTurbulenceDissipationRateAux;

template<>
InputParameters validParams<SpecificTurbulenceDissipationRateAux>();

class SpecificTurbulenceDissipationRateAux : public AuxKernel
{
public:
  SpecificTurbulenceDissipationRateAux(const InputParameters & parameters);

  virtual ~SpecificTurbulenceDissipationRateAux() {}

protected:

  virtual Real computeValue();

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Parameters
  Real _beta_c_star;
};

#endif //SPECIFICTURBULENCEDISSIPATIONRATEAUX_H
