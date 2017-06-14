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

#ifndef KOMEGADYNAMICVISCOSITYAUX_H
#define KOMEGADYNAMICVISCOSITYAUX_H

#include "AuxKernel.h"

class komegaDynamicViscosityAux;

template<>
InputParameters validParams<komegaDynamicViscosityAux>();

class komegaDynamicViscosityAux : public AuxKernel
{
public:
  komegaDynamicViscosityAux(const InputParameters & parameters);

  virtual ~komegaDynamicViscosityAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _k;
  const VariableValue & _omega;

  // Parameters
  Real _gamma_star;
  Real _rho;
  Real _mu_mol;
};

#endif //KOMEGADYNAMICVISCOSITYAUX_H
