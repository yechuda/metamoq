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

#ifndef FARESSCHRODERDYNAMICVISCOSITYAUX_H
#define FARESSCHRODERDYNAMICVISCOSITYAUX_H

#include "AuxKernel.h"

class FaresSchroderDynamicViscosityAux;

template<>
InputParameters validParams<FaresSchroderDynamicViscosityAux>();

class FaresSchroderDynamicViscosityAux : public AuxKernel
{
public:
  FaresSchroderDynamicViscosityAux(const InputParameters & parameters);

  virtual ~FaresSchroderDynamicViscosityAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _nu_tilde;

  // Parameters
  Real _mu_mol;
  Real _rho;

  Real _c_nu_1;
};

#endif //FARESSCHRODERDYNAMICVISCOSITYAUX_H
