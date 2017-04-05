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

#ifndef ABSDENSITYAUX_H
#define ABSDENSITYAUX_H

#include "AuxKernel.h"

class AbsDensityAux;

template<>
InputParameters validParams<AbsDensityAux>();

class AbsDensityAux : public AuxKernel
{
public:
  AbsDensityAux(const InputParameters & parameters);

  virtual ~AbsDensityAux() {}

protected:

  virtual Real computeValue();

  const VariableValue & _density;
};

#endif //ABSDENSITYAUX_H
