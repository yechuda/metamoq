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

#ifndef DENSITYAUX_H
#define DENSITYAUX_H

#include "AuxKernel.h"

class DensityAux;

template<>
InputParameters validParams<DensityAux>();

class DensityAux : public AuxKernel
{
public:
  DensityAux(const InputParameters & parameters);

  virtual ~DensityAux() {}

protected:

  virtual Real computeValue();

  const VariableValue & _log_density;
};

#endif //DENSITYAUX_H
