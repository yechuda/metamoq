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

#ifndef DENSITYOFFSETAUX_H
#define DENSITYOFFSETAUX_H

#include "AuxKernel.h"

class DensityOffsetAux;

template<>
InputParameters validParams<DensityOffsetAux>();

class DensityOffsetAux : public AuxKernel
{
public:
  DensityOffsetAux(const InputParameters & parameters);

  virtual ~DensityOffsetAux() {}

protected:

  virtual Real computeValue();

  const VariableValue & _density;
  const Real _offset;
};

#endif //DENSITYOFFSETAUX_H
