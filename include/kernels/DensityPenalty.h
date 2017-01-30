/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DENSITYPENALTY_H
#define DENSITYPENALTY_H

#include "Kernel.h"
#include "Function.h"

class DensityPenalty;

template<>
InputParameters validParams<DensityPenalty>();

class DensityPenalty : public Kernel
{
public:
  DensityPenalty(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _penalty;
};

#endif //DENSITYPENALTY_H
