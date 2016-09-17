/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DENSITYDIFFUSION_H
#define DENSITYDIFFUSION_H

#include "Kernel.h"
#include "Function.h"

class DensityDiffusion;

template<>
InputParameters validParams<DensityDiffusion>();

class DensityDiffusion : public Kernel
{
public:
  DensityDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _coef;
};

#endif //DENSITYDIFFUSION_H
