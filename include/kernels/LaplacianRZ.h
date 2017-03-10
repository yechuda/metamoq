/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LAPLACIANRZ_H
#define LAPLACIANRZ_H

#include "Kernel.h"
#include "Function.h"

class LaplacianRZ;

template<>
InputParameters validParams<LaplacianRZ>();

class LaplacianRZ : public Kernel
{
public:
  LaplacianRZ(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};

#endif //LAPLACIANRZ_H
