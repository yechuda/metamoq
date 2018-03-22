/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DUMMY_H
#define DUMMY_H

#include "Kernel.h"
#include "Function.h"

class dummy;

template<>
InputParameters validParams<dummy>();

class dummy : public Kernel
{
public:
  dummy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _coef;
  const MaterialProperty<RealVectorValue> & _g;
};

#endif //DUMMY_H
