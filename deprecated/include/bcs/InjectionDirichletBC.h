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

#ifndef INJECTIONDIRICHLETBC_H
#define INJECTIONDIRICHLETBC_H

#include "NodalBC.h"

class InjectionDirichletBC;

template<>
InputParameters validParams<InjectionDirichletBC>();

class InjectionDirichletBC : public NodalBC
{
public:
  InjectionDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const Real & _coef;
  const Real & _E0;
  const VariableValue & _E_mag;
};

#endif /* INJECTIONDIRICHLETBC_H */
