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
#ifndef INJECTIONTRESHOLDBC_H
#define INJECTIONTRESHOLDBC_H

#include "IntegratedBC.h"

class InjectionTresholdBC;
class Function;

template<>
InputParameters validParams<InjectionTresholdBC>();

class InjectionTresholdBC : public IntegratedBC
{
public:

  InjectionTresholdBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  Real _p;
  Real _init_val;
  unsigned int _E_mag_var;
  const VariableValue & _E_mag;
  unsigned int _E0_mag_var;
  const VariableValue & _E0_mag;
  const VariableValue & _u_old;
  Real _E_0;
};

#endif
