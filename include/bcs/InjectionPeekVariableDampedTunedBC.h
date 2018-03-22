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
#ifndef INJECTIONPEEKVARIABLEDAMPEDTUNEDBC_H
#define INJECTIONPEEKVARIABLEDAMPEDTUNEDBC_H

#include "IntegratedBC.h"

class InjectionPeekVariableDampedTunedBC;
class Function;

template<>
InputParameters validParams<InjectionPeekVariableDampedTunedBC>();

class InjectionPeekVariableDampedTunedBC : public IntegratedBC
{
public:

  InjectionPeekVariableDampedTunedBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  Real _p;
  int _r_axis;
  Real _init_val;
  unsigned int _E_x_var;
  const VariableValue & _E_x;
  unsigned int _E_y_var;
  const VariableValue & _E_y;
  unsigned int _E_z_var;
  const VariableValue & _E_z;
  const VariableValue & _u_old;
  int _geometry;
  Real _b;
  Real _coef;
};

#endif