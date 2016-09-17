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

#ifndef COUPLEDSPACECHARGEDENSITY_H
#define COUPLEDSPACECHARGEDENSITY_H

#include "Kernel.h"

class CoupledSpaceChargeDensity;

template<>
InputParameters validParams<CoupledSpaceChargeDensity>();

class CoupledSpaceChargeDensity : public Kernel
{
public:
  CoupledSpaceChargeDensity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const Real _coef;
  unsigned int _v_var;
  const VariableValue & _v;
};

#endif //COUPLEDSPACECHARGEDENSITY_H
