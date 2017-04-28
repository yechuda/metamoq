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

#ifndef FARESSCHRODERSPECIFICTURBULENCEDISSIPATIONRATEBC_H
#define FARESSCHRODERSPECIFICTURBULENCEDISSIPATIONRATEBC_H

#include "IntegratedBC.h"

class FaresSchroderSpecificTurbulenceDissipationRateBC;

template<>
InputParameters validParams<FaresSchroderSpecificTurbulenceDissipationRateBC>();

class FaresSchroderSpecificTurbulenceDissipationRateBC : public IntegratedBC
{
public:

  FaresSchroderSpecificTurbulenceDissipationRateBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Parameters
  Real _p;
  Real _beta_c_star;
};

#endif
