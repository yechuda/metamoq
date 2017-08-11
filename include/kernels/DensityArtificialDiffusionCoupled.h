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

#ifndef DENSITYARTIFICIALDIFFUSIONCOUPLED_H
#define DENSITYARTIFICIALDIFFUSIONCOUPLED_H

#include "Kernel.h"

class DensityArtificialDiffusionCoupled;

template<>
InputParameters validParams<DensityArtificialDiffusionCoupled>();


class DensityArtificialDiffusionCoupled : public Kernel
{
public:
  DensityArtificialDiffusionCoupled(const InputParameters & parameters);
  virtual ~DensityArtificialDiffusionCoupled();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableGradient & _grad_potential;
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  unsigned int _potential_var;
  unsigned int _u_vel_var_number;
  unsigned int _v_vel_var_number;
  unsigned int _w_vel_var_number;
  Real _delta;
  Real _mu;
};


#endif /* DENSITYARTIFICIALDIFFUSIONCOUPLED_H */
