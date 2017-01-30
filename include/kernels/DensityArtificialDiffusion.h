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

#ifndef DENSITYARTIFICIALDIFFUSION_H
#define DENSITYARTIFICIALDIFFUSION_H

#include "Kernel.h"

class DensityArtificialDiffusion;

template<>
InputParameters validParams<DensityArtificialDiffusion>();


class DensityArtificialDiffusion : public Kernel
{
public:
  DensityArtificialDiffusion(const InputParameters & parameters);
  virtual ~DensityArtificialDiffusion();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableGradient & _grad_potential;
  unsigned int _potential_var;
  Real _delta;
  Real _mu;
};


#endif /* DENSITYARTIFICIALDIFFUSION_H */
