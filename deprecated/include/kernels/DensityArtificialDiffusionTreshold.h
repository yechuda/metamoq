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

#ifndef DENSITYARTIFICIALDIFFUSIONTRESHOLD_H
#define DENSITYARTIFICIALDIFFUSIONTRESHOLD_H

#include "Kernel.h"

class DensityArtificialDiffusionTreshold;

template<>
InputParameters validParams<DensityArtificialDiffusionTreshold>();


class DensityArtificialDiffusionTreshold : public Kernel
{
public:
  DensityArtificialDiffusionTreshold(const InputParameters & parameters);
  virtual ~DensityArtificialDiffusionTreshold();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableGradient & _grad_potential;
  unsigned int _potential_var;
  Real _delta;
  Real _mu;
  Real _treshold;
  const VariableValue & _u_old;

};


#endif /* DENSITYARTIFICIALDIFFUSIONTRESHOLD_H */
