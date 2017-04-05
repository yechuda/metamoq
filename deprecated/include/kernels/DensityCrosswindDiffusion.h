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

#ifndef DENSITYCROSSWINDDIFFUSION_H
#define DENSITYCROSSWINDDIFFUSION_H

#include "Kernel.h"

class DensityCrosswindDiffusion;

template<>
InputParameters validParams<DensityCrosswindDiffusion>();

class DensityCrosswindDiffusion : public Kernel
{
public:

  DensityCrosswindDiffusion(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  // Coupled variables
  const VariableValue & _Ex;
  const VariableValue & _Ey;
  const VariableValue & _Ez;

  // Variable numberings
  unsigned int _Ex_var;
  unsigned int _Ey_var;
  unsigned int _Ez_var;

  // Parameters
  const Real _mu;
  const Real _epsilon;
  const Real _scaling;
};

#endif //DENSITYCROSSWINDDIFFUSION_H
