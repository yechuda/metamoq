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

#ifndef DRIFTDIFFUSIONRZ_H
#define DRIFTDIFFUSIONRZ_H

#include "Kernel.h"

class DriftDiffusionRZ;

template<>
InputParameters validParams<DriftDiffusionRZ>();

class DriftDiffusionRZ : public Kernel
{
public:

  DriftDiffusionRZ(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  const Real _mu;
  const Real _D;
  unsigned int _potential_var;
  const VariableGradient & _grad_potential;
};

#endif //DRIFTDIFFUSIONRZ_H
