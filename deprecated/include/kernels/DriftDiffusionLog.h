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

#ifndef DRIFTDIFFUSIONLOG_H
#define DRIFTDIFFUSIONLOG_H

#include "Kernel.h"

class DriftDiffusionLog;

template<>
InputParameters validParams<DriftDiffusionLog>();

class DriftDiffusionLog : public Kernel
{
public:
  DriftDiffusionLog(const InputParameters & parameters);
  virtual ~DriftDiffusionLog();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);


  const Real _mu;
  const Real _diffusivity;
  unsigned int _potential_var;
  const VariableGradient & _grad_potential;
};


#endif /* DRIFTDIFFUSIONLOG_H */
