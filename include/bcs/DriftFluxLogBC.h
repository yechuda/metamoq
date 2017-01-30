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

#ifndef DRIFTFLUXLOGBC_H
#define DRIFTFLUXLOGBC_H

#include "IntegratedBC.h"

class DriftFluxLogBC;

template<>
InputParameters validParams<DriftFluxLogBC>();

class DriftFluxLogBC : public IntegratedBC
{
public:

  DriftFluxLogBC(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  const Real _mu;
  unsigned int _potential_var;
  const VariableGradient & _grad_potential;
};

#endif //DRIFTFLUXLOGBC_H
