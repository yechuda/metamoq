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

#ifndef DRIFTFLUXBC_H
#define DRIFTFLUXBC_H

#include "IntegratedBC.h"

class DriftFluxBC;

template<>
InputParameters validParams<DriftFluxBC>();

class DriftFluxBC : public IntegratedBC
{
public:

  DriftFluxBC(const InputParameters & parameters);

protected:
  
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  
  const Real _coef;
  unsigned int _potential_var;
  const VariableGradient & _grad_potential;
};

#endif //DRIFTFLUXBC_H
