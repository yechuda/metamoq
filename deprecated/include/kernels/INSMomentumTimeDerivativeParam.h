/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTIMEDERIVATIVEPARAM_H
#define INSMOMENTUMTIMEDERIVATIVEPARAM_H

#include "TimeDerivative.h"

// Forward Declarations
class INSMomentumTimeDerivativeParam;

template <>
InputParameters validParams<INSMomentumTimeDerivativeParam>();

/**
 * This class computes the time derivative for the incompressible
 * Navier-Stokes momentum equation.  Could instead use CoefTimeDerivative
 * for this.
 */
class INSMomentumTimeDerivativeParam : public TimeDerivative
{
public:
  INSMomentumTimeDerivativeParam(const InputParameters & parameters);

  virtual ~INSMomentumTimeDerivativeParam() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Parameters
  Real _rho;
};

#endif // INSMOMENTUMTIMEDERIVATIVEPARAM_H
