/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumTimeDerivativeParam.h"

template <>
InputParameters
validParams<INSMomentumTimeDerivativeParam>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addClassDescription("This class computes the time derivative for the incompressible "
                             "Navier-Stokes momentum equation.");
  params.addRequiredParam<Real>("rho", "density");
  return params;
}

INSMomentumTimeDerivativeParam::INSMomentumTimeDerivativeParam(const InputParameters & parameters)
  : TimeDerivative(parameters), _rho(getParam<Real>("rho"))
{
}

Real
INSMomentumTimeDerivativeParam::computeQpResidual()
{
  return _rho * TimeDerivative::computeQpResidual();
}

Real
INSMomentumTimeDerivativeParam::computeQpJacobian()
{
  return _rho * TimeDerivative::computeQpJacobian();
}

Real
INSMomentumTimeDerivativeParam::computeQpOffDiagJacobian(unsigned)
{
  return 0.;
}
