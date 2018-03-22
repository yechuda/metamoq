/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "dummy.h"

template<>
InputParameters validParams<dummy>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("charge_diffusion_coefficient", 0.0, "Charge diffusion coefficient");
  return params;
}

dummy::dummy(const InputParameters & parameters) :
    Kernel(parameters),
    _coef(getParam<Real>("charge_diffusion_coefficient")),
    _g(getMaterialProperty<RealVectorValue>("g"))
{
}

Real
dummy::computeQpResidual()
{
  Real diffusivity = _coef;

  return diffusivity * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real
dummy::computeQpJacobian()
{
  Real diffusivity = _coef;

  return diffusivity * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
