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

#include "DriftDiffusionRZ.h"

template<>
InputParameters validParams<DriftDiffusionRZ>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  return params;
}

DriftDiffusionRZ::DriftDiffusionRZ(const InputParameters & parameters) :
    Kernel(parameters),
    _mu(getParam<Real>("mobility")),
    _D(getParam<Real>("charge_diffusion_coefficient")),
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

Real DriftDiffusionRZ::computeQpResidual()
{
  Real r = _q_point[_qp](0);
  return -_mu * _grad_potential[_qp](0) / r * _u[_qp] * _test[_i][_qp] - _D * _grad_u[_qp](0) / r * _test[_i][_qp];
}

Real DriftDiffusionRZ::computeQpJacobian()
{
  Real r = _q_point[_qp](0);
  return -_mu * _grad_potential[_qp](0) / r * _phi[_j][_qp] * _test[_i][_qp] - _D * _grad_phi[_j][_qp](0) / r * _test[_i][_qp];
}

Real DriftDiffusionRZ::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
    {
      Real r = _q_point[_qp](0);
      return -_mu * _grad_phi[_j][_qp](0) / r * _u[_qp] * _test[_i][_qp];
    }

  else
    return 0.0;
}
