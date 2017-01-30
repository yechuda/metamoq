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

#include "DriftDiffusionLog.h"

template<>
InputParameters validParams<DriftDiffusionLog>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential", "The coupled variable of electric potential");
  params.addParam<Real>("mobility", 0.0, "Ion mobility coefficient");
  params.addParam<Real>("charge_diffusion_coefficient", 0.0, "Charge diffusion coefficient");
  return params;
}

DriftDiffusionLog::DriftDiffusionLog(const InputParameters & parameters) :
    Kernel(parameters),

    _mu(getParam<Real>("mobility")),
    _diffusivity(getParam<Real>("charge_diffusion_coefficient")),

    // Coupled variables
    _potential_var(coupled("potential")),
    _grad_potential(coupledGradient("potential"))
{}

DriftDiffusionLog::~DriftDiffusionLog()
{}

Real DriftDiffusionLog::computeQpResidual()
{
  return _mu * std::exp(_u[_qp]) * _grad_potential[_qp] * _grad_test[_i][_qp]
    + _diffusivity * std::exp(_u[_qp]) * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real DriftDiffusionLog::computeQpJacobian()
{
  return _mu * std::exp(_u[_qp]) * _phi[_j][_qp] * _grad_potential[_qp] * _grad_test[_i][_qp]
    + _diffusivity * (std::exp(_u[_qp]) * _grad_phi[_j][_qp] + std::exp(_u[_qp]) * _phi[_j][_qp] * _grad_u[_qp]) * _grad_test[_i][_qp];
}

Real DriftDiffusionLog::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
  {
    return _mu * std::exp(_u[_qp]) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  }
  else
  {
    return 0.;
  }
}
