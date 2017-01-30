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

#include "CoupledSpaceChargeDensityLog.h"

template<>
InputParameters validParams<CoupledSpaceChargeDensityLog>()
{
  InputParameters params = validParams<Kernel>();

  params.addParam<Real>("permittivity_reciprocal", 0.0, "The reciprocal of the product of free space permittivity and relative permittivity");
  params.addRequiredCoupledVar("log_density", "The natural logarithm of space charge density");

  return params;
}

CoupledSpaceChargeDensityLog::CoupledSpaceChargeDensityLog(const InputParameters & parameters) :
    Kernel(parameters),
    _rho_rec(getParam<Real>("permittivity_reciprocal")),
    _log_density_var(coupled("log_density")),
    _log_density(coupledValue("log_density"))
{
}

Real CoupledSpaceChargeDensityLog::computeQpResidual()
{
  return -_rho_rec * std::exp(_log_density[_qp]) * _test[_i][_qp];
}

Real CoupledSpaceChargeDensityLog::computeQpJacobian()
{
  return 0;
}

Real CoupledSpaceChargeDensityLog::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _log_density_var)
  {
    return -_rho_rec * std::exp(_log_density[_qp]) * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
