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
#include "InjectionTresholdBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionTresholdBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<Real>("initial_value", "Initial value of space charge density on the corona electrode");
  params.addRequiredCoupledVar("E_mag", "Electric field strength magnitude");
  params.addRequiredCoupledVar("E0_mag", "Electric field strength magnitude at corona onset");
  params.addRequiredParam<Real>("E_0", "Treshold value of electric field");

  return params;
}

InjectionTresholdBC::InjectionTresholdBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _init_val(getParam<Real>("initial_value")),
    _E_mag_var(coupled("E_mag")),
    _E_mag(coupledValue("E_mag")),
    _E0_mag_var(coupled("E0_mag")),
    _E0_mag(coupledValue("E0_mag")),
    _u_old(valueOld()),
    _E_0(getParam<Real>("E_0"))
{}

Real
InjectionTresholdBC::computeQpResidual()
{
  if (_E_mag[_qp] > _E_0)
    if (_t_step == 1)
      return _p * _test[_i][_qp] * (_u[_qp] - _init_val);
    else
      return _p * _test[_i][_qp] * (_u[_qp] - (_E_mag[_qp] / _E0_mag[_qp]) * _u_old[_qp]);
  else
    return _p * _test[_i][_qp] * _u[_qp];
}

Real
InjectionTresholdBC::computeQpJacobian()
{
  return _p * _phi[_j][_qp] * _test[_i][_qp];
}

Real
InjectionTresholdBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _E_mag_var)
    if (_E_mag[_qp] > _E_0)
      if (_t_step == 1)
        return 0.0;
      else
        return -_p * _test[_i][_qp] * (_phi[_j][_qp] / _E0_mag[_qp]) * _u_old[_qp];
    else
      return 0.0;
  else if (jvar == _E0_mag_var)
    if (_E_mag[_qp] > _E_0)
      if (_t_step == 1)
        return 0.0;
      else
        return _p * _test[_i][_qp] * (_E_mag[_qp] / std::pow(_E0_mag[_qp], 2.0)) * _phi[_j][_qp] * _u_old[_qp];
    else
      return 0.0;
  else
    return 0.0;
}
