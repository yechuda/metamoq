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
#include "InjectionPenaltyBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionPenaltyBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredParam<Real>("coef", "Space charge density scaling coefficient");
  params.addRequiredParam<Real>("E0", "Electric field strength at corona onset");
  params.addRequiredCoupledVar("E_magnitude", "The magnitude of local electric field strength");

  return params;
}

InjectionPenaltyBC::InjectionPenaltyBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    _coef(getParam<Real>("coef")),
    _E0(getParam<Real>("E0")),
    _E_mag_var(coupled("E_magnitude")),
    _E_mag(coupledValue("E_magnitude"))
{}

Real
InjectionPenaltyBC::computeQpResidual()
{
  Real _difference = _E_mag[_qp] - _E0;

  if (_difference > 0)
    return _p*_test[_i][_qp]*(_coef*(_E0-_E_mag[_qp]) + _u[_qp]);
  else
    return _p*_test[_i][_qp]*_u[_qp];
}

Real
InjectionPenaltyBC::computeQpJacobian()
{
  return _p*_phi[_j][_qp]*_test[_i][_qp];
}

Real
InjectionPenaltyBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real _difference = _E_mag[_qp] - _E0;
  if (jvar == _E_mag_var)
    if (_difference > 0)
      return -_p*_test[_i][_qp]*_coef*_phi[_j][_qp];
    else
      return 0.0;
  return 0.0;
}
