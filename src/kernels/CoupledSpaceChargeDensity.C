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

#include "CoupledSpaceChargeDensity.h"
registerMooseObject("MetamoqApp", CoupledSpaceChargeDensity);

template<>
InputParameters validParams<CoupledSpaceChargeDensity>()
{
  InputParameters params = validParams<Kernel>();

  params.addParam<Real>("permittivity_reciprocal", 0.0, "The reciprocal of the product of free space permittivity and relative permittivity");
  params.addRequiredCoupledVar("space_charge_density", "The coupled variable of space charge density");

  return params;
}

CoupledSpaceChargeDensity::CoupledSpaceChargeDensity(const InputParameters & parameters) :
    Kernel(parameters),
    _coef(getParam<Real>("permittivity_reciprocal")),
    _v_var(coupled("space_charge_density")),
    _v(coupledValue("space_charge_density"))
{
}

Real
CoupledSpaceChargeDensity::computeQpResidual()
{
  Real coefficient = _coef;
  return -coefficient*_v[_qp]*_test[_i][_qp];
}

Real
CoupledSpaceChargeDensity::computeQpJacobian()
{
  return 0;
}

Real
CoupledSpaceChargeDensity::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real coefficient = _coef;
  if (jvar == _v_var)
    return -coefficient*_phi[_j][_qp]*_test[_i][_qp];
  return 0.0;
}
