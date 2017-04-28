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

#include "FaresSchroderDynamicViscosityAux.h"

template<>
InputParameters validParams<FaresSchroderDynamicViscosityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("nu_tilde", "Fares-Schroder variable");

  // Parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");

  params.addParam<Real>("c_nu_1", 9.1, "c_nu_1 parameter");

  return params;
}

FaresSchroderDynamicViscosityAux::FaresSchroderDynamicViscosityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _nu_tilde(coupledValue("nu_tilde")),

    // Parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),

    _c_nu_1(getParam<Real>("c_nu_1"))

{
}

Real FaresSchroderDynamicViscosityAux::computeValue()
{
  Real nu_mol = _mu_mol / _rho;
  Real chi = _nu_tilde[_qp] / nu_mol;
  Real f_nu_1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(_c_nu_1, 3.0));
  return _mu_mol + _rho * _nu_tilde[_qp] * f_nu_1;
}
