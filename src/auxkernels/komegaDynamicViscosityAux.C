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

#include "komegaDynamicViscosityAux.h"

template<>
InputParameters validParams<komegaDynamicViscosityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("k", "turbulence kinetic energy");
  params.addRequiredCoupledVar("omega", "specific turbulence dissipation rate");

  // Parameters
  params.addParam<Real>("gamma_star", 1.0, "gamma_star parameter");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");

  return params;
}

komegaDynamicViscosityAux::komegaDynamicViscosityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _k(coupledValue("k")),
    _omega(coupledValue("omega")),

    // Prameters
    _gamma_star(getParam<Real>("gamma_star")),
    _rho(getParam<Real>("rho")),
    _mu_mol(getParam<Real>("mu_mol"))

{
}

Real komegaDynamicViscosityAux::computeValue()
{
  return _mu_mol + _gamma_star * _rho * _k[_qp] / _omega[_qp];
}
