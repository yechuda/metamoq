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

#include "CurrentDensityMagnitudeAux.h"

template<>
InputParameters validParams<CurrentDensityMagnitudeAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  params.addRequiredCoupledVar("space_charge_density", "The coupled variable of space charge density");
  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  return params;
}

CurrentDensityMagnitudeAux::CurrentDensityMagnitudeAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _grad_potential(coupledGradient("potential")),
    _density(coupledValue("space_charge_density")),
    _grad_density(coupledGradient("space_charge_density")),
    _mu(getParam<Real>("mobility")),
    _D(getParam<Real>("charge_diffusion_coefficient"))
{
}

Real CurrentDensityMagnitudeAux::computeValue()
{
  Real Jx = -_D * _grad_density[_qp](0) - _mu * _grad_potential[_qp](0) * _density[_qp];
  Real Jy = -_D * _grad_density[_qp](1) - _mu * _grad_potential[_qp](1) * _density[_qp];
  Real Jz = -_D * _grad_density[_qp](2) - _mu * _grad_potential[_qp](2) * _density[_qp];

  Real J_mag_squared = std::pow(Jx, 2.0) + std::pow(Jy, 2.0) + std::pow(Jz, 2.0);

  return std::pow(J_mag_squared, 0.5);
}
