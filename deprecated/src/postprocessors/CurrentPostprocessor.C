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

#include "CurrentPostprocessor.h"

template<>
InputParameters validParams<CurrentPostprocessor>()
{
  InputParameters params = validParams<SideIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this postprocessor applies to");
  params.addRequiredCoupledVar("potential", "The coupled variable of electric potential");
  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredParam<Real>("diameter", "Inner diameter of the tube");

  return params;
}

CurrentPostprocessor::CurrentPostprocessor(const InputParameters & parameters) :
    SideIntegralPostprocessor(parameters),
    MooseVariableInterface(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _grad_potential(coupledGradient("potential")),
    _mu(getParam<Real>("mobility")),
    _diffusivity(getParam<Real>("charge_diffusion_coefficient")),
    _diam(getParam<Real>("diameter"))
{
  addMooseVariableDependency(mooseVariable());
}

Real
CurrentPostprocessor::computeQpIntegral()
{
  return (-_diffusivity * _grad_u[_qp] - _mu * _grad_potential[_qp] * _u[_qp]) * M_PI * _diam;
}
