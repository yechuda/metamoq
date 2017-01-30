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

#include "PecletAux.h"

template<>
InputParameters validParams<PecletAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("potential","The potential for calculating the advection velocity.");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredParam<Real>("mobility","Ion mobility coefficient");
  return params;
}

PecletAux::PecletAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _grad_potential(coupledGradient("potential")),
    _D(getParam<Real>("charge_diffusion_coefficient")),
    _mu(getParam<Real>("mobility"))
{
}

Real PecletAux::computeValue()
{
  Real  _beta = _mu * _grad_potential[_qp].norm();
  return (_beta * _current_elem->hmax()) / (2 * _D);
}
