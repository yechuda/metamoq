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

#include "air.h"
registerMooseObject("MetamoqApp", air);

template<>
InputParameters validParams<air>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredCoupledVar("mu", "apparent dynamic viscosity");
  params.addRequiredParam<Real>("rho","density");

  return params;
}

air::air(const InputParameters & parameters) :
    Material(parameters),

    _mu(declareProperty<Real>("mu")),
    _rho(declareProperty<Real>("rho")),

    _mu_in(coupledValue("mu")),
    _rho_in(getParam<Real>("rho"))
{}

void
air::computeQpProperties()
{
  _mu[_qp] = _mu_in[_qp];
  _rho[_qp] = _rho_in;
}
