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

// MOOSE includes
#include "BodyForceFilterAux.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<BodyForceFilterAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("treshold", "Treshold value of body force magnitude");
  params.addRequiredCoupledVar("body_force_x", "x-component of the body force");
  params.addCoupledVar("body_force_y", "y-component of the body force");
  params.addCoupledVar("body_force_z", "z-component of the body force");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on the component of body force to be filtered");

  return params;
}

BodyForceFilterAux::BodyForceFilterAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _treshold(getParam<Real>("treshold")),
    _body_force_x(coupledValue("body_force_x")),
    _body_force_y(_mesh.dimension() >= 2 ? coupledValue("body_force_y") : _zero),
    _body_force_z(_mesh.dimension() >= 3 ? coupledValue("body_force_z") : _zero),
    _component(getParam<unsigned>("component"))
{
}

Real
BodyForceFilterAux::computeValue()
{
  Real _magnitude = std::sqrt((_body_force_x[_qp] * _body_force_x[_qp]) + (_body_force_y[_qp] * _body_force_y[_qp]) + (_body_force_z[_qp] * _body_force_z[_qp]));

  if (_magnitude > _treshold)
    if (_component == 0)
      return _body_force_x[_qp];
    else if (_component == 1)
      return _body_force_y[_qp];
    else
      return _body_force_z[_qp];
  else
    return 0;
}
