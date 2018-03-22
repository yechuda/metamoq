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

#include "LaplaceCurrentPostprocessor.h"

template<>
InputParameters validParams<LaplaceCurrentPostprocessor>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  params.addRequiredParam<Real>("conductivity", "electric conductivity coefficient");
  return params;
}

LaplaceCurrentPostprocessor::LaplaceCurrentPostprocessor(const InputParameters & parameters) :
    SideIntegralVariablePostprocessor(parameters),
    _sigma(getParam<Real>("conductivity"))
{}

Real
LaplaceCurrentPostprocessor::computeQpIntegral()
{
  return -_sigma * _grad_u[_qp] * _normals[_qp];
}
