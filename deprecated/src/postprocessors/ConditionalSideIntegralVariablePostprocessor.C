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

#include "ConditionalSideIntegralVariablePostprocessor.h"

template<>
InputParameters validParams<ConditionalSideIntegralVariablePostprocessor>()
{
  InputParameters params = validParams<SideIntegralPostprocessor>();
  params.addRequiredParam<Real>("E0", "Electric field strength at corona onset");
  params.addRequiredCoupledVar("E_magnitude_Laplace", "The magnitude of local electric field strength from Laplace solution");
  params.addRequiredCoupledVar("variable", "The name of the variable that this boundary condition applies to");
  return params;
}

ConditionalSideIntegralVariablePostprocessor::ConditionalSideIntegralVariablePostprocessor(const InputParameters & parameters) :
    SideIntegralPostprocessor(parameters),
    MooseVariableInterface(this, false),
    _E0(getParam<Real>("E0")),
    _E_mag(coupledValue("E_magnitude_Laplace")),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable"))
{
  addMooseVariableDependency(mooseVariable());
}

Real
ConditionalSideIntegralVariablePostprocessor::computeQpIntegral()
{
  if (_E_mag[_qp] > _E0)
    return _u[_qp];
  else
    return 0;
}
