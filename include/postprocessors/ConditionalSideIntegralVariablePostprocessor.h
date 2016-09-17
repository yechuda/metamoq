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

#ifndef CONDITIONALSIDEINTEGRALVARIABLEPOSTPROCESSOR_H
#define CONDITIONALSIDEINTEGRALVARIABLEPOSTPROCESSOR_H

#include "SideIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

//Forward Declarations
class ConditionalSideIntegralVariablePostprocessor;

template<>
InputParameters validParams<ConditionalSideIntegralVariablePostprocessor>();

class ConditionalSideIntegralVariablePostprocessor :
  public SideIntegralPostprocessor,
  public MooseVariableInterface
{
public:
  ConditionalSideIntegralVariablePostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const Real & _E0;
  const VariableValue & _E_mag;
  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;
};

#endif
