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

#ifndef LAPLACECURRENTPOSTPROCESSOR_H
#define LAPLACECURRENTPOSTPROCESSOR_H

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class LaplaceCurrentPostprocessor;

template<>
InputParameters validParams<LaplaceCurrentPostprocessor>();

class LaplaceCurrentPostprocessor : public SideIntegralVariablePostprocessor
{
public:
  LaplaceCurrentPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const Real _sigma;
};

#endif // LAPLACECURRENTPOSTPROCESSOR_H
