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
#ifndef FARESSCHRODERNOBCBC_H
#define FARESSCHRODERNOBCBC_H

#include "IntegratedBC.h"

class FaresSchroderNoBCBC;

template<>
InputParameters validParams<FaresSchroderNoBCBC>();

class FaresSchroderNoBCBC : public IntegratedBC
{
public:

  FaresSchroderNoBCBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Parameters
  Real _rho;
  Real _mu_mol;
  Real _sigma;

};

#endif
