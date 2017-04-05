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

#ifndef COUPLEDSPACECHARGEDENSITYLOG_H
#define COUPLEDSPACECHARGEDENSITYLOG_H

#include "Kernel.h"

class CoupledSpaceChargeDensityLog;

template<>
InputParameters validParams<CoupledSpaceChargeDensityLog>();

class CoupledSpaceChargeDensityLog : public Kernel
{
public:
  CoupledSpaceChargeDensityLog(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const Real _rho_rec;
  unsigned int _log_density_var;
  const VariableValue & _log_density;
};

#endif //COUPLEDSPACECHARGEDENSITYLOG_H
