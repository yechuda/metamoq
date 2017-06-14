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

#ifndef TURBULENCEKINETICENERGYTRANSPORT_H
#define TURBULENCEKINETICENERGYTRANSPORT_H

#include "Kernel.h"

// Forward Declarations
class TurbulenceKineticEnergyTransport;

template<>
InputParameters validParams<TurbulenceKineticEnergyTransport>();

class TurbulenceKineticEnergyTransport : public Kernel
{
public:
  TurbulenceKineticEnergyTransport(const InputParameters & parameters);

  virtual ~TurbulenceKineticEnergyTransport(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _omega;

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _omega_var_number;

  // Parameters
  Real _rho;
  Real _mu_mol;

  Real _beta_star;
  Real _gamma_star;
  Real _sigma_star;
};

#endif
