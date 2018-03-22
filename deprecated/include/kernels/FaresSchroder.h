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

#ifndef FARESSCHRODER_H
#define FARESSCHRODER_H

#include "Kernel.h"

// Forward Declarations
class FaresSchroder;

template<>
InputParameters validParams<FaresSchroder>();

class FaresSchroder : public Kernel
{
public:
  FaresSchroder(const InputParameters & parameters);

  virtual ~FaresSchroder(){}

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
  const VariableGradient & _grad_omega;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _omega_var_number;

  // Parameters
  Real _rho;
  Real _mu_mol;

  Real _alpha;
  Real _beta_c_star;
  Real _beta_c;
  Real _sigma;

  // Old value
  // const VariableValue & _nu_tilde;

  // Old gradients
  // const VariableGradient & _grad_nu_tilde;
};

#endif
