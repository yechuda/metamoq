/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef FARESSCHRODER88FULLMATERIAL_H
#define FARESSCHRODER88FULLMATERIAL_H

#include "Kernel.h"

// Forward Declarations
class FaresSchroder88FullMaterial;

template <>
InputParameters validParams<FaresSchroder88FullMaterial>();

class FaresSchroder88FullMaterial : public Kernel
{
public:
  FaresSchroder88FullMaterial(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;

  // Parameters
  Real _alpha;
  Real _beta_star;
  Real _beta;
  Real _sigma;

  // Material properties
  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _mu_mol;
  const MaterialProperty<Real> & _omega;
  const MaterialProperty<RealVectorValue> & _gradient_omega;

};

#endif // FARESSCHRODER88FULLMATERIAL_H
