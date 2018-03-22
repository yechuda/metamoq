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

#ifndef AIR_TURBULENT_H
#define AIR_TURBULENT_H

#include "Material.h"

//Forward Declarations
class air_turbulent;

template<>
InputParameters validParams<air_turbulent>();

class air_turbulent : public Material
{
public:
  air_turbulent(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:
  // Material properties declarations
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _mu_mol;
  MaterialProperty<Real> & _omega;
  MaterialProperty<RealVectorValue> & _gradient_omega;

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Second derivative tensors
  const VariableSecond & _second_u_vel;
  const VariableSecond & _second_v_vel;
  const VariableSecond & _second_w_vel;

  // Parameters
  const Real _rho_in;
  const Real _mu_mol_in;
  const Real _beta_c_star;
};

#endif //AIR_TURBULENT_H
