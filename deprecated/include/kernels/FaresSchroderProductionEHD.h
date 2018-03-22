/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FARESSCHRODERPRODUCTIONEHD_H
#define FARESSCHRODERPRODUCTIONEHD_H

#include "Kernel.h"

// Forward Declarations
class FaresSchroderProductionEHD;

template<>
InputParameters validParams<FaresSchroderProductionEHD>();

class FaresSchroderProductionEHD : public Kernel
{
public:
  FaresSchroderProductionEHD(const InputParameters & parameters);

  virtual ~FaresSchroderProductionEHD(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _space_charge_density;
  const VariableValue & _omega;

  // Coupled gradients
  const VariableGradient & _grad_potential;

  // Variable numberings
  unsigned _potential_var;
  unsigned _space_charge_density_var;
  unsigned _omega_var;

  // Required parameters
  Real _c;
  Real _rho;
  Real _mobility;
};

#endif
