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

#ifndef FARESSCHRODERARTIFICIALDIFFUSION_H
#define FARESSCHRODERARTIFICIALDIFFUSION_H

#include "Kernel.h"

class FaresSchroderArtificialDiffusion;

template<>
InputParameters validParams<FaresSchroderArtificialDiffusion>();


class FaresSchroderArtificialDiffusion : public Kernel
{
public:
  FaresSchroderArtificialDiffusion(const InputParameters & parameters);
  virtual ~FaresSchroderArtificialDiffusion();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  Real _delta;
};


#endif /* FARESSCHRODERARTIFICIALDIFFUSION_H */
