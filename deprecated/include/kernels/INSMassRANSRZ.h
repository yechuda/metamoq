/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMASSRANSRZ_H
#define INSMASSRANSRZ_H

#include "INSMassRANS.h"

// Forward Declarations
class INSMassRANSRZ;

template <>
InputParameters validParams<INSMassRANSRZ>();

/**
 * This class computes the mass equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation in RZ coordinates.  Inherits most of its functionality
 * from INSMass, and calls its computeQpXYZ() functions when
 * necessary.
 */
class INSMassRANSRZ : public INSMassRANS
{
public:
  INSMassRANSRZ(const InputParameters & parameters);
  virtual ~INSMassRANSRZ() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled values
  const VariableValue & _u_vel;
};

#endif // INSMASSRANSRZ_H
