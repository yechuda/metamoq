/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMassRANSRZ.h"

template <>
InputParameters
validParams<INSMassRANSRZ>()
{
  InputParameters params = validParams<INSMassRANS>();
  params.addClassDescription("This class computes the mass equation residual and Jacobian "
                             "contributions for the incompressible Navier-Stokes momentum equation "
                             "in RZ coordinates.");
  return params;
}

INSMassRANSRZ::INSMassRANSRZ(const InputParameters & parameters)
  : INSMassRANS(parameters), _u_vel(coupledValue("u"))
{
}

Real
INSMassRANSRZ::computeQpResidual()
{
  // Base class residual contribution
  Real res_base = INSMassRANS::computeQpResidual();

  // The radial coordinate value.
  const Real r = _q_point[_qp](0);

  // The sign of this term is multiplied by -1 here.
  res_base -= _u_vel[_qp] / r * _test[_i][_qp];

  return res_base;
}

Real
INSMassRANSRZ::computeQpOffDiagJacobian(unsigned jvar)
{
  // Base class jacobian contribution
  Real jac_base = INSMassRANS::computeQpOffDiagJacobian(jvar);

  // The radial coordinate value.
  const Real r = _q_point[_qp](0);

  if (jvar == _u_vel_var_number)
    jac_base -= _phi[_j][_qp] / r * _test[_i][_qp];

  return jac_base;
}
