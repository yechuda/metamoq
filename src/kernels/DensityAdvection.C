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

#include "DensityAdvection.h"

template<>
InputParameters validParams<DensityAdvection>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  return params;
}



DensityAdvection::DensityAdvection(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w"))

{
}



Real DensityAdvection::computeQpResidual()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  return U * _grad_u[_qp] * _test[_i][_qp];
}



Real DensityAdvection::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  return U * _grad_phi[_j][_qp] * _test[_i][_qp];
}

Real DensityAdvection::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    return _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];
  }

  else if (jvar == _v_vel_var_number)
  {
    return _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];
  }

  else if (jvar == _w_vel_var_number)
  {
    return _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];
  }

  else
    return 0.0;
}
