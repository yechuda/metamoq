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

#include "SpecificTurbulenceDissipationRateAux.h"

template<>
InputParameters validParams<SpecificTurbulenceDissipationRateAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Parameters
  params.addParam<Real>("beta_c_star", 0.09, "beta_c_star parameter");

  return params;
}

SpecificTurbulenceDissipationRateAux::SpecificTurbulenceDissipationRateAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),

    // Parameters
    _beta_c_star(getParam<Real>("beta_c_star"))

{
}

Real SpecificTurbulenceDissipationRateAux::computeValue()
{
  Real S_mag_squared = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) +
                       2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) +
                       2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2) + (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
  S_mag_squared = std::max(S_mag_squared, std::numeric_limits<Real>::epsilon());

  Real omega = std::pow(S_mag_squared, 0.5) / std::pow(_beta_c_star, 0.5);
  return std::max(omega, std::numeric_limits<Real>::epsilon());
}
