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

#include "DensityArtificialDiffusionCoupled.h"

template<>
InputParameters validParams<DensityArtificialDiffusionCoupled>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential","The potential for calculating the advection velocity.");
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addParam<Real>("delta",0.25,"Tuning parameter.");
  params.addRequiredParam<Real>("mobility","Ion mobility coefficient");
  return params;
}

DensityArtificialDiffusionCoupled::DensityArtificialDiffusionCoupled(const InputParameters & parameters) :
    Kernel(parameters),

    _grad_potential(coupledGradient("potential")),
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _potential_var(coupled("potential")),
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(coupled("v")),
    _w_vel_var_number(coupled("w")),
    _delta(getParam<Real>("delta")),
    _mu(getParam<Real>("mobility"))

{
}

DensityArtificialDiffusionCoupled::~DensityArtificialDiffusionCoupled()
{
}

Real DensityArtificialDiffusionCoupled::computeQpResidual()
{
  Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
  Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
  Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
  Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
  Real _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real DensityArtificialDiffusionCoupled::computeQpJacobian()
{
  Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
  Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
  Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
  Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
  Real _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}

Real DensityArtificialDiffusionCoupled::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
  {
    RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
    Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
    Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
    Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
    Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
    Real _d_beta_d_potential = (_mu * _grad_potential[_qp] + U) * _grad_phi[_j][_qp] / (_beta + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_potential = _d_beta_d_potential * _current_elem->hmax() * _delta;

    return _d_D_art_d_potential * _grad_test[_i][_qp] * _grad_u[_qp];
  }

  else if (jvar == _u_vel_var_number)
  {
    Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
    Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
    Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
    Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
    Real _d_beta_d_u_vel = (_mu * _grad_potential[_qp](0) + _u_vel[_qp]) * _phi[_j][_qp] / (_beta + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_u_vel = _d_beta_d_u_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_u_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }

  else if (jvar == _v_vel_var_number)
  {
    Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
    Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
    Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
    Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
    Real _d_beta_d_v_vel = (_mu * _grad_potential[_qp](1) + _v_vel[_qp]) * _phi[_j][_qp] / (_beta + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_v_vel = _d_beta_d_v_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_v_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }

  else if (jvar == _w_vel_var_number)
  {
    Real _beta_x = _mu * _grad_potential[_qp](0) + _u_vel[_qp];
    Real _beta_y = _mu * _grad_potential[_qp](1) + _v_vel[_qp];
    Real _beta_z = _mu * _grad_potential[_qp](2) + _w_vel[_qp];
    Real _beta = std::pow(std::pow(_beta_x, 2.0) + std::pow(_beta_y, 2.0) + std::pow(_beta_z, 2.0), 0.5);
    Real _d_beta_d_w_vel = (_mu * _grad_potential[_qp](2) + _w_vel[_qp]) * _phi[_j][_qp] / (_beta + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_w_vel = _d_beta_d_w_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_w_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }

  else
    return 0.0;
}
