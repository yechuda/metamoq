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

#include "FaresSchroderArtificialDiffusion.h"

template<>
InputParameters validParams<FaresSchroderArtificialDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addParam<Real>("delta",0.25,"Tuning parameter.");
  return params;
}

FaresSchroderArtificialDiffusion::FaresSchroderArtificialDiffusion(const InputParameters & parameters) :
    Kernel(parameters),

    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(coupled("v")),
    _w_vel_var_number(coupled("w")),
    _delta(getParam<Real>("delta"))

{
}

FaresSchroderArtificialDiffusion::~FaresSchroderArtificialDiffusion()
{
}

Real FaresSchroderArtificialDiffusion::computeQpResidual()
{
  Real  _beta = std::pow(std::pow(_u_vel[_qp], 2) + std::pow(_v_vel[_qp], 2) + std::pow(_w_vel[_qp], 2), 0.5);
  Real  _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real FaresSchroderArtificialDiffusion::computeQpJacobian()
{
  Real  _beta = std::pow(std::pow(_u_vel[_qp], 2) + std::pow(_v_vel[_qp], 2) + std::pow(_w_vel[_qp], 2), 0.5);
  Real  _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}

Real FaresSchroderArtificialDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _u_vel_var_number)
  {
    Real _d_beta_d_u_vel = _u_vel[_qp] * _phi[_j][_qp] / (std::pow(std::pow(_u_vel[_qp], 2) + std::pow(_v_vel[_qp], 2) + std::pow(_w_vel[_qp], 2), 0.5) + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_u_vel = _d_beta_d_u_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_u_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  else if (jvar == _v_vel_var_number)
  {
    Real _d_beta_d_v_vel = _v_vel[_qp] * _phi[_j][_qp] / (std::pow(std::pow(_u_vel[_qp], 2) + std::pow(_v_vel[_qp], 2) + std::pow(_w_vel[_qp], 2), 0.5) + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_v_vel = _d_beta_d_v_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_v_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  else if (jvar == _w_vel_var_number)
  {
    Real _d_beta_d_w_vel = _w_vel[_qp] * _phi[_j][_qp] / (std::pow(std::pow(_u_vel[_qp], 2) + std::pow(_v_vel[_qp], 2) + std::pow(_w_vel[_qp], 2), 0.5) + std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_w_vel = _d_beta_d_w_vel * _current_elem->hmax() * _delta;

    return _d_D_art_d_w_vel * _grad_test[_i][_qp] * _grad_u[_qp];
  }
    return 0.0;
}
