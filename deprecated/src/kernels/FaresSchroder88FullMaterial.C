/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FaresSchroder88FullMaterial.h"

template <>
InputParameters
validParams<FaresSchroder88FullMaterial>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Parameters
  params.addParam<Real>("alpha", 0.29, "alpha parameter");
  params.addParam<Real>("beta_star", 0.09, "beta_star parameter");
  params.addParam<Real>("beta", 0.072, "beta parameter");
  params.addParam<Real>("sigma", 1.2, "sigma parameter");

  return params;
}

FaresSchroder88FullMaterial::FaresSchroder88FullMaterial(const InputParameters & parameters)
  : Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),

  // Coupled gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),

  // Parameters
  _alpha(getParam<Real>("alpha")),
  _beta_star(getParam<Real>("beta_star")),
  _beta(getParam<Real>("beta")),
  _sigma(getParam<Real>("sigma")),

  // Material properties
  _rho(getMaterialProperty<Real>("rho")),
  _mu_mol(getMaterialProperty<Real>("mu_mol")),
  _omega(getMaterialProperty<Real>("omega")),
  _gradient_omega(getMaterialProperty<RealVectorValue>("gradient_omega"))

{
}

Real FaresSchroder88FullMaterial::computeQpResidual()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_u[_qp] * _test[_i][_qp];

  // production part
  RealTensorValue S;
  S(0,0) =                 _grad_u_vel[_qp](0);
  S(0,1) = S(1,0) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
  S(1,1) =                 _grad_v_vel[_qp](1);
  S(0,2) = S(2,0) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));
  S(1,2) = S(2,1) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
  S(2,2) =                 _grad_w_vel[_qp](2);

  Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());

  RealTensorValue grad_U;
  grad_U(0,0) = _grad_u_vel[_qp](0);
  grad_U(0,1) = _grad_u_vel[_qp](1);
  grad_U(0,2) = _grad_u_vel[_qp](2);
  grad_U(1,0) = _grad_v_vel[_qp](0);
  grad_U(1,1) = _grad_v_vel[_qp](1);
  grad_U(1,2) = _grad_v_vel[_qp](2);
  grad_U(2,0) = _grad_w_vel[_qp](0);
  grad_U(2,1) = _grad_w_vel[_qp](1);
  grad_U(2,2) = _grad_w_vel[_qp](2);

  Real production_part = -2.0 * (1.0 - _alpha) / omega * libMesh::TensorTools::inner_product(S, grad_U) * _u[_qp] * _test[_i][_qp];

  // destruction part
  Real destruction_part = (_beta_star - _beta) * omega * _u[_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol[_qp] / _rho[_qp];

  Real diffusion_part = nu_mol * _grad_u[_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = _sigma * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // diffusion destruction part
  Real diffusion_destruction_part = -2.0 * nu_mol / omega * _gradient_omega[_qp] * _grad_u[_qp] * _test[_i][_qp];
  // Real diffusion_destruction_part = 0.0;

  // self diffusion destruction part
  Real self_diffusion_destruction_part = -2.0 * _sigma / omega * _u[_qp] * _gradient_omega[_qp] * _grad_u[_qp] * _test[_i][_qp];
  // Real self_diffusion_destruction_part = 0.0;

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part + diffusion_destruction_part + self_diffusion_destruction_part;
}

Real FaresSchroder88FullMaterial::computeQpJacobian()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // production part
  RealTensorValue S;
  S(0,0) =                 _grad_u_vel[_qp](0);
  S(0,1) = S(1,0) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
  S(1,1) =                 _grad_v_vel[_qp](1);
  S(0,2) = S(2,0) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));
  S(1,2) = S(2,1) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
  S(2,2) =                 _grad_w_vel[_qp](2);

  Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());

  RealTensorValue grad_U;
  grad_U(0,0) = _grad_u_vel[_qp](0);
  grad_U(0,1) = _grad_u_vel[_qp](1);
  grad_U(0,2) = _grad_u_vel[_qp](2);
  grad_U(1,0) = _grad_v_vel[_qp](0);
  grad_U(1,1) = _grad_v_vel[_qp](1);
  grad_U(1,2) = _grad_v_vel[_qp](2);
  grad_U(2,0) = _grad_w_vel[_qp](0);
  grad_U(2,1) = _grad_w_vel[_qp](1);
  grad_U(2,2) = _grad_w_vel[_qp](2);

  Real production_part = -2.0 * (1.0 - _alpha) / omega * libMesh::TensorTools::inner_product(S, grad_U) * _phi[_j][_qp] * _test[_i][_qp];

  // destruction part
  Real destruction_part = (_beta_star - _beta) * omega * _phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol[_qp] / _rho[_qp];

  Real diffusion_part = nu_mol * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = _sigma * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
                             _sigma * _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  // diffusion destruction part
  Real diffusion_destruction_part = -2.0 * nu_mol / omega * _gradient_omega[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp];
  // Real diffusion_destruction_part = 0.0;

  // self diffusion destruction part
  Real self_diffusion_destruction_part = -2.0 * _sigma / omega * _phi[_j][_qp] * _gradient_omega[_qp] * _grad_u[_qp] * _test[_i][_qp] +
                                        (-2.0 * _sigma / omega * _u[_qp] * _gradient_omega[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp]);
  // Real self_diffusion_destruction_part = 0.0;

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part + diffusion_destruction_part + self_diffusion_destruction_part;
}

Real FaresSchroder88FullMaterial::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_u;
    S_u(0) = _grad_u_vel[_qp](0);
    S_u(1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    S_u(2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));
    Real production_part = -2.0 * (1.0 - _alpha) / omega * S_u * _grad_phi[_j][_qp] * _u[_qp] * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_v;
    S_v(0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    S_v(1) = _grad_v_vel[_qp](1);
    S_v(2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
    Real production_part = -2.0 * (1.0 - _alpha) / omega * S_v * _grad_phi[_j][_qp] * _u[_qp] * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_w;
    S_w(0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    S_w(1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    S_w(2) = _grad_w_vel[_qp](2);
    Real production_part = -2.0 * (1.0 - _alpha) / omega * S_w * _grad_phi[_j][_qp] * _u[_qp] * _test[_i][_qp];

    return convection_part + production_part;
  }

  else
    return 0.0;
}
