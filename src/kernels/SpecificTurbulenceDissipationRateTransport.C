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

#include "SpecificTurbulenceDissipationRateTransport.h"

template<>
InputParameters validParams<SpecificTurbulenceDissipationRateTransport>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("k", "turbulence kinetic energy");

  // Parameters
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosiyty");

  params.addParam<Real>("beta", 0.075, "beta parameter");
  params.addParam<Real>("gamma", 0.555556, "gamma parameter");
  params.addParam<Real>("gamma_star", 1.0, "gamma_star parameter");
  params.addParam<Real>("sigma", 0.5, "sigma parameter");

  return params;
}



SpecificTurbulenceDissipationRateTransport::SpecificTurbulenceDissipationRateTransport(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _k(coupledValue("k")),

  // Coupled gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _k_var_number(coupled("k")),

  // Parameters
  _rho(getParam<Real>("rho")),
  _mu_mol(getParam<Real>("mu_mol")),

  _beta(getParam<Real>("beta")),
  _gamma(getParam<Real>("gamma")),
  _gamma_star(getParam<Real>("gamma_star")),
  _sigma(getParam<Real>("sigma"))

{
}



Real SpecificTurbulenceDissipationRateTransport::computeQpResidual()
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

  RealTensorValue KroneckerDelta;
  KroneckerDelta(0,0) = KroneckerDelta(1,1) = KroneckerDelta(2,2) = 1.0;
  KroneckerDelta(0,1) = KroneckerDelta(0,2) = KroneckerDelta(1,0) = KroneckerDelta(1,2) = KroneckerDelta(2,0) = KroneckerDelta(2,1) = 0.0;

  Real production_part = -2.0 * libMesh::TensorTools::inner_product(_gamma * _gamma_star * S - 1.0 / 3.0 * _gamma * _u[_qp] * KroneckerDelta, grad_U) * _test[_i][_qp];

  // destruction part
  Real destruction_part = _beta * _u[_qp] * _u[_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol / _rho;
  Real diffusion_part = nu_mol * _grad_u[_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = _sigma * _gamma_star * _k[_qp] / _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part;
}



Real SpecificTurbulenceDissipationRateTransport::computeQpJacobian()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // production part
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

  RealTensorValue KroneckerDelta;
  KroneckerDelta(0,0) = KroneckerDelta(1,1) = KroneckerDelta(2,2) = 1.0;
  KroneckerDelta(0,1) = KroneckerDelta(0,2) = KroneckerDelta(1,0) = KroneckerDelta(1,2) = KroneckerDelta(2,0) = KroneckerDelta(2,1) = 0.0;

  Real production_part = 2.0 * libMesh::TensorTools::inner_product(1.0 / 3.0 * _gamma * _phi[_j][_qp] * KroneckerDelta, grad_U) * _test[_i][_qp];

  // destruction part
  Real destruction_part = 2.0 * _beta * _u[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol / _rho;
  Real diffusion_part = nu_mol * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = -_sigma * _gamma_star * _k[_qp] / (_u[_qp] * _u[_qp]) * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
                              _sigma * _gamma_star * _k[_qp] / _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part;
}

Real SpecificTurbulenceDissipationRateTransport::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    // production part
    RealVectorValue S_row_0;
    S_row_0(0) = _grad_u_vel[_qp](0);
    S_row_0(1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    S_row_0(2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));
    Real production_part = -2.0 * (2.0 * _gamma * _gamma_star * S_row_0 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _gamma * _u[_qp] * _grad_phi[_j][_qp](0)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    // production part
    RealVectorValue S_row_1;
    S_row_1(0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    S_row_1(1) = _grad_v_vel[_qp](1);
    S_row_1(2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
    Real production_part = -2.0 * (2.0 * _gamma * _gamma_star * S_row_1 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _gamma * _u[_qp] * _grad_phi[_j][_qp](1)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    // production part
    RealVectorValue S_row_2;
    S_row_2(0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    S_row_2(1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    S_row_2(2) = _grad_w_vel[_qp](2);
    Real production_part = -2.0 * (2.0 * _gamma * _gamma_star * S_row_2 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _gamma * _u[_qp] * _grad_phi[_j][_qp](2)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _k_var_number)
  {
    // self diffusion part
    Real self_diffusion_part = _sigma * _gamma_star * _phi[_j][_qp] / _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

    return self_diffusion_part;
  }

  else
    return 0.0;
}
