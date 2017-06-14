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

#include "TurbulenceKineticEnergyTransport.h"

template<>
InputParameters validParams<TurbulenceKineticEnergyTransport>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("omega", "specific turbulence dissipation rate");

  // Parameters
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosiyty");

  params.addParam<Real>("beta_star", 0.09, "beta_star parameter");
  params.addParam<Real>("gamma_star", 1.0, "gamma_star parameter");
  params.addParam<Real>("sigma_star", 0.5, "sigma_star parameter");

  return params;
}



TurbulenceKineticEnergyTransport::TurbulenceKineticEnergyTransport(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _omega(coupledValue("omega")),

  // Coupled gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _omega_var_number(coupled("omega")),

  // Parameters
  _rho(getParam<Real>("rho")),
  _mu_mol(getParam<Real>("mu_mol")),

  _beta_star(getParam<Real>("beta_star")),
  _gamma_star(getParam<Real>("gamma_star")),
  _sigma_star(getParam<Real>("sigma_star"))

{
}



Real TurbulenceKineticEnergyTransport::computeQpResidual()
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

  RealTensorValue KroneckerDelta;
  KroneckerDelta(0,0) = KroneckerDelta(1,1) = KroneckerDelta(2,2) = 1.0;
  KroneckerDelta(0,1) = KroneckerDelta(0,2) = KroneckerDelta(1,0) = KroneckerDelta(1,2) = KroneckerDelta(2,0) = KroneckerDelta(2,1) = 0.0;

  Real production_part = -2.0 * _u[_qp] * libMesh::TensorTools::inner_product(_gamma_star / omega * S - 1.0 / 3.0 * KroneckerDelta, grad_U) * _test[_i][_qp];

  // destruction part
  Real destruction_part = _beta_star * omega * _u[_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol / _rho;
  Real diffusion_part = nu_mol * _grad_u[_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = _sigma_star * _gamma_star / omega * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part;
}



Real TurbulenceKineticEnergyTransport::computeQpJacobian()
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

  RealTensorValue KroneckerDelta;
  KroneckerDelta(0,0) = KroneckerDelta(1,1) = KroneckerDelta(2,2) = 1.0;
  KroneckerDelta(0,1) = KroneckerDelta(0,2) = KroneckerDelta(1,0) = KroneckerDelta(1,2) = KroneckerDelta(2,0) = KroneckerDelta(2,1) = 0.0;

  Real production_part = -2.0 * _phi[_j][_qp] * libMesh::TensorTools::inner_product(_gamma_star / omega * S - 1.0 / 3.0 * KroneckerDelta, grad_U) * _test[_i][_qp];

  // destruction part
  Real destruction_part = _beta_star * omega * _phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real nu_mol = _mu_mol / _rho;
  Real diffusion_part = nu_mol * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  // self diffusion part
  Real self_diffusion_part = _sigma_star * _gamma_star / omega * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
                             _sigma_star * _gamma_star / omega * _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  return convection_part + production_part + destruction_part + diffusion_part + self_diffusion_part;
}

Real TurbulenceKineticEnergyTransport::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    // production part
    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_row_0;
    S_row_0(0) = _grad_u_vel[_qp](0);
    S_row_0(1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    S_row_0(2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));
    Real production_part = -2.0 * _u[_qp] * (2.0 * _gamma_star / omega * S_row_0 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _grad_phi[_j][_qp](0)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    // production part
    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_row_1;
    S_row_1(0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    S_row_1(1) = _grad_v_vel[_qp](1);
    S_row_1(2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));
    Real production_part = -2.0 * _u[_qp] * (2.0 * _gamma_star / omega * S_row_1 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _grad_phi[_j][_qp](1)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    // production part
    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    RealVectorValue S_row_2;
    S_row_2(0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    S_row_2(1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    S_row_2(2) = _grad_w_vel[_qp](2);
    Real production_part = -2.0 * _u[_qp] * (2.0 * _gamma_star / omega * S_row_2 * _grad_phi[_j][_qp] - 1.0 / 3.0 * _grad_phi[_j][_qp](2)) * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _omega_var_number)
  {
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

    Real production_part = 2.0 * _u[_qp] * libMesh::TensorTools::inner_product(_gamma_star / (omega * omega) * _phi[_j][_qp] * S, grad_U) * _test[_i][_qp];

    // destruction part
    Real destruction_part = _beta_star * _phi[_j][_qp] * _u[_qp] * _test[_i][_qp];

    // self diffusion part
    Real self_diffusion_part = -_sigma_star * _gamma_star / (omega * omega) * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

    return production_part + destruction_part + self_diffusion_part;
  }

  else
    return 0.0;
}
