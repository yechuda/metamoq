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

#include "air_turbulent.h"

template<>
InputParameters validParams<air_turbulent>()
{
  InputParameters params = validParams<Material>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Parameters
  params.addRequiredParam<Real>("rho","density");
  params.addRequiredParam<Real>("mu_mol","molecular dynamic viscosity");
  params.addParam<Real>("beta_c_star", 0.09, "beta_c_star parameter");

  return params;
}

air_turbulent::air_turbulent(const InputParameters & parameters) :
    Material(parameters),

    // Material properties declarations
    _rho(declareProperty<Real>("rho")),
    _mu_mol(declareProperty<Real>("mu_mol")),
    _omega(declareProperty<Real>("omega")),
    _gradient_omega(declareProperty<RealVectorValue>("gradient_omega")),

    // Coupled gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),

    // Second derivative tensors
    _second_u_vel(coupledSecond("u")),
    _second_v_vel(coupledSecond("v")),
    _second_w_vel(coupledSecond("w")),

    // Parameters
    _rho_in(getParam<Real>("rho")),
    _mu_mol_in(getParam<Real>("mu_mol")),
    _beta_c_star(getParam<Real>("beta_c_star"))
{}

void
air_turbulent::computeQpProperties()
{
  _mu_mol[_qp] = _mu_mol_in;
  _rho[_qp] = _rho_in;

  RealTensorValue grad_U;
  grad_U(0, 0) = _grad_u_vel[_qp](0);
  grad_U(0, 1) = _grad_u_vel[_qp](1);
  grad_U(0, 2) = _grad_u_vel[_qp](2);

  grad_U(1, 0) = _grad_v_vel[_qp](0);
  grad_U(1, 1) = _grad_v_vel[_qp](1);
  grad_U(1, 2) = _grad_v_vel[_qp](2);

  grad_U(2, 0) = _grad_w_vel[_qp](0);
  grad_U(2, 1) = _grad_w_vel[_qp](1);
  grad_U(2, 2) = _grad_w_vel[_qp](2);

  RealTensorValue Sij = 0.5 * (grad_U + grad_U.transpose());

  Real S = std::pow(2.0 * libMesh::TensorTools::inner_product(Sij, Sij), 0.5);

  _omega[_qp] = S / std::pow(_beta_c_star, 0.5);

  _gradient_omega[_qp](0) = (2.0 / (S * std::pow(_beta_c_star, 0.5))) * (Sij.row(0) * _second_u_vel[_qp].row(0) +
                                                                         Sij.row(1) * _second_v_vel[_qp].row(0) +
                                                                         Sij.row(2) * _second_w_vel[_qp].row(0));

  _gradient_omega[_qp](1) = (2.0 / (S * std::pow(_beta_c_star, 0.5))) * (Sij.row(0) * _second_u_vel[_qp].row(1) +
                                                                         Sij.row(1) * _second_v_vel[_qp].row(1) +
                                                                         Sij.row(2) * _second_w_vel[_qp].row(1));

  _gradient_omega[_qp](2) = (2.0 / (S * std::pow(_beta_c_star, 0.5))) * (Sij.row(0) * _second_u_vel[_qp].row(2) +
                                                                         Sij.row(1) * _second_v_vel[_qp].row(2) +
                                                                         Sij.row(2) * _second_w_vel[_qp].row(2));
}
