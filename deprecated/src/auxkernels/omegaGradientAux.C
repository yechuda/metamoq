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

#include "omegaGradientAux.h"

template<>
InputParameters validParams<omegaGradientAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Parameters
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on desired component of omega gradient");
  params.addParam<Real>("beta_c_star", 0.09, "beta_c_star parameter");

  return params;
}

omegaGradientAux::omegaGradientAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),

    // Second derivative tensors
    _second_u_vel(coupledSecond("u")),
    _second_v_vel(coupledSecond("v")),
    _second_w_vel(coupledSecond("w")),

    // Parameters
    _component(getParam<unsigned>("component")),
    _beta_c_star(getParam<Real>("beta_c_star"))

{
}

Real omegaGradientAux::computeValue()
{
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

  return (2.0 / (S * std::pow(_beta_c_star, 0.5))) * (Sij.row(0) * _second_u_vel[_qp].row(_component) +
                                                      Sij.row(1) * _second_v_vel[_qp].row(_component) +
                                                      Sij.row(2) * _second_w_vel[_qp].row(_component));
}
