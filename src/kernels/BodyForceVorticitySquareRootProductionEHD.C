/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "BodyForceVorticitySquareRootProductionEHD.h"

template<>
InputParameters validParams<BodyForceVorticitySquareRootProductionEHD>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D

  // Required parameters
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("C", "EHD production coefficient");

  return params;
}



BodyForceVorticitySquareRootProductionEHD::BodyForceVorticitySquareRootProductionEHD(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled gradients
  _grad_body_force_x(coupledGradient("body_force_x")),
  _grad_body_force_y(coupledGradient("body_force_y")),
  _grad_body_force_z(coupledGradient("body_force_z")),

  // Variable numberings
  _body_force_x_var(coupled("body_force_x")),
  _body_force_y_var(coupled("body_force_y")),
  _body_force_z_var(coupled("body_force_z")),

  // Required parameters
  _rho(getParam<Real>("rho")),
  _C(getParam<Real>("C"))

{
}



Real BodyForceVorticitySquareRootProductionEHD::computeQpResidual()
{
  Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
  Real V = std::pow(V_squared, 0.5);
  Real W = std::pow(V, 0.5) * std::pow(1.0 / _rho, 0.5);

  return -_C * _u[_qp] * W * _test[_i][_qp];
}



Real BodyForceVorticitySquareRootProductionEHD::computeQpJacobian()
{
  Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
  Real V = std::pow(V_squared, 0.5);
  Real W = std::pow(V, 0.5) * std::pow(1.0 / _rho, 0.5);

  return -_C * _phi[_j][_qp] * W * _test[_i][_qp];
}

Real BodyForceVorticitySquareRootProductionEHD::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _body_force_x_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](0), _grad_body_force_y[_qp](0), _grad_body_force_z[_qp](0));
    Real dW_dBxj = std::pow(1.0 / _rho, 0.5) / (2.0 * std::pow(V_squared, 0.75)) * (_grad_body_force_x[_qp] - grad_B_column) * _grad_phi[_j][_qp];

    return -_C * _u[_qp] * dW_dBxj * _test[_i][_qp];
  }

  else if (jvar == _body_force_y_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](1), _grad_body_force_y[_qp](1), _grad_body_force_z[_qp](1));
    Real dW_dByj = std::pow(1.0 / _rho, 0.5) / (2.0 * std::pow(V_squared, 0.75)) * (_grad_body_force_y[_qp] - grad_B_column) * _grad_phi[_j][_qp];

    return -_C * _u[_qp] * dW_dByj * _test[_i][_qp];
  }

  else if (jvar == _body_force_z_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](2), _grad_body_force_y[_qp](2), _grad_body_force_z[_qp](2));
    Real dW_dBzj = std::pow(1.0 / _rho, 0.5) / (2.0 * std::pow(V_squared, 0.75)) * (_grad_body_force_z[_qp] - grad_B_column) * _grad_phi[_j][_qp];

    return -_C * _u[_qp] * dW_dBzj * _test[_i][_qp];
  }

  else
    return 0.0;
}
