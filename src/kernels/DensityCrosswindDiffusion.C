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

#include "DensityCrosswindDiffusion.h"

template<>
InputParameters validParams<DensityCrosswindDiffusion>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("Ex", "The component of electric field in x direction");
  params.addCoupledVar("Ey", 0, "The component of electric field in y direction"); // only required in 2D and 3D
  params.addCoupledVar("Ez", 0, "The component of electric field in z direction"); // only required in 3D

  // Required parameters
  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredParam<Real>("scaling", "Scaling coefficient");

  return params;
}

DensityCrosswindDiffusion::DensityCrosswindDiffusion(const InputParameters & parameters) :
    Kernel(parameters),

    // Coupled variables
    _Ex(coupledValue("Ex")),
    _Ey(coupledValue("Ey")),
    _Ez(coupledValue("Ez")),

    // Variable numberings
    _Ex_var(coupled("Ex")),
    _Ey_var(coupled("Ey")),
    _Ez_var(coupled("Ez")),

    // Required parameters
    _mu(getParam<Real>("mobility")),
    _epsilon(getParam<Real>("charge_diffusion_coefficient")),
    _scaling(getParam<Real>("scaling"))

{
}

Real DensityCrosswindDiffusion::computeQpResidual()
{
  Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());

  RealVectorValue b;
  b(0) = -_mu * _Ex[_qp];
  b(1) = -_mu * _Ey[_qp];
  b(2) = -_mu * _Ez[_qp];

  Real epsilon_tilde = b.norm() * std::pow(h, 1.5) - _epsilon;

  if (epsilon_tilde > 0.0)
    {
      Real b_norm_2 = std::pow(b.norm(), 2.0);
      Real mu_2 = std::pow(_mu, 2.0);

      RealTensorValue D;
      D(0,0) =          1 - (mu_2 * _Ex[_qp] * _Ex[_qp] / b_norm_2);
      D(0,1) = D(1,0) =   - (mu_2 * _Ex[_qp] * _Ey[_qp] / b_norm_2);
      D(1,1) =          1 - (mu_2 * _Ey[_qp] * _Ey[_qp] / b_norm_2);
      D(0,2) = D(2,0) =   - (mu_2 * _Ex[_qp] * _Ez[_qp] / b_norm_2);
      D(1,2) = D(2,1) =   - (mu_2 * _Ey[_qp] * _Ez[_qp] / b_norm_2);
      D(2,2) =          1 - (mu_2 * _Ez[_qp] * _Ez[_qp] / b_norm_2);

      return _scaling * epsilon_tilde * D * _grad_u[_qp] * _grad_test[_i][_qp];
    }
  else
    return 0.0;
}

Real DensityCrosswindDiffusion::computeQpJacobian()
{
  Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());

  RealVectorValue b;
  b(0) = -_mu * _Ex[_qp];
  b(1) = -_mu * _Ey[_qp];
  b(2) = -_mu * _Ez[_qp];

  Real epsilon_tilde = b.norm() * std::pow(h, 1.5) - _epsilon;

  if (epsilon_tilde > 0.0)
    {
      Real b_norm_2 = std::pow(b.norm(), 2.0);
      Real mu_2 = std::pow(_mu, 2.0);

      RealTensorValue D;
      D(0,0) =          1 - (mu_2 * _Ex[_qp] * _Ex[_qp] / b_norm_2);
      D(0,1) = D(1,0) =   - (mu_2 * _Ex[_qp] * _Ey[_qp] / b_norm_2);
      D(1,1) =          1 - (mu_2 * _Ey[_qp] * _Ey[_qp] / b_norm_2);
      D(0,2) = D(2,0) =   - (mu_2 * _Ex[_qp] * _Ez[_qp] / b_norm_2);
      D(1,2) = D(2,1) =   - (mu_2 * _Ey[_qp] * _Ez[_qp] / b_norm_2);
      D(2,2) =          1 - (mu_2 * _Ez[_qp] * _Ez[_qp] / b_norm_2);

      return _scaling * epsilon_tilde * D * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
    }
  else
    return 0.0;
}

Real DensityCrosswindDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _Ex_var)
    {
      Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());

      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real epsilon_tilde = b.norm() * std::pow(h, 1.5) - _epsilon;

      if (epsilon_tilde > 0.0)
        {
          Real b_norm_2 = std::pow(b.norm(), 2.0);
          Real b_norm_4 = std::pow(b.norm(), 4.0);
          Real mu_2 = std::pow(_mu, 2.0);
          Real mu_4 = std::pow(_mu, 4.0);

          Real d_epsilon_tilde_d_Exj = (mu_2 * _Ex[_qp] * std::pow(h, 1.5) / b.norm()) * _phi[_j][_qp];

          RealTensorValue D;
          D(0,0) =          1 - (mu_2 * _Ex[_qp] * _Ex[_qp] / b_norm_2);
          D(0,1) = D(1,0) =   - (mu_2 * _Ex[_qp] * _Ey[_qp] / b_norm_2);
          D(1,1) =          1 - (mu_2 * _Ey[_qp] * _Ey[_qp] / b_norm_2);
          D(0,2) = D(2,0) =   - (mu_2 * _Ex[_qp] * _Ez[_qp] / b_norm_2);
          D(1,2) = D(2,1) =   - (mu_2 * _Ey[_qp] * _Ez[_qp] / b_norm_2);
          D(2,2) =          1 - (mu_2 * _Ez[_qp] * _Ez[_qp] / b_norm_2);

          RealTensorValue d_D_d_Exj;
          d_D_d_Exj(0,0) =                  (2.0 * mu_4 * _Ex[_qp] * _Ex[_qp] * _Ex[_qp] / b_norm_4 - 2.0 * mu_2 * _Ex[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Exj(0,1) = d_D_d_Exj(1,0) = (2.0 * mu_4 * _Ex[_qp] * _Ex[_qp] * _Ey[_qp] / b_norm_4 - mu_2 * _Ey[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Exj(1,1) =                  (2.0 * mu_4 * _Ey[_qp] * _Ey[_qp] * _Ex[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Exj(0,2) = d_D_d_Exj(2,0) = (2.0 * mu_4 * _Ex[_qp] * _Ex[_qp] * _Ez[_qp] / b_norm_4 - mu_2 * _Ez[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Exj(1,2) = d_D_d_Exj(2,1) = (2.0 * mu_4 * _Ex[_qp] * _Ey[_qp] * _Ez[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Exj(2,2) =                  (2.0 * mu_4 * _Ez[_qp] * _Ez[_qp] * _Ex[_qp] / b_norm_4) * _phi[_j][_qp];

          return _scaling * (d_epsilon_tilde_d_Exj * D + epsilon_tilde * d_D_d_Exj) * _grad_u[_qp] * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else if (jvar == _Ey_var)
    {
      Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());

      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real epsilon_tilde = b.norm() * std::pow(h, 1.5) - _epsilon;

      if (epsilon_tilde > 0.0)
        {
          Real b_norm_2 = std::pow(b.norm(), 2.0);
          Real b_norm_4 = std::pow(b.norm(), 4.0);
          Real mu_2 = std::pow(_mu, 2.0);
          Real mu_4 = std::pow(_mu, 4.0);

          Real d_epsilon_tilde_d_Eyj = (mu_2 * _Ey[_qp] * std::pow(h, 1.5) / b.norm()) * _phi[_j][_qp];

          RealTensorValue D;
          D(0,0) =          1 - (mu_2 * _Ex[_qp] * _Ex[_qp] / b_norm_2);
          D(0,1) = D(1,0) =   - (mu_2 * _Ex[_qp] * _Ey[_qp] / b_norm_2);
          D(1,1) =          1 - (mu_2 * _Ey[_qp] * _Ey[_qp] / b_norm_2);
          D(0,2) = D(2,0) =   - (mu_2 * _Ex[_qp] * _Ez[_qp] / b_norm_2);
          D(1,2) = D(2,1) =   - (mu_2 * _Ey[_qp] * _Ez[_qp] / b_norm_2);
          D(2,2) =          1 - (mu_2 * _Ez[_qp] * _Ez[_qp] / b_norm_2);

          RealTensorValue d_D_d_Eyj;
          d_D_d_Eyj(0,0) =                  (2.0 * mu_4 * _Ex[_qp] * _Ex[_qp] * _Ey[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Eyj(0,1) = d_D_d_Eyj(1,0) = (2.0 * mu_4 * _Ey[_qp] * _Ey[_qp] * _Ex[_qp] / b_norm_4 - mu_2 * _Ex[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Eyj(1,1) =                  (2.0 * mu_4 * _Ey[_qp] * _Ey[_qp] * _Ey[_qp] / b_norm_4 - 2.0 * mu_2 * _Ey[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Eyj(0,2) = d_D_d_Eyj(2,0) = (2.0 * mu_4 * _Ex[_qp] * _Ey[_qp] * _Ez[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Eyj(1,2) = d_D_d_Eyj(2,1) = (2.0 * mu_4 * _Ey[_qp] * _Ey[_qp] * _Ez[_qp] / b_norm_4 - mu_2 * _Ez[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Eyj(2,2) =                  (2.0 * mu_4 * _Ez[_qp] * _Ez[_qp] * _Ey[_qp] / b_norm_4) * _phi[_j][_qp];

          return _scaling * (d_epsilon_tilde_d_Eyj * D + epsilon_tilde * d_D_d_Eyj) * _grad_u[_qp] * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else if (jvar == _Ez_var)
    {
      Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());

      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real epsilon_tilde = b.norm() * std::pow(h, 1.5) - _epsilon;

      if (epsilon_tilde > 0.0)
        {
          Real b_norm_2 = std::pow(b.norm(), 2.0);
          Real b_norm_4 = std::pow(b.norm(), 4.0);
          Real mu_2 = std::pow(_mu, 2.0);
          Real mu_4 = std::pow(_mu, 4.0);

          Real d_epsilon_tilde_d_Ezj = (mu_2 * _Ez[_qp] * std::pow(h, 1.5) / b.norm()) * _phi[_j][_qp];

          RealTensorValue D;
          D(0,0) =          1 - (mu_2 * _Ex[_qp] * _Ex[_qp] / b_norm_2);
          D(0,1) = D(1,0) =   - (mu_2 * _Ex[_qp] * _Ey[_qp] / b_norm_2);
          D(1,1) =          1 - (mu_2 * _Ey[_qp] * _Ey[_qp] / b_norm_2);
          D(0,2) = D(2,0) =   - (mu_2 * _Ex[_qp] * _Ez[_qp] / b_norm_2);
          D(1,2) = D(2,1) =   - (mu_2 * _Ey[_qp] * _Ez[_qp] / b_norm_2);
          D(2,2) =          1 - (mu_2 * _Ez[_qp] * _Ez[_qp] / b_norm_2);

          RealTensorValue d_D_d_Ezj;
          d_D_d_Ezj(0,0) =                  (2.0 * mu_4 * _Ex[_qp] * _Ex[_qp] * _Ez[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Ezj(0,1) = d_D_d_Ezj(1,0) = (2.0 * mu_4 * _Ex[_qp] * _Ey[_qp] * _Ez[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Ezj(1,1) =                  (2.0 * mu_4 * _Ey[_qp] * _Ey[_qp] * _Ez[_qp] / b_norm_4) * _phi[_j][_qp];
          d_D_d_Ezj(0,2) = d_D_d_Ezj(2,0) = (2.0 * mu_4 * _Ex[_qp] * _Ez[_qp] * _Ez[_qp] / b_norm_4 - mu_2 * _Ex[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Ezj(1,2) = d_D_d_Ezj(2,1) = (2.0 * mu_4 * _Ey[_qp] * _Ez[_qp] * _Ez[_qp] / b_norm_4 - mu_2 * _Ey[_qp] / b_norm_2) * _phi[_j][_qp];
          d_D_d_Ezj(2,2) =                  (2.0 * mu_4 * _Ez[_qp] * _Ez[_qp] * _Ez[_qp] / b_norm_4 - 2.0 * mu_2 * _Ez[_qp] / b_norm_2) * _phi[_j][_qp];

          return _scaling * (d_epsilon_tilde_d_Ezj * D + epsilon_tilde * d_D_d_Ezj) * _grad_u[_qp] * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else
    return 0.0;
}
