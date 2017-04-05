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

#include "DensitySUPG.h"

template<>
InputParameters validParams<DensitySUPG>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("Ex", "The component of electric field in x direction");
  params.addCoupledVar("Ey", 0.0, "The component of electric field in y direction"); // only required in 2D and 3D
  params.addCoupledVar("Ez", 0.0, "The component of electric field in z direction"); // only required in 3D

  // Required parameters
  params.addRequiredParam<Real>("mobility", "Ion mobility coefficient");
  params.addRequiredParam<Real>("charge_diffusion_coefficient", "Charge diffusion coefficient");
  params.addRequiredParam<Real>("scaling", "Scaling coefficient");

  return params;
}

DensitySUPG::DensitySUPG(const InputParameters & parameters) :
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
    _scaling(getParam<Real>("scaling")),

    // Second derivative of nonlinear variable
    _second_u(second()),
    _second_phi(secondPhi())

{
}

Real DensitySUPG::computeQpResidual()
{
  RealVectorValue b;
  b(0) = -_mu * _Ex[_qp];
  b(1) = -_mu * _Ey[_qp];
  b(2) = -_mu * _Ez[_qp];

  Real R = -_epsilon * _second_u[_qp].tr() + b * _grad_u[_qp];
  RealVectorValue z = (R / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];

  if ((b.norm() / z.norm() - 1.0) > 0.0)
    {
      Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());
      Real Pe = (b.norm() * h) / (2.0 * _epsilon);
      Real tau = (h / (2.0 * b.norm())) * (1.0 / std::tanh(Pe) - 1.0 / Pe);
      Real sigma = tau * (b.norm() / z.norm() - 1.0);
      return _scaling * R * sigma * z * _grad_test[_i][_qp];
    }
  else
    return 0.0;
}

Real DensitySUPG::computeQpJacobian()
{
  RealVectorValue b;
  b(0) = -_mu * _Ex[_qp];
  b(1) = -_mu * _Ey[_qp];
  b(2) = -_mu * _Ez[_qp];

  Real R = -_epsilon * _second_u[_qp].tr() + b * _grad_u[_qp];
  RealVectorValue z = (R / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];

  if ((b.norm() / z.norm() - 1.0) > 0.0)
    {
      Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());
      Real Pe = (b.norm() * h) / (2.0 * _epsilon);
      Real tau = (h / (2.0 * b.norm())) * (1.0 / std::tanh(Pe) - 1.0 / Pe);
      Real sigma = tau * (b.norm() / z.norm() - 1.0);

      Real dR_duj = -_epsilon * _second_phi[_j][_qp].tr() + b * _grad_phi[_j][_qp];
      RealVectorValue dz_duj = (dR_duj * _grad_u[_qp] / b.norm()) + (_grad_phi[_j][_qp] * R / std::pow(_grad_u[_qp].norm(), 2.0)) - (R * _grad_u[_qp] * (2.0 * _grad_u[_qp] * _grad_phi[_j][_qp] / std::pow(_grad_u[_qp] * _grad_u[_qp], 2.0)));
      Real dsigma_duj = -tau * b.norm() * std::pow(z * z, -1.5) * dz_duj * z;

      return _scaling * (dR_duj * sigma * z + dz_duj * R * sigma + dsigma_duj * R * z) * _grad_test[_i][_qp];
    }
  else
    return 0.0;
}

Real DensitySUPG::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _Ex_var)
    {
      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real R = -_epsilon * _second_u[_qp].tr() + b * _grad_u[_qp];
      RealVectorValue z = (R / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];

      if ((b.norm() / z.norm() - 1.0) > 0.0)
        {
          Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());
          Real Pe = (b.norm() * h) / (2.0 * _epsilon);
          Real tau = (h / (2.0 * b.norm())) * (1.0 / std::tanh(Pe) - 1.0 / Pe);
          Real sigma = tau * (b.norm() / z.norm() - 1.0);

          Real dR_dExj = -_mu * _phi[_j][_qp] * _grad_u[_qp](0);
          RealVectorValue dz_dExj = (dR_dExj / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];
          Real dPe_dExj = (h * _Ex[_qp] * std::pow(_mu, 2.0) / (2.0 * _epsilon * b.norm())) * _phi[_j][_qp];
          Real dtau_dExj = (-h * _Ex[_qp] * std::pow(_mu, 2.0) / (2.0 * std::pow(b.norm(), 3.0))) * (1.0 / std::tanh(Pe) - 1.0 / Pe) * _phi[_j][_qp] + (h / (2.0 * b.norm())) * (1.0 - std::pow(1.0 / std::tanh(Pe), 2.0) + 1.0 / std::pow(Pe, 2.0)) * dPe_dExj;
          Real dbnorm_dExj = (_Ex[_qp] * std::pow(_mu, 2.0) / b.norm()) * _phi[_j][_qp];
          Real dznormrecip_dExj = -std::pow(z * z, -1.5) * z * dz_dExj;
          Real dsigma_dExj = dtau_dExj * b.norm() / z.norm() + dbnorm_dExj * tau / z.norm() + dznormrecip_dExj * tau * b.norm() - dtau_dExj;

          return _scaling * (dR_dExj * sigma * z + dz_dExj * R * sigma + dsigma_dExj * R * z) * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else if (jvar == _Ey_var)
    {
      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real R = -_epsilon * _second_u[_qp].tr() + b * _grad_u[_qp];
      RealVectorValue z = (R / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];

      if ((b.norm() / z.norm() - 1.0) > 0.0)
        {
          Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());
          Real Pe = (b.norm() * h) / (2.0 * _epsilon);
          Real tau = (h / (2.0 * b.norm())) * (1.0 / std::tanh(Pe) - 1.0 / Pe);
          Real sigma = tau * (b.norm() / z.norm() - 1.0);

          Real dR_dEyj = -_mu * _phi[_j][_qp] * _grad_u[_qp](1);
          RealVectorValue dz_dEyj = (dR_dEyj / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];
          Real dPe_dEyj = (h * _Ey[_qp] * std::pow(_mu, 2.0) / (2.0 * _epsilon * b.norm())) * _phi[_j][_qp];
          Real dtau_dEyj = (-h * _Ey[_qp] * std::pow(_mu, 2.0) / (2.0 * std::pow(b.norm(), 3.0))) * (1.0 / std::tanh(Pe) - 1.0 / Pe) * _phi[_j][_qp] + (h / (2.0 * b.norm())) * (1.0 - std::pow(1.0 / std::tanh(Pe), 2.0) + 1.0 / std::pow(Pe, 2.0)) * dPe_dEyj;
          Real dbnorm_dEyj = (_Ey[_qp] * std::pow(_mu, 2.0) / b.norm()) * _phi[_j][_qp];
          Real dznormrecip_dEyj = -std::pow(z * z, -1.5) * z * dz_dEyj;
          Real dsigma_dEyj = dtau_dEyj * b.norm() / z.norm() + dbnorm_dEyj * tau / z.norm() + dznormrecip_dEyj * tau * b.norm() - dtau_dEyj;

          return _scaling * (dR_dEyj * sigma * z + dz_dEyj * R * sigma + dsigma_dEyj * R * z) * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else if (jvar == _Ez_var)
    {
      RealVectorValue b;
      b(0) = -_mu * _Ex[_qp];
      b(1) = -_mu * _Ey[_qp];
      b(2) = -_mu * _Ez[_qp];

      Real R = -_epsilon * _second_u[_qp].tr() + b * _grad_u[_qp];
      RealVectorValue z = (R / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];

      if ((b.norm() / z.norm() - 1.0) > 0.0)
        {
          Real h = 0.5 * (_current_elem->hmin() + _current_elem->hmax());
          Real Pe = (b.norm() * h) / (2.0 * _epsilon);
          Real tau = (h / (2.0 * b.norm())) * (1.0 / std::tanh(Pe) - 1.0 / Pe);
          Real sigma = tau * (b.norm() / z.norm() - 1.0);

          Real dR_dEzj = -_mu * _phi[_j][_qp] * _grad_u[_qp](2);
          RealVectorValue dz_dEzj = (dR_dEzj / std::pow(_grad_u[_qp].norm(), 2.0)) * _grad_u[_qp];
          Real dPe_dEzj = (h * _Ez[_qp] * std::pow(_mu, 2.0) / (2.0 * _epsilon * b.norm())) * _phi[_j][_qp];
          Real dtau_dEzj = (-h * _Ez[_qp] * std::pow(_mu, 2.0) / (2.0 * std::pow(b.norm(), 3.0))) * (1.0 / std::tanh(Pe) - 1.0 / Pe) * _phi[_j][_qp] + (h / (2.0 * b.norm())) * (1.0 - std::pow(1.0 / std::tanh(Pe), 2.0) + 1.0 / std::pow(Pe, 2.0)) * dPe_dEzj;
          Real dbnorm_dEzj = (_Ez[_qp] * std::pow(_mu, 2.0) / b.norm()) * _phi[_j][_qp];
          Real dznormrecip_dEzj = -std::pow(z * z, -1.5) * z * dz_dEzj;
          Real dsigma_dEzj = dtau_dEzj * b.norm() / z.norm() + dbnorm_dEzj * tau / z.norm() + dznormrecip_dEzj * tau * b.norm() - dtau_dEzj;

          return _scaling * (dR_dEzj * sigma * z + dz_dEzj * R * sigma + dsigma_dEzj * R * z) * _grad_test[_i][_qp];
        }
      else
        return 0.0;
    }
  else
    return 0.0;
}
