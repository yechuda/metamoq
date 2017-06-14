/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "FaresSchroderProductionEHD.h"

template<>
InputParameters validParams<FaresSchroderProductionEHD>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("potential", "electric potential");
  params.addRequiredCoupledVar("space_charge_density", "space charge density");
  params.addRequiredCoupledVar("omega", "specific turbulence dissipation rate");

  // Required parameters
  params.addParam<Real>("c", 0.03, "EHD production coefficient");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mobility", "ion mobility coefficient");

  return params;
}



FaresSchroderProductionEHD::FaresSchroderProductionEHD(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _space_charge_density(coupledValue("space_charge_density")),
  _omega(coupledValue("omega")),

  // Coupled gradients
  _grad_potential(coupledGradient("potential")),

  // Variable numberings
  _potential_var(coupled("potential")),
  _space_charge_density_var(coupled("space_charge_density")),
  _omega_var(coupled("omega")),

  // Required parameters
  _c(getParam<Real>("c")),
  _rho(getParam<Real>("rho")),
  _mobility(getParam<Real>("mobility"))

{
}



Real FaresSchroderProductionEHD::computeQpResidual()
{
  if (_t_step == 1)
  {
    return 0.0;
  }
  else
  {
    Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
    Real E_squared = _grad_potential[_qp](0) * _grad_potential[_qp](0) +
                     _grad_potential[_qp](1) * _grad_potential[_qp](1) +
                     _grad_potential[_qp](2) * _grad_potential[_qp](2);
    return -_c / (omega * _rho) * _space_charge_density[_qp] * _mobility * E_squared * _test[_i][_qp];
  }
}



Real FaresSchroderProductionEHD::computeQpJacobian()
{
  return 0.0;
}

Real FaresSchroderProductionEHD::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _potential_var)
  {
    if (_t_step == 1)
    {
      return 0.0;
    }
    else
    {
      Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
      return -_c / (omega * _rho) * _space_charge_density[_qp] * _mobility * 2.0 * _grad_potential[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp];
    }
  }

  else if (jvar == _space_charge_density_var)
  {
    if (_t_step == 1)
    {
      return 0.0;
    }
    else
    {
      Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
      Real E_squared = _grad_potential[_qp](0) * _grad_potential[_qp](0) +
                       _grad_potential[_qp](1) * _grad_potential[_qp](1) +
                       _grad_potential[_qp](2) * _grad_potential[_qp](2);
      return -_c / (omega * _rho) * _phi[_j][_qp] * _mobility * E_squared * _test[_i][_qp];
    }
  }

  else if (jvar == _omega_var)
  {
    if (_t_step == 1)
    {
      return 0.0;
    }
    else
    {
      Real omega = std::max(_omega[_qp], std::numeric_limits<Real>::epsilon());
      Real E_squared = _grad_potential[_qp](0) * _grad_potential[_qp](0) +
                       _grad_potential[_qp](1) * _grad_potential[_qp](1) +
                       _grad_potential[_qp](2) * _grad_potential[_qp](2);
      return _c / (omega * omega * _rho) * _phi[_j][_qp] * _space_charge_density[_qp] * _mobility * E_squared * _test[_i][_qp];
    }
  }

  else
    return 0.0;
}
