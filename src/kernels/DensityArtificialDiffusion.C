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

#include "DensityArtificialDiffusion.h"

template<>
InputParameters validParams<DensityArtificialDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("potential","The potential for calculating the advection velocity.");
  params.addParam<Real>("delta",0.25,"Tuning parameter.");
  params.addRequiredParam<Real>("mobility","Ion mobility coefficient");
  return params;
}

DensityArtificialDiffusion::DensityArtificialDiffusion(const InputParameters & parameters) :
    Kernel(parameters),

    _grad_potential(coupledGradient("potential")),
    _potential_var(coupled("potential")),
    _delta(getParam<Real>("delta")),
    _mu(getParam<Real>("mobility"))

{
}

DensityArtificialDiffusion::~DensityArtificialDiffusion()
{
}

Real DensityArtificialDiffusion::computeQpResidual()
{
  Real  _beta = _mu * _grad_potential[_qp].norm();
  Real  _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real DensityArtificialDiffusion::computeQpJacobian()
{
  Real  _beta = _mu * _grad_potential[_qp].norm();
  Real  _D_art = _beta * _current_elem->hmax() * _delta;

  return _D_art * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}

Real DensityArtificialDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_var)
  {
    Real  _d_beta_d_potential = _mu * _grad_potential[_qp] * _grad_phi[_j][_qp] / (_grad_potential[_qp].norm()+std::numeric_limits<double>::epsilon());
    Real _d_D_art_d_potential = _d_beta_d_potential * _current_elem->hmax() * _delta;

    return _d_D_art_d_potential * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  else
    return 0.0;
}
