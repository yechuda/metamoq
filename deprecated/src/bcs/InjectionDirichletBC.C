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

#include "InjectionDirichletBC.h"

template<>
InputParameters validParams<InjectionDirichletBC>()
{
  InputParameters p = validParams<NodalBC>();
  p.addRequiredParam<Real>("coef", "Space charge density scaling coefficient");
  p.addRequiredParam<Real>("E0", "Electric field strength at corona onset");
  p.addRequiredCoupledVar("E_magnitude_Laplace", "The magnitude of local electric field strength from Laplace solution");
  return p;
}


InjectionDirichletBC::InjectionDirichletBC(const InputParameters & parameters) :
  NodalBC(parameters),
  _coef(getParam<Real>("coef")),
  _E0(getParam<Real>("E0")),
  _E_mag(coupledValue("E_magnitude_Laplace"))
{}

Real
InjectionDirichletBC::computeQpResidual()
{
  if (_E_mag[_qp] > _E0)
    return _u[_qp] - _coef*(_E_mag[_qp]-_E0);
  else
    return _u[_qp];
}
