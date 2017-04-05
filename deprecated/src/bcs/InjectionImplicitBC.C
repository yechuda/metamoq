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
#include "InjectionImplicitBC.h"
#include "Function.h"

template<>
InputParameters validParams<InjectionImplicitBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("E0", "Electric field strength at corona onset");
  params.addRequiredCoupledVar("E_magnitude", "The magnitude of local electric field strength");

  return params;
}

InjectionImplicitBC::InjectionImplicitBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _E0(getParam<Real>("E0")),
    _E_mag_var(coupled("E_magnitude")),
    _E_mag(coupledValue("E_magnitude"))
{}

Real
InjectionImplicitBC::computeQpResidual()
{
  return _test[_i][_qp]*_u[_qp]*(_E_mag[_qp]-_E0);
}

Real
InjectionImplicitBC::computeQpJacobian()
{
  return _test[_i][_qp]*_phi[_j][_qp]*(_E_mag[_qp]-_E0);
}

Real
InjectionImplicitBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _E_mag_var)
      return _test[_i][_qp]*_u[_qp]*_phi[_j][_qp];
  else
      return 0.0;
}
