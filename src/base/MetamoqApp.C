#include "MetamoqApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Kernels
#include "InverseWallDistance.h"
#include "CoupledSpaceChargeDensity.h"
#include "DensityDiffusion.h"
#include "DensityArtificialDiffusion.h"
#include "CoupledPotentialGradient.h"
#include "BodyForceComponent.h"
#include "INSMomentumTractionFormRANS.h"
#include "INSMomentumTractionFormRANSRZ.h"
#include "ElectricFieldBodyForceExplicit.h"
#include "ElectricFieldBodyForceExplicitRamp.h"
#include "ElectricFieldBodyForceExplicitTimeRamp.h"
#include "NeeKovasznay.h"
#include "NeeKovasznayProductionEHD.h"
#include "FaresSchroderSpecificTurbulenceDissipationRate.h"
#include "FaresSchroder.h"
#include "FaresSchroder88.h"
#include "FaresSchroderProductionEHD.h"
#include "BodyForceVorticitySquareRootProductionEHD.h"
#include "TurbulenceKineticEnergyTransport.h"
#include "SpecificTurbulenceDissipationRateTransport.h"
#include "FaresSchroder88Full.h"
#include "FaresSchroderArtificialDiffusion.h"
#include "DensityAdvection.h"
#include "DensityArtificialDiffusionCoupled.h"
#include "INSMomentumTimeDerivativeParam.h"
#include "INSMassRANS.h"
#include "INSMassRANSRZ.h"

// Auxkernels
#include "ApparentDynamicViscosityAux.h"
#include "WallDistanceAux.h"
#include "FaresSchroderDynamicViscosityAux.h"
#include "komegaDynamicViscosityAux.h"
#include "SpecificTurbulenceDissipationRateAux.h"

// BCs
#include "OnePointBoundedInverseDistanceDirichletBC.h"
#include "TwoPointsMinInverseDistanceDirichletBC.h"
#include "InjectionPeekConstantDampedBC.h"
#include "InjectionPeekVariableDampedBC.h"
#include "DriftFluxBC.h"
#include "BodyForceBC.h"
#include "NeeKovasznayNoBCBC.h"
#include "FaresSchroderSpecificTurbulenceDissipationRateBC.h"
#include "CathodeFluxBC.h"

// Postprocessors
#include "CurrentPostprocessor.h"

template<>
InputParameters validParams<MetamoqApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

MetamoqApp::MetamoqApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  MetamoqApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  MetamoqApp::associateSyntax(_syntax, _action_factory);
}

MetamoqApp::~MetamoqApp()
{
}

// External entry point for dynamic application loading
extern "C" void MetamoqApp__registerApps() { MetamoqApp::registerApps(); }
void
MetamoqApp::registerApps()
{
  registerApp(MetamoqApp);
}

// External entry point for dynamic object registration
extern "C" void MetamoqApp__registerObjects(Factory & factory) { MetamoqApp::registerObjects(factory); }
void
MetamoqApp::registerObjects(Factory & factory)
{
  // Kernels
  registerKernel(InverseWallDistance);
  registerKernel(CoupledSpaceChargeDensity);
  registerKernel(DensityDiffusion);
  registerKernel(DensityArtificialDiffusion);
  registerKernel(CoupledPotentialGradient);
  registerKernel(BodyForceComponent);
  registerKernel(INSMomentumTractionFormRANS);
  registerKernel(INSMomentumTractionFormRANSRZ);
  registerKernel(ElectricFieldBodyForceExplicit);
  registerKernel(ElectricFieldBodyForceExplicitRamp);
  registerKernel(ElectricFieldBodyForceExplicitTimeRamp);
  registerKernel(NeeKovasznay);
  registerKernel(NeeKovasznayProductionEHD);
  registerKernel(FaresSchroderSpecificTurbulenceDissipationRate);
  registerKernel(FaresSchroder);
  registerKernel(FaresSchroder88);
  registerKernel(FaresSchroderProductionEHD);
  registerKernel(BodyForceVorticitySquareRootProductionEHD);
  registerKernel(TurbulenceKineticEnergyTransport);
  registerKernel(SpecificTurbulenceDissipationRateTransport);
  registerKernel(FaresSchroder88Full);
  registerKernel(FaresSchroderArtificialDiffusion);
  registerKernel(DensityAdvection);
  registerKernel(DensityArtificialDiffusionCoupled);
  registerKernel(INSMomentumTimeDerivativeParam);
  registerKernel(INSMassRANS);
  registerKernel(INSMassRANSRZ);

  // Auxkernels
  registerAux(ApparentDynamicViscosityAux);
  registerAux(WallDistanceAux);
  registerAux(FaresSchroderDynamicViscosityAux);
  registerAux(komegaDynamicViscosityAux);
  registerAux(SpecificTurbulenceDissipationRateAux);

  // BCs
  registerBoundaryCondition(OnePointBoundedInverseDistanceDirichletBC);
  registerBoundaryCondition(TwoPointsMinInverseDistanceDirichletBC);
  registerBoundaryCondition(InjectionPeekConstantDampedBC);
  registerBoundaryCondition(InjectionPeekVariableDampedBC);
  registerBoundaryCondition(DriftFluxBC);
  registerBoundaryCondition(BodyForceBC);
  registerBoundaryCondition(NeeKovasznayNoBCBC);
  registerBoundaryCondition(FaresSchroderSpecificTurbulenceDissipationRateBC);
  registerBoundaryCondition(CathodeFluxBC);

  // Postprocessors
  registerPostprocessor(CurrentPostprocessor);
}

// External entry point for dynamic syntax association
extern "C" void MetamoqApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { MetamoqApp::associateSyntax(syntax, action_factory); }
void
MetamoqApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
