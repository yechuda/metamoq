#include "MetamoqApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "CoupledSpaceChargeDensity.h"
#include "CoupledPotentialGradient.h"
#include "DensityDiffusion.h"
#include "DriftFluxBC.h"
#include "BodyForceComponentAux.h"
#include "InjectionDirichletBC.h"
#include "ConditionalSideIntegralVariablePostprocessor.h"
#include "ConditionalSideAverageValue.h"
#include "ElectricFieldBodyForce.h"
#include "InjectionPenaltyBC.h"
#include "InjectionImplicitBC.h"
#include "InjectionPenaltyNormalBC.h"
#include "InjectionRatioPenaltyNormalBC.h"
#include "InjectionPeekConstantBC.h"
#include "InjectionPeekVariableBC.h"
#include "BodyForceFilterAux.h"
#include "DriftDiffusionLog.h"
#include "CoupledSpaceChargeDensityLog.h"
#include "DriftFluxLogBC.h"
#include "InjectionPeekConstantLogBC.h"
#include "InjectionPeekVariableLogBC.h"
#include "DensityAux.h"
#include "DensityArtificialDiffusion.h"
#include "PecletAux.h"
#include "DensityArtificialDiffusionTreshold.h"
#include "ValueJumpIndicator.h"
#include "LogStabilization.h"
#include "InjectionPeekConstantDampedBC.h"
#include "InjectionPeekVariableDampedBC.h"
#include "DensityCrosswindDiffusion.h"
#include "DensityPenalty.h"
#include "AbsDensityAux.h"
#include "PoissonOffset.h"
#include "DriftDiffusionOffset.h"
#include "DriftFluxOffsetBC.h"
#include "DensityOffsetAux.h"
#include "InjectionPeekConstantOffsetBC.h"
#include "InjectionPeekVariableOffsetBC.h"
#include "DensitySUPG.h"
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
  registerKernel(CoupledSpaceChargeDensity);
  registerKernel(CoupledPotentialGradient);
  registerKernel(DensityDiffusion);
  registerBoundaryCondition(DriftFluxBC);
  registerAux(BodyForceComponentAux);
  registerBoundaryCondition(InjectionDirichletBC);
  registerPostprocessor(ConditionalSideIntegralVariablePostprocessor);
  registerPostprocessor(ConditionalSideAverageValue);
  registerKernel(ElectricFieldBodyForce);
  registerBoundaryCondition(InjectionPenaltyBC);
  registerBoundaryCondition(InjectionImplicitBC);
  registerBoundaryCondition(InjectionPenaltyNormalBC);
  registerBoundaryCondition(InjectionRatioPenaltyNormalBC);
  registerBoundaryCondition(InjectionPeekConstantBC);
  registerBoundaryCondition(InjectionPeekVariableBC);
  registerAux(BodyForceFilterAux);
  registerKernel(DriftDiffusionLog);
  registerKernel(CoupledSpaceChargeDensityLog);
  registerBoundaryCondition(DriftFluxLogBC);
  registerBoundaryCondition(InjectionPeekConstantLogBC);
  registerBoundaryCondition(InjectionPeekVariableLogBC);
  registerAux(DensityAux);
  registerKernel(DensityArtificialDiffusion);
  registerAux(PecletAux);
  registerKernel(DensityArtificialDiffusionTreshold);
  registerIndicator(ValueJumpIndicator);
  registerKernel(LogStabilization);
  registerBoundaryCondition(InjectionPeekConstantDampedBC);
  registerBoundaryCondition(InjectionPeekVariableDampedBC);
  registerKernel(DensityCrosswindDiffusion);
  registerKernel(DensityPenalty);
  registerAux(AbsDensityAux);
  registerKernel(PoissonOffset);
  registerKernel(DriftDiffusionOffset);
  registerBoundaryCondition(DriftFluxOffsetBC);
  registerAux(DensityOffsetAux);
  registerBoundaryCondition(InjectionPeekConstantOffsetBC);
  registerBoundaryCondition(InjectionPeekVariableOffsetBC);
  registerKernel(DensitySUPG);
  registerPostprocessor(CurrentPostprocessor);
}

// External entry point for dynamic syntax association
extern "C" void MetamoqApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { MetamoqApp::associateSyntax(syntax, action_factory); }
void
MetamoqApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
