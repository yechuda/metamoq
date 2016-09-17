#include "MetamoqApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "CoupledSpaceChargeDensity.h"
#include "CoupledPotentialGradient.h"
#include "DensityDiffusion.h"
#include "DriftFluxBC.h"
#include "BodyForceComponent.h"
#include "InjectionDirichletBC.h"
#include "ConditionalSideIntegralVariablePostprocessor.h"
#include "ConditionalSideAverageValue.h"

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
  registerAux(BodyForceComponent);
  registerBoundaryCondition(InjectionDirichletBC);
  registerPostprocessor(ConditionalSideIntegralVariablePostprocessor);
  registerPostprocessor(ConditionalSideAverageValue);
}

// External entry point for dynamic syntax association
extern "C" void MetamoqApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { MetamoqApp::associateSyntax(syntax, action_factory); }
void
MetamoqApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
