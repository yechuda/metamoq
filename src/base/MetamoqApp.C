#include "MetamoqApp.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

#include "NavierStokesApp.h"

// #include "Moose.h"
// #include "ModulesApp.h"

// Kernels
#include "CoupledSpaceChargeDensity.h"
#include "DensityDiffusion.h"
#include "DensityArtificialDiffusion.h"
#include "CoupledPotentialGradient.h"
#include "BodyForceComponent.h"
#include "ElectricFieldBodyForceExplicit.h"

// Materials
#include "air.h"

// BCs
#include "DriftFluxBC.h"
#include "BodyForceBC.h"
#include "InjectionPeekConstantDampedTunedBC.h"
#include "InjectionPeekVariableDampedTunedBC.h"
#include "InjectionTresholdBC.h"

// Postprocessors
#include "CurrentPostprocessor.h"

InputParameters
MetamoqApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  params.set<bool>("automatic_automatic_scaling") = false;

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  // params.set<bool>("use_legacy_dirichlet_bc") = false;
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

registerKnownLabel("MetamoqApp");

MetamoqApp::MetamoqApp(const InputParameters & parameters) : MooseApp(parameters)
{
  MetamoqApp::registerAll(_factory, _action_factory, _syntax);
}

MetamoqApp::~MetamoqApp()
{
}

void
MetamoqApp::registerApps()
{
  registerApp(MetamoqApp);
}

void
MetamoqApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"MetamoqApp"});
  Registry::registerActionsTo(af, {"MetamoqApp"});

  NavierStokesApp::registerAll(f, af, s);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
MetamoqApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MetamoqApp::registerAll(f, af, s);
}

extern "C" void
MetamoqApp__registerApps()
{
  MetamoqApp::registerApps();
}
