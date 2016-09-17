#ifndef METAMOQAPP_H
#define METAMOQAPP_H

#include "MooseApp.h"

class MetamoqApp;

template<>
InputParameters validParams<MetamoqApp>();

class MetamoqApp : public MooseApp
{
public:
  MetamoqApp(InputParameters parameters);
  virtual ~MetamoqApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* METAMOQAPP_H */
