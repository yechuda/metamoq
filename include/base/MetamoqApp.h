#ifndef METAMOQAPP_H
#define METAMOQAPP_H

#include "MooseApp.h"

class MetamoqApp : public MooseApp
{
public:
  static InputParameters validParams();

  MetamoqApp(const InputParameters & parameters);
  virtual ~MetamoqApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* METAMOQAPP_H */
