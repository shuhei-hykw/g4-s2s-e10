// -*- C++ -*-

#include "VisActionInitialization.hh"
#include "VisPrimaryGeneratorAction.hh"

//_____________________________________________________________________________
VisActionInitialization::VisActionInitialization()
{
}

//_____________________________________________________________________________
VisActionInitialization::~VisActionInitialization()
{}

//_____________________________________________________________________________
void
VisActionInitialization::BuildForMaster() const
{
}

//_____________________________________________________________________________
void
VisActionInitialization::Build() const
{
  SetUserAction(new VisPrimaryGeneratorAction);
}
