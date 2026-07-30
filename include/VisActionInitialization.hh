// -*- C++ -*-
//
// Minimal action initialization for the visualization-only binary
// (G4S2SVis): only a trivial primary generator, no run/event/stepping
// actions -- deliberately avoids S2SRunAction/S2SAnaManager, which
// are not safe to invoke twice (master + worker) under an MT run
// manager (see analysis-note.md).

#ifndef VIS_ACTION_INITIALIZATION_HH
#define VIS_ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

//_____________________________________________________________________________
class VisActionInitialization : public G4VUserActionInitialization
{
public:
  VisActionInitialization();
  ~VisActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;
};

#endif
