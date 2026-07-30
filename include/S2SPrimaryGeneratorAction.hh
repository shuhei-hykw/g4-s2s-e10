// -*- C++ -*-

#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <globals.hh>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;

//_____________________________________________________________________________
class S2SPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  static G4String ClassName();
  S2SPrimaryGeneratorAction();
  ~S2SPrimaryGeneratorAction();

private:
  G4ParticleGun* m_particleGun;
  G4int          m_generator;

public:
  TFile *profileK18;
  TTree *k18track;
  G4double p_3rd[500], xtgt[500], ytgt[500], utgt[500], vtgt[500];
  G4int trigflag[32];
  G4int ntK18;
  G4double chisqrK18[500], CBtof0[500];
  TGraph *gr;
  G4double Sum;
  G4int inum;
  
public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  void GenerateDemo(G4Event* anEvent);
  void GenerateUniform0(G4Event* anEvent);
  void GenerateFocusCheck(G4Event* anEvent);
  void GenerateMonochromeBeam(G4Event* anEvent);
  void GenerateUniformSpherical(G4Event* anEvent);
  void GenerateBeam(G4Event* anEvent);
  void GenerateBeamThrough(G4Event* anEvent);
  void GenerateBeamGausProfile(G4Event* anEvent);
  void GenerateBeamFixSeed(G4Event* anEvent);
  void GenerateScatParticles(G4Event* anEvent);
  void GenerateDefocusBeam(G4Event* anEvent);
  void Generate12XiBeryllium(G4Event* anEvent);
  void GenerateElementaryXiMinus(G4Event* anEvent);
  void GenerateElementarySigmaMinus(G4Event* anEvent);
  void GenerateElementarySigmaPlus(G4Event* anEvent);
  void Generate7XiH(G4Event* anEvent);
  void GenerateRich7XiHSpectrum(G4Event* anEvent);
  void GenerateKH7XiHSpectrum(G4Event* anEvent);
  void GenerateE63_7LambdaLi(G4Event* anEvent, G4int MassNum);
  void GenerateProton(G4Event* anEvent);
  void GenerateSigmaNCusp(G4Event* anEvent);
  void GenerateQFLambda(G4Event* anEvent);
  void GenerateQFSigmaZ(G4Event* anEvent);
  void GenerateQFSigmaP(G4Event* anEvent);
  void GenerateK0Production(G4Event* anEvent);
  void GenerateK0CoherentProduction(G4Event* anEvent);
  void GenerateE63_LiProfileMonoGamma(G4Event* anEvent);
  void GenerateE63_4LHGamma(G4Event* anEvent);
};

//_____________________________________________________________________________
inline G4String
S2SPrimaryGeneratorAction::ClassName()
{
  static G4String s_name("S2SPrimaryGeneratorAction");
  return s_name;
}

#endif
