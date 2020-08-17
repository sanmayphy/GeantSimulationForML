#ifndef H02FullTrajectoryInfo_h
#define H02FullTrajectoryInfo_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include "G4TrajectoryContainer.hh"


struct FullTrajectoryInfo {
G4int fTrackID;
G4int fParentID;
G4int fPDGCode;
G4int fprocessId;
G4double fPDGCharge;
G4ThreeVector fMomentum;
G4ThreeVector fMomentumDir;
G4double    fEnergy ;
G4ThreeVector fVertexPosition;
G4double fGlobalTime;

std::vector<G4int> vTrackID;
std::vector<G4int> vParentID;
std::vector<G4ThreeVector> vTrackMomentumDir;
std::vector<G4ThreeVector> vPreStepTrackPos; 
std::vector<G4ThreeVector> vPostStepTrackPos;
std::vector<G4ThreeVector> vTrackPos;
std::vector<G4double> vTrackTime;
std::vector<G4double> vTrackEnergy;
std::vector<G4double> vTrackPDGID;

}; 


#endif
