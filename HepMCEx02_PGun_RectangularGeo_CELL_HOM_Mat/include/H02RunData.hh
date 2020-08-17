//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file H02RunData.hh
/// \brief Definition of the H02RunData class

#ifndef H02RunData_h
#define H02RunData_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "G4Event.hh"

#include "H02TrackInformation.hh"
#include "H02FullTrajectoryInfo.hh"
#include "GlobalVariables.hh"
#include <array>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int kAbs = 0;
const G4int kGap = 1;
const G4int kDim = 2;

// const G4int kLayer(10); 
// const G4int K_NETA(20) ;
// const G4int K_NPHI(20) ;

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
/// 
/// In order to reduce the number of data members a 2-dimensions array 
/// is introduced for each quantity:
/// - fEdep[], fTrackLength[].
///
/// The data are collected step by step in H02SteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in B4EventAction.

class H02RunData : public G4Run
{
public:
  H02RunData();
  virtual ~H02RunData();

  void Add(G4int id, G4double de, G4double dl);
  void FillPerEvent();

  void SetTrackParentPerEvent(const G4Event* event);

  void AddCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu) ;
  void AddTotalEnergy(G4double ECh, G4double ENu);

  void Reset();

  // Get methods
  G4String  GetVolumeName(G4int id) const;
  G4double  GetEdep(G4int id) const;
  G4double  GetTrackLength(G4int id) const; 
  //std::vector<G4int> GetTrackParents() const ;
  std::vector<FullTrajectoryInfo> fAllTrajectoryInfo;

  std::vector<double> fTrackStepNumb, fTrackPDGCode, fTrackPDGCharge, fTrackEnergy;
  std::vector<double> fTrackPx, fTrackPy, fTrackPz ;
  std::vector<std::vector<double>> fTrackStepID, fTrackStepPDGID, fTrackStepEnergy, fTrackStepTime;
  std::vector<std::vector<double>> fTrackStepInitPosX, fTrackStepInitPosY, fTrackStepInitPosZ;
  std::vector<std::vector<double>> fTrackStepFinalPosX, fTrackStepFinalPosY, fTrackStepFinalPosZ;

  std::vector<double> fTruthParticleE, fTruthParticlePx, fTruthParticlePy, fTruthParticlePz, fTruthParticleM, fTruthParticlePDGID;  
  std::vector<double> fTruthParticleMother1, fTruthParticleMother2 ;

  double fCell_Energy[kLayer][K_NETA][K_NPHI] ;
  double fCh_Cell_Energy[kLayer][K_NETA][K_NPHI] ;
  double fNu_Cell_Energy[kLayer][K_NETA][K_NPHI] ;

  double fTotal_Ch_Energy, fTotal_Nu_Energy;
  

  float fVariable;

private:
  std::array<G4String, kDim>  fVolumeNames;
  std::array<G4double, kDim>  fEdep;
  std::array<G4double, kDim>  fTrackLength; 

  //std::vector<G4int> fTrackParents;
};

// inline functions

// inline std::vector<G4int> H02RunData::GetTrackParents() const {
//   return fTrackParents;
// }

inline void H02RunData::Add(G4int id, G4double de, G4double dl) {
  fEdep[id] += de; 
  fTrackLength[id] += dl;
}

inline void H02RunData::AddTotalEnergy(G4double ECh, G4double ENu) {

  fTotal_Ch_Energy +=  ECh ;
  fTotal_Nu_Energy  += ENu ;
}

inline void H02RunData::AddCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu) {
  fCell_Energy[iLayer][iCellEta][iCellPhi] += ETot ;
  fCh_Cell_Energy[iLayer][iCellEta][iCellPhi] += ECh ;
  fNu_Cell_Energy[iLayer][iCellEta][iCellPhi] += ENu ;
}

inline G4String  H02RunData::GetVolumeName(G4int id) const {
  return fVolumeNames[id];
}

inline G4double  H02RunData::GetEdep(G4int id) const {
  return fEdep[id];
}   

inline G4double  H02RunData::GetTrackLength(G4int id) const {
  return fTrackLength[id];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

