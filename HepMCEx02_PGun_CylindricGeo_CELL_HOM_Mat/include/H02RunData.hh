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
#include "H02RunAction.hh"
#include "H02TrackInformation.hh"
#include "H02FullTrajectoryInfo.hh"
#include "GlobalVariables.hh"
#include <array>
// #include "PCell.hh"

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

  void AddCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu);//,G4double EHadCh, G4double EHadNu, G4double EEM) ;
  // void AddPCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu);
  void AddTotalEnergy(G4double ECh, G4double ENu, G4double Eesc, G4double Enotrc);//,G4double EHadCh, G4double EHadNu, G4double EEM);
  void AddTrueEnergy(G4double ECh, G4double ENu);
  void AddPCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu, MParticle cell);
  void Reset();
  inline void Mass( G4double dn);
  inline void SumNbRadLength( G4double dn);

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

  PCell fPArray[6][kMaxPixel][kMaxPixel];

  double fCell_Energy[kLayer][K_NETA][K_NPHI] ;
  double fCh_Cell_Energy[kLayer][K_NETA][K_NPHI] ;
  double fNu_Cell_Energy[kLayer][K_NETA][K_NPHI] ;
  double fEnergy_in_gap;
  double fTrue_Ch_Energy, fTrue_Nu_Energy;

  double fTotal_Ch_Energy, fTotal_Nu_Energy, fTotal_Esc_Energy, fTotal_Enotrc_Energy,fKin_Ch_Energy, fKin_Nu_Energy;
  double fTotal_Had_Ch_Energy, fTotal_Had_Nu_Energy, fTotal_EM_Energy;
  double fnbRadLen;

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
inline void H02RunData::SumNbRadLength( G4double dn)//Geantino
{
  fnbRadLen += dn;
}
inline void H02RunData::Mass( G4double dn)//Geantino
{
  fVariable += dn;
}

inline void H02RunData::AddTotalEnergy(G4double ECh, G4double ENu, G4double Eesc,G4double Enotrc){//, G4double EHadCh, G4double EHadNu, G4double EEM) {

  fTotal_Ch_Energy +=  ECh ;
  fTotal_Nu_Energy  += ENu ;
  // fTotal_Esc_Energy  += Eesc ;
  // fTotal_Enotrc_Energy += Enotrc;
  // fTotal_Had_Ch_Energy +=  EHadCh ;
  // fTotal_Had_Nu_Energy  += EHadNu ;
  // fTotal_EM_Energy  += EEM ;
}

inline void H02RunData::AddTrueEnergy(G4double ECh, G4double ENu){//, G4double EHadCh, G4double EHadNu, G4double EEM) {

  fTrue_Ch_Energy +=  ECh ;
  fTrue_Nu_Energy += ENu ;
  // fTotal_Had_Ch_Energy +=  EHadCh ;
  // fTotal_Had_Nu_Energy  += EHadNu ;
  // fTotal_EM_Energy  += EEM ;
}
inline void H02RunData::AddPCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu, MParticle ptrc) 
{
    fPArray[iLayer][iCellEta][iCellPhi].chEnergy+=ECh;
    fPArray[iLayer][iCellEta][iCellPhi].nuEnergy+=ENu;
    int sizeptr = fPArray[iLayer][iCellEta][iCellPhi].prtcl.size();
    bool check = true;
    for (int part = 0; part < sizeptr; part++)
    {
      if (fPArray[iLayer][iCellEta][iCellPhi].prtcl[part].PosInList == ptrc.PosInList)
      {
        fPArray[iLayer][iCellEta][iCellPhi].prtcl[part].Energy+=ETot;
        check = false;
      }
    }
    if (fPArray[iLayer][iCellEta][iCellPhi].prtcl.empty())
    {
      fPArray[iLayer][iCellEta][iCellPhi].prtcl.push_back(ptrc);
    }
}


inline void H02RunData::AddCell(G4int iLayer, G4int iCellEta, G4int iCellPhi, G4double ETot, G4double ECh, G4double ENu){//, G4double EHadCh, G4double EHadNu, G4double EEM) {
  fCell_Energy   [iLayer][iCellEta][iCellPhi] += ETot ;
  fCh_Cell_Energy[iLayer][iCellEta][iCellPhi] += ECh ;
  fNu_Cell_Energy[iLayer][iCellEta][iCellPhi] += ENu ;
  // fPArray        [iLayer][iCellEta][iCellPhi].chEnergy+=ECh;
  // fPArray        [iLayer][iCellEta][iCellPhi].nuEnergy+=ENu;
  // for (int i = 0; i < kMaxPixel; i++)
  // {
  //   for (int j = 0; j < kMaxPixel; j++)
  //   {
  //     fCell_Energy   [4][i][j] = 1000 ;
  //     fCh_Cell_Energy[4][i][j] = 1000 ;
  //     fNu_Cell_Energy[4][i][j] = 1000 ;
  //   }
  // }
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

