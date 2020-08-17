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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02SteppingAction.cc
/// \brief Implementation of the H02SteppingAction class
//
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "H02RunData.hh"
#include "H02SteppingAction.hh"
#include "H02FullTrajectoryInfo.hh"
#include <math.h>
#include "H02RunAction.hh"
// #include "CaloRConstants.hh"

//#include "GlobalVariables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02SteppingAction::H02SteppingAction()
 : G4UserSteppingAction()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02SteppingAction::~H02SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int * H02SteppingAction::CellIndex(double XPos, double YPos, double ZPos) {

  int R_Bin(-1), Eta_Bin(-1), Phi_Bin(-1);



   
  double r = pow(XPos, 2) + pow(YPos, 2);

  if(r >= pow(r_inn_ECAL1, 2) && r < pow(r_out_ECAL1, 2)) R_Bin = 0 ;
  if(r >= pow(r_inn_ECAL2, 2) && r < pow(r_out_ECAL2, 2)) R_Bin = 1 ;
  if(r >= pow(r_inn_ECAL3, 2) && r < pow(r_out_ECAL3, 2)) R_Bin = 2 ;
  if(r >= pow(r_inn_HCAL1, 2) && r < pow(r_out_HCAL1, 2)) R_Bin = 3 ;
  if(r >= pow(r_inn_HCAL2, 2) && r < pow(r_out_HCAL2, 2)) R_Bin = 4 ;
  if(r >= pow(r_inn_HCAL3, 2) && r < pow(r_out_HCAL3, 2)) R_Bin = 5 ;

  double EtaPos = -1*log(tan(0.5*acos(ZPos/pow(r+pow(ZPos, 2), 0.5))));
  Eta_Bin = (int) floor((eta_max+EtaPos)/d_eta);

  double PhiPos = atan2(YPos, XPos);
  if (PhiPos < 0) PhiPos += tube_dPhi;
  Phi_Bin = (int) floor(PhiPos/divided_tube_dPhi);

  static int ZXYBin[3] = {-1};
  ZXYBin[0] = R_Bin; ZXYBin[1] = Eta_Bin; ZXYBin[2] = Phi_Bin;


  return ZXYBin;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02SteppingAction::UserSteppingAction(const G4Step* astep)
{
  // G4cout<<"B0"<<G4endl;
  G4Track* aTrack= astep-> GetTrack();
  auto kine = aTrack->GetDynamicParticle()->GetKineticEnergy();
  auto edep = astep->GetTotalEnergyDeposit();

  auto runData = static_cast<H02RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  auto event
    = static_cast<const G4Event*>(
        G4RunManager::GetRunManager()->GetCurrentEvent());

  G4ThreeVector PreStepPoint = astep->GetPreStepPoint()->GetPosition() ;
  G4ThreeVector PostStepPoint = astep->GetPostStepPoint()->GetPosition() ;
  G4TouchableHandle touch1 = astep->GetPreStepPoint()->GetTouchableHandle();
  G4LogicalVolume* lvol = touch1->GetVolume()->GetLogicalVolume();
  G4int ParentID = aTrack->GetParentID() ;


  if(aTrack-> GetTrackStatus() != fAlive && edep==0.0 ) return;

  if((int)runData->fAllTrajectoryInfo.size()==0) return;
  // G4cout<<"B1"<<G4endl;
  //!Geantino
  
  G4ParticleDefinition* particle0 = aTrack->GetDefinition();
  G4double mass = 0;
  if (particle0 == G4ChargedGeantino::ChargedGeantino()&&lvol->GetName()!="EXP_HALL_LV")
  {
    // mass+=lvol->GetMass();
    G4double radl  = lvol->GetMaterial()->GetRadlen()/cm;
    G4double density  = lvol->GetMaterial()->GetDensity()/(g/cm3);
    G4double stepl = astep->GetStepLength()/cm;
    // runData->SumNbRadLength((lvol->GetMass()/g)*density*(stepl)/radl);
    runData->SumNbRadLength((stepl)/radl);
    runData->Mass(lvol->GetMass()/g);
    G4cout<<"PreStepPoint.x() "<<PreStepPoint.x()<<" PreStepPoint.y() "<<PreStepPoint.y()<<" PreStepPoint.z() "<<PreStepPoint.z()<<G4endl;
  }
  
 
  //!Geantino end
  // G4cout<<"B2"<<G4endl;
  // ----- assigning the trajectory to step ---- // 
  bool foundTraj(false);
  int mTraj(-1), mParent(-1);
      // ------ search over all the trajectories ---------- //
  for(int iTraj = (int)runData->fAllTrajectoryInfo.size()-1; iTraj >= 0; iTraj--) 
  {
      // ----- search for all the existing trackID in a given trajectory ---- //
    for(int iParent = (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size()-1; iParent >= 0; iParent--  ) 
    {

      if(   (runData->fAllTrajectoryInfo[iTraj].vTrackID[iParent] == ParentID)  
            // && (runData->fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iParent] == PreStepPoint)
            )
      {
           foundTraj = true; 
           mTraj = iTraj ;
           mParent = iParent ;
           break ;

      } // if( runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] == ParentID )

    } // for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  )

    if(foundTraj) break;

  } // for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++)
  // G4cout<<"B3"<<G4endl;
  if(foundTraj &&  edep > 0.) 
  {
    int *Bin = CellIndex(PreStepPoint.x(), PreStepPoint.y() , PreStepPoint.z() );
    int *Bin_End = CellIndex(PostStepPoint.x(), PostStepPoint.y() , PostStepPoint.z() );
    G4double Charge = runData->fAllTrajectoryInfo[mTraj].fPDGCharge ;

    G4double Etot(0), Ech(0.), Enu(0.);//, EHadCh(0.), EHadNu(0.), EEM(0.);
    
    if(Charge == 0) 
    {
      Enu = edep ;
      // G4cout<<"NurunData->fAllTrajectoryInfo "<< runData->fAllTrajectoryInfo[mTraj].fParentID<<G4endl;
      // G4cout<<"NurunData->fAllTrajectoryInfo "<< runData->fAllTrajectoryInfo[mTraj].fEnergy<<G4endl;
    }
    else            
    {
      Ech = edep ;
      // G4cout<<"ChrunData->fAllTrajectoryInfo "<< runData->fAllTrajectoryInfo[mTraj].fParentID<<G4endl;
    }

    Etot = Ech + Enu ;

    G4int ChTraj(mTraj);

    auto runAction
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();

    if((*Bin >= 0 && *Bin < 256 ) && (*(Bin+1) >= 0 && *(Bin+1) < 256)&& (*(Bin+2) >= 0 && *(Bin+2) < 256) ) 
    {
      runData->AddCell(*Bin, *(Bin+1), *(Bin+2), Etot, Ech, Enu) ;
      runData->AddTotalEnergy(Ech, Enu, 0, 0);
      // MParticle ptrc;
      // ptrc.Energy = Etot;
      // ptrc.PosInList = runData->fAllTrajectoryInfo[mTraj].fParentID;
      // runData->AddPCell(*Bin, *(Bin+1), *(Bin+2), Etot, Ech, Enu, ptrc);
      
    } // if(*Bin >= 0 && *(Bin+1) >= 0 && *(Bin+2) >= 0 )
  } // if(foundTraj && edep > 0.)

  G4ParticleDefinition* particleType= aTrack-> GetDefinition();
  if((particleType == G4MuonPlus::MuonPlusDefinition())
     || (particleType == G4MuonMinus::MuonMinusDefinition())) return;

}
