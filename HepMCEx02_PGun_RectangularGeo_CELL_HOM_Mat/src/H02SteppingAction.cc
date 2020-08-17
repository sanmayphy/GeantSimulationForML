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

#include "CaloRConstants.hh"

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

  int Z_Bin(-1), X_Bin(-1), Y_Bin(-1);

  

  if(ZPos >= -1137.33 && ZPos < -1020.42) Z_Bin = 0 ;
  if(ZPos >= -1020.42 && ZPos < -396.901) Z_Bin = 1 ;
  if(ZPos >= -396.901 && ZPos < -163.081) Z_Bin = 2 ;
  if(ZPos >= -153.081 && ZPos < 108.489 ) Z_Bin = 3 ;
  if(ZPos >= 108.489  && ZPos < 823.447 ) Z_Bin = 4 ;
  if(ZPos >= 823.447  && ZPos < 1137.33 ) Z_Bin = 5 ;

  int N_ETA = CELL_NUMBER[Z_Bin]; 
  int N_PHI = CELL_NUMBER[Z_Bin];

  double EtaBinWidth = (ETA_MAX - ETA_MIN)/N_ETA ;
  for(G4int iEtacell = 0; iEtacell < N_ETA; iEtacell++) {
    if( ( XPos >=  ETA_MIN + EtaBinWidth * iEtacell  ) &&
        ( XPos <  ETA_MIN + EtaBinWidth * (iEtacell + 1)  )
      )
      X_Bin = iEtacell ;
  }


  double PhiBinWidth = (PHI_MAX - PHI_MIN)/N_PHI ;
  for(G4int iPhicell = 0; iPhicell < N_PHI; iPhicell++) {
    if( ( YPos >=  PHI_MIN + PhiBinWidth * iPhicell  ) &&
        ( YPos <  PHI_MIN + PhiBinWidth * (iPhicell + 1)  )
      )
      Y_Bin = iPhicell ;
  }


  //G4cout<<" Final ZXY bin = "<< zbin <<" : " << xbin  << " : " << " : " << ybin << G4endl ;

  static int ZXYBin[3] = {-1};
  ZXYBin[0] = Z_Bin; ZXYBin[1] = X_Bin; ZXYBin[2] = Y_Bin;
  

  return ZXYBin;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02SteppingAction::UserSteppingAction(const G4Step* astep)
{
  G4Track* aTrack= astep-> GetTrack();

  if(aTrack-> GetTrackStatus() != fAlive) return;
  //if(atrack-> GetParentID() == 0) return;

  auto runData = static_cast<H02RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  auto event
    = static_cast<const G4Event*>(
        G4RunManager::GetRunManager()->GetCurrentEvent());


  G4ThreeVector PreStepPoint = astep->GetPreStepPoint()->GetPosition() ;
  G4ThreeVector PostStepPoint = astep->GetPostStepPoint()->GetPosition() ;

  G4double Eta = PreStepPoint.eta();

  G4int ParentID = aTrack->GetParentID() ;

  // ----- assigning the trajectory to step ---- // 
  bool foundTraj(false);
      int mTraj(-1), mParent(-1);
      //G4cout << "Scanning over existing traj = " << runData->fAllTrajectoryInfo.size() << G4endl ;
      // ------ search over all the trajectories ---------- //
      for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++) {
      // ----- search for all the existing trackID in a given trajectory ---- //
       for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  ) {

         if(   (runData->fAllTrajectoryInfo[iTraj].vTrackID[iParent] == ParentID)  
            // && (runData->fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iParent] == PreStepPoint)
            ) {
           
           foundTraj = true; 
           mTraj = iTraj ;
           mParent = iParent ;
           break ;

         } // if( runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] == ParentID )

       } // for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  )

       if(foundTraj) break;

    } // for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++)

   G4double Charge = runData->fAllTrajectoryInfo[mTraj].fPDGCharge ;
   G4int    PDGId  = runData->fAllTrajectoryInfo[mTraj].fPDGCode ;
   auto edep = astep->GetTotalEnergyDeposit();

   
   //auto edep = astep->GetDeltaEnergy();

   //edep -= astep->GetNonIonizingEnergyDeposit();

 //  if(edep<=0) return;

   if(foundTraj &&  edep > 0.) {
     // G4cout << " Step Energy : " << edep 
     //     << " Origin charge : " << Charge  
     //     << "PDG ID : " << runData->fAllTrajectoryInfo[mTraj].fPDGCode
     //      << " Eta : " << PreStepPoint.eta()  
     //      << " Phi : " << PreStepPoint.phi()
     //      << ", Trajectory Number = " << mTraj
     //      << ", Step Number = " << mParent
     //      << ", StepX = " << PreStepPoint.x()
     //      << ", StepY = " << PreStepPoint.y()
     //      << ", StepZ = " << PreStepPoint.z()
     //      << ", Energy = " << edep
     //      << ", Time = "  << runData->fAllTrajectoryInfo[mTraj].vTrackTime[mParent]
     //      << G4endl ;
   
  

   G4double Etot(0), Ech(0.), Enu(0.) ;
   if(Charge == 0) Enu = edep ;
   else            Ech = edep ;

   Etot = Ech + Enu ;

   G4int ChTraj(mTraj);

   // if( fabs(Charge) > 0 )   ChTraj =  mTraj ;
   // else { ChTraj =  fabs(mTraj - 1) ; }
   
   G4double chParticleEta = runData->fAllTrajectoryInfo[ChTraj].fMomentumDir.eta();
   G4double chParticlePhi = runData->fAllTrajectoryInfo[ChTraj].fMomentumDir.phi();

   int *Bin = CellIndex(PreStepPoint.x(), PreStepPoint.y() , PreStepPoint.z() );

   int *Bin_End = CellIndex(PostStepPoint.x(), PostStepPoint.y() , PostStepPoint.z() );
   //Etot = Ech + Enu ;

   if(*Bin >= 0 && *(Bin+1) >= 0 && *(Bin+2) >= 0 ) {
   // G4cout << " Z Bin : "<< *Bin 
   //        << " X Bin : "<< *(Bin+1)
   //        << " Y Bin : "<< *(Bin+2)
   //        << " Z Bin_F : "<< *Bin_End 
   //        << " X Bin_F : "<< *(Bin_End+1)
   //        << " Y Bin_F : "<< *(Bin_End+2)

   // << G4endl ;
   
  runData->AddCell(*Bin, *(Bin+1), *(Bin+2), Etot, Ech, Enu) ;
  runData->AddTotalEnergy(Ech, Enu);
  } // if(*Bin >= 0)

  } // if(foundTraj && edep > 0.)

  G4ParticleDefinition* particleType= aTrack-> GetDefinition();
  if((particleType == G4MuonPlus::MuonPlusDefinition())
     || (particleType == G4MuonMinus::MuonMinusDefinition())) return;

  G4StepPoint* prestep= astep-> GetPreStepPoint();
  G4VPhysicalVolume* pv= prestep-> GetPhysicalVolume();
  G4String pvname= pv-> GetName();
  if(pvname=="BARREL_CAL_PV" || pvname=="ENDCAP_CAL_PV" ) {
    aTrack-> SetTrackStatus(fKillTrackAndSecondaries);
  }
}
