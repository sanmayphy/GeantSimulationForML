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
/// \file runAndEvent/H02/src/H02TrackingAction.cc
/// \brief Implementation of the H02TrackingAction class
//
//

#include "H02TrackingAction.hh"
//#include "H02Trajectory.hh"
#include "H02TrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4RunManager.hh"
#include "H02RunData.hh"
#include "H02SteppingAction.hh"
#include "H02FullTrajectoryInfo.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
H02TrackingAction::H02TrackingAction()
:G4UserTrackingAction()
{;}

void H02TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // // Create trajectory only for track in tracking region
  // H02TrackInformation* trackInfo = 
  //   (H02TrackInformation*)(aTrack->GetUserInformation());

  // if(trackInfo->GetTrackingStatus() > 0)
  // {
  //   fpTrackingManager->SetStoreTrajectory(true);
  //   fpTrackingManager->SetTrajectory(new H02Trajectory(aTrack));
  // }
  // else
  // { fpTrackingManager->SetStoreTrajectory(false); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void H02TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  // if(secondaries)
  // {
  //   H02TrackInformation* info = 
  //     (H02TrackInformation*)(aTrack->GetUserInformation());
  //   size_t nSeco = secondaries->size();
  //   if(nSeco>0)
  //   {
  //     for(size_t i=0;i<nSeco;i++)
  //     { 
  //       H02TrackInformation* infoNew = new H02TrackInformation(info);
  //       (*secondaries)[i]->SetUserInformation(infoNew);
  //     }
  //   }
  // }

  auto runData = static_cast<H02RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());


  auto event
    = static_cast<const G4Event*>(
        G4RunManager::GetRunManager()->GetCurrentEvent());

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();


      if(aTrack->GetParentID() == 0) {

      FullTrajectoryInfo trjInfo;
      trjInfo.fParentID  = aTrack->GetParentID()        ;
      trjInfo.fTrackID  = aTrack->GetTrackID()         ;
      trjInfo.fPDGCharge = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetCharge() ;
      trjInfo.fPDGCode   = aTrack->GetDefinition()->GetPDGEncoding();
      trjInfo.fMomentum  = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentum() ;
      trjInfo.fMomentumDir  = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentumDirection() ;
      trjInfo.fEnergy      = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
      //trjInfo.fEnergy      = aTrack->GetDynamicParticle()->GetKineticEnergy();
      trjInfo.fVertexPosition = aTrack->GetVertexPosition();
      trjInfo.fGlobalTime = aTrack->GetGlobalTime();

      trjInfo.vTrackMomentumDir.push_back( aTrack->GetMomentum()  );
      trjInfo.vTrackID.push_back(  aTrack->GetTrackID() ) ;
      trjInfo.vParentID.push_back( aTrack->GetParentID() ) ;
      trjInfo.vPreStepTrackPos.push_back( aTrack->GetStep()->GetPreStepPoint()->GetPosition() ) ;
      trjInfo.vPostStepTrackPos.push_back( aTrack->GetStep()->GetPostStepPoint()->GetPosition() ) ;
      trjInfo.vTrackPos.push_back( aTrack->GetPosition() ) ;
      trjInfo.vTrackTime.push_back( aTrack->GetGlobalTime() );
      trjInfo.vTrackEnergy.push_back( aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy() ) ;
      //trjInfo.vTrackEnergy.push_back( aTrack->GetMomentum().mag()  );
      trjInfo.vTrackPDGID.push_back( aTrack->GetDefinition()->GetPDGEncoding() ) ;

      runData->fAllTrajectoryInfo.push_back(trjInfo) ;

      //break;
    } // if(trj->GetParentID() == 0)

    else {
      
      G4int ParentID = aTrack->GetParentID() ;
      G4ThreeVector PreStepPoint = aTrack->GetStep()->GetPreStepPoint()->GetPosition() ;
      G4ThreeVector PostStepPoint = aTrack->GetStep()->GetPostStepPoint()->GetPosition() ;

      bool foundTraj(false);
      int mTraj(-1), mParent(-1);
      //G4cout << "Scanning over existing traj = " << runData->fAllTrajectoryInfo.size() << G4endl ;
      // ------ search over all the trajectories ---------- //
      for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++) {
      //G4cout << "Scanning over existing step = " << runData->fAllTrajectoryInfo[iTraj].vParentID.size() << G4endl ;
      // ----- search for all the existing trackID in a given trajectory ---- //
       for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  ) {

         if(   (runData->fAllTrajectoryInfo[iTraj].vTrackID[iParent] == ParentID)  
           // && (runData->fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iParent].isNear(PreStepPoint) ) 
            ) {
           
           foundTraj = true; 
           mTraj = iTraj ;
           mParent = iParent ;
           // G4cout << "Found Secondary Track with ParentID = " 
           // << runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] 
           // << G4endl;
           break ;

         } // if( runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] == ParentID )

       } // for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  )

       if(foundTraj) break;

      } // for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++)

      if(foundTraj) {

      runData->fAllTrajectoryInfo[mTraj].vTrackMomentumDir.push_back( aTrack->GetMomentum()  );
      runData->fAllTrajectoryInfo[mTraj].vParentID.push_back( aTrack->GetParentID() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackID.push_back(  aTrack->GetTrackID() ) ;
      runData->fAllTrajectoryInfo[mTraj].vPreStepTrackPos.push_back( PreStepPoint ) ;
      runData->fAllTrajectoryInfo[mTraj].vPostStepTrackPos.push_back( PostStepPoint ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackPos.push_back( aTrack->GetPosition() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackTime.push_back( aTrack->GetGlobalTime() );
      runData->fAllTrajectoryInfo[mTraj].vTrackEnergy.push_back( aTrack->GetStep()->GetTotalEnergyDeposit() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackPDGID.push_back( aTrack->GetDefinition()->GetPDGEncoding() ) ;
     } //  if(foundTraj)

     // G4cout << "Secondary Track with ParentID = " <<  aTrack->GetParentID() 
     // << " Trajectory index = " << mTraj
     // << G4endl;

    } //if(aTrack->GetParentID() != 0)
    

}


