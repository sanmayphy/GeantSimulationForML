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
#include "H02RunAction.hh"
#include "H02RunAction.hh"

#include "H02SteppingAction.hh"
#include "H02FullTrajectoryInfo.hh"
#include "H02PrimaryGeneratorAction.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
H02TrackingAction::H02TrackingAction()
:G4UserTrackingAction()
{;}

void H02TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // G4cout<<"preTR "<<aTrack->GetStep()->GetTotalEnergyDeposit()<<G4endl;
  // // Create trajectory only for track in tracking region
  // H02TrackInformation* trackInfo = 
    // (H02TrackInformation*)(aTrack->GetUserInformation());
  // G4cout << "GetTrackingStatus = " << trackInfo->GetTrackingStatus() << G4endl ;
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
  auto runGeneratorAction = static_cast<const H02PrimaryGeneratorAction*>
                 (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String GeneratorName =  runGeneratorAction->GetGeneratorName();

  auto runData = static_cast<H02RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());


  auto event
    = static_cast<const G4Event*>(
        G4RunManager::GetRunManager()->GetCurrentEvent());

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  auto runAction
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  G4int ParentID = aTrack->GetParentID() ;
  float Ech = 0;
  float Enu = 0;
  // G4cout<<"C0"<<G4endl;

  if(aTrack->GetParentID() == 0) 
  {

    FullTrajectoryInfo trjInfo;
    trjInfo.fParentID  = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTrackID();//GetParentID()        ;
    trjInfo.fTrackID  = aTrack->GetTrackID()         ;
    trjInfo.fPDGCharge = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetCharge() ;
    trjInfo.fPDGCode   = aTrack->GetDefinition()->GetPDGEncoding();
    trjInfo.fMomentum  = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentum() ;
    trjInfo.fMomentumDir  = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentumDirection() ;
    trjInfo.fEnergy      = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
    //trjInfo.fEnergy      = aTrack->GetDynamicParticle()->GetKineticEnergy();

    trjInfo.fVertexPosition = aTrack->GetVertexPosition();
    trjInfo.fGlobalTime = aTrack->GetGlobalTime();
    // G4cout << "ID = " 
    // << aTrack->GetDefinition()->GetPDGEncoding()
    // << " rjInfo.fEnergy = " 
    // << aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy() 
    // << " m = "
    // << aTrack->GetDynamicParticle()->GetPrimaryParticle()-> GetMass()
    // << G4endl ;

    trjInfo.vTrackMomentumDir.push_back( aTrack->GetMomentum()  );
    trjInfo.vTrackID.push_back(  aTrack->GetTrackID() ) ;
    trjInfo.vParentID.push_back( aTrack->GetParentID() ) ;
    trjInfo.vPreStepTrackPos.push_back( aTrack->GetStep()->GetPreStepPoint()->GetPosition() ) ;
    trjInfo.vPostStepTrackPos.push_back( aTrack->GetStep()->GetPostStepPoint()->GetPosition() ) ;
    trjInfo.vTrackPos.push_back( aTrack->GetPosition() ) ;
    trjInfo.vTrackTime.push_back( aTrack->GetGlobalTime() );
    trjInfo.vTrackEnergy.push_back( aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy() ) ;//GetTotalEnergy
    //trjInfo.vTrackEnergy.push_back( aTrack->GetMomentum().mag()  );
    trjInfo.vTrackPDGID.push_back( aTrack->GetDefinition()->GetPDGEncoding() ) ;

    runData->fAllTrajectoryInfo.push_back(trjInfo) ;
    //!Pythia8
    if (GeneratorName=="pythia8")
    {
      if(aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetCharge() == 0)
      {
        Enu = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
      }
      else
      {
        Ech = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
      }
      runData->AddTrueEnergy(Ech, Enu);
    }
    //!Pythia8 END
  //break;
  } // if(trj->GetParentID() == 0)
  else 
  {
    //!Particel gun
    // G4int ParentID = aTrack->GetParentID() ;
    //!Particel gun end
    
    G4ThreeVector PreStepPoint = aTrack->GetStep()->GetPreStepPoint()->GetPosition();
    G4ThreeVector PostStepPoint = aTrack->GetStep()->GetPostStepPoint()->GetPosition();

    bool foundTraj(false);
    int mTraj(-1), mParent(-1);
    // G4cout << "Scanning over existing traj = " << runData->fAllTrajectoryInfo.size() << G4endl ;
    // ------ search over all the trajectories ---------- //
    // G4cout<<""<<G4endl;
    // G4cout << "Nymber of Trj = " << (int)runData->fAllTrajectoryInfo.size() << G4endl ;
    for(int iTraj = (int)runData->fAllTrajectoryInfo.size()-1; iTraj >= 0; iTraj--) 
    {
      // G4cout << "Scanning over existing step = " << (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size() << G4endl ;
      // G4cout << "iTraj = " << iTraj << G4endl ;
    // ----- search for all the existing trackID in a given trajectory ---- //
      for(int iParent = (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size()-1; iParent >=0 ; iParent--  ) 
      {

        if(   (runData->fAllTrajectoryInfo[iTraj].vTrackID[iParent] == ParentID)  
          // && (runData->fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iParent].isNear(PreStepPoint) ) 
          ) {
          
          foundTraj = true; 
          mTraj = iTraj ;
          mParent = iParent ;
          // G4cout << "Found Secondary Track with ParentID = " 
          // << runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] 
          // << " mTraj = "<< mTraj
          // << G4endl;
          break ;

        } // if( runData->fAllTrajectoryInfo[iTraj].vParentID[iParent] == ParentID )

      } // for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vParentID.size(); iParent++  )

      if(foundTraj) break;

    } // for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++)
    // G4cout<<"C2"<<G4endl;
    if(foundTraj) 
    {
      // G4cout<<"Check1"<<G4endl;
      runData->fAllTrajectoryInfo[mTraj].vTrackMomentumDir.push_back( aTrack->GetMomentum()  );
      runData->fAllTrajectoryInfo[mTraj].vParentID.push_back( aTrack->GetParentID() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackID.push_back(  aTrack->GetTrackID() ) ;
      runData->fAllTrajectoryInfo[mTraj].vPreStepTrackPos.push_back( PreStepPoint ) ;
      runData->fAllTrajectoryInfo[mTraj].vPostStepTrackPos.push_back( PostStepPoint ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackPos.push_back( aTrack->GetPosition() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackTime.push_back( aTrack->GetGlobalTime() );
      runData->fAllTrajectoryInfo[mTraj].vTrackEnergy.push_back( aTrack->GetStep()->GetTotalEnergyDeposit() ) ;
      runData->fAllTrajectoryInfo[mTraj].vTrackPDGID.push_back( aTrack->GetDefinition()->GetPDGEncoding() ) ;
        // G4cout<<"Check2"<<G4endl;
      // G4cout << "Traking TrackID "<<aTrack->GetTrackID()<<" Energ Dep "<< aTrack->GetStep()->GetTotalEnergyDeposit() <<G4endl;
      // G4cout << G4endl;
    } //  if(foundTraj)

    // G4cout << "Secondary Track with ParentID = " <<  aTrack->GetParentID() 
    // << " Trajectory index = " << mTraj
    // << G4endl;

  } //if(aTrack->GetParentID() != 0)
  // cout<<"bb2"<<endl;

}


