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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02EventAction.cc
/// \brief Implementation of the H02EventAction class
//
//

#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "H02EventAction.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "HepMCG4Pythia8Interface.hh"
#include "H02MuonSD.hh"
#include "H02RunData.hh"
#include "H02RunAction.hh"
#include "Pythia8/Pythia.h"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02EventAction::H02EventAction()
 : G4UserEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02EventAction::~H02EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02EventAction::BeginOfEventAction(const G4Event* anEvent)
{
  const G4Event* ev = anEvent; ev=0;
#ifdef DEBUG_HEPMC
  // printout primary information
  G4cout << "Print out primary information" << G4endl;
  G4int nVtx= anEvent-> GetNumberOfPrimaryVertex();
  G4int i;
  for(i=0; i< nVtx; i++) {
    const G4PrimaryVertex* primaryVertex= anEvent-> GetPrimaryVertex(i);
    primaryVertex-> Print();
  }
#endif

auto runData
    = static_cast<H02RunData*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << " Print out hit information" << G4endl;
  // G4SDManager* SDManager= G4SDManager::GetSDMpointer();
  // H02MuonSD* muonSD=
  //   (H02MuonSD*)SDManager-> FindSensitiveDetector("/mydet/muon");
  // muonSD-> PrintAll();
  // muonSD-> DrawAll();
  // G4cout << G4endl;

  auto runGeneratorAction = static_cast<const H02PrimaryGeneratorAction*>
                 (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String GeneratorName =  runGeneratorAction->GetGeneratorName();

  // G4cout << "Current Generator Name : " << GeneratorName << G4endl ;

  // G4VPrimaryGenerator* primaryGen = runGeneratorAction->GetGenerator();
  // HepMCG4Pythia8Interface *py8Gen = (HepMCG4Pythia8Interface*)primaryGen ;

  // Pythia8::Event pythiaEvent = py8Gen->GetPythiaObject();
                                   

  // G4int NTruth = 0;
  // G4double ETruth = 0. ;

  // for(int ip = 0; ip < (int)pythiaEvent.size(); ip++) {
  //   if(pythiaEvent[ip].isFinal()) {
  //     cout <<"Particle : " << ip+1 <<", ID = "<< pythiaEvent[ip].id() 
  //     << " , Energy val = " << pythiaEvent[ip].e()
  //     << ", Mother1 = " << pythiaEvent[ip].mother1()
  //     << ", Mother2 = " << pythiaEvent[ip].mother2()
  //     << ", Production time = " << pythiaEvent[ip].tProd()
  //     << endl;

  //     NTruth ++ ;
  //     ETruth += pythiaEvent[ip].e() ;
  //    }
  //   }

  //   G4cout << "Truth Particle in event : " << NTruth << G4endl ;

   // ===========  Fill the ntuples ============ //
   auto runData
    = static_cast<H02RunData*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
   runData->FillPerEvent();

   auto runAction 
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();

   runAction->outTree->Fill();   

   // =================================================== //
   // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // extract the trajectories and print them out
  G4cout << G4endl;
  G4cout << "Trajectories in tracker "<< n_trajectories << G4endl;

  G4int nUniqueTraj(0), nSecTraj(0);
  for(G4int i=0; i<n_trajectories; i++)
    {
      G4Trajectory* trj =
        (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      //trj->ShowTrajectory();
      if(trj && trj->GetParentID() == 0) nUniqueTraj ++ ;
      // if(trj && trj->GetParentID() == 2)
      //   {
          nSecTraj ++;
          //G4cout << "Traj " << nUniqueTraj <<" : Particle = " << trj->GetParticleName() << G4endl ;
          //G4cout << "Traj " << nSecTraj <<" : TrackID = " << trj->GetTrackID() << G4endl ;
         // G4cout << i+1<<" th traj corresponds to step, parentID = " << trj->GetParentID() << ", 1st Truth Charge = " <<
         // trj0->GetCharge() /*<< " , 2nd Truth Chg = " << trj1->GetCharge() */<< G4endl ;
         //break;
        // }
    }
    
    int TotalSecondary(0);
    double TotalInitialE(0);

    for(int iTraj = 0; iTraj < (int)runData->fAllTrajectoryInfo.size(); iTraj++){
      
      for(int iParent = 0; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vTrackID.size(); iParent++ ) {
          TotalSecondary ++ ;
          
          TotalInitialE += runData->fAllTrajectoryInfo[iTraj].vTrackEnergy.at(iParent);
          // G4cout << "Step " << iParent +1 
          // << ", Time = " << runData->fAllTrajectoryInfo[iTraj].vTrackTime[iParent]
          // << G4endl ;
      }


      G4cout << " Traj : " << iTraj +1 << ", Energy = " 
      << runData->fAllTrajectoryInfo[iTraj].fEnergy 
      << ", Px = "
      << runData->fAllTrajectoryInfo[iTraj].fMomentum.x()
      << ", Py = "
      << runData->fAllTrajectoryInfo[iTraj].fMomentum.y()
      << ", Pz = "
      << runData->fAllTrajectoryInfo[iTraj].fMomentum.z()
      << ", TotalStepE = "
      << TotalInitialE
      << G4endl;
      //TotalInitialE += runData->fAllTrajectoryInfo[iTraj].fEnergy;


    	// G4cout << "Traj : " << iTraj +1 << ", Particle = " 
    	// << runData->fAllTrajectoryInfo[iTraj].fPDGCode 
    	// << ", Momentum : " 
    	// << runData->fAllTrajectoryInfo[iTraj].fEnergy / (1000 * CLHEP::GeV )
    	// << ", Vertex X = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].x() 
    	// << ", Vertex Y = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].y()
    	// << ", Vertex Z = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].z()
    	// << ", Charge = " << runData->fAllTrajectoryInfo[iTraj].fPDGCharge
     //  << ", Secondary showers = " << runData->fAllTrajectoryInfo[iTraj].vTrackID.size()
    	// << G4endl;
    }

    G4cout << "Unique Traj " << nUniqueTraj <<", Secondary Traj =  " << nSecTraj 
    << ", N Unique Traj RunData = " << runData->fAllTrajectoryInfo.size()
    << ", N Secondary Traj RunData = " <<  TotalSecondary
    // << ", Truth-E = "<< ETruth
    // << ", Traj-E = " << TotalInitialE * 1e-6

    << G4endl ;

    G4cout << " Total Neutral energy " << runAction->Total_Nu_Energy << G4endl;
    G4cout << " Total Charged energy " << runAction->Total_Ch_Energy << G4endl;
  // ====================================================== //
  // G4int n_vertex = evt->GetNumberOfPrimaryVertex();
  // for(G4int iv=0;iv<n_vertex;iv++)
  // {
  //   G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
  //   G4cout << G4endl;
  //   G4cout << "Primary vertex "
  //          << G4ThreeVector(pv->GetX0(),pv->GetY0(),pv->GetZ0())
  //          << "   at t = " << (pv->GetT0())/CLHEP::ns << " [ns]" 
  //          <<", PDG ID = " << pv->GetPrimary()->GetPDGcode()
  //          <<", Energy = " << pv->GetPrimary()->GetTotalEnergy()
  //          << G4endl;
  //   // if(fpEventManager->GetVerboseLevel()>0)
  //   // {
  //   //   G4PrimaryParticle* pp = pv->GetPrimary();
  //   //   while(pp)
  //   //   {
  //   //     //PrintPrimary(pp,0);
  //   //     G4cout << " G4Track ID " << pp->GetTrackID() << G4endl;
  //   //     pp = pp->GetNext();
  //   //   }
  //   // }
  // }



}
