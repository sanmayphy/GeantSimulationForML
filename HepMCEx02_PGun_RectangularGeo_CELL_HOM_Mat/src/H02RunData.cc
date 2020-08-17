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
/// \file H02RunData.cc
/// \brief Implementation of the H02RunData class

#include "H02RunData.hh"
#include "H02RunAction.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "CaloRConstants.hh"

#include "G4TrajectoryContainer.hh"
#include "H02PrimaryGeneratorAction.hh"
#include "HepMCG4Pythia8Interface.hh"
#include "Pythia8/Pythia.h"
//#include "RE01Trajectory.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

H02RunData::H02RunData() 
 : G4Run(),
   fVolumeNames{ { "Absorber", "Gap" } }
{
  for ( auto& edep : fEdep ) { 
    edep = 0.; 
  }
  for ( auto& trackLength : fTrackLength ) {
    trackLength = 0.; 
  }

  for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
    for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
      for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
        fCell_Energy[iLayer][iEtacell][iPhicell] = 0. ;
        fCh_Cell_Energy[iLayer][iEtacell][iPhicell] = 0. ;
        fNu_Cell_Energy[iLayer][iEtacell][iPhicell] = 0. ;

      }
    }
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

H02RunData::~H02RunData()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02RunData::FillPerEvent()
{
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // accumulate statistic
  // in the order od the histograms, ntuple columns declarations
  // G4int counter = 0;
  // for ( auto edep : fEdep ) {
  //   analysisManager->FillH1(counter, edep);
  //   analysisManager->FillNtupleDColumn(counter++, edep);
  // }
  // for ( auto trackLength : fTrackLength ) {
  //   analysisManager->FillH1(counter, trackLength);
  //   analysisManager->FillNtupleDColumn(counter++, trackLength);
  // }  

  // for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
  //   for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
  //     for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
  //       analysisManager->FillNtupleDColumn(counter++, fCell_Energy[iLayer][iEtacell][iPhicell]);
  //       analysisManager->FillNtupleDColumn(counter++, fCh_Cell_Energy[iLayer][iEtacell][iPhicell]);
  //       analysisManager->FillNtupleDColumn(counter++, fNu_Cell_Energy[iLayer][iEtacell][iPhicell]);
  //     }
  //   }
  // }


  auto runGeneratorAction = static_cast<const H02PrimaryGeneratorAction*>
                 (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String GeneratorName =  runGeneratorAction->GetGeneratorName();

  G4cout << "Current Generator Name : " << GeneratorName << G4endl ;

  G4VPrimaryGenerator* primaryGen = runGeneratorAction->GetGenerator();

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

  //     fTruthParticleE.push_back( pythiaEvent[ip].e() ) ;
  //     fTruthParticlePx.push_back( pythiaEvent[ip].px() ) ;
  //     fTruthParticlePy.push_back( pythiaEvent[ip].py() ) ;
  //     fTruthParticlePz.push_back( pythiaEvent[ip].pz() ) ;
  //     fTruthParticleM.push_back( pythiaEvent[ip].m() ) ;
  //     fTruthParticleMother1.push_back( pythiaEvent[ip].mother1() ) ;
  //     fTruthParticleMother2.push_back( pythiaEvent[ip].mother2() ) ;
      
  //    }
  //   }

  //   G4cout << "Truth Particle in event : " << NTruth << G4endl ;

  //G4cout<<"Cell Energy : " << fCell_Energy[0][10][10] << G4endl;

  auto runAction 
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();

   runAction->cellE = fCell_Energy[0][10][10];

   runAction->Total_Ch_Energy = fTotal_Ch_Energy;
   runAction->Total_Nu_Energy = fTotal_Nu_Energy;


   for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
    for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
      for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
        
        runAction->Cell_Energy[iLayer][iEtacell][iPhicell] = fCell_Energy[iLayer][iEtacell][iPhicell] ;
        runAction->Ch_Cell_Energy[iLayer][iEtacell][iPhicell] = fCh_Cell_Energy[iLayer][iEtacell][iPhicell] ;
        runAction->Nu_Cell_Energy[iLayer][iEtacell][iPhicell] = fNu_Cell_Energy[iLayer][iEtacell][iPhicell] ;

        
      }
    }
  }
  // analysisManager->FillNtupleDColumn(counter++, fAllTrajectoryInfo.size());

  G4int NTracj = fAllTrajectoryInfo.size();

  for(G4int iTraj = 0; iTraj < NTracj; iTraj++){
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].vTrackID.size()) : (-99) );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fPDGCode) : (-99) );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fPDGCharge) : (-99)  );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fEnergy) : (-99) );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fMomentum.x()) : (-99) );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fMomentum.y()) : (-99) );
    // analysisManager->FillNtupleDColumn(counter++, (iTraj < NTracj) ? (fAllTrajectoryInfo[iTraj].fMomentum.z()) : (-99) );

    G4cout << "Source particleID = " << fAllTrajectoryInfo[iTraj].fPDGCode << G4endl;

    G4int NTrackSteps = fAllTrajectoryInfo[iTraj].vTrackID.size() ;
    fTrackStepNumb.push_back( double(NTrackSteps) );
    fTrackPDGCode.push_back( fAllTrajectoryInfo[iTraj].fPDGCode );
    fTrackPDGCharge.push_back( fAllTrajectoryInfo[iTraj].fPDGCharge );
    fTrackEnergy.push_back( fAllTrajectoryInfo[iTraj].fEnergy );
    fTrackPx.push_back( fAllTrajectoryInfo[iTraj].fMomentum.x() );
    fTrackPy.push_back( fAllTrajectoryInfo[iTraj].fMomentum.y() );
    fTrackPz.push_back( fAllTrajectoryInfo[iTraj].fMomentum.z() );

    std::vector<double> StepID, StepPDGID, StepEnergy, StepTime ;
    std::vector<double> StepInitPosX, StepInitPosY, StepInitPosZ;
    std::vector<double> StepFinalPosX, StepFinalPosY, StepFinalPosZ ;

    for(G4int iStep = 0; iStep < NTrackSteps; iStep++) {

      // G4cout << "Source particleID = " << fAllTrajectoryInfo[iTraj].fPDGCode 
      // << ", Shower : " << iStep+1 << ", PDG ID = " << fAllTrajectoryInfo[iTraj].vTrackPDGID[iStep]
      // << ", Radial Pos : " << sqrt( pow(fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].x(),2) + pow(fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].y(),2) )
      // << G4endl;

      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackID[iStep]) : (-99) ); 
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackPDGID[iStep]) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackEnergy[iStep]) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackTime[iStep]) : (-99) );

      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackPos[iStep].x()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackPos[iStep].y()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vTrackPos[iStep].z()) : (-99) );

      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].x()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].y()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].z()) : (-99) );

      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].x()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].y()) : (-99) );
      // analysisManager->FillNtupleDColumn(counter++, (iStep < std::min(NTrackSteps, N_STORE_STEP) ) ? (fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].z()) : (-99) );
  

      StepID.push_back(fAllTrajectoryInfo[iTraj].vTrackID[iStep]);   
      StepPDGID.push_back(fAllTrajectoryInfo[iTraj].vTrackPDGID[iStep]);
      StepEnergy.push_back(fAllTrajectoryInfo[iTraj].vTrackEnergy[iStep]);
      StepTime.push_back(fAllTrajectoryInfo[iTraj].vTrackTime[iStep]);

      StepInitPosX.push_back(fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].x());
      StepInitPosY.push_back(fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].y());
      StepInitPosZ.push_back(fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].z());

      StepFinalPosX.push_back(fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].x());
      StepFinalPosY.push_back(fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].y());
      StepFinalPosZ.push_back(fAllTrajectoryInfo[iTraj].vPostStepTrackPos[iStep].z());


     /*  
      if(fAllTrajectoryInfo[iTraj].fPDGCode == 111) {
        if(iStep < 20)
          G4cout << " Step = " << iStep + 1 << ", PDG = " << fAllTrajectoryInfo[iTraj].vTrackPDGID[iStep] 
                 << ", z pos = " << fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].z()
                 << ", x pos = " << fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].x() 
                 << ", y pos = " << fAllTrajectoryInfo[iTraj].vPreStepTrackPos[iStep].y()
                 << ", t pos = " << fAllTrajectoryInfo[iTraj].vTrackTime[iStep]
                 << ", E = " << fAllTrajectoryInfo[iTraj].vTrackEnergy[iStep]
                 << G4endl;
      }
     */ 

    } // for(G4int iStep = 0; iStep < N_STORE_STEP; iStep++)


    if(fAllTrajectoryInfo[iTraj].fPDGCode == 111) {
    runAction->Photon1_E = StepEnergy[1];
    runAction->Photon2_E = StepEnergy[3];

    G4cout << "Photon1 energy = " << StepEnergy[1] << G4endl;
    G4cout << "Photon2 energy = " << StepEnergy[3] << G4endl;

    runAction->Photon1_Theta = fAllTrajectoryInfo[iTraj].vTrackMomentumDir[1].theta();
    runAction->Photon2_Theta = fAllTrajectoryInfo[iTraj].vTrackMomentumDir[3].theta();

    runAction->Photon1_Phi = fAllTrajectoryInfo[iTraj].vTrackMomentumDir[1].phi();
    runAction->Photon2_Phi = fAllTrajectoryInfo[iTraj].vTrackMomentumDir[3].phi();


    } // if(fAllTrajectoryInfo[iTraj].fPDGCode == 111) storing photon from pi0


    fTrackStepID.push_back(StepID);
    fTrackStepPDGID.push_back(StepPDGID);
    fTrackStepEnergy.push_back(StepEnergy);
    fTrackStepTime.push_back(StepTime);

    fTrackStepInitPosX.push_back(StepInitPosX);
    fTrackStepInitPosY.push_back(StepInitPosY);
    fTrackStepInitPosZ.push_back(StepInitPosZ);

    fTrackStepFinalPosX.push_back(StepFinalPosX);
    fTrackStepFinalPosY.push_back(StepFinalPosY);
    fTrackStepFinalPosZ.push_back(StepFinalPosZ);


  } // for(G4int iTraj = 0; iTraj < N_STORE_TRAJ; iTraj++)

  // analysisManager->AddNtupleRow();  

  
  auto evt
    = static_cast<const G4Event*>(
        G4RunManager::GetRunManager()->GetCurrentEvent());

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // extract the trajectories and print them out
  G4cout << G4endl;
  G4cout << "Trajectories in tracker "<< n_trajectories << G4endl;

  G4cout << ", N Unique Traj RunData = " << fAllTrajectoryInfo.size() << G4endl;
    


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void H02RunData::SetTrackParentPerEvent(const G4Event* event) {

//   G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
//   G4int n_trajectories = 0;
//   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

//   //G4cout <<"N trajectories in the event : " << n_trajectories << G4endl ;

//   // extract the trajectories and print them out
  
//   for(G4int i=0; i<n_trajectories; i++) 
//     {
//       RE01Trajectory* trj = 
//         (RE01Trajectory*)((*(event->GetTrajectoryContainer()))[i]);
//       //trj->ShowTrajectory();
//       if(trj)
//       fTrackParents.push_back( trj->GetParentID() ) ;
//     }

// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void H02RunData::Reset()
{ 
  for ( auto& edep : fEdep ) { 
    edep = 0.; 
  }
  for ( auto& trackLength : fTrackLength ) {
    trackLength = 0.; 
  }

  for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
   for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
    for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
        fCell_Energy[iLayer][iEtacell][iPhicell] = 0. ;
        fCh_Cell_Energy[iLayer][iEtacell][iPhicell] = 0. ;
        fNu_Cell_Energy[iLayer][iEtacell][iPhicell] = 0. ;

      }
    }
  }

  fTotal_Ch_Energy = 0;
  fTotal_Nu_Energy = 0;

  //fTrackParents.clear();
  fAllTrajectoryInfo.clear();

  fTrackStepNumb.clear();
  fTrackPDGCode.clear();
  fTrackPDGCharge.clear();
  fTrackEnergy.clear();
  fTrackPx.clear();
  fTrackPy.clear();
  fTrackPz.clear();

  fTrackStepID.clear();
  fTrackStepPDGID.clear();
  fTrackStepEnergy.clear();
  fTrackStepTime.clear();
  fTrackStepInitPosX.clear();
  fTrackStepInitPosY.clear();
  fTrackStepInitPosZ.clear();
  fTrackStepFinalPosX.clear();
  fTrackStepFinalPosY.clear();
  fTrackStepFinalPosZ.clear();

  fTruthParticleE.clear();
  fTruthParticlePx.clear();
  fTruthParticlePy.clear();
  fTruthParticlePz.clear();
  fTruthParticleM.clear(); 
  fTruthParticlePDGID.clear();
  fTruthParticleMother1.clear();
  fTruthParticleMother2.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
