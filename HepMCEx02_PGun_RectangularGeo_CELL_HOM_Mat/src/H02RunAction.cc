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
/// \file H02RunAction.cc
/// \brief Implementation of the H02RunAction class

#include "H02RunAction.hh"
#include "H02RunData.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "CaloRConstants.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

H02RunAction::H02RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  
  // // Set up sizes. (HEIGHT x WIDTH)
  // Cell_Energy.resize(kNLayers);
  // Ch_Cell_Energy.resize(kNLayers);
  // Nu_Cell_Energy.resize(kNLayers);

  // for (int i_layer = 0; i_layer < kNLayers; i_layer++) {
  //   Cell_Energy[i_layer].resize(CELL_NUMBER[i_layer]);
  //   Ch_Cell_Energy[i_layer].resize(CELL_NUMBER[i_layer]);
  //   Nu_Cell_Energy[i_layer].resize(CELL_NUMBER[i_layer]);

  //   for (int x_cell = 0; x_cell < CELL_NUMBER[i_layer]; x_cell++) {
  //     Cell_Energy[i_layer][x_cell].resize(CELL_NUMBER[i_layer]);
  //     Ch_Cell_Energy[i_layer][x_cell].resize(CELL_NUMBER[i_layer]); 
  //     Nu_Cell_Energy[i_layer][x_cell].resize(CELL_NUMBER[i_layer]);
  //   }
  // }

  // Book histograms, ntuple
  //
  
  // // Creating histograms
  // analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
  // analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
  // analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  // analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

  // // Creating ntuple
  // //
  // analysisManager->CreateNtuple("B4", "Edep and TrackL");
  // analysisManager->CreateNtupleDColumn("Eabs");
  // analysisManager->CreateNtupleDColumn("Egap");
  // analysisManager->CreateNtupleDColumn("Labs");
  // analysisManager->CreateNtupleDColumn("Lgap");


  // // int total_bins = 10 * 100 * 100;  // 3 overflow bins for the three calo layers

  // // for (int i = 0; i < total_bins; ++i) {

  // //   std::stringstream out;
  // //   out << i;
  // //   analysisManager->CreateNtupleDColumn("cell_" + out.str());
  // // }
  // // analysisManager->CreateNtupleDColumn("TotalEnergy");


  // char name[500]; 

  // for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
  //   for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
  //     for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
        
  //       sprintf (name, "cell_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       analysisManager->CreateNtupleDColumn(name);
  //       sprintf (name, "cellCh_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       analysisManager->CreateNtupleDColumn(name);
  //       sprintf (name, "cellNu_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       analysisManager->CreateNtupleDColumn(name);

  //     }
  //   }
  // }

  // analysisManager->CreateNtupleDColumn("NTrajectories");

  // for(G4int iTraj = 0; iTraj < N_STORE_TRAJ; iTraj++){

  //     sprintf (name, "Trajectory%i_NSTEP", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_PDGID", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_PDGCharge", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Energy", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Px", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Py", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Pz", iTraj+1);
  //     analysisManager->CreateNtupleDColumn(name);

  //   for(G4int iStep = 0; iStep < N_STORE_STEP; iStep++){

  //     sprintf (name, "Trajectory%i_Step%i_ID", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_PDGID", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_E", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_Time", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);

  //     sprintf (name, "Trajectory%i_Step%i_PosX", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_PosY", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_PosZ", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);

  //     sprintf (name, "Trajectory%i_Step%i_InitPosX", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_InitPosY", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_InitPosZ", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);

  //     sprintf (name, "Trajectory%i_Step%i_FinalPosX", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_FinalPosY", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);
  //     sprintf (name, "Trajectory%i_Step%i_FinalPosZ", iTraj+1, iStep+1);
  //     analysisManager->CreateNtupleDColumn(name);


  //   } // for(G4int iStep = 0; iStep < N_STORE_STEP; iStep++)
  // } // for(G4int iTraj = 0; iTraj < N_STORE_TRAJ; iTraj++)
  

  // analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

H02RunAction::~H02RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* H02RunAction::GenerateRun()
{
  return (new H02RunData);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02RunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "PFlowNtuple";
  //analysisManager->OpenFile(fileName);

  outf = new TFile("PFlowNtupleFile_QCD.root", "RECREATE");
  outTree = new TTree("EventTree", "EventTree");

  auto runData = static_cast<H02RunData*>
    (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  Float_t new_v  ;
  char name[500];

   
  // outTree->Branch("Trajectory_NStep","vector<double>", &runData->fTrackStepNumb);
  // outTree->Branch("Trajectory_PDGID","vector<double>", &runData->fTrackPDGCode);
  // outTree->Branch("Trajectory_PDGCharge","vector<double>", &runData->fTrackPDGCharge);
  // outTree->Branch("Trajectory_Energy","vector<double>", &runData->fTrackEnergy);
  // outTree->Branch("Trajectory_Px","vector<double>", &runData->fTrackPx);
  // outTree->Branch("Trajectory_Py","vector<double>", &runData->fTrackPy);
  // outTree->Branch("Trajectory_Pz","vector<double>", &runData->fTrackPz);

  // outTree->Branch("Trajectory_StepID","vector<vector<double>>", &runData->fTrackStepID);
  // outTree->Branch("Trajectory_StepPDGID","vector<vector<double>>", &runData->fTrackStepPDGID);
  // outTree->Branch("Trajectory_StepEnergy","vector<vector<double>>", &runData->fTrackStepEnergy);
  // outTree->Branch("Trajectory_StepTime","vector<vector<double>>", &runData->fTrackStepTime);

  // outTree->Branch("Trajectory_StepInitPosX","vector<vector<double>>", &runData->fTrackStepInitPosX);
  // outTree->Branch("Trajectory_StepInitPosY","vector<vector<double>>", &runData->fTrackStepInitPosY);
  // outTree->Branch("Trajectory_StepInitPosZ","vector<vector<double>>", &runData->fTrackStepInitPosZ);

  // outTree->Branch("Trajectory_StepFinalPosX","vector<vector<double>>", &runData->fTrackStepFinalPosX);
  // outTree->Branch("Trajectory_StepFinalPosY","vector<vector<double>>", &runData->fTrackStepFinalPosY);
  // outTree->Branch("Trajectory_StepFinalPosZ","vector<vector<double>>", &runData->fTrackStepFinalPosZ);

  // outTree->Branch("TruthParticleE","vector<double>", &runData->fTruthParticleE);
  // outTree->Branch("TruthParticlePx","vector<double>", &runData->fTruthParticlePx);
  // outTree->Branch("TruthParticlePy","vector<double>", &runData->fTruthParticlePy);
  // outTree->Branch("TruthParticlePz","vector<double>", &runData->fTruthParticlePz);
  // outTree->Branch("TruthParticleM","vector<double>", &runData->fTruthParticleM);
  // outTree->Branch("TruthParticleMother1","vector<double>", &runData->fTruthParticleMother1);
  // outTree->Branch("TruthParticleMother2","vector<double>", &runData->fTruthParticleMother2);

  //outTree->Branch("Event_Cell_Energy", "Cell_Energy[kLayer][K_NETA][K_NPHI]/F", &runData->fVariable);  
  
  
  outTree->Branch("new_v", &new_v, "new_v/F");
  outTree->Branch("Cell_E", &cellE, "Cell_E/F");

  outTree->Branch("Total_Ch_Energy", &Total_Ch_Energy, "Total_Ch_Energy/F");
  outTree->Branch("Total_Nu_Energy", &Total_Nu_Energy, "Total_Nu_Energy/F");

  outTree->Branch("True_Ch_Energy", &True_Ch_Energy, "True_Ch_Energy/F");
  outTree->Branch("True_Nu_Energy", &True_Nu_Energy, "True_Nu_Energy/F");
  outTree->Branch("Smeared_Ch_Energy", &Smeared_Ch_Energy, "Smeared_Ch_Energy/F");

  outTree->Branch("Trk_X_indx", &Trk_X_indx, "Trk_X_indx/I");
  outTree->Branch("Trk_Y_indx", &Trk_Y_indx, "Trk_Y_indx/I");

  outTree->Branch("Trk_X_pos", &Trk_X_pos, "Trk_X_pos/F");
  outTree->Branch("Trk_Y_pos", &Trk_Y_pos, "Trk_Y_pos/F");

  outTree->Branch("Trk_Theta", &Trk_Theta, "Trk_Theta/F");
  outTree->Branch("Trk_Phi", &Trk_Phi, "Trk_Phi/F");

  outTree->Branch("Pi0_X_pos", &Pi0_X_pos, "Pi0_X_pos/F");
  outTree->Branch("Pi0_Y_pos", &Pi0_Y_pos, "Pi0_Y_pos/F");

  outTree->Branch("Pi0_Theta", &Pi0_Theta, "Pi0_Theta/F");
  outTree->Branch("Pi0_Phi", &Pi0_Phi, "Pi0_Phi/F");

  outTree->Branch("Photon1_E", &Photon1_E, "Photon1_E/F");
  outTree->Branch("Photon2_E", &Photon2_E, "Photon2_E/F");

  outTree->Branch("Photon1_Theta", &Photon1_Theta, "Photon1_Theta/F");
  outTree->Branch("Photon2_Theta", &Photon2_Theta, "Photon2_Theta/F");

  outTree->Branch("Photon1_Phi", &Photon1_Phi, "Photon1_Phi/F");
  outTree->Branch("Photon2_Phi", &Photon2_Phi, "Photon2_Phi/F");

  outTree->Branch("cell_Energy", &Cell_Energy, "Cell_Energy[6][128][128]/F");
  outTree->Branch("cellCh_Energy", &Ch_Cell_Energy, "Ch_Cell_Energy[6][128][128]/F");
  outTree->Branch("cellNu_Energy", &Nu_Cell_Energy, "Nu_Cell_Energy[6][128][128]/F");
  outTree->Branch("Noise_cell_Energy", &Noise_Cell_Energy, "Noise_Cell_Energy[6][128][128]/F");
  // outTree->Branch("cell_Energy", "std::vector< std::vector< std::vector< Float_t > > >", &Cell_Energy);
  // outTree->Branch("cellCh_Energy", "std::vector< std::vector< std::vector< Float_t > > >", &Ch_Cell_Energy);
  // outTree->Branch("cellNu_Energy", "std::vector< std::vector< std::vector< Float_t > > >", &Nu_Cell_Energy);


  // for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
  //   for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
  //     for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){
        
  //       sprintf (name, "cell_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       outTree->Branch(name, &Cell_Energy[iLayer][iEtacell][iPhicell], name);
  //       sprintf (name, "cellCh_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       outTree->Branch(name, &Ch_Cell_Energy[iLayer][iEtacell][iPhicell], name);
  //       sprintf (name, "cellNu_Layer%i_EtaCell%i_PhiCell%i", iLayer, iEtacell, iPhicell);
  //       outTree->Branch(name, &Nu_Cell_Energy[iLayer][iEtacell][iPhicell], name);

  //     }
  //   }
  // }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  // if ( analysisManager->GetH1(1) ) {
  //   G4cout << G4endl << " ----> print histograms statistic ";
  //   if(isMaster) {
  //     G4cout << "for the entire run " << G4endl << G4endl; 
  //   }
  //   else {
  //     G4cout << "for the local thread " << G4endl << G4endl; 
  //   }
    
  //   G4cout << " EAbs : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    
  //   G4cout << " EGap : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
  //   G4cout << " LAbs : mean = " 
  //     << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length") 
  //     << " rms = " 
  //     << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

  //   G4cout << " LGap : mean = " 
  //     << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
  //     << " rms = " 
  //     << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
  // }

  // // save histograms & ntuple
  // //
  // analysisManager->Write();
  // analysisManager->CloseFile();

  outf->cd();
  outTree->Write();
  outf->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
