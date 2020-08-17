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

#include "HepMCG4Interface.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// #include "CaloRConstants.hh"
// #include <vector>

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
  // int xRange =256;
  // int yRange = 500
  // vector< vector<int> > vec(xRange, vector<int>(yRange, initialValue));
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
  //* Sum of charge energy across all cells  deposits; How it is look like after all events generated Total_Ch_Energy[nevent];
  outTree->Branch("Total_Ch_Energy", &Total_Ch_Energy, "Total_Ch_Energy/F");

  //* Sum of neutral energy across all cells  deposits; How it is look like after all events generated Total_Nu_Energy[nevent];
  outTree->Branch("Total_Nu_Energy", &Total_Nu_Energy, "Total_Nu_Energy/F");

  //* Sum of charge energy across all true particles; How it is look like after all events generated True_Ch_Energy[nevent];
  outTree->Branch("True_Ch_Energy", &True_Ch_Energy, "True_Ch_Energy/F");

  //* Sum of neutral energy across all true particles; How it is look like after all events generated True_Nu_Energy[nevent];
  outTree->Branch("True_Nu_Energy", &True_Nu_Energy, "True_Nu_Energy/F");

  // outTree->Branch("Smeared_Ch_Energy", &Smeared_Ch_Energy, "Smeared_Ch_Energy/F");

  //* Array with cells with total energy(Charge+Neutral+Noise), which are in topoclusters
  outTree->Branch("cell_TopoEnergy", &Cell_TopoEnergy, "cell_TopoEnergy[6][256][256]/F");

  //* Array with all cells with charge energy(Charge)
  outTree->Branch("cellCh_Energy", &Ch_Cell_Energy, "Ch_Cell_Energy[6][256][256]/F");

  //* Array with all cells with charge energy(Neutral)
  outTree->Branch("cellNu_Energy", &Nu_Cell_Energy, "Nu_Cell_Energy[6][256][256]/F");

  //* Array with cells with charge energy(Charge), which are in topoclusters
  outTree->Branch("cellCh_TopoEnergy", &Ch_Cell_TopoEnergy, "Ch_Cell_TopoEnergy[6][256][256]/F");
  

  //* Array with cells with neutral energy(Neutral), which are in topoclusters
  outTree->Branch("cellNu_TopoEnergy", &Nu_Cell_TopoEnergy, "Nu_Cell_TopoEnergy[6][256][256]/F");

  //* Array with cells with noise, which are in topoclusters
  outTree->Branch("Noise_Cell_TopoEnergy", &Noise_Cell_TopoEnergy, "Noise_Cell_TopoEnergy[6][256][256]/F");

    //* Array with cells with total energy(Charge+Neutral+Noise), which are in topoclusters
  outTree->Branch("cell_SuperTopoEnergy", &Cell_SuperTopoEnergy, "cell_SuperTopoEnergy[6][256][256]/F");

  //* Array with cells with charge energy(Charge), which are in topoclusters
  outTree->Branch("cellCh_SuperTopoEnergy", &Ch_Cell_SuperTopoEnergy, "Ch_Cell_SuperTopoEnergy[6][256][256]/F");
  

  //* Array with cells with neutral energy(Neutral), which are in topoclusters
  outTree->Branch("cellNu_SuperTopoEnergy", &Nu_Cell_SuperTopoEnergy, "Nu_Cell_SuperTopoEnergy[6][256][256]/F");

  //* Array with cells with noise, which are in topoclusters
  outTree->Branch("Noise_Cell_SuperTopoEnergy", &Noise_Cell_SuperTopoEnergy, "Noise_Cell_SuperTopoEnergy[6][256][256]/F");

  //* Array with cells with noise
  // outTree->Branch("Noise_Cell_Energy", &Noise_Cell_Energy, "Noise_Cell_Energy[6][256][256]/F");

  //* Array with cells with topological clusters labels
  outTree->Branch("TopoCluster_Cell", &TopoCluster_Cell, "TopoCluster_Cell[6][256][256]/F");
  //* Array with cells with topological clusters labels
  outTree->Branch("SuperTopoCluster_Cell", &SuperTopoCluster_Cell, "SuperTopoCluster_Cell[6][256][256]/F");
  outTree->Branch("Ch_TopoCluster_Cell", &Ch_TopoCluster_Cell, "Ch_TopoCluster_Cell[6][256][256]/F");
  outTree->Branch("Nu_TopoCluster_Cell", &Nu_TopoCluster_Cell, "Nu_TopoCluster_Cell[6][256][256]/F");

  //* Array (image) of tracks
  outTree->Branch("Track_Cell", &Track_Cell, "Track_Cell[256][256]/F");

  //* Particle Flow Array with neutral energy
  // outTree->Branch("PFlowArray", &PFlowArray, "PFlowArray[6][256][256]/F");

  //* Information about real particles that was generated in Pythia8.
  outTree->Branch("RealParticalInfo", "vector<vector<double> >", &RealParticalInfo);

  //* Information about trajectory of charge particles in magnetic field.[PDG, Enrgy, momentum, Eta, Phi ,indexEta, indExphi, \Delta R(distance to cluster), Associated Cluster ]
  //* eta[-1.5,1.5] indexEta[0,256]; phi[0,2\pi] phi[0,256]
  outTree->Branch("ChargeParticlTraj", "vector<vector<double> >", &ChargeParticlTraj);

  //* List of topoclusters [Energy, eta, phi, z, sigma_eta, sigma_phi Energy_ch, Energy_nu, Noise]; label of cluster is its serial.
  outTree->Branch("ListTopo", "vector<vector<float> >", &ListTopo);

  //* List of topoclusters [Energy, eta, phi, z, sigma_eta, sigma_phi Energy_ch, Energy_nu, Noise]; label of cluster is its serial.
  outTree->Branch("ChListTopo", "vector<vector<float> >", &ChListTopo);

//* List of topoclusters [Energy, eta, phi, z, sigma_eta, sigma_phi Energy_ch, Energy_nu, Noise]; label of cluster is its serial.
  outTree->Branch("NuListTopo", "vector<vector<float> >", &NuListTopo);

  //* List of topoclusters [Energy, eta, phi, z, sigma_eta, sigma_phi Energy_ch, Energy_nu, Noise]; label of cluster is its serial.
  outTree->Branch("NoiseListTopo", "vector<vector<float> >", &NoiseListTopo);

  //* RandomSeed for making recovery event easy.
  outTree->Branch("Seed", &Seed, "Seed/F");

  //* Some variables.
  outTree->Branch("SomeAddlVar", &ro0, "SomeAddlVar/F");

  //* E_ref for PFlow
  outTree->Branch("Eref", &Eref, "Eref/F");




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void H02RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{
  cout<<"End"<<endl;
  auto analysisManager = G4AnalysisManager::Instance();

  outf->cd();
  outTree->Write();
  outf->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
