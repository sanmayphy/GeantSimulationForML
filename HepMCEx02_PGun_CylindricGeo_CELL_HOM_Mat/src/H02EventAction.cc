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
#include "TopoCluster.hh"
#include "H02TrajField.hh"
#include "H02PFlow.hh"
#include "DownResolution.hh"
// #include "Pythia8Plugins/FastJet3.h"

// #include "Cell.hh"




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
  // ! This is not active
  // printout primary information
  G4cout << "Print out primary information" << G4endl;
  G4int nVtx= anEvent-> GetNumberOfPrimaryVertex();
  G4int i;
  for(i=0; i< nVtx; i++) {
    const G4PrimaryVertex* primaryVertex= anEvent-> GetPrimaryVertex(i);
    primaryVertex-> Print();
  }
#endif
// ! Reset variables
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
  G4double ETruth = 0. ;
  auto runData
    = static_cast<H02RunData*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
   runData->FillPerEvent();

  auto runAction
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  //! Pythia8
  if (GeneratorName=="pythia8")
  {
    G4cout << "Current Generator Name : " << GeneratorName << G4endl ;

    runAction->True_Nu_Energy = 0;
    runAction->True_Ch_Energy = 0;
    int SizeRealParticalInfo = runAction->RealParticalInfo.size();
    for (int nParticle = 0; nParticle < SizeRealParticalInfo; nParticle++)
    {
      vector <double> &Particle = runAction->RealParticalInfo[nParticle];
      auto DefParticle =  G4ParticleTable::GetParticleTable()->FindParticle(runAction->RealParticalInfo[nParticle][0]);
      double charge = DefParticle->GetPDGCharge();
      if (charge==0)
      {
        runAction->True_Nu_Energy+=Particle[4];
      }
      else
      {
        runAction->True_Ch_Energy+=Particle[4];
      }
    }
    G4VPrimaryGenerator* primaryGen = runGeneratorAction->GetGenerator();
    HepMCG4Pythia8Interface *py8Gen = (HepMCG4Pythia8Interface*)primaryGen ;

    Pythia8::Event pythiaEvent = py8Gen->GetPythiaObject();


    G4int NTruth = 0;

    for(int ip = 0; ip < (int)pythiaEvent.size(); ip++)
    {
      if(pythiaEvent[ip].isFinal()) 
      {
        cout <<"Particle : " << ip+1 <<", ID = "<< pythiaEvent[ip].id()
        << " , Energy val = " << pythiaEvent[ip].e()
        << ", Mother1 = " << pythiaEvent[ip].mother1()
        << ", Mother2 = " << pythiaEvent[ip].mother2()
        << ", Px = "
        << pythiaEvent[ip].px()
        << ", Py = "
        << pythiaEvent[ip].py()
        << ", Pz = "
        << pythiaEvent[ip].pz()
        // << ", Production time = " << pythiaEvent[ip].tProd()
        << endl;

        NTruth ++ ;
        ETruth += pythiaEvent[ip].e() ;
      }
    }

      G4cout << "Truth Particle in event : " << NTruth << G4endl ;
    //? ChargeTraj calculation
    G4cout << "Check1 " << G4endl ;
    runAction->ChargeParticlTraj.clear();
    ChargeParticlTraj *ChTraj = new ChargeParticlTraj(runAction->RealParticalInfo, runAction->ChargeParticlTraj, runAction->Track_Cell);
    G4cout << "Check2 " << G4endl ;
    delete (ChTraj);
    G4cout << "Check!!! " << NTruth << G4endl ;
    // double px_1 = runAction->RealParticalInfo[1][1];
    // double py_1 = runAction->RealParticalInfo[1][2];
    // double pz_1 = runAction->RealParticalInfo[1][3];
    // double px_2 = runAction->RealParticalInfo[2][1];
    // double py_2 = runAction->RealParticalInfo[2][2];
    // double pz_2 = runAction->RealParticalInfo[2][3];
    // double theta = (px_1*px_2 + py_1*py_2 + pz_1*pz_2)/((px_1*px_1 + py_1*py_1 + pz_1*pz_1) * (px_2*px_2 + py_2*py_2 + pz_2*pz_2));
    // runAction->ro0 = theta*10e+5;
    //? ChargeTraj calculation
  }
  //! Pythia8 END
   // ===========  Fill the ntuples ============ //

    for (int l = 0; l < 6; l++)
    {
      for (int i = 0; i < 256; i++)
      {
        for (int j = 0; j < 256; j++)
        {
          runAction->Noise_Cell_TopoEnergy[l][i][j] = 0;
        }
      }
    }
    //? Down Resoulution
    DownResolution *Down_ptr = new DownResolution(runAction->Cell_Energy, runAction->Noise_Cell_Energy,
                                                  runAction->Ch_Cell_Energy,runAction->Nu_Cell_Energy);
    Down_ptr->FillTotEn(runAction->Cell_TopoEnergy);
    Down_ptr->FillChEn(runAction->Ch_Cell_TopoEnergy);
    Down_ptr->FillNuEn(runAction->Nu_Cell_TopoEnergy);
    Down_ptr->FillNoise(runAction->Noise_Cell_TopoEnergy);
    delete (Down_ptr);
    //? Down Resoulution end

    


    //? TopoClisters maker
    try
    {
      TopoCluster *ptr_1 = new TopoCluster(runAction->Cell_SuperTopoEnergy, runAction->Noise_Cell_TopoEnergy,
                                          runAction->Ch_Cell_TopoEnergy,runAction->Nu_Cell_TopoEnergy, LayersPix);
      ptr_1->Run();
      ptr_1->FillNoise(runAction->Cell_SuperTopoEnergy);
      runAction->NoiseListTopo = ptr_1->FinalListClust;
      delete (ptr_1);
      TopoCluster *ptr_2 = new TopoCluster(runAction->Ch_Cell_TopoEnergy, runAction->Noise_Cell_TopoEnergy,
                                          runAction->Ch_Cell_TopoEnergy,runAction->Nu_Cell_TopoEnergy, LayersPix);
      ptr_2->Run();
      runAction->ChListTopo = ptr_2->FinalListClust;
      delete (ptr_2);
      TopoCluster *ptr_3 = new TopoCluster(runAction->Nu_Cell_TopoEnergy, runAction->Noise_Cell_TopoEnergy,
                                          runAction->Ch_Cell_TopoEnergy,runAction->Nu_Cell_TopoEnergy, LayersPix);
      ptr_3->Run();
      runAction->NuListTopo = ptr_3->FinalListClust;
      delete (ptr_3);
      
      G4cout<<"runAction->NoiseListTopo "<<runAction->NoiseListTopo.size()<<G4endl;
      G4cout<<"runAction->ChListTopo "<<runAction->ChListTopo.size()<<G4endl;
      G4cout<<"runAction->NuListTopo "<<runAction->NuListTopo.size()<<G4endl;
      // ptr_1->FillTopo(runAction->SuperTopoCluster_Cell   );
      // G4cout<<"A4"<<G4endl;
      // ptr_1->FillTotEn(runAction->Cell_SuperTopoEnergy  );
      // // G4cout<<"A5"<<G4endl;
      // ptr_1->FillChEn(runAction->Ch_Cell_SuperTopoEnergy);
      // // G4cout<<"A6"<<G4endl;
      // ptr_1->FillNuEn(runAction->Nu_Cell_SuperTopoEnergy);
      // // G4cout<<"A7"<<G4endl;
      // ptr_1->FillNoise(runAction->Noise_Cell_SuperTopoEnergy);
      // G4cout<<"A8"<<G4endl;
      
    }
    catch (...)// (const std::bad_alloc& e) 
    {
        G4cout << "AntonC Allocation failed" << G4endl;
    }
    try
    {
      // G4cout<<"A0"<<G4endl;
      TopoCluster *ptr_1 = new TopoCluster(runAction->Cell_TopoEnergy, runAction->Noise_Cell_TopoEnergy,
                                          runAction->Ch_Cell_TopoEnergy,runAction->Nu_Cell_TopoEnergy, LayersPix);
      // G4cout<<"A1"<<G4endl;
      ptr_1->Run();
      // G4cout<<"A2"<<G4endl;
      runAction->ListTopo = ptr_1->FinalListClust;
      // G4cout<<"A3"<<G4endl;
      ptr_1->FillTopo(runAction->TopoCluster_Cell   );
      // G4cout<<"A4"<<G4endl;
      ptr_1->FillTotEn(runAction->Cell_TopoEnergy  );
      // G4cout<<"A5"<<G4endl;
      ptr_1->FillChEn(runAction->Ch_Cell_TopoEnergy);
      // G4cout<<"A6"<<G4endl;
      ptr_1->FillNuEn(runAction->Nu_Cell_TopoEnergy);
      // G4cout<<"A7"<<G4endl;
      // ptr_1->FillNoise(runAction->Noise_Cell_TopoEnergy);
      // G4cout<<"A8"<<G4endl;
      delete (ptr_1);
    }
    catch (...)// (const std::bad_alloc& e) 
    {
        G4cout << "AntonC Allocation failed" << G4endl;
    }


    
    //? TopoClisters maker end
    float E_ref = 0;
    int numberOfTopocLust = runAction->ListTopo.size();
    for (int track = 0; track < runAction->ChargeParticlTraj.size(); track++)
    {
      for (int clust = 0;clust < numberOfTopocLust; clust++ )
      {
        float deta = runAction->ChargeParticlTraj[track][3] - runAction->ListTopo[clust][1];
            // G4cout<<"A9.3"<<G4endl;
        float dphi = abs(runAction->ChargeParticlTraj[track][4] - runAction->ListTopo[clust][2]);
        // G4cout<<"A9.4"<<G4endl;
        if (dphi > M_PI)
        {
            dphi -= 2*M_PI;
        }
        float R = sqrt(sqr(dphi)+sqr(deta));
        if (R<0.4)
        {
          E_ref += runAction->ListTopo[clust][0];
        }
      }
    }
    runAction->Eref = E_ref;

  

    // //? PF
    PFlow *pflow = new PFlow(runAction->ChargeParticlTraj, runAction->ListTopo, runAction->Cell_TopoEnergy, runAction->TopoCluster_Cell); 
    pflow->Fill(runAction->PFlowArray);
    // //? PF
    float nuetralTopo = 0;
    float pflowneutral = 0;
    float Chenerg = 0 ;
    for (int i = 0; i < 6; i++)
    {
        
        for (int j = 0; j < LayersPix[i]; j++)
        {
            for (int k = 0; k < LayersPix[i]; k++)
            {
                nuetralTopo +=  runAction->Nu_Cell_TopoEnergy[i][j][k];
                pflowneutral += runAction->PFlowArray[i][j][k];
                Chenerg += runAction->Ch_Cell_TopoEnergy[i][j][k];
            }
        }
    }
    // runAction->ro0 = (pflowneutral-nuetralTopo)/nuetralTopo;
    runAction->outTree->Fill();
    // G4cout<<"Chenerg: "<<Chenerg<<G4endl;
    // G4cout<<"nuetralTopo: "<<nuetralTopo<<G4endl;
    // G4cout<<"pflowneutral: "<<pflowneutral<<G4endl;
    // G4cout<<"Relative Residiual: "<<(pflowneutral-nuetralTopo)/nuetralTopo<<G4endl;
    // vector < vector <float> > TrackVSTopoCLust(numberOfTrack, vector <float> (numberOfTopocLust,0));
    // for (int i = 0; i < 256; i++)
    // {
    //   for (int j = 0; j < 256; j++)
    //   {
    //     G4cout<<"CellCoord: "<<i<<" "<<j<<G4endl;
    //     G4cout<<"CellChen1:  "<<runData->fPArray[0][i][j].chEnergy<<G4endl;
    //     G4cout<<"CellNuen1:  "<<runData->fPArray[0][i][j].nuEnergy<<G4endl;
    //     G4cout<<"CellChen2:  "<<runAction->Ch_Cell_Energy[0][i][j]<<G4endl;
    //     G4cout<<"CellChen2:  "<<runAction->Nu_Cell_Energy[0][i][j]<<G4endl;
    //     int sizep = runData->fPArray[0][i][j].prtcl.size() ;
    //     for (int k = 0; k<sizep;k++)
    //     {
    //       G4cout<<"CellPar:  "<<runData->fPArray[0][i][j].prtcl[k].PosInList <<G4endl;
    //       G4cout<<"CellEn:   "<<runData->fPArray[0][i][j].prtcl[k].Energy <<G4endl;
    //     }
        
    //   }
    // }


    


   // =================================================== //
   // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // extract the trajectories and print them out
  G4cout << G4endl;
  G4cout << "Trajectories in tracker "<< n_trajectories << G4endl;

  G4int nUniqueTraj(0), nSecTraj(0);
  for(G4int i=0; i<n_trajectories; i++){
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
    double TotalInitialE_B(0);
    double TotalInitialE_A(0);
    // G4cout << "fAllTrajectoryInfo.size() "<< (int)runData->fAllTrajectoryInfo.size() << G4endl;
  for(int iTraj = (int)runData->fAllTrajectoryInfo.size()-1; iTraj >=0; iTraj--)
  {
    TotalInitialE = 0;
    for(int iParent = 1; iParent < (int)runData->fAllTrajectoryInfo[iTraj].vTrackID.size(); iParent++ )
    {
      TotalSecondary ++ ;
      TotalInitialE += runData->fAllTrajectoryInfo[iTraj].vTrackEnergy.at(iParent);
          // G4cout << "Step " << iParent
          // << ", Energy = " << runData->fAllTrajectoryInfo[iTraj].vTrackEnergy.at(iParent)
          // << G4endl ;
    }


      // G4cout << " Traj : " << iTraj +1
      // << " ID = "
      // <<runData->fAllTrajectoryInfo[iTraj].fPDGCode
      // << ", Energy = "
      // << runData->fAllTrajectoryInfo[iTraj].fEnergy
      // << ", Px = "
      // << runData->fAllTrajectoryInfo[iTraj].fMomentum.x()
      // << ", Py = "
      // << runData->fAllTrajectoryInfo[iTraj].fMomentum.y()
      // << ", Pz = "
      // << runData->fAllTrajectoryInfo[iTraj].fMomentum.z()
      // << ", TotalStepE = "
      // << TotalInitialE
      // << G4endl;
      TotalInitialE_A += runData->fAllTrajectoryInfo[iTraj].fEnergy;
      TotalInitialE_B += TotalInitialE;


    // 	// G4cout << "Traj : " << iTraj +1 << ", Particle = " 
    // 	// << runData->fAllTrajectoryInfo[iTraj].fPDGCode 
    // 	// << ", Momentum : " 
    // 	<< runData->fAllTrajectoryInfo[iTraj].fEnergy / (1000 * CLHEP::GeV )
    // 	// << ", Vertex X = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].x() 
    // 	// << ", Vertex Y = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].y()
    // 	// << ", Vertex Z = " << runData->fAllTrajectoryInfo[iTraj].vPreStepTrackPos[1].z()
    // 	// << ", Charge = " << runData->fAllTrajectoryInfo[iTraj].fPDGCharge
    //  //  << ", Secondary showers = " << runData->fAllTrajectoryInfo[iTraj].vTrackID.size()
    // 	// << G4endl;
  }
    if (GeneratorName=="pythia8")
    {
      G4cout << "Unique Traj " << nUniqueTraj <<", Secondary Traj =  " << nSecTraj
      << ", N Unique Traj RunData = " << runData->fAllTrajectoryInfo.size()
      << ", N Secondary Traj RunData = " <<  TotalSecondary
      << ", Truth-E = "<< ETruth
      << ", Traj-E_B = " << TotalInitialE_B * 1e-3
      << ", Traj-E = " << TotalInitialE_A * 1e-3//* 1e-6

      << G4endl ;
    }

    G4cout << " Total Neutral energy " << runAction->Total_Nu_Energy << G4endl;
    G4cout << " Total Charged energy " << runAction->Total_Ch_Energy << G4endl;
    G4cout << " Charged and nuetral  energy " << runAction->Total_Ch_Energy + runAction->Total_Nu_Energy<< G4endl;
    // G4cout << " Charged and nuetral  energy " << runAction->Total_Esc_Energy+runAction->Total_Ch_Energy +  runAction->Total_Nu_Energy<< G4endl;
    // G4cout << " Total Dep trk energy " << runAction->Total_Enotrc_Energy << G4endl;
    // G4cout << " strange " << runAction->Total_Esc_Energy << G4endl;
    // G4cout << " Total not trk energy " << runAction->Total_Enotrc_Energy << G4endl;
    // G4cout << " Total with not trk energy " << runAction->Total_Enotrc_Energy+
    // runAction->Total_Esc_Energy+
    // runAction->Total_Ch_Energy +  runAction->Total_Nu_Energy<< G4endl;
    // fTotal_Enotrc_Energy
    // G4cout << " Total energy " << runData->fAllTrajectoryInfo << G4endl;
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
  // fCell_Energy;



}
