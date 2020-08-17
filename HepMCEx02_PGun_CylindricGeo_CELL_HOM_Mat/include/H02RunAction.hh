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
/// \file H02RunAction.hh
/// \brief Definition of the H02RunAction class

#ifndef H02RunAction_h
#define H02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "GlobalVariables.hh"
#include "PCell.hh"

#include "CaloRConstants.hh"

#include "TFile.h"
#include "TTree.h"
#include <vector>

class G4Run;

using namespace std;

//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<vector<float> > >+;
//#endif


/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in B4Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

class H02RunAction : public G4UserRunAction
{
  public:
    H02RunAction();
    virtual ~H02RunAction();

    virtual G4Run* GenerateRun();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    TFile *outf;
    TTree *outTree;

    Float_t cellE;


    Int_t   NTrajectories;
    Float_t Cell_Energy[6][kMaxPixel][kMaxPixel] ;
    Float_t Ch_Cell_Energy[6][kMaxPixel][kMaxPixel] ;
    Float_t Nu_Cell_Energy[6][kMaxPixel][kMaxPixel] ;
    Float_t Cell_TopoEnergy[6][kMaxPixel][kMaxPixel] ;
    Float_t Ch_Cell_TopoEnergy[6][kMaxPixel][kMaxPixel] ;
    Float_t Nu_Cell_TopoEnergy[6][kMaxPixel][kMaxPixel] ;
    Float_t Noise_Cell_Energy[6][kMaxPixel][kMaxPixel] ;
    Float_t Noise_Cell_TopoEnergy[6][kMaxPixel][kMaxPixel] ;
    Float_t TopoCluster_Cell[6][kMaxPixel][kMaxPixel];
    Float_t PFlowArray[6][kMaxPixel][kMaxPixel];
    Float_t Cell_SuperTopoEnergy[6][kMaxPixel][kMaxPixel];
    Float_t Ch_Cell_SuperTopoEnergy[6][kMaxPixel][kMaxPixel];
    Float_t Nu_Cell_SuperTopoEnergy[6][kMaxPixel][kMaxPixel];
    Float_t Noise_Cell_SuperTopoEnergy[6][kMaxPixel][kMaxPixel];
    Float_t SuperTopoCluster_Cell[6][kMaxPixel][kMaxPixel];
    PCell PArray[6][kMaxPixel][kMaxPixel];
    
    Float_t Ch_TopoCluster_Cell[6][kMaxPixel][kMaxPixel];
    Float_t Nu_TopoCluster_Cell[6][kMaxPixel][kMaxPixel];
    Float_t Track_Cell[kMaxPixel][kMaxPixel] = {};
    std::vector<vector <double> > RealParticalInfo;
    std::vector<vector <float> > ListTopo;
    std::vector<vector <float> > ChListTopo;
    std::vector<vector <float> > NuListTopo;
    std::vector<vector <float> > NoiseListTopo;
    std::vector<vector <double> > ChargeParticlTraj;
    
    

    Float_t Total_Ch_Energy, Total_Nu_Energy;
    Float_t Total_Had_Ch_Energy, Total_Had_Nu_Energy, Total_EM_Energy;
    Float_t True_Ch_Energy, True_Nu_Energy ;
    Float_t Smeared_Ch_Energy;


    Float_t Trk_Theta_pi_p, Trk_Phi_pi_p;
    Float_t Trk_Theta_pi_0, Trk_Phi_pi_0;

    Float_t Photon1_E, Photon2_E;
    Float_t Photon1_Theta, Photon2_Theta;
    Float_t Photon1_Phi, Photon2_Phi;
    Float_t Seed;
    Float_t Kin_Ch_Energy, Kin_Nu_Energy;
    Float_t nbRadLen;//Geantino
    Float_t Eref;//Geantino

    Float_t ro0;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

