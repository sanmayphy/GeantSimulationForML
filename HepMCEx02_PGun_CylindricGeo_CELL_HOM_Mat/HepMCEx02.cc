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
/// \file eventgenerator/HepMC/HepMCEx02/HepMCEx02.cc
/// \brief Main program of the eventgenerator/HepMC/HepMCEx02 example
//
//
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example of HepMC-interface
// --------------------------------------------------------------

#include "G4Types.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "H02DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "H02PrimaryGeneratorAction.hh"
#include "H02EventAction.hh"
#include "H02SteppingAction.hh"
#include "H02RunData.hh"
#include "H02RunAction.hh"
#include "H02TrackingAction.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include <cstdlib>
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

int main(int argc, char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  int seed;
  G4cout << "argc: " << argc<< G4endl;
  G4cout << "argv: " << argv<< G4endl;
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  else if ( argc == 2 ) {
    G4long seed;
    time_t systime = time(NULL);
    seed = (long) systime;
  }
  else{
    stringstream str;
    str << argv[2];
    str >> seed;
    G4cout << argv[2] << G4endl;
    G4cout << seed << G4endl;
  }
  
   

    
  // }

  //choose the Random engine

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  // G4Random::setTheSeed(seed);
  // G4cout << "G4Random: " << G4Random::getTheSeed()<< G4endl;
  G4cout << "seeds: " << seed<< G4endl;
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << "!!!!!!!Seed: " << CLHEP::HepRandom::getTheSeed() << G4endl;


  G4RunManager* runManager= new G4RunManager;

  // User Initialization classes (mandatory)
  //
  G4VUserDetectorConstruction* detector = new H02DetectorConstruction;
  runManager-> SetUserInitialization(detector);

  
  //
  G4VUserPhysicsList* physics = new FTFP_BERT;
  runManager-> SetUserInitialization(physics);

  runManager-> Initialize();


  G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* fParticleDef = fParticleTable->FindParticle("pi0");
  G4VDecayChannel* fMode =
                  new G4PhaseSpaceDecayChannel("pi0",1,2,"gamma","gamma");


  // G4ParticleDefinition* fParticleDef1 = fParticleTable->FindParticle("pi+");
  // G4VDecayChannel* fMode1 =
  //                 new G4PhaseSpaceDecayChannel("pi+",1.00,2,"mu+","nu_mu");

  G4DecayTable* fTable = new G4DecayTable();
  fTable->Insert(fMode);
  //fTable->Insert(fMode1);
      // fMode = new G4PionRadiativeDecayChannel("pi+",0.000017);
      //  fTable->Insert(fMode);
  fParticleDef->SetDecayTable(fTable);
  //fParticleDef1->SetDecayTable(fTable);
  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new H02PrimaryGeneratorAction;
  runManager-> SetUserAction(gen_action);

  runManager-> SetUserAction(new H02RunAction);
  //
  G4UserEventAction* event_action = new H02EventAction;
  runManager-> SetUserAction(event_action);
  //

  G4UserTrackingAction* track_action = new H02TrackingAction;
  runManager-> SetUserAction(track_action);

  G4UserSteppingAction* stepping_action = new H02SteppingAction;
  runManager-> SetUserAction(stepping_action);

  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout << G4endl;

 //get the pointer to the User Interface manager   
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (!ui) { // batch mode
    visManager-> SetVerboseLevel("quiet");
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager-> ApplyCommand(command+fileName);    
  } else {  // interactive mode : define UI session

    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }

    ui-> SessionStart();
    delete ui;
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;
  delete runManager;
}

