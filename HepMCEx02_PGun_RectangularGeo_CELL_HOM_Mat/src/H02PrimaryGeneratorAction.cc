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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02PrimaryGeneratorAction.cc
/// \brief Implementation of the H02PrimaryGeneratorAction class
//
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"
#include "HepMCG4Pythia8Interface.hh"
#include "H02PrimaryGeneratorAction.hh"
#include "H02PrimaryGeneratorMessenger.hh"
#include "G4GeneralParticleSource.hh"

#include "H02RunData.hh"
#include "H02RunAction.hh"
#include "G4RunManager.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorAction::H02PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  // default generator is particle gun.
  fCurrentGenerator= fParticleGun= new G4ParticleGun();
  //fCurrentGenerator= fParticleGun= new G4GeneralParticleSource();
  fCurrentGeneratorName= "particleGun";
  fHepmcAscii= new HepMCG4AsciiReader();
//#ifdef G4LIB_USE_PYTHIA
  fPythiaGen= new HepMCG4Pythia8Interface();
//#else
//  fPythiaGen= 0;
//#endif

  fParticleGun_ = new G4ParticleGun(1);

  fGentypeMap["particleGun"]= fParticleGun;
  fGentypeMap["hepmcAscii"]= fHepmcAscii;
  fGentypeMap["pythia8"]= fPythiaGen;

  fMessenger= new H02PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorAction::~H02PrimaryGeneratorAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fCurrentGenerator) {

  	//fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    //fParticleGun-> GeneratePrimaryVertex(anEvent);


  // G4double Theta1 = CLHEP::pi/2. - CLHEP::pi/20. + CLHEP::pi/10. *G4UniformRand();
  // G4double Phi1 = CLHEP::twopi*G4UniformRand() ;

  // G4double Theta2 = Theta1 - CLHEP::pi/16. + CLHEP::pi/8. *G4UniformRand();
  // G4double Phi2 = Phi1 - CLHEP::pi/16. + CLHEP::pi/8. *G4UniformRand();

  G4double Theta1 = 0;
  G4double Phi1 = CLHEP::twopi*G4UniformRand() ;

  G4double Theta2 = Theta1 + CLHEP::pi/100. + CLHEP::pi/200. *G4UniformRand();
  //G4double Phi2 = CLHEP::twopi*G4UniformRand();
  G4double Phi2 = Phi1 - CLHEP::pi/30. + CLHEP::pi/60. *G4UniformRand();


  G4int sign = 0;
  if(G4UniformRand() > 0.5) sign = 1;
  else                      sign = -1;
/*
  G4double Theta3 = Theta2 + sign * ( CLHEP::pi/30 + G4UniformRand() * CLHEP::pi/30  ) ;
  G4double Phi3   = Phi2   + sign * ( CLHEP::pi/30 + G4UniformRand() * CLHEP::pi/30  ) ;*/

 G4double Theta3 = Theta2 + sign * ( CLHEP::pi/200 + G4UniformRand() * CLHEP::pi/200  ) ;
  G4double Phi3   = Phi2   + sign * ( CLHEP::pi/200 + G4UniformRand() * CLHEP::pi/200  ) ; 

  fParticleGun_ = nullptr ;
  fParticleGun_ = new G4ParticleGun(1);

  double L_lim = 17.5;
  double U_lim = 20.0;

  double ChEnergy = (L_lim + (U_lim - L_lim) * G4UniformRand() ) * GeV;
 // double NuEnergy = (15 + 15 * G4UniformRand() ) * GeV;

   //double ChEnergy = (20  ) * GeV;
   //double NuEnergy = ChEnergy ;
   double NuEnergy = ( (L_lim + U_lim)/2.  ) * GeV;
  
  auto runAction
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();

  G4double xSource, ySource ;
  G4double Off_Source = 10 * cm ;
  
  G4double posRand = G4UniformRand();

  if( posRand < 0.25 ) { xSource = Off_Source; ySource = Off_Source; runAction->Trk_X_indx = 5; runAction->Trk_Y_indx = 5;}
  if( posRand > 0.25 && posRand < 0.50 ) { xSource = Off_Source; ySource = -1 * Off_Source; runAction->Trk_X_indx = 5; runAction->Trk_Y_indx = -5;}
  if( posRand > 0.50 && posRand < 0.75 ) { xSource = -1 * Off_Source; ySource = Off_Source; runAction->Trk_X_indx = -5; runAction->Trk_Y_indx = 5;}
  if( posRand > 0.75 && posRand < 1.00 ) { xSource = -1 * Off_Source; ySource = -1 * Off_Source; runAction->Trk_X_indx = -5; runAction->Trk_Y_indx = -5;}

  runAction->Trk_X_pos = xSource + Z_Source * tan(Theta2) * cos(Phi2);
  runAction->Trk_Y_pos = ySource + Z_Source * tan(Theta2) * sin(Phi2);

  runAction->Trk_Theta = Theta2; 
  runAction->Trk_Phi =   Phi2;

  runAction->Pi0_X_pos = xSource + Z_Source * tan(Theta3) * cos(Phi3);
  runAction->Pi0_Y_pos = ySource + Z_Source * tan(Theta3) * sin(Phi3);

  runAction->Pi0_Theta = Theta3;
  runAction->Pi0_Phi =   Phi3;  
  
  auto particleDefinition1
    = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun_->SetParticleDefinition(particleDefinition1);
  fParticleGun_->SetParticleMomentumDirection(  G4ThreeVector( sin(Theta2) * cos(Phi2), sin(Theta2) * sin(Phi2), cos(Theta2) )  );
  //fParticleGun_->SetParticleMomentumDirection(  G4ThreeVector( 0, 0, 1 )  );
  fParticleGun_->SetParticleEnergy(ChEnergy);
  /*fParticleGun_
    ->SetParticlePosition(G4ThreeVector(xSource, ySource, -1 * (227.466 * cm/2 + 250 * cm) ));
  */  
  fParticleGun_
    ->SetParticlePosition(G4ThreeVector(xSource, ySource, -1 *  (227.466 * cm/2 + Z_Source) ) );
  fParticleGun_->GeneratePrimaryVertex(anEvent);
 

  fParticleGun_ = nullptr ;
  fParticleGun_ = new G4ParticleGun(1);
  auto particleDefinition2
    = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
  fParticleGun_->SetParticleDefinition(particleDefinition2);


  // G4double Theta = CLHEP::pi/2.;
  // G4double Phi = CLHEP::twopi*G4UniformRand() ;
  fParticleGun_->SetParticleMomentumDirection(G4ThreeVector( sin(Theta3) * cos(Phi3), sin(Theta3) * sin(Phi3), cos(Theta3)));
  //fParticleGun_->SetParticleMomentumDirection(G4ThreeVector( 0., 0., 1.));
  fParticleGun_->SetParticleEnergy(NuEnergy);
  /*fParticleGun_
    ->SetParticlePosition(G4ThreeVector(xSource, ySource, -1 *  (227.466 * cm/2 + 250 * cm) ) );
  */ 
  fParticleGun_
    ->SetParticlePosition(G4ThreeVector(xSource, ySource, -1 *  (227.466 * cm/2 + Z_Source) ) );
  fParticleGun_->GeneratePrimaryVertex(anEvent);

 /* auto runAction 
    = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
*/


  runAction->True_Ch_Energy = ChEnergy;
  runAction->True_Nu_Energy = NuEnergy;

  runAction->Smeared_Ch_Energy = gRandom->Gaus( ChEnergy/GeV, 5 * 0.0001 * (ChEnergy/GeV) * (ChEnergy/GeV)  ) * GeV; 

  //G4double NoiseValues[6] = {10, 20, 15, 20, 20, 20};
  //G4double NoiseValues[6] = {13, 34, 17, 54, 33, 54};
  G4double NoiseValues[6] = {13.,  34., 17., 14.,  8., 14.};

  for(G4int iLayer = 0; iLayer < kLayer; iLayer++){
    for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++){
      for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++){

        runAction->Noise_Cell_Energy[iLayer][iEtacell][iPhicell] = gRandom->Gaus(0., NoiseValues[iLayer] ) * MeV  ;


      }
    }
  } 

  }
  else
    G4Exception("H02PrimaryGeneratorAction::GeneratePrimaries",
                "InvalidSetup", FatalException,
                "Generator is not instanciated.");
}
