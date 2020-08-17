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
/// \file eventgenerator/HepMC/HepMCEx03/src/HepMCG4Pythia8Interface.cc
/// \brief Implementation of the HepMCG4PythiaInterface class for Pythia8
//

#ifdef G4LIB_USE_PYTHIA8
// ! Defenetly active when we use Pythia
#include "HepMCG4Pythia8Interface.hh"
#include "HepMCG4Pythia8Messenger.hh"

#include "HepMC/GenEvent.h"
#include "Randomize.hh"


using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Pythia8Interface::HepMCG4Pythia8Interface()
   : verbose(0)
{
  messenger= new HepMCG4Pythia8Messenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Pythia8Interface::~HepMCG4Pythia8Interface()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::CallPythiaReadString(G4String par)
{
   pythia.readString(par);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::CallPythiaInit(G4int beam,
                                        G4int target, G4double eCM)
{
   //pythia.init(beam,target, eCM);
   pythia.init();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::CallPythiaStat()
{
   pythia.stat();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Pythia8::Event HepMCG4Pythia8Interface::GetPythiaObject()
{
   // Pythia8::Event _Event = sum_events;//pythia.event;
   return sum_events;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::SetRandomSeed(G4int iseed)
{
   
   pythia.readString("Random:setSeed = on");
   ostringstream Seed;
   Seed<<"Random:seed = "<<iseed;
   pythia.readString(Seed.str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::PrintRandomStatus(std::ostream& ostr) const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::SetUserParameters()
{
  G4cout << "set user parameters of PYTHIA common." << G4endl
         << "nothing to be done in default."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMC::GenEvent* HepMCG4Pythia8Interface::GenerateHepMCEvent()
{  
   //? Pile-up
   // double nPileupAvg = 2.5;
   // int nPileup = 0;//poisson(nPileupAvg, pythia.rndm); 
   // G4cout << "Number of Pile-up events "<< nPileup <<G4endl;
   // HepMC::GenEvent* hepmcevt;
   // hepmcevt = new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM);
   // //signall event 
   // pythia.next();
   // sum_events = pythia.event;
   // ToHepMC.fill_next_event( pythia, hepmcevt,-1, false );
   // if(verbose>0)
   //    {
   //       hepmcevt-> print();
   //    } 
   // //pile-up
   // for (int iPileup = 0; iPileup < nPileup; ++iPileup) 
   // {
   //    pythia.next();
   //    sum_events += pythia.event;
   //    ToHepMC.fill_next_event( pythia, hepmcevt,-1, false );
   //    if(verbose>0)
   //    {
   //       hepmcevt-> print();
   //    } 


   // }
   // G4cout << " particles_size "<< hepmcevt-> particles_size()<< G4endl;
   // return hepmcevt;
   //?  Pile-up end
   //?  Quark and Gluon
   Pythia8::Event& event  = pythia.event;
   pythia.readString("ProcessLevel:all = off");
   pythia.readString("HardQCD:all = on");
   SetRandomSeed( CLHEP::RandFlat::shootInt(900000000));
   pythia.init();

   int id = messenger->GQparticle;
   double minenergy = messenger->MinEnergy;
   double maxenergy = messenger->MaxEnergy;
   double mineta    = messenger->MinEta;
   double maxeta    = messenger->MaxEta;


   double ee  = CLHEP::RandFlat::shoot(minenergy,maxenergy);//100;
   G4cout<<"Energy: "<<ee<<G4endl;
   double eta = CLHEP::RandFlat::shoot(mineta,maxeta);
   G4cout<<"Eta: "<<eta<<G4endl;
   double phi = CLHEP::RandFlat::shoot(0.,2*M_PI);
   Pythia8::ParticleData& pdt = pythia.particleData;
   event.reset();
   
   if (id==21)//? Gluon
   {
      G4cout<<"Gluon!"<<G4endl;
      double pt = ee/cosh(eta);
      double pz = pt*sinh(eta);
      double px = pt*cos(phi);
      double py = pt*sin(phi);
      event.append( id, 23, 101, 102,  px,  py,  pz, ee);
      event.append( id, 23, 102, 101, -px, -py, -pz, ee);
   }
   else if (id==2)//?  Quark
   {
      G4cout<<"Quark!"<<G4endl;
      double mm = pdt.m0(id);
      double pp = Pythia8::sqrtpos(ee*ee - mm*mm);
      double pt = pp/cosh(eta);
      double pz = pt*sinh(eta);
      double px = pt*cos(phi);
      double py = pt*sin(phi);
      event.append(  id, 23, 101,   0, px, py,  pz, ee, mm);
      event.append( -id, 23,   0, 101, -px, -py,  -pz, ee, mm);
   }
   else if (id==211)
   {
      double delta_r = CLHEP::RandFlat::shoot(M_PI/100, 3*M_PI/100 );
      double a = eta - delta_r;
      double b = eta + delta_r;
      double eta0 = CLHEP::RandFlat::shoot(a, b);
      int sign = 0;
      if (CLHEP::RandFlat::shoot(0.0, 1.0)>0.5) sign = 1;
      else sign = -1;
      double phi0 = sign*sqrt(sqr(delta_r)-sqr(eta - eta0))+phi;
      G4cout<<"PiPlus PiZero!"<<G4endl;
      double mm = pdt.m0(id);
      double pp = Pythia8::sqrtpos(ee*ee - mm*mm);
      double pt = pp/cosh(eta);
      double pz = pt*sinh(eta);
      double px = pt*cos(phi);
      double py = pt*sin(phi);
      event.append(  id, 23, 0,   0, px, py,  pz, ee, mm);
      double ee0 = CLHEP::RandFlat::shoot(minenergy,maxenergy);
      double mm0 = pdt.m0(111);
      double pp0 = Pythia8::sqrtpos(ee0*ee0 - mm0*mm0);
      double pt0 = pp0/cosh(eta0);
      double pz0 = pt0*sinh(eta0);
      double px0 = pt0*cos(phi0);
      double py0 = pt0*sin(phi0);
      event.append(  111, 23, 0,   0, px0, py0,  pz0, ee0, mm0);
   }
   double scale = ee;
   event[1].scale( scale);
   event[2].scale( scale);
   cout<<"AntonC! "<< event[2].canDecay()<<endl;
   cout<<"AntonC! "<< event[1].canDecay()<<endl;
   cout<<"AntonC! "<< event[2].mayDecay()<<endl;
   cout<<"AntonC! "<< event[1].mayDecay()<<endl;
   pythia.forceTimeShower( 1, 2, scale, 0);
   if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
    }
   HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM);
   ToHepMC.fill_next_event( pythia, hepmcevt);
   if(verbose>0){
      hepmcevt-> print();
   } 
   G4cout << " particles_size "<< hepmcevt-> particles_size()<< G4endl;
   sum_events = event;
   return hepmcevt;
   //? Quark end
   //? Usual code
   // pythia.next();

   // HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM);
   // ToHepMC.fill_next_event( pythia, hepmcevt );
   // if(verbose>0)
   // {
   //    hepmcevt-> print();
   // } 
   // sum_events = pythia.event;
   // G4cout << " particles_size "<< hepmcevt-> particles_size()<< G4endl;
   // return hepmcevt;
   //? Usual code end
}








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Pythia8Interface::Print() const
{
  G4cout << "Pythia8Interface::Print()..." << G4endl;
}



// !ToDo  for pile-up
int HepMCG4Pythia8Interface::poisson(double nAvg, Pythia8::Rndm& rndm) 
{

  // Set maximum to avoid overflow.
  const int NMAX = 100;

  // Random number.
  double rPoisson = rndm.flat() * exp(nAvg);

  // Initialize.
  double rSum  = 0.;
  double rTerm = 1.;

  // Add to sum and check whether done.
  for (int i = 0; i < NMAX; ) {
    rSum += rTerm;
    if (rSum > rPoisson) return i;

    // Evaluate next term. 
    ++i;
    rTerm *= nAvg / i;
  }
}

#endif
