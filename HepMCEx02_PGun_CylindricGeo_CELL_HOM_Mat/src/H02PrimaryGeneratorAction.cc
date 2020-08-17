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
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
H02PrimaryGeneratorAction::H02PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  // default generator is particle gun.
  fCurrentGenerator = fParticleGun = new G4ParticleGun();
  //fCurrentGenerator= fParticleGun= new G4GeneralParticleSource();
  fCurrentGeneratorName = "particleGun";
  fHepmcAscii = new HepMCG4AsciiReader();
  //#ifdef G4LIB_USE_PYTHIA
  fPythiaGen = new HepMCG4Pythia8Interface();
  //#else
  //  fPythiaGen= 0;
  //#endif

  fParticleGun_ = new G4ParticleGun(1);

  fGentypeMap["particleGun"] = fParticleGun;
  fGentypeMap["hepmcAscii"] = fHepmcAscii;
  fGentypeMap["pythia8"] = fPythiaGen;

  fMessenger = new H02PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02PrimaryGeneratorAction::~H02PrimaryGeneratorAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  auto runAction
  = (H02RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  if (fCurrentGenerator)
  {
    for (int i = 0; i < kMaxPixel; i++)
    {
      for (int j = 0; j < kMaxPixel; j++)
      {
        runAction->Track_Cell[i][j] = 0;
      }
    }
    //! particleGun
    if (fCurrentGeneratorName=="particleGun")
    {
      // G4cout<<"AAAA0"<<G4endl;

      G4double delta_r = CLHEP::pi/100 + (3*CLHEP::pi/100  - CLHEP::pi/100) * G4UniformRand();
      // G4double delta_r_sq = pow(delta_r,2);
      G4int sign = 0;
      G4double eta_fr = -1;
      G4double eta_to = 1;
      G4double eta0 = eta_fr + (eta_to - eta_fr) * G4UniformRand();
      G4double Theta0 = 2*atan(exp(-1*eta0));//1.27519765;//CLHEP::pi;//0;
      G4double phi_fr = 0;//CLHEP::pi;
      G4double phi_to = 2*CLHEP::pi;
      G4double Phi0 = phi_fr + (phi_to - phi_fr) * G4UniformRand();
      G4double a = eta0 - delta_r;
      G4double b = delta_r + eta0;
      G4double eta1 = eta0;//a + (b - a) * G4UniformRand();
      G4double eta2 = a + (b - a) * G4UniformRand();

      G4double Theta1 = 2*atan(exp(-1*eta1));
      G4double Theta2 = 2*atan(exp(-1*eta2));


      if(G4UniformRand() > 0.5) sign = 1;
      else                      sign = -1;
      G4double Phi1 = sign*sqrt(pow(delta_r,2)-pow((eta0-eta1),2))+Phi0;
      if(G4UniformRand() > 0.5) sign = 1;
      else                      sign = -1;
      G4double Phi2 = sign*sqrt(pow(delta_r,2)-pow((eta0-eta2),2))+Phi0;




      fParticleGun_ = nullptr ;
      fParticleGun_ = new G4ParticleGun(1);

      // double EnPi_P = (5 + 5 * G4UniformRand())  * GeV;
      // double EnPi_0 = (5 + 5 * G4UniformRand())  * GeV;
      double PtPi_P = (5 + 5 * G4UniformRand())  * GeV;
      double PtPi_0 = (5 + 5 * G4UniformRand())  * GeV;
      double pxPi_P_init = PtPi_P * cos(Phi0);
      double pyPi_P_init = PtPi_P * sin(Phi0);
      double pzPi_P_init = PtPi_P * sinh(eta1);

      double pxPi_0_init = PtPi_0 * cos(Phi1);
      double pyPi_0_init = PtPi_0 * sin(Phi1);
      double pzPi_0_init = PtPi_0 * sinh(eta1);



      fParticleGun_ = nullptr ;
      fParticleGun_ = new G4ParticleGun(1);
      auto particleDefinition1 = G4ParticleTable::GetParticleTable()->FindParticle("pi+");//"pi0"geantino#chargedgeantino
      fParticleGun_->SetParticleDefinition(particleDefinition1);
      fParticleGun_->SetParticleMomentum(G4ThreeVector(pxPi_P_init,pyPi_P_init,pzPi_P_init));
      // fParticleGun_->SetParticleMomentumDirection(
      // G4ThreeVector( sin(Theta1) * cos(Phi0), sin(Theta1) * sin(Phi0), cos(Theta1)));
      // fParticleGun_->SetParticleEnergy(EnPi_P);
      fParticleGun_->SetParticlePosition(G4ThreeVector(0., 0., 0. ) );
      fParticleGun_->GeneratePrimaryVertex(anEvent);


      fParticleGun_ = nullptr ;
      fParticleGun_ = new G4ParticleGun(1);
      auto particleDefinition2 = G4ParticleTable::GetParticleTable()->FindParticle("pi0");//"pi0"geantino#chargedgeantino
      fParticleGun_->SetParticleDefinition(particleDefinition2);
      fParticleGun_->SetParticleMomentum(G4ThreeVector(pxPi_0_init,pyPi_0_init,pzPi_0_init));
      // fParticleGun_->SetParticleMomentumDirection(
      // G4ThreeVector( sin(Theta2) * cos(Phi0), sin(Theta2) * sin(Phi0), cos(Theta2)));
      // fParticleGun_->SetParticleEnergy(EnPi_0);
      fParticleGun_->SetParticlePosition(G4ThreeVector(0., 0., 0. ) );
      fParticleGun_->GeneratePrimaryVertex(anEvent);



      


      double Hmagnetic = fieldValue;//3.8*tesla;//* tesla
      double mass   = particleDefinition1->GetPDGMass();//* Mev
      double charge = particleDefinition1->GetPDGCharge();//* 1 elementary charge
      runAction->ChargeParticlTraj.clear();
      // G4cout<<"AAAA1"<<G4endl;
      double px = pxPi_P_init;
      double py = pyPi_P_init;
      double pz = pzPi_P_init;
      double pt = sqrt(sqr(px)+sqr(py));
      double LighSpeed = 299792458 * m/s;//* mm/ ns
      double momentum = sqrt(sqr(px)+sqr(py)+sqr(pz));
      double energy = sqrt(sqr(momentum)+pow(mass,2))-mass;
      double ChEnergy = energy;
      // double NuEnergy = EnPi_0;
      runAction->True_Ch_Energy = ChEnergy+mass;
      // runAction->True_Nu_Energy = NuEnergy+particleDefinition2->GetPDGMass();
      // runAction->Seed = CLHEP::HepRandom::getTheSeed();
      if (charge!=0)
      {
        if (Hmagnetic!=0)
        {
            G4cout<<Hmagnetic<<" "<<charge<<G4endl;
            // double px = px_init;
            // double py = py_init;
            // double pz = pz_init;
            // double LighSpeed = 299792458 * m/s;//* mm/ ns
            // float momentum = sqrt(sqr(px)+sqr(py)+sqr(pz));
            // double energy = sqrt(sqr(momentum)+pow(mass,2))-mass; //EnPi_P;//Particle[4];//* It is kinetic energy
            // double mom = pow(pow((EnPi_P+mass),2)-pow(mass,2),0.5);
            // double px = mom*sin(Theta1) * cos(Phi0);//Particle[1];//* Mev
            // double py = mom*sin(Theta1) * sin(Phi0); //Particle[2];//* Mev
            // double pz = mom*cos(Theta1);//Particle[3];//* Mev
            double x0 = (py)/(charge*Hmagnetic*LighSpeed);//* mm
            double y0 = -1*(px)/(charge*Hmagnetic*LighSpeed);//* mm
            double alpha = sgn(charge*Hmagnetic) * r_inn_ECAL2/(2*pow(pow(x0,2)+pow(y0,2),0.5));
            if (abs(alpha)<1)
            {
              double px_f = (px*(1-2*sqr(alpha))+2*py*alpha*sqrt(1-sqr(alpha)));
              double py_f = (py*(1-2*sqr(alpha))-2*px*alpha*sqrt(1-sqr(alpha)));

              int IndexEta = -1;
              int IndexPhi = -1;

              vector <float> phi(6,0);
              vector <float> eta(6,0);
              vector <float> IndexPhi_v(6,-1);
              vector <float> IndexEta_v(6,-1);
              float R_inn[6] = {r_inn_ECAL1,r_inn_ECAL2,r_inn_ECAL3,r_inn_HCAL1,r_inn_HCAL2,r_inn_HCAL3};

              double x = 2*(x0*pow(alpha,2)-y0*alpha*pow((1-pow(alpha,2)),0.5));
              double y = 2*(y0*pow(alpha,2)+x0*alpha*pow((1-pow(alpha,2)),0.5));
              double z = (pz)/(abs(charge*Hmagnetic)*LighSpeed)*(2*asin(alpha));//*mm
              double r2 = pow(x,2)+pow(y,2);

              double a = sqr(pt/mass);
              double b = 2*(px_f*x+py_f*y)/mass; 


              for (int lay = 0; lay < 6; lay++)
              {
                  // G4cout<<"G2"<<G4endl;
                  float c = sqr(R_inn[lay])-r2;
                  float t = (-b+sqrt(4*a*c+sqr(b)))/(2*a);
                  double x_f = x+(px/mass)*t;
                  double y_f = y+(py/mass)*t;
                  double z_f = z+(pz/mass)*t;
                  double phi_ = atan2(y_f,x_f);
                  // G4cout<<"G2"<<G4endl;
                  if (phi_<0) phi_+=tube_dPhi;
                  phi[lay]=phi_;
                  double eta_ = -1*log(tan(0.5*acos(z_f/sqrt(sqr(R_inn[lay])+sqr(z_f)))));
                  // G4cout<<"G3"<<G4endl;
                  eta[lay]=eta_;
                  float deta = d_eta*(kMaxPixel/LayersPix[lay]);
                  float dphi = divided_tube_dPhi*(kMaxPixel/LayersPix[lay]);
                  IndexPhi = (float)floor(phi_/dphi);
                  IndexEta = (float)floor((eta_max+eta_)/deta);
                  if ((IndexPhi>LayersPix[lay]) || (IndexEta>LayersPix[lay]) ) 
                  {
                      IndexPhi=-1;
                      IndexEta=-1;
                  }
                  IndexPhi_v[lay] = IndexPhi;
                  IndexEta_v[lay] = IndexEta;
                  // G4cout<<"G4"<<G4endl;
              }
              //* [abseta,absphi,celleta,cellphi]
              if ((IndexEta!=-1 && IndexPhi!=-1) && sqrt(sqr(px)+sqr(py))<40e+03)
              {
                vector <double> TrajOfParticl;
                TrajOfParticl.push_back(particleDefinition1->GetPDGEncoding());
                TrajOfParticl.push_back(energy);
                TrajOfParticl.push_back(momentum);
                TrajOfParticl.push_back(pt);
                for (int lay = 0; lay<6;lay++)
                {
                    TrajOfParticl.push_back(eta[lay]);
                    TrajOfParticl.push_back(phi[lay]);
                    TrajOfParticl.push_back(IndexEta_v[lay]);
                    TrajOfParticl.push_back(IndexPhi_v[lay]);
                }
                TrajOfParticl.push_back(100000);
                TrajOfParticl.push_back(0);

                runAction->ChargeParticlTraj.push_back(TrajOfParticl);
                runAction->Track_Cell[IndexEta][IndexPhi] += energy;
              }
            }
            // G4cout<<"AAAA"<<G4endl;
        }
        else 
        {
          int IndexEta = -1;
          int IndexPhi = -1;
          
          for(int iEtacell = 0; iEtacell < kMaxPixel; iEtacell++)
          {
              if( ( eta1 >=  -1 * eta_max + d_eta * iEtacell ) &&
                  ( eta1 <  -1 * eta_max + d_eta * (iEtacell + 1)  ))
              {
                  IndexEta = iEtacell ;
                  break;
              }
          }
          for(int iPhicell = 0; iPhicell < kMaxPixel; iPhicell++)
          {
            if( ( Phi0 >=   divided_tube_dPhi * iPhicell  ) &&
                ( Phi0 <  divided_tube_dPhi * (iPhicell + 1)  ))
            {
                IndexPhi = iPhicell ;
                break;
            }
          }
          if ((IndexEta!=-1 && IndexPhi!=-1) && sqrt(sqr(px)+sqr(py))<40e+03)
            {
              vector <double> TrajOfParticl;
              vector <float> phi(6,Phi0);
              vector <float> eta(6,eta1);
              vector <float> IndexPhi_v(6,IndexPhi);
              vector <float> IndexEta_v(6,IndexEta);
              TrajOfParticl.push_back(particleDefinition1->GetPDGEncoding());
              TrajOfParticl.push_back(energy);
              TrajOfParticl.push_back(momentum);
              TrajOfParticl.push_back(pt);
              if (Phi0<0) Phi0 += tube_dPhi;
              // G4cout<<Phi0<<G4endl;
              
              for (int lay = 0; lay<6;lay++)
              {
                  TrajOfParticl.push_back(eta[lay]);
                  TrajOfParticl.push_back(phi[lay]);
                  TrajOfParticl.push_back(IndexEta_v[lay]);
                  TrajOfParticl.push_back(IndexPhi_v[lay]);
              }
              // TrajOfParticl.push_back(eta1);
              // TrajOfParticl.push_back(Phi0);
              // TrajOfParticl.push_back(IndexEta);
              // TrajOfParticl.push_back(IndexPhi);
              TrajOfParticl.push_back(100000);
              TrajOfParticl.push_back(0);
              

              runAction->ChargeParticlTraj.push_back(TrajOfParticl);
              runAction->Track_Cell[IndexEta][IndexPhi] += energy;
            }
          // runAction->Track_Cell[IndexEta][IndexPhi] += EnPi_P;
        }
      }
      



      G4double NoiseValues[6] = {13.,  34., 17., 54.,  33., 54.};

      for(G4int iLayer = 0; iLayer < kLayer; iLayer++)
      {
        for(G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++)
        {
          for(G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++)
          {
            runAction->Noise_Cell_Energy[iLayer][iEtacell][iPhicell] = gRandom->Gaus(0., NoiseValues[iLayer] ) * MeV  ;
          }
        }
      }
    }
    //! particleGun end
    //! pythia8
    else if (fCurrentGeneratorName=="pythia8")
    {
      G4cout << " !!!1!!! " << G4endl;
      fGentypeMap["pythia8"]->GeneratePrimaryVertex(anEvent);
      G4cout << " !!!2!!! " << G4endl;
      auto runAction = (H02RunAction *)G4RunManager::GetRunManager()->GetUserRunAction();
    
      runAction->Seed = CLHEP::HepRandom::getTheSeed();

      G4double NoiseValues[6] = {13., 34., 41., 75., 50., 25.};//{13., 34., 17., 54., 33., 54.}

      for (G4int iLayer = 0; iLayer < kLayer; iLayer++)
      {
        for (G4int iEtacell = 0; iEtacell < K_NETA; iEtacell++)
        {
          for (G4int iPhicell = 0; iPhicell < K_NPHI; iPhicell++)
          {

            runAction->Noise_Cell_Energy[iLayer][iEtacell][iPhicell] = gRandom->Gaus(0., NoiseValues[iLayer]) * MeV;
          }
        }
      }
    }
    //! Pythia8 end
  }
  else
    G4Exception("H02PrimaryGeneratorAction::GeneratePrimaries",
                "InvalidSetup", FatalException,
                "Generator is not instanciated.");
}
