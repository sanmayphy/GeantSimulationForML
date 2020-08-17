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
/// \file CaloR/src/CaloRCellParameterisation.cc
/// \brief Implementation of the CaloRCellParameterisation class

#include "CaloRCellParameterisation.hh"
// #include "CaloRConstants.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4CutTubs.hh"
#include "G4SystemOfUnits.hh"

Em1CellParameterisation::Em1CellParameterisation() : G4VPVParameterisation() {
  // for (auto copyNo=0; copyNo<kNofEm1Cells; copyNo++) {
  for (auto copyNo=0; copyNo<kNofEm1Rows; copyNo++) {
    // auto column = copyNo / kNofEm1Rows ;
    auto column = copyNo;
    // auto row = copyNo % kNofEm1Rows;
    fXEm1Cell[copyNo] = (column+0.5)*ECAL1_dz - calorSizeXY/2;
    // G4cout << (column+0.5)*ECAL1_dz - calorSizeXY/2 << G4endl;
    // fXEm1Cell[copyNo] = (column+0.5)*ECAL1_dx - calorSizeXY/2;
    // fYEm1Cell[copyNo] = (row+0.5)*ECAL1_dy - calorSizeXY/2;
  }
}
Em1CellParameterisation::~Em1CellParameterisation() {}
void Em1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  // physVol->SetTranslation(G4ThreeVector(fXEm1Cell[copyNo],fYEm1Cell[copyNo],0.));}
  physVol->SetTranslation(G4ThreeVector(0.,0.,fXEm1Cell[copyNo]));}
 
  // physVol->SetTranslation(G4ThreeVector(0.,0.,100));}
// void Em1CellParameterisation::ComputeDimensions(G4CutTubs& Calorim,const G4int copyNo, G4VPhysicalVolume* physVol )const{
//   G4double tube_dPhi = 2.* M_PI * rad;  
//   G4double X0_ECAL = 3.897 * cm; //1./X0_ECAL_Inv ;
//   // Calorim = G4CutTubs("ECAL1cell", 50. * cm,50. * cm + 3 * X0_ECAL, ECAL1_dz, 0., tube_dPhi, G4ThreeVector(0, 0., -1 ), G4ThreeVector(0., 0.,1));

// }


Em2CellParameterisation::Em2CellParameterisation() : G4VPVParameterisation() {
  // for (auto copyNo=0; copyNo<kNofEm2Cells; copyNo++) {
  for (auto copyNo=0; copyNo<kNofEm2Rows; copyNo++) {
    // auto column = copyNo / kNofEm2Rows ;
    // auto row = copyNo % kNofEm2Rows;
    auto column = copyNo;
    fXEm2Cell[copyNo] = (column+0.5)*ECAL2_dz - calorSizeXY/2;
    // fYEm2Cell[copyNo] = (row+0.5)*ECAL2_dy - calorSizeXY/2;
  }
}
Em2CellParameterisation::~Em2CellParameterisation() {}
void Em2CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  // physVol->SetTranslation(G4ThreeVector(fXEm2Cell[copyNo],fYEm2Cell[copyNo],0.));}
  physVol->SetTranslation(G4ThreeVector(0.,0.,fXEm2Cell[copyNo]));}

Em3CellParameterisation::Em3CellParameterisation() : G4VPVParameterisation() {
  // for (auto copyNo=0; copyNo<kNofEm3Cells; copyNo++) {
  for (auto copyNo=0; copyNo<kNofEm3Rows; copyNo++) {
    // auto column = copyNo / kNofEm3Rows ;
    // auto row = copyNo % kNofEm3Rows;
    auto column = copyNo;
    fXEm3Cell[copyNo] = (column+0.5)*ECAL3_dz - calorSizeXY/2;
    // fYEm3Cell[copyNo] = (row+0.5)*ECAL3_dy - calorSizeXY/2;
  }
}
Em3CellParameterisation::~Em3CellParameterisation() {}
void Em3CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
  // physVol->SetTranslation(G4ThreeVector(fXEm3Cell[copyNo],fYEm3Cell[copyNo],0.));}
  physVol->SetTranslation(G4ThreeVector(0.,0.,fXEm3Cell[copyNo]));}

// Had1CellParameterisation::Had1CellParameterisation() : G4VPVParameterisation() {
//   for (auto copyNo=0; copyNo<kNofHad1Cells; copyNo++) {
//     auto column = copyNo / kNofHad1Rows ;
//     auto row = copyNo % kNofHad1Rows;
//     fXHad1Cell[copyNo] = (column+0.5)*HCAL1_dx - calorSizeXY/2;
//     fYHad1Cell[copyNo] = (row+0.5)*HCAL1_dy - calorSizeXY/2;
//   }
// }
// Had1CellParameterisation::~Had1CellParameterisation() {}
// void Had1CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
//   physVol->SetTranslation(G4ThreeVector(fXHad1Cell[copyNo],fYHad1Cell[copyNo],0.));}


// Had2CellParameterisation::Had2CellParameterisation() : G4VPVParameterisation() {
//   for (auto copyNo=0; copyNo<kNofHad2Cells; copyNo++) {
//     auto column = copyNo / kNofHad2Rows ;
//     auto row = copyNo % kNofHad2Rows;
//     fXHad2Cell[copyNo] = (column+0.5)*HCAL2_dx - calorSizeXY/2;
//     fYHad2Cell[copyNo] = (row+0.5)*HCAL2_dy - calorSizeXY/2;
//   }
// }
// Had2CellParameterisation::~Had2CellParameterisation() {}
// void Had2CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
//   physVol->SetTranslation(G4ThreeVector(fXHad2Cell[copyNo],fYHad2Cell[copyNo],0.));}

// Had3CellParameterisation::Had3CellParameterisation() : G4VPVParameterisation() {
//   for (auto copyNo=0; copyNo<kNofHad3Cells; copyNo++) {
//     auto column = copyNo / kNofHad3Rows ;
//     auto row = copyNo % kNofHad3Rows;
//     fXHad3Cell[copyNo] = (column+0.5)*HCAL3_dx - calorSizeXY/2;
//     fYHad3Cell[copyNo] = (row+0.5)*HCAL3_dy - calorSizeXY/2;
//   }
// }
// Had3CellParameterisation::~Had3CellParameterisation() {}
// void Had3CellParameterisation::ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const{
//   physVol->SetTranslation(G4ThreeVector(fXHad3Cell[copyNo],fYHad3Cell[copyNo],0.));}
