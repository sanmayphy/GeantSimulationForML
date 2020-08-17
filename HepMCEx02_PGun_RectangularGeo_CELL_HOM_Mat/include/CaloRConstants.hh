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
/// \file HepMCEx02_PGun_NewRectangularGeo/include/CaloRConstants.hh
/// \brief Definition of detector constants used in Calo_RectangularGeo project.

#ifndef CaloRConstants_h
#define CaloRConstants_h 1

#include "G4SystemOfUnits.hh"

#include "globals.hh" // used for G4int, G4double

enum class CaloIdx { ECAL1, ECAL2, ECAL3, HCAL1, HCAL2, HCAL3 };

// Geometry parameters cell size of dR=0.1 in mm
constexpr G4double dR_01 = 125 * cm;

constexpr G4double Z_Source = 150 * cm;

constexpr G4int kNLayers = 6;

constexpr G4int kMaxPixel = 128;

constexpr G4int kNoEM1_Xcell(kMaxPixel), kNoEM1_Ycell(kMaxPixel) ;
constexpr G4int kNoEM2_Xcell(kMaxPixel), kNoEM2_Ycell(kMaxPixel) ;
constexpr G4int kNoEM3_Xcell(kMaxPixel), kNoEM3_Ycell(kMaxPixel) ;

constexpr G4int kNoHCAL1_Xcell(kMaxPixel), kNoHCAL1_Ycell(kMaxPixel) ;
constexpr G4int kNoHCAL2_Xcell(kMaxPixel), kNoHCAL2_Ycell(kMaxPixel) ;
constexpr G4int kNoHCAL3_Xcell(kMaxPixel), kNoHCAL3_Ycell(kMaxPixel) ;

constexpr G4int CELL_NUMBER[kNLayers] = {kNoEM1_Xcell, kNoEM2_Xcell, kNoEM3_Xcell, kNoHCAL1_Xcell, kNoHCAL2_Xcell, kNoHCAL3_Xcell}; 

//ECAL (use Pb-LAr mix with 1.5:2.1 ratio)
constexpr G4double  ECAL1_X0 = 4, ECAL2_X0 = 16, ECAL3_X0 = 3, X0_PbWO4 = 18.74 * mm, Lint_PbWO4 = 26.065 * cm;
constexpr G4double  ECAL1_dx = dR_01/kNoEM1_Xcell, ECAL1_dy = dR_01/kNoEM1_Ycell, ECAL1_dz = ECAL1_X0*X0_PbWO4;
constexpr G4double  ECAL2_dx = dR_01/kNoEM2_Xcell, ECAL2_dy = dR_01/kNoEM2_Ycell, ECAL2_dz = ECAL2_X0*X0_PbWO4;
constexpr G4double  ECAL3_dx = dR_01/kNoEM3_Xcell, ECAL3_dy = dR_01/kNoEM3_Ycell, ECAL3_dz = ECAL3_X0*X0_PbWO4;

//HCAL (use boxes of Fe+Scint with 4.7:1 ratio)
constexpr G4double  HCAL1_Lint = 1.5, HCAL2_Lint = 4.1, HCAL3_Lint = 1.8, Lint_Fe = 17.486 * cm;
constexpr G4double  HCAL1_dx = dR_01/kNoHCAL1_Xcell, HCAL1_dy = dR_01/kNoHCAL1_Ycell, HCAL1_dz = HCAL1_Lint*Lint_Fe;
constexpr G4double  HCAL2_dx = dR_01/kNoHCAL2_Xcell, HCAL2_dy = dR_01/kNoHCAL2_Ycell, HCAL2_dz = HCAL2_Lint*Lint_Fe;
constexpr G4double  HCAL3_dx = dR_01/kNoHCAL3_Xcell, HCAL3_dy = dR_01/kNoHCAL3_Ycell, HCAL3_dz = HCAL3_Lint*Lint_Fe;

// Detector geometry
constexpr G4double ECALThickness = ECAL1_dz + ECAL2_dz + ECAL3_dz;
constexpr G4double HCALThickness = HCAL1_dz + HCAL2_dz + HCAL3_dz;

constexpr G4double calorSizeXY = dR_01; // to cover dR  = 0.8 
constexpr G4double calorThickness = ECALThickness + HCALThickness;

// Calorimeter cell geometry
constexpr G4int kNofEm1Columns = (calorSizeXY/ECAL1_dx);
constexpr G4int kNofEm1Rows    = (calorSizeXY/ECAL1_dy);
constexpr G4int kNofEm1Cells   = kNofEm1Columns * kNofEm1Rows;
constexpr G4int kNofEm2Columns = (calorSizeXY/ECAL2_dx);
constexpr G4int kNofEm2Rows    = (calorSizeXY/ECAL2_dy);
constexpr G4int kNofEm2Cells   = kNofEm2Columns * kNofEm2Rows;
constexpr G4int kNofEm3Columns = (calorSizeXY/ECAL3_dx);
constexpr G4int kNofEm3Rows    = (calorSizeXY/ECAL3_dy);
constexpr G4int kNofEm3Cells   = kNofEm3Columns * kNofEm3Rows;

constexpr G4int kNofHad1Columns = (calorSizeXY/HCAL1_dx);
constexpr G4int kNofHad1Rows    = (calorSizeXY/HCAL1_dy);
constexpr G4int kNofHad1Cells   = kNofHad1Columns * kNofHad1Rows;
constexpr G4int kNofHad2Columns = (calorSizeXY/HCAL2_dx);
constexpr G4int kNofHad2Rows    = (calorSizeXY/HCAL2_dy);
constexpr G4int kNofHad2Cells   = kNofHad2Columns * kNofHad2Rows;
constexpr G4int kNofHad3Columns = (calorSizeXY/HCAL3_dx);
constexpr G4int kNofHad3Rows    = (calorSizeXY/HCAL3_dy);
constexpr G4int kNofHad3Cells   = kNofHad3Columns * kNofHad3Rows;

#endif
