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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02DetectorConstruction.cc
/// \brief Implementation of the H02DetectorConstruction class
//

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "H02Field.hh"
#include "H02DetectorConstruction.hh"
#include "H02MuonSD.hh"
#include "GlobalVariables.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "CaloRCellParameterisation.hh"
#include "CaloRConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// constants (detector parameters)
// [experimental hall]
// static const G4double R_EXPHALL= 5.*m;
// static const G4double DZ_EXPHALL= 10.*m;

// // [calorimeter]
// static const G4double RIN_BARREL_CAL= 2.*m;
// static const G4double ROUT_BARREL_CAL= 3.*m;
// static const G4double DZ_BARREL_CAL= 5.*m;

// static const G4double RIN_ENDCAP_CAL= 1.*m;
// static const G4double ROUT_ENDCAP_CAL= 3.*m;
// static const G4double DZ_ENDCAP_CAL= 0.5*m;

// // [muon system]
// static const G4double RIN_BARREL_MUON= 4.3*m;
// // static const G4double ROUT_BARREL_MUON= 4.5*m;
// static const G4double DX_BARREL_MUON= RIN_BARREL_MUON*std::cos(67.5*deg)-5.*cm;
// static const G4double DY_BARREL_MUON= 10.*cm;
// static const G4double DZ_BARREL_MUON= 7.*m;

// static const G4double RIN_ENDCAP_MUON=  1.*m;
// static const G4double ROUT_ENDCAP_MUON= 4.5*m;
// static const G4double DZ_ENDCAP_MUON= 10.*cm;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02DetectorConstruction::H02DetectorConstruction()
 : G4VUserDetectorConstruction() ,
 fAbsorberPV(nullptr),
 fGapPV(nullptr),
 fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H02DetectorConstruction::~H02DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* H02DetectorConstruction::Construct()
{
  // ==============================================================
  // Materials
  // ==============================================================

  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* lead = nistManager->FindOrBuildMaterial("G4_Pb");
  G4Material* iron = nistManager->FindOrBuildMaterial("G4_Fe");
  G4Material* plastic = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  

  // Argon gas
  G4double a, z, density, density_liq;
  //a= 39.95*g/mole;
  density= 1.782e-03*g/cm3;
  density_liq = 1.390*g/cm3;
  G4Material* ar= new G4Material("ArgonGas", z=18., a = 39.95*g/mole, density);
  G4Material* ar_liq= new G4Material("LiquidArgon", z = 18., a = 39.95*g/mole, density_liq);

  G4Material* defaultMaterial = new G4Material("Galactic", z=1., a=1.01*g/mole,density= CLHEP::universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);


  // G4Material *absorberMaterial_ECAL = lead;
  // G4Material *gapMaterial = ar_liq;

  // ====== Composite material construction ====== //
  G4double rho_Pb = 11.35 * g/cm3;
  G4double rho_LAr = 1.40 * g/cm3; 

  G4double rho_Fe = 7.874 * g/cm3;
  G4double rho_Sc = 1.032 * g/cm3;

  G4double ECAL_Pb_Prop(1.2), ECAL_LAr_Prop(4.5) ;
  G4double HCAL_Fe_Prop(4.7), HCAL_Sc_Prop(1.0) ;

  //G4double ECAL_Pb_Prop(0.), ECAL_LAr_Prop(8.2) ;
  //G4double HCAL_Fe_Prop(0.), HCAL_Sc_Prop(1.0) ;


  G4double rho_ECAL = (rho_Pb * ECAL_Pb_Prop + rho_LAr * ECAL_LAr_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) ;
  G4double rho_HCAL = (rho_Fe * HCAL_Fe_Prop + rho_Sc * HCAL_Sc_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) ;

  G4double fracMass;

  G4Material *Material_ECAL = new G4Material("Material_ECAL", rho_ECAL, 2);
  Material_ECAL->AddMaterial(lead, fracMass = (ECAL_Pb_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) * 100 * perCent  );
  Material_ECAL->AddMaterial(ar_liq, fracMass = (ECAL_LAr_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) * 100 * perCent  ); 

  G4Material *Material_HCAL = new G4Material("Material_HCAL", rho_HCAL, 2);
  Material_HCAL->AddMaterial(iron, fracMass = (HCAL_Fe_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) * 100 * perCent  );
  Material_HCAL->AddMaterial(plastic, fracMass = (HCAL_Sc_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) * 100 * perCent  );

  G4double X0_Pb = 0.56 * cm;
  G4double X0_LAr = 14.065 * cm;
  G4double Lambda_int_Fe = 16.99 * cm;
  G4double Lambda_int_Sc = 69.97 * cm;

  G4double X0_ECAL_Inv = 1./(ECAL_Pb_Prop + ECAL_LAr_Prop) * ( ECAL_Pb_Prop/X0_Pb + ECAL_LAr_Prop/X0_LAr ) ;
  G4double Lambda_int_HCAL_Inv = 1./(HCAL_Fe_Prop + HCAL_Sc_Prop) * ( HCAL_Fe_Prop/Lambda_int_Fe + HCAL_Sc_Prop/Lambda_int_Sc ) ;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // ====== Define the sizes ======== //
  G4int nofECAL_Layers = 3;
  G4int nofHCAL_Layers = 3;

  G4double X0_ECAL = 3.897 * cm; //1./X0_ECAL_Inv ;
  //G4double X0_ECAL = 3.897 * cm; 
  G4double Lambda_int_HCAL = 17.438 * cm; //1./Lambda_int_HCAL_Inv ;

  G4double Total_ECAL_Length = 3 * X0_ECAL + 16 * X0_ECAL + 6 * X0_ECAL ;
  G4double Total_HCAL_Length = 1.5 * Lambda_int_HCAL + 4.1 * Lambda_int_HCAL + 1.8 * Lambda_int_HCAL ;

  G4double Total_Calo_Length = Total_ECAL_Length + Total_HCAL_Length + 1.0 * cm; // -- the 1.0cm is the assumed gap between ECAL and HCAL
  
  G4double worldSizeXY = 1.2 * Total_Calo_Length;
  G4double worldSizeZ  = 1.5 * Total_Calo_Length;

  G4double Total_Calo_Length_XY = calorSizeXY ;


  G4cout << "Combined Radiation length : " << X0_ECAL << G4endl;
  G4cout << "Combined Interaction length : " << Lambda_int_HCAL << G4endl;
  G4cout << "Total Calo Length : " << Total_Calo_Length << G4endl;
  // ==============================================================
  // Experimental Hall (world)
  // ==============================================================
  G4Box* expHallSolid=
    new G4Box("EXP_HALL", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2 );

  G4LogicalVolume* expHallLV=
    new G4LogicalVolume(expHallSolid, defaultMaterial, "EXP_HALL_LV");

  // visualization attributes
  G4VisAttributes* expHallVisAtt=
    new G4VisAttributes(true, G4Colour(1., 1., 1.));
  //expHallVisAtt-> SetForceWireframe(TRUE);
  expHallLV-> SetVisAttributes(expHallVisAtt);

  G4PVPlacement* expHall= new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 expHallLV,          // its logical volume                         
                 "Exp_HALL",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );
  //            ...                                    MV, MANY, copy#

  // =========== ECAL_LAYER-1 ============== //
  G4Box *ECAL1_Solid = new G4Box( "ECAL_Absorber_LAYER-1_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 3 * X0_ECAL /2 );
  G4LogicalVolume *ECAL1_LV = new G4LogicalVolume(ECAL1_Solid, Material_ECAL, "ECAL_LAYER-1_LV");

  G4double zpos_ECAL1 = -1 * Total_Calo_Length/2 + 3 * X0_ECAL/2;
  G4PVPlacement *ECAL1 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_ECAL1 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 ECAL1_LV,          // its logical volume                         
                 "ECAL1",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* ECAL1_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 255, 0)); 
  ECAL1_LV-> SetVisAttributes(ECAL1_VisAtt);

  // ECAL Cells
  //
  auto ECAL1cellSolid = new G4Box("ECAL1cell",ECAL1_dx/2,ECAL1_dy/2,3 * X0_ECAL/2);
  auto ECAL1cellLV = new G4LogicalVolume(ECAL1cellSolid,Material_ECAL,"ECAL1cellLogical");
  new G4PVParameterised("ECAL1cellPhysical",ECAL1cellLV,ECAL1_LV,
                        kXAxis,kNofEm1Cells,new Em1CellParameterisation());
  ECAL1cellLV-> SetVisAttributes(ECAL1_VisAtt);

  G4cout <<"Layer1_zI : " << (zpos_ECAL1 - 3 * X0_ECAL/2) << ", Layer1_zF : " 
         << (zpos_ECAL1 + 3 * X0_ECAL/2) 
         << "Layer1_Width : " << 3 * X0_ECAL
         << G4endl;

  // // ---------------------------------------------------------- //
  // G4Box *ECAL1_Scintillator_Solid = new G4Box( "ECAL_Scintillator_LAYER-1_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 3 * 14.0 * cm/2 );
  // G4LogicalVolume *ECAL1_Scintillator_LV = new G4LogicalVolume(ECAL1_Scintillator_Solid, ar_liq, "ECAL_Scintillator_LAYER-1_LV");

  // G4double zpos_ECAL1_Scintillator_Solid = zpos_ECAL1_Absorber + 3 * 0.5 * cm/2 + 3 * 14.0 * cm/2;
  // G4PVPlacement *ECAL1_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_ECAL1_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                ECAL1_Scintillator_LV,          // its logical volume                         
  //                "ECAL1_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps
  //   );

  // G4VisAttributes* ECAL1_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(128, 0, 0));
  // ECAL1_Scintillator_LV-> SetVisAttributes(ECAL1_Scintillator_VisAtt);  

  // // Cells 
  // auto ECAL1cellSolid_Scn = new G4Box("ECAL1cell_Scn",ECAL1_dx/2,ECAL1_dy/2,3 * 14.0 * cm/2);
  // auto ECAL1cellLV_Scn = new G4LogicalVolume(ECAL1cellSolid_Scn,ar_liq,"ECAL1cellLogical_Scn");
  // new G4PVParameterised("ECAL1cellPhysical_Scn",ECAL1cellLV_Scn,ECAL1_Scintillator_LV,
  //                       kXAxis,kNofEm1Cells,new Em1CellParameterisation());
  // ECAL1cellLV_Scn-> SetVisAttributes(ECAL1_Scintillator_VisAtt);

  // G4cout <<"Layer1_zI : " << (-1 * Total_Calo_Length/2 - 3 * 0.5 * cm/2) << ", Layer1_zF : " 
  //        << (zpos_ECAL1_Scintillator_Solid + 3 * 14.0 * cm/2) << G4endl;

  // =========== ECAL_LAYER-2 ============== //
  G4Box *ECAL2_Solid = new G4Box( "ECAL_LAYER-2_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 16 * X0_ECAL/2 );
  G4LogicalVolume *ECAL2_LV = new G4LogicalVolume(ECAL2_Solid, Material_ECAL, "ECAL_LAYER-2_LV");

  G4double zpos_ECAL2 = zpos_ECAL1 + 3 * X0_ECAL/2 + 16  * X0_ECAL/2;
  G4PVPlacement *ECAL2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_ECAL2 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 ECAL2_LV,          // its logical volume                         
                 "ECAL2_Absorber",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* ECAL2_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 255, 0)); 
  ECAL2_LV-> SetVisAttributes(ECAL2_VisAtt);

  // ECAL Cells
  //
  auto ECAL2cellSolid = new G4Box("ECAL2cell",ECAL2_dx/2,ECAL2_dy/2,16 * X0_ECAL/2);
  auto ECAL2cellLV = new G4LogicalVolume(ECAL2cellSolid,Material_ECAL,"ECAL2cellLogical");
  new G4PVParameterised("ECAL2cellPhysical",ECAL2cellLV,ECAL2_LV,
                        kXAxis,kNofEm2Cells,new Em2CellParameterisation());
  ECAL2cellLV-> SetVisAttributes(ECAL2_VisAtt);

  G4cout <<"Layer2_zI : " << (zpos_ECAL2 - 16 * X0_ECAL/2) << ", Layer2_zF : " 
         << (zpos_ECAL2 + 16 * X0_ECAL/2) 
         << "Layer2_Width : " << 16 * X0_ECAL
         << G4endl;
  // // ---------------------------------------------------------- //
  // G4Box *ECAL2_Scintillator_Solid = new G4Box( "ECAL_Scintillator_LAYER-2_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 16 * 14.0 * cm/2 );
  // G4LogicalVolume *ECAL2_Scintillator_LV = new G4LogicalVolume(ECAL2_Scintillator_Solid, ar_liq, "ECAL_Scintillator_LAYER-2_LV");

  // G4double zpos_ECAL2_Scintillator_Solid = zpos_ECAL2_Absorber + 16 * 0.5 * cm/2 + 16 * 14.0 * cm/2;
  // G4PVPlacement *ECAL2_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_ECAL2_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                ECAL2_Scintillator_LV,          // its logical volume                         
  //                "ECAL2_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps
  //   );

  // G4VisAttributes* ECAL2_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(128, 0, 0));
  // ECAL2_Scintillator_LV-> SetVisAttributes(ECAL2_Scintillator_VisAtt);

  // // Cells 
  // auto ECAL2cellSolid_Scn = new G4Box("ECAL2cell_Scn",ECAL2_dx/2,ECAL2_dy/2,16 * 14.0 * cm/2);
  // auto ECAL2cellLV_Scn = new G4LogicalVolume(ECAL2cellSolid_Scn,ar_liq,"ECAL2cellLogical_Scn");
  // new G4PVParameterised("ECAL2cellPhysical_Scn",ECAL2cellLV_Scn,ECAL2_Scintillator_LV,
  //                       kXAxis,kNofEm2Cells,new Em2CellParameterisation());
  // ECAL2cellLV_Scn-> SetVisAttributes(ECAL2_Scintillator_VisAtt);

  // G4cout <<"Layer2_zI : " << (zpos_ECAL2_Absorber - 16 * 0.5 * cm/2) << ", Layer2_zF : " 
  //        << (zpos_ECAL2_Scintillator_Solid + 16 * 14.0 * cm/2) << G4endl;

  // =========== ECAL_LAYER-3 ============== //
  G4Box *ECAL3_Solid = new G4Box( "ECAL_LAYER-3_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 6 * X0_ECAL/2 );
  G4LogicalVolume *ECAL3_LV = new G4LogicalVolume(ECAL3_Solid, Material_ECAL, "ECAL_LAYER-3_LV");

  G4double zpos_ECAL3 = zpos_ECAL2 + 16 * X0_ECAL/2 + 6 * X0_ECAL/2;
  G4PVPlacement *ECAL3 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_ECAL3 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 ECAL3_LV,          // its logical volume                         
                 "ECAL3_Absorber",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* ECAL3_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 255, 0)); 
  ECAL3_LV-> SetVisAttributes(ECAL3_VisAtt);

  // ECAL Cells
  //
  auto ECAL3cellSolid = new G4Box("ECAL3cell",ECAL3_dx/2,ECAL3_dy/2,6 * X0_ECAL/2);
  auto ECAL3cellLV = new G4LogicalVolume(ECAL3cellSolid,Material_ECAL,"ECAL3cellLogical");
  new G4PVParameterised("ECAL3cellPhysical",ECAL3cellLV,ECAL3_LV,
                        kXAxis,kNofEm3Cells,new Em3CellParameterisation());
  ECAL3cellLV-> SetVisAttributes(ECAL3_VisAtt);

  G4cout <<"Layer3_zI : " << (zpos_ECAL3 - 6 * X0_ECAL/2) << ", Layer3_zF : " 
         << (zpos_ECAL3 + 6 * X0_ECAL/2) 
         << "Layer3_Width : " << 6 * X0_ECAL
         << G4endl;
  // // ---------------------------------------------------------- //
  // G4Box *ECAL3_Scintillator_Solid = new G4Box( "ECAL_Scintillator_LAYER-3_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 6 * 14.0 * cm/2 );
  // G4LogicalVolume *ECAL3_Scintillator_LV = new G4LogicalVolume(ECAL3_Scintillator_Solid, ar_liq, "ECAL_Scintillator_LAYER-3_LV");

  // G4double zpos_ECAL3_Scintillator_Solid = zpos_ECAL3_Absorber + 6 * 0.5 * cm/2 + 6 * 14.0 * cm/2;
  // G4PVPlacement *ECAL3_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_ECAL3_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                ECAL3_Scintillator_LV,          // its logical volume                         
  //                "ECAL3_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number         
  //               fCheckOverlaps
  //   );

  // G4VisAttributes* ECAL3_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(128, 0, 0));
  // ECAL3_Scintillator_LV-> SetVisAttributes(ECAL3_Scintillator_VisAtt);

  // // Cells 
  // auto ECAL3cellSolid_Scn = new G4Box("ECAL3cell_Scn",ECAL3_dx/2,ECAL3_dy/2,6 * 14.0 * cm/2);
  // auto ECAL3cellLV_Scn = new G4LogicalVolume(ECAL3cellSolid_Scn,ar_liq,"ECAL2cellLogical_Scn");
  // new G4PVParameterised("ECAL3cellPhysical_Scn",ECAL3cellLV_Scn,ECAL3_Scintillator_LV,
  //                       kXAxis,kNofEm3Cells,new Em3CellParameterisation());
  // ECAL3cellLV_Scn-> SetVisAttributes(ECAL3_Scintillator_VisAtt);

  // G4cout <<"Layer3_zI : " << (zpos_ECAL3_Absorber - 6 * 0.5 * cm/2) << ", Layer3_zF : " 
  //        << (zpos_ECAL3_Scintillator_Solid + 6 * 14.0 * cm/2) << G4endl;

  // =========== ECAL - HCAL GAP ============== //
  G4Box *GAP_Solid = new G4Box( "GAP_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 1 * cm/2 );
  G4LogicalVolume *GAP_LV = new G4LogicalVolume(GAP_Solid, defaultMaterial, "GAP_LV");

  G4double zpos_GAP = zpos_ECAL3 + 6 * X0_ECAL/2 + 1 * cm/2;
  G4PVPlacement *GAP = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_GAP ),  // at (0,0,-1 * Total_Calo_Length/2)
                 GAP_LV,          // its logical volume                         
                 "GAP",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* GAP_VisAtt=
    new G4VisAttributes(true, G4Colour(255, 0, 255)); 
  GAP_LV-> SetVisAttributes(GAP_VisAtt);

  // =========== HCAL_LAYER-1 ============== //
  G4Box *HCAL1_Solid = new G4Box( "HCAL_LAYER-1_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 1.5  * Lambda_int_HCAL/2 );
  G4LogicalVolume *HCAL1_LV = new G4LogicalVolume(HCAL1_Solid, iron, "HCAL_Absorber_LAYER-1_LV");

  G4double zpos_HCAL1 = zpos_GAP + 1 * cm/2 + 1.5 * Lambda_int_HCAL/2;
  G4PVPlacement *HCAL1 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_HCAL1 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 HCAL1_LV,          // its logical volume                         
                 "HCAL1_Absorber",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* HCAL1_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 0, 255, 0.8)); 
  HCAL1_LV-> SetVisAttributes(HCAL1_VisAtt);

  // HCAL Cells
  //
  auto HCAL1cellSolid = new G4Box("HCAL1cell",HCAL1_dx/2,HCAL1_dy/2,1.5 * Lambda_int_HCAL/2);
  auto HCAL1cellLV = new G4LogicalVolume(HCAL1cellSolid,Material_HCAL,"HCAL1cellLogical");
  new G4PVParameterised("HCAL1cellPhysical",HCAL1cellLV,HCAL1_LV,
                        kXAxis,kNofHad1Cells,new Had1CellParameterisation());
  HCAL1cellLV-> SetVisAttributes(HCAL1_VisAtt);

  G4cout <<"Layer4_zI : " << (zpos_HCAL1 - 1.5  * Lambda_int_HCAL/2) << ", Layer4_zF : " 
         << (zpos_HCAL1 + 1.5  * Lambda_int_HCAL/2) 
         << "Layer4_Width : " << 1.5  * Lambda_int_HCAL
         << G4endl;
  // // ---------------------------------------------------------- //
  // G4Box *HCAL1_Scintillator_Solid = new G4Box( "HCAL_Scintillator_LAYER-1_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 1.5 * 79.4 * cm/2 );
  // G4LogicalVolume *HCAL1_Scintillator_LV = new G4LogicalVolume(HCAL1_Scintillator_Solid, plastic, "HCAL_Scintillator_LAYER-1_LV");

  // G4double zpos_HCAL1_Scintillator_Solid = zpos_HCAL1_Absorber + 1.5 * 16.8 * cm/2 + 1.5 * 79.4 * cm/2;
  // G4PVPlacement *HCAL1_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_HCAL1_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                HCAL1_Scintillator_LV,          // its logical volume                         
  //                "HCAL1_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps
  //   );

  // G4VisAttributes* HCAL1_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(255, 255, 0));
  // HCAL1_Scintillator_LV-> SetVisAttributes(HCAL1_Scintillator_VisAtt);

  // // Cells 
  // auto HCAL1cellSolid_Scn = new G4Box("HCAL1cell_Scn",HCAL1_dx/2,HCAL1_dy/2, 1.5 * 79.4 * cm/2);
  // auto HCAL1cellLV_Scn = new G4LogicalVolume(HCAL1cellSolid_Scn,plastic,"HCAL1cellLogical_Scn");
  // new G4PVParameterised("HCAL1cellPhysical_Scn",HCAL1cellLV_Scn,HCAL1_Scintillator_LV,
  //                       kXAxis,kNofHad1Cells,new Had1CellParameterisation());
  // HCAL1cellLV_Scn-> SetVisAttributes(HCAL1_Scintillator_VisAtt);

  // G4cout <<"Layer4_zI : " << (zpos_HCAL1_Absorber - 1.5 * 16.8 * cm/2) << ", Layer4_zF : " 
  //        << (zpos_HCAL1_Scintillator_Solid + 1.5 * 79.4 * cm/2) << G4endl;


  // =========== HCAL_LAYER-2 ============== //
  G4Box *HCAL2_Solid = new G4Box( "HCAL_LAYER-2_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 4.1 * Lambda_int_HCAL/2 );
  G4LogicalVolume *HCAL2_LV = new G4LogicalVolume(HCAL2_Solid, iron, "HCAL_LAYER-2_LV");

  G4double zpos_HCAL2 = zpos_HCAL1 + 1.5 * Lambda_int_HCAL/2 + 4.1 * Lambda_int_HCAL/2;
  G4PVPlacement *HCAL2 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_HCAL2 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 HCAL2_LV,          // its logical volume                         
                 "HCAL2",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* HCAL2_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 0, 255, 0.8)); 
  HCAL2_LV-> SetVisAttributes(HCAL2_VisAtt); 

  // HCAL Cells
  //
  auto HCAL2cellSolid = new G4Box("HCAL2cell",HCAL2_dx/2,HCAL2_dy/2, 4.1 * Lambda_int_HCAL/2);
  auto HCAL2cellLV = new G4LogicalVolume(HCAL2cellSolid,Material_HCAL,"HCAL2cellLogical");
  new G4PVParameterised("HCAL2cellPhysical",HCAL2cellLV,HCAL2_LV,
                        kXAxis,kNofHad2Cells,new Had2CellParameterisation());
  HCAL2cellLV-> SetVisAttributes(HCAL2_VisAtt); 

  G4cout <<"Layer5_zI : " << (zpos_HCAL2 - 4.1  * Lambda_int_HCAL/2) << ", Layer5_zF : " 
         << (zpos_HCAL2 + 4.1  * Lambda_int_HCAL/2) 
         << "Layer5_Width : " << 4.1  * Lambda_int_HCAL
         << G4endl;
  // // ---------------------------------------------------------- //
  // G4Box *HCAL2_Scintillator_Solid = new G4Box( "HCAL_Scintillator_LAYER-2_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 4.1 * 79.4 * cm/2 );
  // G4LogicalVolume *HCAL2_Scintillator_LV = new G4LogicalVolume(HCAL2_Scintillator_Solid, plastic, "HCAL_Scintillator_LAYER-2_LV");

  // G4double zpos_HCAL2_Scintillator_Solid = zpos_HCAL2_Absorber + 4.1 * 16.8 * cm/2 + 4.1 * 79.4 * cm/2;
  // G4PVPlacement *HCAL2_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_HCAL2_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                HCAL2_Scintillator_LV,          // its logical volume                         
  //                "HCAL2_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps
  //   );

  // G4VisAttributes* HCAL2_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(255, 255, 0));
  // HCAL2_Scintillator_LV-> SetVisAttributes(HCAL2_Scintillator_VisAtt);

  // // Cells 
  // auto HCAL2cellSolid_Scn = new G4Box("HCAL2cell_Scn",HCAL2_dx/2,HCAL2_dy/2, 4.1 * 79.4 * cm/2);
  // auto HCAL2cellLV_Scn = new G4LogicalVolume(HCAL2cellSolid_Scn,plastic,"HCAL2cellLogical_Scn");
  // new G4PVParameterised("HCAL2cellPhysical_Scn",HCAL2cellLV_Scn,HCAL2_Scintillator_LV,
  //                       kXAxis,kNofHad2Cells,new Had2CellParameterisation());
  // HCAL2cellLV_Scn-> SetVisAttributes(HCAL2_Scintillator_VisAtt);

  // G4cout <<"Layer5_zI : " << (zpos_HCAL2_Absorber - 4.1 * 16.8 * cm/2) << ", Layer5_zF : " 
  //        << (zpos_HCAL2_Scintillator_Solid + 4.1 * 79.4 * cm/2) << G4endl;

  // =========== HCAL_LAYER-3 ============== //
  G4Box *HCAL3_Solid = new G4Box( "HCAL_LAYER-3_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 1.8 * Lambda_int_HCAL/2 );
  G4LogicalVolume *HCAL3_LV = new G4LogicalVolume(HCAL3_Solid, Material_HCAL, "HCAL_LAYER-3_LV");

  G4double zpos_HCAL3 = zpos_HCAL2 + 4.1 * Lambda_int_HCAL/2 + 1.8 * Lambda_int_HCAL/2;
  G4PVPlacement *HCAL3 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., zpos_HCAL3 ),  // at (0,0,-1 * Total_Calo_Length/2)
                 HCAL3_LV,          // its logical volume                         
                 "HCAL3_Absorber",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps
    );

  G4VisAttributes* HCAL3_VisAtt=
    new G4VisAttributes(true, G4Colour(0, 0, 255, 0.8)); 
  HCAL3_LV-> SetVisAttributes(HCAL3_VisAtt);

  // HCAL Cells
  //
  auto HCAL3cellSolid = new G4Box("HCAL3cell",HCAL3_dx/2,HCAL3_dy/2, 1.8 * Lambda_int_HCAL/2);
  auto HCAL3cellLV = new G4LogicalVolume(HCAL3cellSolid,Material_HCAL,"HCAL3cellLogical");
  new G4PVParameterised("HCAL3cellPhysical",HCAL3cellLV,HCAL3_LV,
                        kXAxis,kNofHad3Cells,new Had3CellParameterisation());
  HCAL3cellLV-> SetVisAttributes(HCAL3_VisAtt);

  G4cout <<"Layer6_zI : " << (zpos_HCAL3 - 1.8  * Lambda_int_HCAL/2) << ", Layer6_zF : " 
         << (zpos_HCAL3 + 1.8  * Lambda_int_HCAL/2) 
         << "Layer6_Width : " << 1.8  * Lambda_int_HCAL
         << G4endl;
  // // ---------------------------------------------------------- //
  // G4Box *HCAL3_Scintillator_Solid = new G4Box( "HCAL_Scintillator_LAYER-3_Solid", Total_Calo_Length_XY/2, Total_Calo_Length_XY/2, 1.8 * 79.4 * cm/2 );
  // G4LogicalVolume *HCAL3_Scintillator_LV = new G4LogicalVolume(HCAL3_Scintillator_Solid, plastic, "HCAL_Scintillator_LAYER-3_LV");

  // G4double zpos_HCAL3_Scintillator_Solid = zpos_HCAL3_Absorber + 1.8 * 16.8 * cm/2 + 1.8 * 79.4 * cm/2;
  // G4PVPlacement *HCAL3_Scintillator = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., zpos_HCAL3_Scintillator_Solid ),  // at (0,0,-1 * Total_Calo_Length/2)
  //                HCAL3_Scintillator_LV,          // its logical volume                         
  //                "HCAL3_Scintillator",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //                fCheckOverlaps
  //   );

  // G4VisAttributes* HCAL3_Scintillator_VisAtt=
  //   new G4VisAttributes(true, G4Colour(255, 255, 0));
  // HCAL3_Scintillator_LV-> SetVisAttributes(HCAL3_Scintillator_VisAtt);

  // // Cells 
  // auto HCAL3cellSolid_Scn = new G4Box("HCAL3cell_Scn",HCAL3_dx/2,HCAL3_dy/2, 1.8 * 79.4 * cm/2);
  // auto HCAL3cellLV_Scn = new G4LogicalVolume(HCAL3cellSolid_Scn,plastic,"HCAL3cellLogical_Scn");
  // new G4PVParameterised("HCAL3cellPhysical_Scn",HCAL3cellLV_Scn,HCAL3_Scintillator_LV,
  //                       kXAxis,kNofHad3Cells,new Had3CellParameterisation());
  // HCAL3cellLV_Scn-> SetVisAttributes(HCAL3_Scintillator_VisAtt);

  // G4cout <<"Layer6_zI : " << (zpos_HCAL3_Absorber - 1.8 * 16.8 * cm/2) << ", Layer6_zF : " 
  //        << (zpos_HCAL3_Scintillator_Solid + 1.8 * 79.4 * cm/2) << G4endl;


  return expHall;
}
