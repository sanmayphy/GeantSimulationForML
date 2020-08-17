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
/// \file H02DetectorConstruction.cc
/// \brief Implementation of the H02DetectorConstruction class
#include "H02DetectorConstruction.hh"

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
#include "G4CutTubs.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVParameterised.hh"
#include "CaloRCellParameterisation.hh"
#include "CaloRConstants.hh"
#include "G4ReflectionFactory.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "H02RunAction.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4ThreadLocal
// G4GlobalMagFieldMessenger* H02DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

H02DetectorConstruction::H02DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fCheckOverlaps(false)
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
  long double a, z, density, density_liq;
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
  long double rho_Pb = 11.35 * g/cm3;
  long double rho_LAr = 1.40 * g/cm3;

  long double rho_Fe = 7.874 * g/cm3;
  long double rho_Sc = 1.032 * g/cm3;

  long double ECAL_Pb_Prop(1.2), ECAL_LAr_Prop(4.5) ;
  long double HCAL_Fe_Prop(4.7), HCAL_Sc_Prop(1.0) ;

  //long double ECAL_Pb_Prop(0.), ECAL_LAr_Prop(8.2) ;
  //long double HCAL_Fe_Prop(0.), HCAL_Sc_Prop(1.0) ;


  long double rho_ECAL = (rho_Pb * ECAL_Pb_Prop + rho_LAr * ECAL_LAr_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) ;
  long double rho_HCAL = (rho_Fe * HCAL_Fe_Prop + rho_Sc * HCAL_Sc_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) ;

  long double fracMass;

  G4Material *Material_ECAL = new G4Material("Material_ECAL", rho_ECAL, 2);
  Material_ECAL->AddMaterial(lead, fracMass = (ECAL_Pb_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) * 100 * perCent  );
  Material_ECAL->AddMaterial(ar_liq, fracMass = (ECAL_LAr_Prop)/(ECAL_Pb_Prop + ECAL_LAr_Prop) * 100 * perCent  );

  G4Material *Material_HCAL = new G4Material("Material_HCAL", rho_HCAL, 2);
  Material_HCAL->AddMaterial(iron, fracMass = (HCAL_Fe_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) * 100 * perCent  );
  Material_HCAL->AddMaterial(plastic, fracMass = (HCAL_Sc_Prop)/(HCAL_Fe_Prop + HCAL_Sc_Prop) * 100 * perCent  );

  long double X0_Pb = 0.56 * cm;
  long double X0_LAr = 14.065 * cm;
  long double Lambda_int_Fe = 16.99 * cm;
  long double Lambda_int_Sc = 69.97 * cm;

  long double X0_ECAL_Inv = 1./(ECAL_Pb_Prop + ECAL_LAr_Prop) * ( ECAL_Pb_Prop/X0_Pb + ECAL_LAr_Prop/X0_LAr ) ;
  long double Lambda_int_HCAL_Inv = 1./(HCAL_Fe_Prop + HCAL_Sc_Prop) * ( HCAL_Fe_Prop/Lambda_int_Fe + HCAL_Sc_Prop/Lambda_int_Sc ) ;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // ====== Define the sizes ======== //
  G4int nofECAL_Layers = 3;
  G4int nofHCAL_Layers = 3;

  long double X0_ECAL = 3.897 * cm; //1./X0_ECAL_Inv ;
  //long double X0_ECAL = 3.897 * cm;
  long double Lambda_int_HCAL = 17.438 * cm; //1./Lambda_int_HCAL_Inv ;

  long double Total_ECAL_Length = 3 * X0_ECAL + 16 * X0_ECAL + 6 * X0_ECAL ;
  long double Total_HCAL_Length = 1.5 * Lambda_int_HCAL + 4.1 * Lambda_int_HCAL + 1.8 * Lambda_int_HCAL ;

  long double Total_Calo_Length = Total_ECAL_Length + Total_HCAL_Length + 1.0 * cm; // -- the 1.0cm is the assumed gap between ECAL and HCAL

  long double theta_min = 2*atan(exp(-1*eta_max));
  long double d_theta;
  // long double half_lengthZ = r_inn*(1.0/(tan(theta_min)));
  long double l =  r_inn_HCAL3/tan(theta_min);
  long double vx,vy,vz,vx_o,vz_o;
  long double worldSizeXY = 4. * r_out_HCAL3;//#1.2 * Total_Calo_Length;
  long double worldSizeZ  = 4. * l; //1.5 * Total_Calo_Length;

  long double Total_Calo_Length_XY = calorSizeXY ;



  G4cout << "Combined Radiation length : " << X0_ECAL << G4endl;
  G4cout << "Combined Interaction length : " << Lambda_int_HCAL << G4endl;
  G4cout << "Total Calo Length : " << Total_Calo_Length << G4endl;
  G4cout << "Tan : " << tan(theta_min) << G4endl;
  G4cout << "theta_min : " << theta_min << G4endl;
  // G4cout << "Length : " << half_lengthZ << G4endl;
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
                 0,             // copy number
                 fCheckOverlaps
    );
  //            ...                                    MV, MANY, copy#



  // G4Tubs *Detector = new G4Tubs("Detector", 0,r_out_HCAL3,2.5*l,0,tube_dPhi);
  // G4LogicalVolume *Detector_LV = new G4LogicalVolume(Detector, defaultMaterial, "Detector_LV");
  // G4PVPlacement *Detector_PL = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(),
  //                Detector_LV,          // its logical volume
  //                "Detector_PL",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,            // copy number
  //               fCheckOverlaps
  //   );

  //*Skeleton
  G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
  zRot->rotateZ(0);
  G4ThreeVector zTrans(0, 0, -l/2);
  G4ReflectZ3D  reflection(0);
  G4Transform3D rotation = G4Rotate3D(*zRot);
  G4CutTubs *Cut_ECAL1_buff = new G4CutTubs("Cut_ECAL1_buf", r_inn_ECAL1 - 10*cm, r_out_ECAL1 + 20*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_ECAL1;
  G4LogicalVolume *ECAL1_LV;
  G4PVPlacement *ECAL1;
  G4CutTubs *Cut_ECAL2_buff = new G4CutTubs("Cut_ECAL2_buf", r_inn_ECAL2 - 10*cm, r_out_ECAL2 + 20*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_ECAL2;
  G4LogicalVolume *ECAL2_LV;
  G4PVPlacement *ECAL2;
  G4CutTubs *Cut_ECAL3_buff = new G4CutTubs("Cut_ECAL3_buf", r_inn_ECAL3 - 10*cm, r_out_ECAL3 + 20*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_ECAL3;
  G4LogicalVolume *ECAL3_LV;
  G4PVPlacement *ECAL3;
  G4CutTubs *Cut_HCAL1_buff = new G4CutTubs("Cut_HCAL1_buf", r_inn_HCAL1 - 10*cm, r_out_HCAL1 + 20*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_HCAL1;
  G4LogicalVolume *HCAL1_LV;
  G4PVPlacement *HCAL1;
  G4CutTubs *Cut_HCAL2_buff = new G4CutTubs("Cut_HCAL2_buf", r_inn_HCAL2 - 20*cm, r_out_HCAL2 + 30*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_HCAL2;
  G4LogicalVolume *HCAL2_LV;
  G4PVPlacement *HCAL2;
  G4CutTubs *Cut_HCAL3_buff = new G4CutTubs("Cut_HCAL3_buf", r_inn_HCAL3 - 10*cm, r_out_HCAL3 + 20*cm, l, 0, tube_dPhi, G4ThreeVector (0,0,-1), G4ThreeVector (0,0,1));
  G4CutTubs *Cut_HCAL3;
  G4LogicalVolume *HCAL3_LV;
  G4PVPlacement *HCAL3;
  G4SubtractionSolid *subtraction;
  //* Skeleton
  G4Translate3D translation(0., 0., l/2);

  //*Preporation for visualization
  G4VisAttributes* ECAL1_VisAtt=
  new G4VisAttributes(true, G4Colour(100, 0, 255));
  G4VisAttributes* ECAL2_VisAtt=
  new G4VisAttributes(true, G4Colour(30, 143, 255));
  G4VisAttributes* ECAL3_VisAtt=
  new G4VisAttributes(true, G4Colour(0, 0, 255));
  G4VisAttributes* HCAL1_VisAtt=
  new G4VisAttributes(true, G4Colour(0, 255, 0));
  G4VisAttributes* HCAL2_VisAtt=
  new G4VisAttributes(true, G4Colour(255, 255, 0));
  G4VisAttributes* HCAL3_VisAtt=
  new G4VisAttributes(true, G4Colour(255, 0, 0));
  //*Preporation for visualization end


  G4int index_ECAL1 = 0;
  G4int index_ECAL2 = 0;
  G4int index_ECAL3 = 0;
  G4int index_HCAL1 = 0;
  G4int index_HCAL2 = 0;
  G4int index_HCAL3 = 0;
  G4int index_Gap = 0;

  // There are no way to make paramtrisation for G4CuTubs in the current Geant4 version that is way we have so many loops.
  for (auto copyNo=1; copyNo<=(kMaxPixel/2); copyNo++)//*loop that creates detector pixeles in positive z direction
  {

    d_theta = 2*atan(exp(-1*copyNo*d_eta));//* transformation from eta to theta
    vy = -cos( d_theta);
    vz = sin( d_theta);

    Cut_ECAL1 = new G4CutTubs("Cut", r_inn_ECAL1, r_out_ECAL1, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL1, Cut_ECAL1_buff, 0,zTrans);
    ECAL1_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL1_LV");

    Cut_ECAL2 = new G4CutTubs("Cut", r_inn_ECAL2, r_out_ECAL2, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL2, Cut_ECAL2_buff, 0,zTrans);
    ECAL2_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL2_LV");

    Cut_ECAL3 = new G4CutTubs("Cut", r_inn_ECAL3, r_out_ECAL3, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL3, Cut_ECAL3_buff, 0,zTrans);
    ECAL3_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL3_LV");

    Cut_HCAL1 = new G4CutTubs("Cut", r_inn_HCAL1, r_out_HCAL1, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL1, Cut_HCAL1_buff, 0,zTrans);
    HCAL1_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL1_LV");

    Cut_HCAL2 = new G4CutTubs("Cut", r_inn_HCAL2, r_out_HCAL2, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL2, Cut_HCAL2_buff, 0,zTrans);
    HCAL2_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL2_LV");

    Cut_HCAL3 = new G4CutTubs("Cut", r_inn_HCAL3, r_out_HCAL3, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,0,-1),G4ThreeVector (0,vy,vz));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL3, Cut_HCAL3_buff, 0,zTrans);
    HCAL3_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL3_LV");

    for (auto angle=0; angle<kMaxPixel; angle++)//* loop that creates detector pixels in phi axis
    {
      zRot = new G4RotationMatrix;
      zRot->rotateZ(angle*divided_tube_dPhi);
      G4Transform3D rotation = G4Rotate3D(*zRot);//* Rotate

      ECAL1 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL1_LV,          //* its logical volume
                    "ECAL1",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL1,                //* copy number
                    fCheckOverlaps
        );
      index_ECAL1++;


      ECAL2 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL2_LV,          //* its logical volume
                    "ECAL2",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL2,              //* copy number
                    fCheckOverlaps
        );
      index_ECAL2++;

      ECAL3 = new G4PVPlacement(
                    zRot,               //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL3_LV,          //* its logical volume
                    "ECAL3",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL3,             //* copy number
                    fCheckOverlaps
        );
      index_ECAL3++;


      HCAL1 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL1_LV,          //* its logical volume
                    "HCAL1",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL1 ,              //* copy number
                    fCheckOverlaps
        );
      index_HCAL1++;


      HCAL2 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL2_LV,          //* its logical volume
                    "HCAL2",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL2,              //* copy number
                    fCheckOverlaps
        );
      index_HCAL2++;


      HCAL3 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,-l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL3_LV,          //* its logical volume
                    "HCAL3",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL3,              //* copy number
                    fCheckOverlaps
        );
      index_HCAL3++;

      ECAL1_LV->SetVisAttributes(ECAL1_VisAtt);
      ECAL2_LV->SetVisAttributes(ECAL2_VisAtt);
      ECAL3_LV->SetVisAttributes(ECAL3_VisAtt);
      HCAL1_LV->SetVisAttributes(HCAL1_VisAtt);
      HCAL2_LV->SetVisAttributes(HCAL2_VisAtt);
      HCAL3_LV->SetVisAttributes(HCAL3_VisAtt);
      // G4cout<<"New circle element"<<G4endl;
    }
  // G4cout<<"New \eta element"<<G4endl;
  }
  zTrans = G4ThreeVector (0, 0, l/2);
  index_ECAL1 = 0;
  index_ECAL2 = 0;
  index_ECAL3 = 0;
  index_HCAL1 = 0;
  index_HCAL2 = 0;
  index_HCAL3 = 0;
  for (auto copyNo=1; copyNo<=(kMaxPixel/2); copyNo++) //*loop that creates detector pixeles in negative z direction
  {

    d_theta = 2*atan(exp(-1*copyNo*d_eta));
    vy = -cos( d_theta);
    vz = -sin( d_theta);


    Cut_ECAL1 = new G4CutTubs("Cut", r_inn_ECAL1, r_out_ECAL1, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL1, Cut_ECAL1_buff, 0,zTrans);
    ECAL1_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL1_LV");

    Cut_ECAL2 = new G4CutTubs("Cut", r_inn_ECAL2, r_out_ECAL2, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL2, Cut_ECAL2_buff, 0,zTrans);
    ECAL2_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL2_LV");

    Cut_ECAL3 = new G4CutTubs("Cut", r_inn_ECAL3, r_out_ECAL3, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_ECAL3, Cut_ECAL3_buff, 0,zTrans);
    ECAL3_LV = new G4LogicalVolume(subtraction, Material_ECAL, "ECAL3_LV");

    Cut_HCAL1 = new G4CutTubs("Cut", r_inn_HCAL1, r_out_HCAL1, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL1, Cut_HCAL1_buff, 0,zTrans);
    HCAL1_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL1_LV");

    Cut_HCAL2 = new G4CutTubs("Cut", r_inn_HCAL2, r_out_HCAL2, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL2, Cut_HCAL2_buff, 0,zTrans);
    HCAL2_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL2_LV");

    Cut_HCAL3 = new G4CutTubs("Cut", r_inn_HCAL3, r_out_HCAL3, l/2, (M_PI-divided_tube_dPhi)/2, divided_tube_dPhi ,G4ThreeVector (0,vy,vz),G4ThreeVector (0,0,1));
    subtraction = new G4SubtractionSolid("Sub", Cut_HCAL3, Cut_HCAL3_buff, 0,zTrans);
    HCAL3_LV = new G4LogicalVolume(subtraction, Material_HCAL, "HCAL3_LV");

    for (auto angle=0; angle<kMaxPixel; angle++)//* loop that creates detector pixels in phi axis
    {
      zRot = new G4RotationMatrix;
      zRot->rotateZ(1*angle*divided_tube_dPhi);

      ECAL1 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL1_LV,          //* its logical volume
                    "ECAL1_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL1,             //* copy number
                    fCheckOverlaps
        );


      index_ECAL1++;
      ECAL2 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL2_LV,          //* its logical volume
                    "ECAL2_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL2 ,               //* copy number
                    fCheckOverlaps
        );
      index_ECAL2++;
      ECAL3 = new G4PVPlacement(
                    zRot,                //*  rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    ECAL3_LV,          //* its logical volume
                    "ECAL3_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_ECAL3 ,             //* copy number
                    fCheckOverlaps
        );
      index_ECAL3++;

      HCAL1 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL1_LV,          //* its logical volume
                    "HCAL1_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL1 ,               //* copy number
                    fCheckOverlaps
        );
      index_HCAL1++;
      HCAL2 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL2_LV,          //* its logical volume
                    "HCAL2_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL2 ,              //* copy number
                    fCheckOverlaps
        );
      index_HCAL2++;
      HCAL3 = new G4PVPlacement(
                    zRot,                //* rotation
                    G4ThreeVector(0,0,l/2),  //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
                    HCAL3_LV,          //* its logical volume
                    "HCAL3_back",          //* its name
                    expHallLV,                //* its mother  volume
                    false,            //* no boolean operation
                    index_HCAL3  ,            //* copy number
                    fCheckOverlaps
        );
      index_HCAL3++;
    }


    ECAL1_LV->SetVisAttributes(ECAL1_VisAtt);
    ECAL2_LV->SetVisAttributes(ECAL2_VisAtt);
    ECAL3_LV->SetVisAttributes(ECAL3_VisAtt);
    HCAL1_LV->SetVisAttributes(HCAL1_VisAtt);
    HCAL2_LV->SetVisAttributes(HCAL2_VisAtt);
    HCAL3_LV->SetVisAttributes(HCAL3_VisAtt);
  }

  //* Make gap between ECAL3 and HCAL1
  // G4Tubs *Gap = new G4Tubs("Ineer_part", r_inn_GAP ,r_out_GAP, l, 0, tube_dPhi);
  // G4LogicalVolume *Gap_LV = new G4LogicalVolume(Gap, defaultMaterial, "Gap_LV");
  // G4PVPlacement *Gap_PL = new G4PVPlacement(
  //                0,                // no rotation
  //                G4ThreeVector(0., 0., 0. ),
  //                Gap_LV,          // its logical volume
  //                "Gap_PL",          // its name
  //                expHallLV,                // its mother  volume
  //                false,            // no boolean operation
  //                0,                // copy number
  //               fCheckOverlaps
  //   );


  //? Inner Tracker valuem
  long double eta = 1.5;
  long double thet = 2*atan(exp(-1*eta));
  long double widthiron_add = 350 * um;
  long double l_Pix0 =  r_out_trkPix0/tan(thet);
  long double l_Pix0_iron =  (r_out_trkPix0+widthiron_add)/tan(thet);
  long double l_Pix1 =  r_out_trkPix1/tan(thet);
  long double l_Pix1_iron =   (r_out_trkPix1+widthiron_add)/tan(thet);
  long double l_Pix2 =  r_out_trkPix2/tan(thet);
  long double l_Pix2_iron =   (r_out_trkPix2+widthiron_add)/tan(thet);
  long double l_Pix3 =  r_out_trkPix3/tan(thet);
  long double l_Pix3_iron =   (r_out_trkPix3+widthiron_add)/tan(thet);
  long double l_Pix4 =  r_out_trkPix4/tan(thet);
  long double l_Pix4_iron =   (r_out_trkPix4+widthiron_add)/tan(thet);

  long double l_Str0 =  r_out_trkStr0/tan(thet);
  long double l_Str0_iron =  (r_out_trkStr0+widthiron_add)/tan(thet);
  long double l_Str1 =  r_out_trkStr1/tan(thet);
  long double l_Str1_iron =  (r_out_trkStr1+widthiron_add)/tan(thet);
  long double l_Str2 =  r_out_trkStr2/tan(thet);
  long double l_Str2_iron =  (r_out_trkStr2+widthiron_add)/tan(thet);
  long double l_Str3 =  r_out_trkStr3/tan(thet);
  long double l_Str3_iron =  (r_out_trkStr3+widthiron_add)/tan(thet);


  G4Material* elSi = nistManager->FindOrBuildMaterial("G4_Si");
  //* Inner created to have magnetic field only in the inner part of detector
  G4Tubs *Inner= new G4Tubs("Inner_trkPix0",0,r_inn_ECAL1,l,0,tube_dPhi);
  G4LogicalVolume *Inner_LV = new G4LogicalVolume(Inner, defaultMaterial, "Inner_LV");
  G4PVPlacement *Inner_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_LV,          // its logical volume
                 "Inner_PL",          // its name
                 expHallLV,                // its mother  volume
                 false,            // no boolean operation
                 0  ,              // copy number
                fCheckOverlaps
    );


  G4Tubs *Inner_trkPix0 = new G4Tubs("Inner_trkPix0",r_inn_trkPix0,r_out_trkPix0,l_Pix0,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix0_LV = new G4LogicalVolume(Inner_trkPix0, elSi, "Inner_trkPix0_LV");
  G4PVPlacement *Inner_trkPix0_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix0_LV,          // its logical volume
                 "Inner_trkPix0_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0 ,               // copy number
                fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkPix0_iron = new G4Tubs("Inner_trkPix0_iron", r_out_trkPix0, r_out_trkPix0+widthiron_add, l_Pix0_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix0_LV_iron = new G4LogicalVolume(Inner_trkPix0_iron, iron, "Inner_trkPix0_LV_iron");
  G4PVPlacement *Inner_trkPix0_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix0_LV_iron,          // its logical volume
                 "Inner_trkPix0_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0 ,              // copy number
                 fCheckOverlaps
    );


  G4Tubs *Inner_trkPix1 = new G4Tubs("Inner_trkPix1",r_inn_trkPix1,r_out_trkPix1,l_Pix1,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix1_LV = new G4LogicalVolume(Inner_trkPix1, elSi, "Inner_trkPix1_LV");
  G4PVPlacement *Inner_trkPix1_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix1_LV,          // its logical volume
                 "Inner_trkPix1_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkPix1_iron = new G4Tubs("Inner_trkPix1_iron", r_out_trkPix1, r_out_trkPix1+widthiron_add, l_Pix1_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix1_LV_iron = new G4LogicalVolume(Inner_trkPix1_iron, iron, "Inner_trkPix1_LV_iron");
  G4PVPlacement *Inner_trkPix1_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix1_LV_iron,          // its logical volume
                 "Inner_trkPix1_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                fCheckOverlaps
    );



  G4Tubs *Inner_trkPix2 = new G4Tubs("Inner_trkPix2", r_inn_trkPix2, r_out_trkPix2, l_Pix2, 0, tube_dPhi);
  G4LogicalVolume *Inner_trkPix2_LV = new G4LogicalVolume(Inner_trkPix2, elSi, "Inner_trkPix2_LV");
  G4PVPlacement *Inner_trkPix2_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix2_LV,          // its logical volume
                 "Inner_trkPix2_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0      ,          // copy number
                fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkPix2_iron = new G4Tubs("Inner_trkPix2_iron", r_out_trkPix2, r_out_trkPix2+widthiron_add, l_Pix2_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix2_LV_iron = new G4LogicalVolume(Inner_trkPix2_iron, iron, "Inner_trkPix2_LV_iron");
  G4PVPlacement *Inner_trkPix2_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix2_LV_iron,          // its logical volume
                 "Inner_trkPix2_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0    ,           // copy number
                 fCheckOverlaps
    );


  G4Tubs *Inner_trkPix3 = new G4Tubs("Inner_trkPix3", r_inn_trkPix3, r_out_trkPix3, l_Pix3, 0, tube_dPhi);
  G4LogicalVolume *Inner_trkPix3_LV = new G4LogicalVolume(Inner_trkPix3, elSi, "Inner_trkPix3_LV");
  G4PVPlacement *Inner_trkPix3_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix3_LV,          // its logical volume
                 "Inner_trkPix3_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkPix3_iron = new G4Tubs("Inner_trkPix3_iron",r_out_trkPix3, r_out_trkPix3+widthiron_add, l_Pix3_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix3_LV_iron = new G4LogicalVolume(Inner_trkPix3_iron, iron, "Inner_trkPix3_LV_iron");
  G4PVPlacement *Inner_trkPix3_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix3_LV_iron,          // its logical volume
                 "Inner_trkPix3_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                 fCheckOverlaps
    );

  G4Tubs *Inner_trkPix4 = new G4Tubs("Inner_trkPix4", r_inn_trkPix4, r_out_trkPix4, l_Pix4, 0, tube_dPhi);
  G4LogicalVolume *Inner_trkPix4_LV = new G4LogicalVolume(Inner_trkPix4, elSi, "Inner_trkPix4_LV");
  G4PVPlacement *Inner_trkPix4_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix4_LV,          // its logical volume
                 "Inner_trkPix4_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkPix4_iron = new G4Tubs("Inner_trkPix4_iron",r_out_trkPix4,r_out_trkPix4+widthiron_add,l_Pix4_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkPix4_LV_iron = new G4LogicalVolume(Inner_trkPix4_iron, iron, "Inner_trkPix4_LV_iron");
  G4PVPlacement *Inner_trkPix4_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkPix4_LV_iron,          // its logical volume
                 "Inner_trkPix4_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                 fCheckOverlaps
    );

  G4Tubs *Inner_trkStr0 = new G4Tubs("Inner_trkStr0",r_inn_trkStr0,r_out_trkStr0,l_Str0,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr0_LV = new G4LogicalVolume(Inner_trkStr0, elSi, "Inner_trkStr0_LV");
  G4PVPlacement *Inner_trkStr0_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr0_LV,          // its logical volume
                 "Inner_trkStr0_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                  fCheckOverlaps
    );

  //additional iron layer to have 0.2 radiation length withou final iron layer
  G4Tubs *Inner_trkStr0_iron = new G4Tubs("Inner_trkStr0_iron",r_out_trkStr0,r_out_trkStr0+widthiron_add,l_Str0_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr0_LV_iron = new G4LogicalVolume(Inner_trkStr0_iron, iron, "Inner_trkStr0_LV_iron");
  G4PVPlacement *Inner_trkStr0_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr0_LV_iron,          // its logical volume
                 "Inner_trkStr0_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0  ,             // copy number
                 fCheckOverlaps
    );

  G4Tubs *Inner_trkStr1 = new G4Tubs("Inner_trkStr1",r_inn_trkStr1,r_out_trkStr1,l_Str1,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr1_LV = new G4LogicalVolume(Inner_trkStr1, elSi, "Inner_trkStr1_LV");
  G4PVPlacement *Inner_trkStr1_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr1_LV,          // its logical volume
                 "Inner_trkStr1_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                 fCheckOverlaps
    );


  G4Tubs *Inner_trkStr1_iron = new G4Tubs("Inner_trkStr1_iron",r_out_trkStr1,r_out_trkStr1+widthiron_add,l_Str1_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr1_LV_iron = new G4LogicalVolume(Inner_trkStr1_iron, iron, "Inner_trkStr1_LV_iron");
  G4PVPlacement *Inner_trkStr1_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr1_LV_iron,          // its logical volume
                 "Inner_trkStr1_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0    ,           // copy number
                 fCheckOverlaps
    );

  G4Tubs *Inner_trkStr2 = new G4Tubs("Inner_trkStr2",r_inn_trkStr2,r_out_trkStr2,l_Str2,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr2_LV = new G4LogicalVolume(Inner_trkStr2, elSi, "Inner_trkStr2_LV");
  G4PVPlacement *Inner_trkStr2_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr2_LV,          // its logical volume
                 "Inner_trkStr2_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                fCheckOverlaps
    );

    G4Tubs *Inner_trkStr2_iron = new G4Tubs("Inner_trkStr2_iron",r_out_trkStr2,r_out_trkStr2+widthiron_add,l_Str2_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr2_LV_iron = new G4LogicalVolume(Inner_trkStr2_iron, iron, "Inner_trkStr2_LV_iron");
  G4PVPlacement *Inner_trkStr2_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr2_LV_iron,          // its logical volume
                 "Inner_trkStr2_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0    ,           // copy number
                fCheckOverlaps
    );


  G4Tubs *Inner_trkStr3 = new G4Tubs("Inner_trkStr3",r_inn_trkStr3,r_out_trkStr3,l_Str3,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr3_LV = new G4LogicalVolume(Inner_trkStr3, elSi, "Inner_trkStr3_LV");
  G4PVPlacement *Inner_trkStr3_PL = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr3_LV,          // its logical volume
                 "Inner_trkStr3_PL",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0     ,          // copy number
                fCheckOverlaps
    );

  G4Tubs *Inner_trkStr3_iron = new G4Tubs("Inner_trkStr3_iron",r_out_trkStr3,r_out_trkStr3+widthiron_add,l_Str3_iron,0,tube_dPhi);
  G4LogicalVolume *Inner_trkStr3_LV_iron = new G4LogicalVolume(Inner_trkStr3_iron, iron, "Inner_trkStr3_LV_iron");
  G4PVPlacement *Inner_trkStr3_PL_iron = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Inner_trkStr3_LV_iron,          // its logical volume
                 "Inner_trkStr3_PL_iron",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,               // copy number
                fCheckOverlaps
    );



  long double iron_r_inn = 125 * cm;
  long double iron_width = 2*cm;
  long double l_Iron_Layer =  (iron_r_inn+iron_width)/tan(thet);
  G4Tubs *Iron_Layer = new G4Tubs("Inner_Iron_layer",iron_r_inn,iron_r_inn+iron_width,l_Iron_Layer,0,tube_dPhi);
  G4LogicalVolume *Iron_Layer_LV = new G4LogicalVolume(Iron_Layer, iron, "Inner_Iron_layer_LV");
  G4PVPlacement *Iron_Layer_PV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0. ),
                 Iron_Layer_LV,          // its logical volume
                 "Iron_gap",          // its name
                 Inner_LV,                // its mother  volume
                 false,            // no boolean operation
                 0               // copy number
                //  fCheckOverlaps
    );
  Inner_trkPix0_LV->SetVisAttributes(HCAL3_VisAtt);
  Inner_trkPix1_LV->SetVisAttributes(HCAL3_VisAtt);
  Inner_trkPix2_LV->SetVisAttributes(HCAL3_VisAtt);
  Inner_trkPix3_LV->SetVisAttributes(HCAL3_VisAtt);
  Inner_trkPix4_LV->SetVisAttributes(HCAL3_VisAtt);
  Inner_trkStr0_LV->SetVisAttributes(ECAL3_VisAtt);
  Inner_trkStr1_LV->SetVisAttributes(ECAL3_VisAtt);
  Inner_trkStr2_LV->SetVisAttributes(ECAL3_VisAtt);
  Inner_trkStr3_LV->SetVisAttributes(ECAL3_VisAtt);
  //? Inner Tracker valuem end
  //? EndCaps draft
  // G4VisAttributes* Cap_VisAtt =
  // new G4VisAttributes(true, G4Colour(30, 143, 255));
  // Cap_VisAtt-> SetForceSolid(true);
  // auto Cap_left = new G4Tubs("Cap_g", 0., r_out_HCAL3, r_out_HCAL3, 0., tube_dPhi);
  // G4LogicalVolume* Cap_left_LV = new G4LogicalVolume(Cap_left, Material_HCAL, "ECAL_LAYER-1_LV");
  // auto Cap_left_PL = new G4PVPlacement(0, G4ThreeVector(0,0,-l-r_out_HCAL3/2),Cap_left_LV,"Cap", Detector_LV, true,0);
  // auto Cap_right_PL = new G4PVPlacement(0, G4ThreeVector(0,0,l+r_out_HCAL3/2),Cap_left_LV,"Cap", Detector_LV, true,0);
  // Cap_left_LV->SetVisAttributes(Cap_VisAtt);
  //? EndCaps draft end
  //? Magnetic Field
  // long double fieldValue = 0;//3.8 * tesla;
  if (fieldValue!=0)
  {
    G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()
    ->GetFieldManager();
    Inner_LV->SetFieldManager(fieldMgr,true);
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    G4cout<<"Chord "<<fieldMgr->GetChordFinder()->GetDeltaChord()<<G4endl;
    fieldMgr->GetChordFinder()->SetDeltaChord(0.01);
    G4cout<<"Chord "<<fieldMgr->GetChordFinder()->GetDeltaChord()<<G4endl;
    
  }
  // //? Magnetic Field end
    // G4cout<<"Material_ECAL "<<Material_ECAL->GetRadlen()/cm<<G4endl;
    // G4cout<<"Material_HCAL "<<Material_HCAL->GetRadlen()/cm<<G4endl;
    return expHall;
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
