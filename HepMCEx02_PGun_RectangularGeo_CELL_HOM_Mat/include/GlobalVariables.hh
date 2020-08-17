
#ifndef GlobalVariables_h
#define GlobalVariables_h 1

#include "G4SystemOfUnits.hh"

#include "CaloRConstants.hh"

// constants (detector parameters)
// [experimental hall]
static const G4double R_EXPHALL= 5.*m;
static const G4double DZ_EXPHALL= 10.*m;

// [calorimeter]
static const G4double RIN_BARREL_ECAL= 2.*m;
static const G4double ROUT_BARREL_ECAL= 3.*m;
static const G4double DZ_BARREL_ECAL= 5.*m;
static const G4int    NLAYER_BARREL_CAL = 3 ; 
static const G4int NLAYER_ECAL = 10;

static const G4double RIN_BARREL_CAL= 2.*m;
static const G4double ROUT_BARREL_CAL= 3.*m;
static const G4double DZ_BARREL_CAL= 5.*m;

static const G4double RIN_BARREL_HCAL= 3.*m;
static const G4double ROUT_BARREL_HCAL= 4.*m;
static const G4double DZ_BARREL_HCAL= 5.*m;
static const G4int    NLAYER_BARREL_HCAL = 3 ;

static const G4double RIN_ENDCAP_CAL= 1.*m;
static const G4double ROUT_ENDCAP_CAL= 3.*m;
static const G4double DZ_ENDCAP_CAL= 0.5*m;

const G4int kLayer(6); 
const G4int K_NETA(kMaxPixel) ;
const G4int K_NPHI(kMaxPixel) ;
static const G4double ETA_MIN = -1 * calorSizeXY/2 ; 
static const G4double ETA_MAX =  calorSizeXY/2 ;
static const G4double PHI_MIN = -1 * calorSizeXY/2 ;
static const G4double PHI_MAX =  calorSizeXY/2 ;

const G4int N_STORE_TRAJ = 100; 
const G4int N_STORE_STEP = 500;


// [muon system]
static const G4double RIN_BARREL_MUON= 4.3*m;
// static const G4double ROUT_BARREL_MUON= 4.5*m;
static const G4double DX_BARREL_MUON= RIN_BARREL_MUON*std::cos(67.5*deg)-5.*cm;
static const G4double DY_BARREL_MUON= 10.*cm;
static const G4double DZ_BARREL_MUON= 7.*m;

static const G4double RIN_ENDCAP_MUON=  1.*m;
static const G4double ROUT_ENDCAP_MUON= 4.5*m;
static const G4double DZ_ENDCAP_MUON= 10.*cm;

#endif
