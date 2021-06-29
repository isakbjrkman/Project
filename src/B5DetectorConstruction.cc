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
/// \file B5DetectorConstruction.cc
/// \brief Implementation of the B5DetectorConstruction class

#include "B5DetectorConstruction.hh"
#include "B5MagneticField.hh"
#include "B5CellParameterisation.hh"
#include "B5HodoscopeSD.hh"
#include "B5DriftChamberSD.hh"
#include "B5EmCalorimeterSD.hh"
#include "B5HadCalorimeterSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal B5MagneticField* B5DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* B5DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::B5DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(nullptr),
  fHodoscope1Logical(nullptr), fHodoscope2Logical(nullptr),
  fWirePlane1Logical(nullptr), fWirePlane2Logical(nullptr),
  fCellLogical(nullptr), fHadCalScintiLogical(nullptr),
  fMagneticLogical(nullptr),
  fVisAttributes(),
  fArmAngle(30.*deg), fArmRotation(nullptr), fSecondArmPhys(nullptr)

{
  fArmRotation = new G4RotationMatrix();
  fArmRotation->rotateY(fArmAngle);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::~B5DetectorConstruction()
{
  delete fArmRotation;
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*cm, env_sizeZ = 10*cm;
  auto air = G4Material::GetMaterial("G4_AIR");
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",1.5*env_sizeXY,0.5*env_sizeXY,1.5*env_sizeZ);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  //colors
 
  G4Color blue(0.537, 0.812, 0.941);
  G4VisAttributes* blueVis = new G4VisAttributes(blue);
  
  G4Color purple(0.5, 0.0, 0.5);
  G4VisAttributes* purpleVis = new G4VisAttributes(purple);
  
  G4Color yellow(1.0, 1.0, 0.0);
  G4VisAttributes* yellowVis = new G4VisAttributes(yellow);
  
  G4Color grey(0.5, 0.5, 0.5);
  G4VisAttributes* greyVis = new G4VisAttributes(grey);
  
  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.4*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        air,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    worldLogical,            //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  
  
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4ThreeVector pos1 = G4ThreeVector(-13.255*mm, 13.255*mm, 0*mm);
        
  // Conical section shape       

  G4Box* solidShape1 =    
    new G4Box("Shape1", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  logicShape1->SetVisAttributes(blueVis);    
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //     
  // Shape 2
  //
  G4ThreeVector pos2 = G4ThreeVector(-13.255*mm, -13.255*mm, 0*mm);

  // Trapezoid shape       
  
  G4Box* solidShape2 =    
    new G4Box("Shape2",                      //its name
              0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm); //its size
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape1_mat,          //its material
                        "Shape2");           //its name
                        
  logicShape2->SetVisAttributes(blueVis);                      
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
                    
  G4ThreeVector pos3 = G4ThreeVector(13.255*mm, 13.255*mm, 0*mm);
        
  // Conical section shape       

  G4Box* solidShape3 =    
    new G4Box("Shape3", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicShape3 =                         
    new G4LogicalVolume(solidShape3,         //its solid
                        shape1_mat,          //its material
                        "Shape3");           //its name
      
   logicShape3->SetVisAttributes(blueVis);            
  new G4PVPlacement(0,                       //no rotation
                    pos3,                    //at position
                    logicShape3,             //its logical volume
                    "Shape3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                  
  
   G4ThreeVector pos4 = G4ThreeVector(13.255*mm, -13.255*mm, 0*mm);
        
  // Conical section shape       

  G4Box* solidShape4 =    
    new G4Box("Shape4", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicShape4 =                         
    new G4LogicalVolume(solidShape4,         //its solid
                        shape1_mat,          //its material
                        "Shape4");           //its name
   logicShape4->SetVisAttributes(blueVis);            
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape4,             //its logical volume
                    "Shape4",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                  
     G4ThreeVector pos5 = G4ThreeVector(0*mm, 0*mm, 11*mm);
        
  // Quartz window (monolithic)      

  G4Box* solidShape5 =    
    new G4Box("Shape5", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*2.0*mm);
                      
  G4LogicalVolume* logicShape5 =                         
    new G4LogicalVolume(solidShape5,         //its solid
                        shape1_mat,          //its material
                        "Shape5");           //its name
               
  logicShape5->SetVisAttributes(purpleVis);           
               
  new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logicShape5,             //its logical volume
                    "Shape5",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                                  
                                       
  
    G4Material* shape6_mat = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");
  G4ThreeVector pos7 = G4ThreeVector(-13.255*mm, 13.255*mm, 12.015*mm);
        
  // Conical section shape       
/*
  G4Box* solidShape7 =    
    new G4Box("Shape7", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
                      
  G4LogicalVolume* logicShape7 =                         
    new G4LogicalVolume(solidShape7,         //its solid
                        shape6_mat,          //its material
                        "Shape7");           //its name
        
  logicShape7->SetVisAttributes(yellowVis);      
               
  new G4PVPlacement(0,                       //no rotation
                    pos7,                    //at position
                    logicShape7,             //its logical volume
                    "Shape7",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  //     
  // Shape 2
  //

  G4ThreeVector pos8 = G4ThreeVector(-13.255*mm, -13.255*mm, 12.015*mm);

  // Trapezoid shape       
  
  G4Box* solidShape8=    
    new G4Box("Shape8",                      //its name
              0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm); //its size
                
  G4LogicalVolume* logicShape8 =                         
    new G4LogicalVolume(solidShape8,         //its solid
                        shape6_mat,          //its material
                        "Shape8");           //its name
                        
  logicShape8->SetVisAttributes(yellowVis);  
               
  new G4PVPlacement(0,                       //no rotation
                    pos8,                    //at position
                    logicShape8,             //its logical volume
                    "Shape8",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    

  G4ThreeVector pos9 = G4ThreeVector(13.255*mm, 13.255*mm, 12.015*mm);
        
  // Conical section shape       

  G4Box* solidShape9 =    
    new G4Box("Shape9", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
                      
  G4LogicalVolume* logicShape9 =                         
    new G4LogicalVolume(solidShape9,         //its solid
                        shape6_mat,          //its material
                        "Shape9");           //its name
                        
  logicShape9->SetVisAttributes(yellowVis);                        
               
  new G4PVPlacement(0,                       //no rotation
                    pos9,                    //at position
                    logicShape9,             //its logical volume
                    "Shape9",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                  
  
   G4ThreeVector pos10 = G4ThreeVector(13.255*mm, -13.255*mm, 12.015*mm);
        
  // Conical section shape       

  G4Box* solidShape10 =    
    new G4Box("Shape10", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
                      
  G4LogicalVolume* logicShape10 =                         
    new G4LogicalVolume(solidShape10,         //its solid
                        shape6_mat,          //its material
                        "Shape10");           //its name
                        
  logicShape10->SetVisAttributes(yellowVis);                        
               
  new G4PVPlacement(0,                       //no rotation
                    pos10,                    //at position
                    logicShape10,             //its logical volume
                    "Shape10",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                  
   */
   
  //Photocathode
  G4Material* shape11_mat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4ThreeVector pos11 = G4ThreeVector(0.0*mm, 0.0*mm, 21.03*mm);
           
  G4Box* solidShape11 =    
    new G4Box("Shape11", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*18.0*mm);
                      
  G4LogicalVolume* logicShape11 =                         
    new G4LogicalVolume(solidShape11,         //its solid
                        shape11_mat,          //its material
                        "Shape11");           //its name
                        
  logicShape11->SetVisAttributes(greyVis);                        
               
  new G4PVPlacement(0,                       //no rotation
                    pos11,                    //at position
                    logicShape11,             //its logical volume
                    "Shape11",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking                 
                                                          
                    
                    
                    
  /*
  // set step limit in tube with magnetic field  
  G4UserLimits* userLimits = new G4UserLimits(1*m);
  fMagneticLogical->SetUserLimits(userLimits);
  
  // first arm
  */
  auto firstArmSolid 
    = new G4Box("firstArmBox",1*27.5*mm, 1*27.5*mm, 1*0.01*mm);
  auto firstArmLogical
    = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,12.015*mm),firstArmLogical,
                    "firstArmPhysical",logicEnv,
                    false,0,checkOverlaps);
  
  // second arm
  /*auto secondArmSolid 
    = new G4Box("secondArmBox",1*27.5*mm, 1*27.5*mm, 1*0.01*mm);
  auto secondArmLogical
    = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
  //auto x = -5.*m * std::sin(fArmAngle);
  //auto z = 5.*m * std::cos(fArmAngle);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,35.015*mm),secondArmLogical,
                        "secondArmPhys",logicEnv,
                        false,0,checkOverlaps);
  */
  // hodoscopes in first arm
  
  auto chamber1Solid
    = new G4Box("chamber1Box",0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
  auto chamber1Logical
    = new G4LogicalVolume(chamber1Solid,shape6_mat,"chamber1Logical");

  for (auto i=0;i<2;i++) {
      G4double x1 = -13.255*mm+i*2*13.255*mm;
      for (auto j=0;j<2; j++){
      G4double y1 = -13.255*mm+j*2*13.255*mm;
      new G4PVPlacement(0,G4ThreeVector(x1,y1,0.0*mm),chamber1Logical,
                        "chamber1Physical",firstArmLogical,
                        false,i+j+1,checkOverlaps);
  	}
  }
    chamber1Logical->SetVisAttributes(yellowVis);  
    
/*    
    auto chamber2Solid
    = new G4Box("chamber2Box",0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
  auto chamber2Logical
    = new G4LogicalVolume(chamber2Solid,shape6_mat,"chamber2Logical");

  for (auto i=0;i<2;i++) {
      G4double x1 = -13.255*mm+i*2*13.255*mm;
      for (auto j=0;j<2; j++){
      G4double y1 = -13.255*mm+j*2*13.255*mm;
      new G4PVPlacement(0,G4ThreeVector(x1,y1,0.0*mm),chamber2Logical,
                        "chamber2Physical",secondArmLogical,
                        false,i+j+1,checkOverlaps);
  	}
  }
    chamber2Logical->SetVisAttributes(yellowVis);  
    
  */  
    
   
    
  /*
  // drift chambers in first arm
  auto chamber1Solid 
    = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
  auto chamber1Logical
    = new G4LogicalVolume(chamber1Solid,argonGas,"chamber1Logical");

  for (auto i=0;i<kNofChambers;i++) {
    G4double z1 = (i-kNofChambers/2)*0.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                      "chamber1Physical",firstArmLogical,
                      false,i,checkOverlaps);
  }
  */
  //"virtual" wire plane
  auto wirePlane1Solid 
    = new G4Box("wirePlane1Box", 0.5*26.5*mm, 0.5*26.5*mm, 0.25*0.01*mm);
  fWirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,air,"wirePlane1Logical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0),fWirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,
                    false,0,checkOverlaps);
                    
  
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  auto hodoscope1 = new B5HodoscopeSD(SDname="/hodoscope1");
  sdManager->AddNewDetector(hodoscope1);
  fWirePlane1Logical->SetSensitiveDetector(hodoscope1);

  auto hodoscope2 = new B5HodoscopeSD(SDname="/hodoscope2");
  sdManager->AddNewDetector(hodoscope2);
  fWirePlane1Logical->SetSensitiveDetector(hodoscope2);
  
  auto chamber1 = new B5DriftChamberSD(SDname="/chamber1");
  sdManager->AddNewDetector(chamber1);
  fWirePlane1Logical->SetSensitiveDetector(chamber1);

  auto chamber2 = new B5DriftChamberSD(SDname="/chamber2");
  sdManager->AddNewDetector(chamber2);
  fWirePlane1Logical->SetSensitiveDetector(chamber2);
  
  auto emCalorimeter = new B5EmCalorimeterSD(SDname="/EMcalorimeter");
  sdManager->AddNewDetector(emCalorimeter);
  fWirePlane1Logical->SetSensitiveDetector(emCalorimeter);
  
  auto hadCalorimeter = new B5HadCalorimeterSD(SDname="/HadCalorimeter");
  sdManager->AddNewDetector(hadCalorimeter);
  fWirePlane1Logical->SetSensitiveDetector(hadCalorimeter);

  // magnetic field ----------------------------------------------------------

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");
  // With a density different from the one defined in NIST
  // G4double density = 1.782e-03*g/cm3; 
  // nistManager->BuildMaterialWithNewDensity("B5_Ar","G4_Ar",density);
  // !! cases segmentation fault

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  
  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Vacuum "Galactic"
  // nistManager->FindOrBuildMaterial("G4_Galactic");

  // Vacuum "Air with low density"
  // auto air = G4Material::GetMaterial("G4_AIR");
  // G4double density = 1.0e-5*air->GetDensity();
  // nistManager
  //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::SetArmAngle(G4double val)
{
  if (!fSecondArmPhys) {
      G4cerr << "Detector has not yet been constructed." << G4endl;
      return;
  }
  
  fArmAngle = val;
  *fArmRotation = G4RotationMatrix();  // make it unit vector
  fArmRotation->rotateY(fArmAngle);
  auto x = -5.*m * std::sin(fArmAngle);
  auto z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));
  
  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::DefineCommands()
{
  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/B5/detector/", 
                                      "Detector control");

  // armAngle command

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
