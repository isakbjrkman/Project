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
/// \file FT0DetectorConstruction.cc
/// \brief Implementation of the FT0DetectorConstruction class
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <array>

#include "FT0Constants.hh"
#include "FT0DetectorConstruction.hh"
#include "FT0HadCalorimeterSD.hh"


#include "G4OpBoundaryProcess.hh"
#include "G4Cerenkov.hh"
#include "G4OpticalPhoton.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
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
#include "G4LogicalBorderSurface.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0DetectorConstruction::FT0DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(nullptr),
  fDetectorPlaneLogical(nullptr),
  fVisAttributes()

{
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0DetectorConstruction::~FT0DetectorConstruction()
{
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* FT0DetectorConstruction::Construct()
{
   
 G4double mPhotonEnergyD[354];
 G4double mEfficMet[354];  //variable used in had class, i.e. unnecessary here
 G4double mRindexMet[354];
 G4double mAbs[354];
 G4double mRefractiveIndexCathode[354];
 G4double mAbsCathode[354];
 G4double mReflQuartz[354];
 G4double mReflBlackPaper[354];
 G4double mAbsBlackPaper[354]; 
 G4double mRindexBlackPaper[354];
 G4double mRindexQuartz[354];

 G4double refl[354];
 G4double transmission[354];

 G4double mReflMirror[354];
 G4double mRindexMirror[354];
 
 G4double mReflCathode[354];
 
 int nBins = 354; //sizeof(mPhotonEnergyD)/sizeof(mPhotonEnergyD[0]);
 
   
 //Read quartz properties data from .txt file.

 std::ifstream myfile("quartz.txt");
 std::string line;
 if (!myfile) {
 std::cout << "File read failed" << std::endl;
 }
 std::getline(myfile, line);
 
 int num = 0;
 

 while(std::getline(myfile, line))  { 
 
        double energy;
        double abs;
        double ref;
        double eff;
std::istringstream ss(line);

ss >> energy >> abs >> ref >> eff;
 mPhotonEnergyD[num] = energy*eV;
 mAbs[num] = abs*cm;
 mRindexMet[num] = ref;
 mEfficMet[num] = eff;  
 std::cout << mPhotonEnergyD[num] << " " << mAbs[num] << " " << mRindexMet[num] << " " << mEfficMet[num] << "\n";
 ++num;
     } 

 for (auto i = 0; i < nBins; i++) {
    mReflBlackPaper[i] = 0.0;
    mRindexBlackPaper[i] = 1.;                  
     
    mRefractiveIndexCathode[i] = 1.0;  //3.3 
    mAbsCathode[i] = 1*pow(10,-4)*cm; //1000000.*cm;

   transmission[i] = 0.0; 
     
   mReflMirror[i] = 1.0;
   mRindexMirror[i] = 1.0;
  }
 
  G4double deltaY = 7.43*cm;
  G4double deltaZ = 335.5*cm;
  
  
  // Construct materials
  ConstructMaterials();
  auto nist = G4NistManager::Instance();
  
  //Vacuum properties
  //--------------------------------------------------------------
    G4double zet      = 1.0;
    G4double amass    = 1.01*g/mole;
    G4double density  = universe_mean_density;
    G4double pressure = 3.e-18*pascal;
    G4double tempture = 2.73*kelvin;
    G4Material* vacuum = new G4Material("Vacuum", zet, amass, density,
                            kStateGas, tempture, pressure);
                           

  // Envelope parameters
  //
  G4double env_sizeXY = 400*cm, env_sizeZ = 350*cm; 
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  
  auto worldSolid 
    = new G4Box("worldBox",0.6*env_sizeXY,0.5*env_sizeXY,1.1*env_sizeZ);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");
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
  
  G4Color black(0., 0., 0.);
  G4VisAttributes* blackVis = new G4VisAttributes(black);
    
    
  G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.0855*g/mole);
  G4Element* elO = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Material* SiO2 = new G4Material("Silicon_dioxide", 2.533*g/cm3, 2);
  SiO2->AddElement(elO, 2);
  SiO2->AddElement(elSi, 1);
  
 
  
  G4MaterialPropertiesTable *MPTquartz = new G4MaterialPropertiesTable();
  MPTquartz->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins)->SetSpline(true);      
  MPTquartz->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbs, nBins)->SetSpline(true);
  
  G4cout << "Quartz G4MaterialPropertiesTable:" << G4endl;  
  SiO2->SetMaterialPropertiesTable(MPTquartz);

 
    
  //     
  // Quartz radiators
  //   --------------------------------------------------------------
 

  G4Box* solidQuartz1 =    
    new G4Box("QuartzRadiator", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicQuartz1 =                         
    new G4LogicalVolume(solidQuartz1,         //its solid
                        SiO2,          //its material
                        "QuartzRadiator");           //its name                  

   G4VPhysicalVolume* logicVolumeQuartz1 =  new G4PVPlacement(0,G4ThreeVector(-(13.25+dist/2)*mm,-(13.25+dist/2)*mm+deltaY,-12.005*mm+deltaZ),logicQuartz1,"QuartzRadiator1",worldLogical,false,1,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz2 =  new G4PVPlacement(0,G4ThreeVector(-(13.25+dist/2)*mm,(13.25+dist/2)*mm+deltaY,-12.005*mm+deltaZ),logicQuartz1,"QuartzRadiator2",worldLogical,false,2,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz3 =  new G4PVPlacement(0,G4ThreeVector((13.25+dist/2)*mm,-(13.25+dist/2)*mm+deltaY,-12.005*mm+deltaZ),logicQuartz1,"QuartzRadiator3",worldLogical,false,3,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz4 =  new G4PVPlacement(0,G4ThreeVector((13.25+dist/2)*mm,(13.25+dist/2)*mm+deltaY,-12.005*mm+deltaZ),logicQuartz1,"QuartzRadiator4",worldLogical,false,4,checkOverlaps);                     
   
   
         
  
  G4MaterialPropertiesTable *MPTquartzMirror = new G4MaterialPropertiesTable();
  MPTquartzMirror->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins)->SetSpline(true);      
  MPTquartzMirror->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMirror, nBins)->SetSpline(true);

  //mirror-like properties facing the vacuum
  G4OpticalSurface* opQuartzSurface = new G4OpticalSurface("QuartzSurfaceMirror");
  opQuartzSurface->SetType(dielectric_dielectric);     //dielectric_metal
  opQuartzSurface->SetFinish(polishedbackpainted);  //polishedbackpainted   
  opQuartzSurface->SetModel(unified);    //unified													
  opQuartzSurface->SetMaterialPropertiesTable(MPTquartzMirror);
  
    new G4LogicalBorderSurface(
    "QuartzSurface1", logicVolumeQuartz1, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface2", logicVolumeQuartz2, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface3", logicVolumeQuartz3, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface4", logicVolumeQuartz4, worldPhysical, opQuartzSurface);  
    
  
  logicQuartz1->SetVisAttributes(blueVis);              
   
  
  //Black paper on top of radiator    
  //--------------------------------------------------------------

  
  G4Element* elC = new G4Element("Carbon", "C", 6., 12.0107*g/mole); 
  G4Material* C = new G4Material("Carbon Material", 3.52*g/cm3, 1); 
  C->AddElement(elC, 1); 
   
    
  G4MaterialPropertiesTable *MPTtopPaper = new G4MaterialPropertiesTable();
  MPTtopPaper->AddProperty("RINDEX", mPhotonEnergyD, mRindexBlackPaper, nBins); //commented recently
  MPTtopPaper->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflBlackPaper, nBins);
  MPTtopPaper->AddProperty("TRANSMISSION", mPhotonEnergyD, transmission, nBins);
  
  
  C->SetMaterialPropertiesTable(MPTtopPaper);
   
  
   G4double paperThick = 0.1;     //CHANGED FROM 0.1
        
   G4Box* solidTopPaper =    
    new G4Box("TopBlackPaper", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*paperThick*mm);
                          
  G4LogicalVolume* logicTopPaper1 =                         
    new G4LogicalVolume(solidTopPaper,         //its solid
                        C,          //its material
                        "TopBlackPaper");           //its name                  
  
  //absorbing non-reflective properties
   
  G4OpticalSurface* opTopPaperSurface = new G4OpticalSurface("TopBlackPaperSurface");
  opTopPaperSurface->SetType(dielectric_dielectric);     //dielectric_metal
  opTopPaperSurface->SetFinish(ground);  //polishedbackpainted   
  opTopPaperSurface->SetModel(unified);  //unified													
  opTopPaperSurface->SetMaterialPropertiesTable(MPTtopPaper);

      // paper->vacuum  
    new G4LogicalSkinSurface(
    "TopBlackPaperSurface1", logicTopPaper1, opTopPaperSurface);
    
    new G4LogicalSkinSurface(
    "TopBlackPaperSurface2",logicTopPaper1, opTopPaperSurface);
    
    new G4LogicalSkinSurface(
    "TopBlackPaperSurface3",  logicTopPaper1, opTopPaperSurface);
    
    new G4LogicalSkinSurface(
    "TopBlackPaperSurface4", logicTopPaper1, opTopPaperSurface);
    
   
   G4VPhysicalVolume* logicVolumeTopPaper1 =  new G4PVPlacement(0,G4ThreeVector(-(13.25+dist/2)*mm,-(13.25+dist/2)*mm+deltaY,deltaZ-22.005*mm-paperThick/2*mm),logicTopPaper1,"TopBlackPaper1",worldLogical,false,1,checkOverlaps);
   G4VPhysicalVolume* logicVolumeTopPaper2 =  new G4PVPlacement(0,G4ThreeVector(-(13.25+dist/2)*mm,(13.25+dist/2)*mm+deltaY,deltaZ-22.005*mm-paperThick/2*mm),logicTopPaper1,"TopBlackPaper2",worldLogical,false,2,checkOverlaps);
   G4VPhysicalVolume* logicVolumeTopPaper3 =  new G4PVPlacement(0,G4ThreeVector((13.25+dist/2)*mm,-(13.25+dist/2)*mm+deltaY,deltaZ-22.005*mm-paperThick/2*mm),logicTopPaper1,"TopBlackPaper3",worldLogical,false,3,checkOverlaps);
   G4VPhysicalVolume* logicVolumeTopPaper4 =  new G4PVPlacement(0,G4ThreeVector((13.25+dist/2)*mm,(13.25+dist/2)*mm+deltaY,deltaZ-22.005*mm-paperThick/2*mm),logicTopPaper1,"TopBlackPaper4",worldLogical,false,4,checkOverlaps);   
   
       
                   
   logicTopPaper1->SetVisAttributes(blackVis);     
        
        
  // Quartz window (monolithic)      
  //--------------------------------------------------------------
  
  G4MaterialPropertiesTable *MPTquartzWindow = new G4MaterialPropertiesTable();
  MPTquartzWindow->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins)->SetSpline(true);      
  MPTquartzWindow->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbs, nBins)->SetSpline(true);
    
    
  G4ThreeVector posQuartzWindow = G4ThreeVector(0*mm, 0*mm+deltaY, -1.005*mm+deltaZ);
  
  G4Box* solidQuartz5 =    
    new G4Box("QuartzWindow", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*2.0*mm);
                      
  G4LogicalVolume* logicQuartz5 =                         
    new G4LogicalVolume(solidQuartz5,         //its solid
                        SiO2,          //its material
                        "QuartzWindow");           //its name
               
  logicQuartz5->SetVisAttributes(purpleVis);           
               
  
  G4VPhysicalVolume* logicQuartzWindow = new G4PVPlacement(0,                      
                    posQuartzWindow,             //at position
                    logicQuartz5,             //its logical volume
                    "QuartzWindow",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


 //Photocathode                    
 //-------------------------------------------------------------- 
  G4Element* elGa = new G4Element("Gallium", "Ga", 31., 69.723*g/mole);
  G4Element* elAr = new G4Element("Arsenic", "Ar", 33., 74.9216*g/mole);
  G4Material* GaAr = new G4Material("Galliumarsenide", 5.32*g/cm3, 2);
  GaAr->AddElement(elGa, 1);
  GaAr->AddElement(elAr, 1);

  G4MaterialPropertiesTable* MPTgal = new G4MaterialPropertiesTable();
  MPTgal->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins); 
  MPTgal->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbsCathode, nBins);
  
  G4cout << "GaAr G4MaterialPropertiesTable:" << G4endl;

  GaAr->SetMaterialPropertiesTable(MPTgal); 

  G4double cathodeX = 26.5*mm;
  G4double cathodeY = 26.5*mm;
  G4double cathodeZ = 0.01*mm;  
 
  auto chamber1Solid
    = new G4Box("chamber1Box",0.5*cathodeX, 0.5*cathodeY, 0.5*cathodeZ);
  fDetectorPlaneLogical									
    = new G4LogicalVolume(chamber1Solid,GaAr,"fDetectorPlaneLogical");
  G4double cathodePosX = (13.25+dist/2)*mm;  
  G4double cathodePosY = (13.25+dist/2)*mm;
  G4double cathodePosZ = 0*mm;  
         
      
      
       G4VPhysicalVolume* logicCathode1 = new G4PVPlacement(0,G4ThreeVector(-cathodePosX,-cathodePosY+deltaY,cathodePosZ+deltaZ),fDetectorPlaneLogical,   
                        "chamber1Physical",worldLogical,
                        false,1,checkOverlaps);
    G4VPhysicalVolume* logicCathode2 = new G4PVPlacement(0,G4ThreeVector(-cathodePosX,cathodePosY+deltaY,cathodePosZ+deltaZ),fDetectorPlaneLogical,   
                        "chamber1Physical",worldLogical,
                        false,2,checkOverlaps);
    G4VPhysicalVolume* logicCathode3 = new G4PVPlacement(0,G4ThreeVector(cathodePosX,-cathodePosY+deltaY,cathodePosZ+deltaZ),fDetectorPlaneLogical,   
                        "chamber1Physical",worldLogical,
                        false,3,checkOverlaps);
    G4VPhysicalVolume* logicCathode4 = new G4PVPlacement(0,G4ThreeVector(cathodePosX,cathodePosY+deltaY,cathodePosZ+deltaZ),fDetectorPlaneLogical,   
                        "chamber1Physical",worldLogical,
                        false,4,checkOverlaps);   
    
      
  //Print geometry data             
  G4cout << "---Detector Positions in World---" << G4endl;                                                          
  G4cout << "For DetectorID = " << logicCathode1->GetCopyNo() << "  x = " <<  logicCathode1->GetTranslation()(0)-0.5*cathodeX  << " - " << logicCathode1->GetTranslation()(0)+0.5*cathodeX  
  << " mm, y = " <<  logicCathode1->GetTranslation()(1)-0.5*cathodeY  << " - " << logicCathode1->GetTranslation()(1)+0.5*cathodeY 
  << " mm, z = " <<  logicCathode1->GetTranslation()(2)-0.5*cathodeZ  << " - " << logicCathode1->GetTranslation()(2)+0.5*cathodeZ << " mm" << G4endl;
  
  G4cout << "For DetectorID = " << logicCathode2->GetCopyNo() << "  x = " <<  logicCathode2->GetTranslation()(0)-0.5*cathodeX  << " - " << logicCathode2->GetTranslation()(0)+0.5*cathodeX 
  << " mm, y = " <<  logicCathode2->GetTranslation()(1)-0.5*cathodeY  << " - " << logicCathode2->GetTranslation()(1)+0.5*cathodeY 
  << " mm, z = " <<  logicCathode2->GetTranslation()(2)-0.5*cathodeZ  << " - " << logicCathode2->GetTranslation()(2)+0.5*cathodeZ << " mm" << G4endl;
  
  G4cout << "For DetectorID = " << logicCathode3->GetCopyNo() << "  x = " <<  logicCathode3->GetTranslation()(0)-0.5*cathodeX  << " - " << logicCathode3->GetTranslation()(0)+0.5*cathodeX 
  << " mm, y = " <<  logicCathode3->GetTranslation()(1)-0.5*cathodeY  << " - " << logicCathode3->GetTranslation()(1)+0.5*cathodeY
  << " mm, z = " <<  logicCathode3->GetTranslation()(2)-0.5*cathodeZ  << " - " << logicCathode3->GetTranslation()(2)+0.5*cathodeZ << " mm" << G4endl;
  
  G4cout << "For DetectorID = " << logicCathode4->GetCopyNo() << "  x = " <<  logicCathode4->GetTranslation()(0)-0.5*cathodeX  << " - " << logicCathode4->GetTranslation()(0)+0.5*cathodeX  
  << " mm, y = " <<  logicCathode4->GetTranslation()(1)-0.5*cathodeY  << " - " <<  logicCathode4->GetTranslation()(1)+0.5*cathodeY 
  << " mm, z = " <<  logicCathode4->GetTranslation()(2)-0.5*cathodeZ  << " - " <<  logicCathode4->GetTranslation()(2)+0.5*cathodeZ << " mm" << G4endl;
  
    
  fDetectorPlaneLogical->SetVisAttributes(yellowVis);
 
 
  //MCP-PMT (ceramic)
  //--------------------------------------------------------------
  
  G4Element* elH = new G4Element("Hydrogen", "H", 1., 1.00784*g/mole);
  G4Element* elAl = new G4Element("Aluminium", "Al", 13., 26.981539*g/mole);
  
  G4Material* Ceramic = new G4Material("Ceramic", 2.65*g/cm3, 4); //kaolinite
  Ceramic->AddElement(elAl, 2);
  Ceramic->AddElement(elSi, 2);
  Ceramic->AddElement(elO, 9);
  Ceramic->AddElement(elH, 4);
    
  
  G4ThreeVector posMCP = G4ThreeVector(0.0*mm, 0.0*mm+deltaY, 9.005*mm+deltaZ);    
   
           
  G4Box* solidMCP =    
    new G4Box("MCP", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*18.0*mm);
                      
  G4LogicalVolume* logicMCP =                         
    new G4LogicalVolume(solidMCP,         //its solid
                        Ceramic,          //its material
                        "MCP");           //its name
                                                             
  G4VPhysicalVolume* logicPhysMCP = new G4PVPlacement(0,                       //no rotation
                    posMCP,                    //at position
                    logicMCP,             //its logical volume
                    "MCP",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking                 
       
                                                                                        
  logicMCP->SetVisAttributes(greyVis);  
  
 
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;  //worldPhysical
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  auto hadCalorimeter = new FT0HadCalorimeterSD(SDname="/HadCalorimeter");
  sdManager->AddNewDetector(hadCalorimeter);
  fDetectorPlaneLogical->SetSensitiveDetector(hadCalorimeter);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0DetectorConstruction::ConstructMaterials()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void FT0DetectorConstruction::DefineCommands()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
