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
 G4double mEfficMet[354];
 G4double mRindexMet[354];
 G4double mAbs[354];
 G4double mRefractiveIndexAir[354];
 G4double mAbsAir[354];
 G4double mRefractiveIndexCathode[354];
 G4double mAbsCathode[354];
 G4double mReflMet[354];
 G4double mReflMCP[354];

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
    mRefractiveIndexAir[i] = 1.0;       
    mAbsAir[i] = 0.3;
    mRefractiveIndexCathode[i] = 1.; 
    mAbsCathode[i] = 1.;
    mReflMet[i] = 0.9;  
    mReflMCP[i] = 0.5;	
  }
 

  
  // Construct materials
  ConstructMaterials();
  auto nist = G4NistManager::Instance();
  
  //Air properties
  //--------------------------------------------------------------
    G4double zet      = 1.0;
    G4double amass    = 1.01*g/mole;
    G4double density  = universe_mean_density;
    G4double pressure = 3.e-18*pascal;
    G4double tempture = 2.73*kelvin;
    G4Material* air = new G4Material("Vacuum", zet, amass, density,
                            kStateGas, tempture, pressure);
                           
  
  G4MaterialPropertiesTable* MPTair = new G4MaterialPropertiesTable();
  MPTair->AddProperty("RINDEX", mPhotonEnergyD, mRefractiveIndexAir, nBins);
  MPTair->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbsAir, nBins);
  
  G4cout << "Air G4MaterialPropertiesTable:" << G4endl;

  air->SetMaterialPropertiesTable(MPTair);
  
  // Envelope parameters
  //
  G4double env_sizeXY = 400*cm, env_sizeZ = 350*cm; 
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  
  auto worldSolid 
    = new G4Box("worldBox",0.6*env_sizeXY,0.5*env_sizeXY,1.1*env_sizeZ);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  G4OpticalSurface* opAirSurface = new G4OpticalSurface("AirSurface");
  opAirSurface->SetType(dielectric_dielectric);  
  opAirSurface->SetFinish(polished);  
  opAirSurface->SetModel(glisur); 
  opAirSurface->SetMaterialPropertiesTable(MPTair);	//recently added												
  
  new G4LogicalBorderSurface(
    "AirSurface", worldPhysical, worldPhysical, opAirSurface);


  //colors
 
  G4Color blue(0.537, 0.812, 0.941);
  G4VisAttributes* blueVis = new G4VisAttributes(blue);
  
  G4Color purple(0.5, 0.0, 0.5);
  G4VisAttributes* purpleVis = new G4VisAttributes(purple);
  
  G4Color yellow(1.0, 1.0, 0.0);
  G4VisAttributes* yellowVis = new G4VisAttributes(yellow);
  
  G4Color grey(0.5, 0.5, 0.5);
  G4VisAttributes* greyVis = new G4VisAttributes(grey);
  

  //     
  // Quartz radiator 1
  //   --------------------------------------------------------------
  G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.0855*g/mole);
  G4Element* elO = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Material* SiO2 = new G4Material("Silicon_dioxide", 2.533*g/cm3, 2);
  SiO2->AddElement(elO, 2);
  SiO2->AddElement(elSi, 1);
  
  
  
  G4double deltaY = 7.43*cm;
  G4double deltaZ = 335.5*cm-12.005*mm;
  
  
  G4MaterialPropertiesTable *MPTquartz = new G4MaterialPropertiesTable();
  MPTquartz->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins)->SetSpline(true);      
  MPTquartz->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbs, nBins)->SetSpline(true);
  MPTquartz->AddProperty("EFFICIENCY", mPhotonEnergyD, mEfficMet, nBins)->SetSpline(true);   //has no effect. Registered in FT0HadCalorimeterSD instead
  //MPTquartz->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMet, nBins);  //seems to be the issue. When enabled it decreases the number of registered photons by a factor >2
  
  
  G4cout << "Quartz G4MaterialPropertiesTable:" << G4endl;  
  SiO2->SetMaterialPropertiesTable(MPTquartz);
 

  G4Box* solidQuartz1 =    
    new G4Box("QuartzRadiator", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicQuartz1 =                         
    new G4LogicalVolume(solidQuartz1,         //its solid
                        SiO2,          //its material
                        "QuartzRadiator");           //its name                  

   G4VPhysicalVolume* logicVolumeQuartz1 =  new G4PVPlacement(0,G4ThreeVector(-13.255*mm,-13.255*mm+deltaY,0.0*mm+deltaZ),logicQuartz1,"QuartzRadiator",worldLogical,false,1,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz2 =  new G4PVPlacement(0,G4ThreeVector(-13.255*mm,13.255*mm+deltaY,0.0*mm+deltaZ),logicQuartz1,"QuartzRadiator",worldLogical,false,2,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz3 =  new G4PVPlacement(0,G4ThreeVector(13.255*mm,-13.255*mm+deltaY,0.0*mm+deltaZ),logicQuartz1,"QuartzRadiator",worldLogical,false,3,checkOverlaps);
   G4VPhysicalVolume* logicVolumeQuartz4 =  new G4PVPlacement(0,G4ThreeVector(13.255*mm,13.255*mm+deltaY,0.0*mm+deltaZ),logicQuartz1,"QuartzRadiator",worldLogical,false,4,checkOverlaps);                     
                                               
  
  
  G4OpticalSurface* opQuartzSurface = new G4OpticalSurface("QuartzSurface");
  opQuartzSurface->SetType(dielectric_metal);     //dielectric_metal
  opQuartzSurface->SetFinish(polishedbackpainted);  //polishedbackpainted   
  opQuartzSurface->SetModel(unified);  //unified													
  opQuartzSurface->SetMaterialPropertiesTable(MPTquartz);
  
    new G4LogicalBorderSurface(
    "QuartzSurface1", logicVolumeQuartz1, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface2", logicVolumeQuartz2, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface3", logicVolumeQuartz3, worldPhysical, opQuartzSurface);
    
    new G4LogicalBorderSurface(
    "QuartzSurface4", logicVolumeQuartz4, worldPhysical, opQuartzSurface);   
 
  
  logicQuartz1->SetVisAttributes(blueVis);              
   
        
  // Quartz window (monolithic)      
  //--------------------------------------------------------------
  
  G4ThreeVector posQuartzWindow = G4ThreeVector(0*mm, 0*mm+deltaY, 11*mm+deltaZ);
  
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
                
               
  G4OpticalSurface* opQuartzSurface5 = new G4OpticalSurface("QuartzSurface5");
  opQuartzSurface5->SetType(dielectric_metal);     
  opQuartzSurface5->SetFinish(polishedbackpainted);  
  opQuartzSurface5->SetModel(unified);
  opQuartzSurface5->SetMaterialPropertiesTable(MPTquartz); 
  
   new G4LogicalBorderSurface(
    "QuartzSurface5", logicQuartzWindow, worldPhysical, opQuartzSurface5);

         
                                      
 //Photocathode                    
 //--------------------------------------------------------------
  
  
  G4Material* GaAr = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");

  G4MaterialPropertiesTable* MPTgal = new G4MaterialPropertiesTable();
  MPTgal->AddProperty("RINDEX", mPhotonEnergyD, mRefractiveIndexCathode, nBins);
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
  G4double cathodePosX = 13.255*mm;  
  G4double cathodePosY = 13.255*mm;
  G4double cathodePosZ = 12.005*mm;  
    
    
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
 
  
  G4OpticalSurface* opCathodeSurface = new G4OpticalSurface("CathodeSurface");
  opCathodeSurface->SetType(dielectric_dielectric);     
  opCathodeSurface->SetFinish(polishedbackpainted);  
  opCathodeSurface->SetModel(unified);   
  opCathodeSurface->SetMaterialPropertiesTable(MPTgal);
  												
  
   new G4LogicalBorderSurface(
    "CathodeSurface1", logicCathode1, worldPhysical, opCathodeSurface);    
    
   new G4LogicalBorderSurface(
    "CathodeSurface2", logicCathode2, worldPhysical, opCathodeSurface);    
    
   new G4LogicalBorderSurface(
    "CathodeSurface3", logicCathode3, worldPhysical, opCathodeSurface);     
    
   new G4LogicalBorderSurface(
    "CathodeSurface4", logicCathode4, worldPhysical, opCathodeSurface);     
         
    
  fDetectorPlaneLogical->SetVisAttributes(yellowVis);
 
 
  //MCP-PMT (ceramic)
  //--------------------------------------------------------------
  G4Material* GaOx = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  
  G4MaterialPropertiesTable* MPTalu = new G4MaterialPropertiesTable();
  MPTalu->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMCP, nBins);

  G4cout << "GaOx G4MaterialPropertiesTable:" << G4endl;

  GaOx->SetMaterialPropertiesTable(MPTalu);
    
  G4ThreeVector posMCP = G4ThreeVector(0.0*mm, 0.0*mm+deltaY, 21.01*mm+deltaZ);    
           
  G4Box* solidMCP =    
    new G4Box("MCP", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*18.0*mm);
                      
  G4LogicalVolume* logicMCP =                         
    new G4LogicalVolume(solidMCP,         //its solid
                        GaOx,          //its material
                        "MCP");           //its name
                        
  logicMCP->SetVisAttributes(greyVis);                        
               
  G4VPhysicalVolume* logicPhysMCP = new G4PVPlacement(0,                       //no rotation
                    posMCP,                    //at position
                    logicMCP,             //its logical volume
                    "MCP",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking                 
  
  G4OpticalSurface* opMCPSurface = new G4OpticalSurface("MCPSurface");
  opMCPSurface->SetType(dielectric_dielectric);       
  opMCPSurface->SetFinish(polishedbackpainted);		
  opMCPSurface->SetModel(unified);			
  opMCPSurface->SetMaterialPropertiesTable(MPTalu);													
  
  new G4LogicalBorderSurface(
    "MCPSurface", logicPhysMCP, worldPhysical, opMCPSurface);
                                                                                        
 
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
