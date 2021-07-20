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
  fWirePlane1Logical(nullptr),
  fHadCalScintiLogical(nullptr),
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
 int nBins = sizeof(mPhotonEnergyD)/sizeof(mPhotonEnergyD[0]);
 
 
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
  
  
    G4double zet      = 1.0;
    G4double amass    = 1.01*g/mole;
    G4double density  = universe_mean_density;
    G4double pressure = 3.e-18*pascal;
    G4double tempture = 2.73*kelvin;
    G4Material* air = new G4Material("Vacuum", zet, amass, density,
                            kStateGas, tempture, pressure);
                           
  //Air properties
  //--------------------------------------------------------------
  G4MaterialPropertiesTable* MPT2 = new G4MaterialPropertiesTable();
  MPT2->AddProperty("RINDEX", mPhotonEnergyD, mRefractiveIndexAir, nBins);
  MPT2->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbsAir, nBins);
  
  G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
  MPT2->DumpTable();

  air->SetMaterialPropertiesTable(MPT2);
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*cm, env_sizeZ = 10*cm; 
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  
  auto worldSolid 
    = new G4Box("worldBox",2.0*env_sizeXY,0.8*env_sizeXY,2.0*env_sizeZ);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  G4OpticalSurface* opAirSurface = new G4OpticalSurface("AirSurface");
  opAirSurface->SetType(dielectric_dielectric);  
  opAirSurface->SetFinish(polished);  
  opAirSurface->SetModel(glisur); 
  													
  
   G4LogicalBorderSurface* airSurface = new G4LogicalBorderSurface(
    "AirSurface", worldPhysical, worldPhysical, opAirSurface);

  G4OpticalSurface* opticalSurface4 = dynamic_cast<G4OpticalSurface*>(
    airSurface->GetSurface(worldPhysical, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface4)
    opticalSurface4->DumpInfo();

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
  
  
  G4MaterialPropertiesTable *MPT = new G4MaterialPropertiesTable();
  MPT->AddProperty("RINDEX", mPhotonEnergyD, mRindexMet, nBins)->SetSpline(true);      
  MPT->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbs, nBins)->SetSpline(true);
  MPT->AddProperty("EFFICIENCY", mPhotonEnergyD, mEfficMet, nBins)->SetSpline(true);
  MPT->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMet, nBins);
  
  
  G4cout << "Quartz G4MaterialPropertiesTable:" << G4endl;
  MPT->DumpTable();  
  SiO2->SetMaterialPropertiesTable(MPT);
 

  G4Box* solidShape1 =    
    new G4Box("Shape1", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        SiO2,          //its material
                        "Shape1");           //its name                  

   G4VPhysicalVolume* logicQuartz1 =  new G4PVPlacement(0,G4ThreeVector(-13.255*mm,-13.255*mm,0.0*mm),logicShape1,"Shape1",worldLogical,false,1,checkOverlaps);
   G4VPhysicalVolume* logicQuartz2 =  new G4PVPlacement(0,G4ThreeVector(-13.255*mm,13.255*mm,0.0*mm),logicShape1,"Shape1",worldLogical,false,2,checkOverlaps);
   G4VPhysicalVolume* logicQuartz3 =  new G4PVPlacement(0,G4ThreeVector(13.255*mm,-13.255*mm,0.0*mm),logicShape1,"Shape1",worldLogical,false,3,checkOverlaps);
   G4VPhysicalVolume* logicQuartz4 =  new G4PVPlacement(0,G4ThreeVector(13.255*mm,13.255*mm,0.0*mm),logicShape1,"Shape1",worldLogical,false,4,checkOverlaps);                     
                                               
  
  
  G4OpticalSurface* opQuartzSurface = new G4OpticalSurface("QuartzSurface");
  opQuartzSurface->SetType(dielectric_metal);     //dielectric_metal
  opQuartzSurface->SetFinish(Polished_LUT);  //polishedbackpainted   works: Rough_LUT
  opQuartzSurface->SetModel(unified);  //unified													
  
   G4LogicalBorderSurface* quartzSurface1 = new G4LogicalBorderSurface(
    "QuartzSurface1", logicQuartz1, worldPhysical, opQuartzSurface);
    
   G4LogicalBorderSurface* quartzSurface2 = new G4LogicalBorderSurface(
    "QuartzSurface2", logicQuartz2, worldPhysical, opQuartzSurface);
    
   G4LogicalBorderSurface* quartzSurface3 = new G4LogicalBorderSurface(
    "QuartzSurface3", logicQuartz3, worldPhysical, opQuartzSurface);
    
   G4LogicalBorderSurface* quartzSurface4 = new G4LogicalBorderSurface(
    "QuartzSurface4", logicQuartz4, worldPhysical, opQuartzSurface);   

  G4OpticalSurface* opticalSurface11 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface1->GetSurface(logicQuartz1, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface11)
    opticalSurface11->DumpInfo();
    
  G4OpticalSurface* opticalSurface12 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface2->GetSurface(logicQuartz2, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface12)
    opticalSurface12->DumpInfo();
    
  G4OpticalSurface* opticalSurface13 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface3->GetSurface(logicQuartz3, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface13)
    opticalSurface13->DumpInfo();
    
  G4OpticalSurface* opticalSurface14 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface4->GetSurface(logicQuartz4, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface14)
    opticalSurface14->DumpInfo();      
  
  
  logicShape1->SetVisAttributes(blueVis);              
 
  G4ThreeVector pos5 = G4ThreeVector(0*mm, 0*mm, 11*mm);
      
        
  // Quartz window (monolithic)      
   //--------------------------------------------------------------
  G4Box* solidShape5 =    
    new G4Box("Shape5", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*2.0*mm);
                      
  G4LogicalVolume* logicShape5 =                         
    new G4LogicalVolume(solidShape5,         //its solid
                        SiO2,          //its material
                        "Shape5");           //its name
               
  logicShape5->SetVisAttributes(purpleVis);           
               
  
  G4VPhysicalVolume* logicPhys2 = new G4PVPlacement(0,                      
                    pos5,                    //at position
                    logicShape5,             //its logical volume
                    "Shape5",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
               
  G4OpticalSurface* opQuartzSurface5 = new G4OpticalSurface("QuartzSurface5");
  opQuartzSurface5->SetType(dielectric_metal);     
  opQuartzSurface5->SetFinish(Polished_LUT);  
  opQuartzSurface5->SetModel(unified); 
  
   G4LogicalBorderSurface* quartzSurface5 = new G4LogicalBorderSurface(
    "QuartzSurface5", logicPhys2, worldPhysical, opQuartzSurface5);

  G4OpticalSurface* opticalSurface5 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface5->GetSurface(logicPhys2, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface5)
    opticalSurface5->DumpInfo();         
         
                                      
 //Photocathode                    
 //--------------------------------------------------------------
  
  
  G4Material* GaAr = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");

  G4MaterialPropertiesTable* MPT3 = new G4MaterialPropertiesTable();
  MPT3->AddProperty("RINDEX", mPhotonEnergyD, mRefractiveIndexCathode, nBins);
  MPT3->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbsCathode, nBins);
  
  G4cout << "GaAr G4MaterialPropertiesTable:" << G4endl;
  MPT3->DumpTable();

  GaAr->SetMaterialPropertiesTable(MPT3);
 
  G4double cathodeX = 26.5*mm;
  G4double cathodeY = 26.5*mm;
  G4double cathodeZ = 0.01*mm;
 
  auto chamber1Solid
    = new G4Box("chamber1Box",0.5*cathodeX, 0.5*cathodeY, 0.5*cathodeZ);
  fWirePlane1Logical									
    = new G4LogicalVolume(chamber1Solid,GaAr,"fWirePlane1Logical");
  G4double cathodePosX = 13.255*mm;  
  G4double cathodePosY = 13.255*mm;
  G4double cathodePosZ = 12.005*mm;  //12.015*mm correct?  if no space between::12.005
    
    
    G4VPhysicalVolume* logicCathode1 = new G4PVPlacement(0,G4ThreeVector(-cathodePosX,-cathodePosY,cathodePosZ),fWirePlane1Logical,   
                        "chamber1Physical",worldLogical,
                        false,1,checkOverlaps);
    G4VPhysicalVolume* logicCathode2 = new G4PVPlacement(0,G4ThreeVector(-cathodePosX,cathodePosY,cathodePosZ),fWirePlane1Logical,   
                        "chamber1Physical",worldLogical,
                        false,2,checkOverlaps);
    G4VPhysicalVolume* logicCathode3 = new G4PVPlacement(0,G4ThreeVector(cathodePosX,-cathodePosY,cathodePosZ),fWirePlane1Logical,   
                        "chamber1Physical",worldLogical,
                        false,3,checkOverlaps);
    G4VPhysicalVolume* logicCathode4 = new G4PVPlacement(0,G4ThreeVector(cathodePosX,cathodePosY,cathodePosZ),fWirePlane1Logical,   
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
  opCathodeSurface->SetFinish(Polished_LUT);  
  opCathodeSurface->SetModel(unified);   
  												
  
  G4LogicalBorderSurface* cathodeSurface1 = new G4LogicalBorderSurface(
    "CathodeSurface1", logicCathode1, worldPhysical, opCathodeSurface);    
    
  G4LogicalBorderSurface* cathodeSurface2 = new G4LogicalBorderSurface(
    "CathodeSurface2", logicCathode2, worldPhysical, opCathodeSurface);    
    
  G4LogicalBorderSurface* cathodeSurface3 = new G4LogicalBorderSurface(
    "CathodeSurface3", logicCathode3, worldPhysical, opCathodeSurface);     
    
  G4LogicalBorderSurface* cathodeSurface4 = new G4LogicalBorderSurface(
    "CathodeSurface4", logicCathode4, worldPhysical, opCathodeSurface);     
    
     
  G4OpticalSurface* opticalSurfaceCathode1 = dynamic_cast<G4OpticalSurface*>(
    cathodeSurface1->GetSurface(logicCathode1, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurfaceCathode1) {
    opticalSurfaceCathode1->DumpInfo();
  } 
  
  
  G4OpticalSurface* opticalSurfaceCathode2 = dynamic_cast<G4OpticalSurface*>(
    cathodeSurface2->GetSurface(logicCathode2, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurfaceCathode2) {
    opticalSurfaceCathode2->DumpInfo();
  } 
  
  
  G4OpticalSurface* opticalSurfaceCathode3 = dynamic_cast<G4OpticalSurface*>(
    cathodeSurface3->GetSurface(logicCathode3, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurfaceCathode3) {
    opticalSurfaceCathode3->DumpInfo();
  }  


  G4OpticalSurface* opticalSurfaceCathode4 = dynamic_cast<G4OpticalSurface*>(
    cathodeSurface4->GetSurface(logicCathode4, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurfaceCathode4) {
    opticalSurfaceCathode4->DumpInfo();
  } 
    
  fWirePlane1Logical->SetVisAttributes(yellowVis);
 
  //MCP-PMT (ceramic)
  //--------------------------------------------------------------
  G4Material* GaOx = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  
  G4MaterialPropertiesTable* MPT4 = new G4MaterialPropertiesTable();
  MPT4->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMCP, nBins);

  G4cout << "GaOx G4MaterialPropertiesTable:" << G4endl;
  MPT4->DumpTable();

  GaOx->SetMaterialPropertiesTable(MPT4);
  
  
  G4ThreeVector pos7 = G4ThreeVector(0.0*mm, 0.0*mm, 21.01*mm);    //correct:21.03?  if no space between:21.01
           
  G4Box* solidShape7 =    
    new G4Box("Shape7", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*18.0*mm);
                      
  G4LogicalVolume* logicShape7 =                         
    new G4LogicalVolume(solidShape7,         //its solid
                        GaOx,          //its material
                        "Shape7");           //its name
                        
  logicShape7->SetVisAttributes(greyVis);                        
               
  G4VPhysicalVolume* logicPhys5 = new G4PVPlacement(0,                       //no rotation
                    pos7,                    //at position
                    logicShape7,             //its logical volume
                    "Shape7",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking                 
  
  G4OpticalSurface* opMCPSurface = new G4OpticalSurface("MCPSurface");
  opMCPSurface->SetType(dielectric_dielectric);       
  opMCPSurface->SetFinish(Polished_LUT);		
  opMCPSurface->SetModel(unified);			
  													
  
   G4LogicalBorderSurface* mcpSurface = new G4LogicalBorderSurface(
    "MCPSurface", logicPhys5, worldPhysical, opMCPSurface);

  G4OpticalSurface* opticalSurfaceMCP = dynamic_cast<G4OpticalSurface*>(
    mcpSurface->GetSurface(logicPhys5, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurfaceMCP) {
    opticalSurfaceMCP->DumpInfo();
  }
  
  
                                                                                           
 
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
  fWirePlane1Logical->SetSensitiveDetector(hadCalorimeter);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0DetectorConstruction::ConstructMaterials()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void FT0DetectorConstruction::DefineCommands()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
