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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <array>


#include "B5DetectorConstruction.hh"
#include "B5HadCalorimeterSD.hh"

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

B5DetectorConstruction::B5DetectorConstruction()
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

B5DetectorConstruction::~B5DetectorConstruction()
{
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
  

  
 G4double mPhotonEnergyD[354];
 G4double mEfficMet[354];
 G4double mRindexMet[354];
 G4double mAbs[354];
 G4double mRefractiveIndexAir[354];
 G4double mAbsAir[354];
 G4double mReflMCP[354];
 G4double mEffMCP[354];
 G4double mAbsMCP[354];
 G4double mRefractiveIndexCathode[354];
 G4double mAbsCathode[354];
 G4double mReflMet[354];
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
    mReflMCP[i] = 0.9;   //test corr: 0.
    mEffMCP[i] = 0.;
    mAbsMCP[i] = 1.;
    mRefractiveIndexCathode[i] = 1.;
    mAbsCathode[i] = 1.;
    mReflMet[i] = 0.9;
    //mEfficMet[i] = 0.;
  }
 
  
  // Construct materials
  ConstructMaterials();
  auto nist = G4NistManager::Instance();
  
  /*
    G4double zet      = 1.0;
    G4double amass    = 1.01*g/mole;
    G4double density  = universe_mean_density;
    G4double pressure = 3.e-18*pascal;
    G4double tempture = 2.73*kelvin;
    G4Material* air = new G4Material("Vacuum", zet, amass, density,
                            kStateGas, tempture, pressure);
   
  */ 
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");                        
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
  /*
  G4OpticalSurface* opAirSurface = new G4OpticalSurface("AirSurface");
  opAirSurface->SetType(dielectric_dielectric);  //dielectric_dielectric
  opAirSurface->SetFinish(ground);  //ground
  opAirSurface->SetModel(unified); //unified
  													//needed??
  
   G4LogicalBorderSurface* airSurface = new G4LogicalBorderSurface(
    "AirSurface", worldPhysical, worldPhysical, opAirSurface);

  G4OpticalSurface* opticalSurface4 = dynamic_cast<G4OpticalSurface*>(
    airSurface->GetSurface(worldPhysical, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface4)
    opticalSurface4->DumpInfo();
*/
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
  opQuartzSurface->SetFinish(Rough_LUT);  //polishedbackpainted
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
  opQuartzSurface5->SetType(dielectric_metal);     //dielectric_metal
  opQuartzSurface5->SetFinish(Rough_LUT);  //polishedbackpainted
  opQuartzSurface5->SetModel(unified);  //unified
  
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

  G4int o = 0;  
  auto chamber1Solid
    = new G4Box("chamber1Box",0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
  fWirePlane1Logical									
    = new G4LogicalVolume(chamber1Solid,GaAr,"fWirePlane1Logical");
    
  for (auto i=0;i<2;i++) {
      G4double x1 = -13.255*mm+i*2*13.255*mm;
      for (auto j=0;j<2; j++){
      o++;   
      G4double y1 = -13.255*mm+j*2*13.255*mm;
      new G4PVPlacement(0,G4ThreeVector(x1,y1,12.005*mm),fWirePlane1Logical,    //12.015*mm correct   test:12.005
                        "chamber1Physical",worldLogical,
                        false,o,checkOverlaps);
  	}
  } 
 /*   
  G4OpticalSurface* opCathodeSurface = new G4OpticalSurface("CathodeSurface");
  opCathodeSurface->SetType(dielectric_metal);     //dielectric_metal
  opCathodeSurface->SetFinish(polishedbackpainted);  //polishedbackpainted
  opCathodeSurface->SetModel(unified);   //unified
  												
  
   G4LogicalBorderSurface* cathodeSurface = new G4LogicalBorderSurface(
    "CathodeSurface", logicPhys3, worldPhysical, opCathodeSurface);    //name has changed MODIFY reminder

  G4OpticalSurface* opticalSurface3 = dynamic_cast<G4OpticalSurface*>(
    cathodeSurface->GetSurface(logicPhys3, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface3) {
    opticalSurface3->DumpInfo();
  } */
    fWirePlane1Logical->SetVisAttributes(yellowVis);
 
  //MCP-PMT (ceramic)
  //--------------------------------------------------------------
  G4Material* GaOx = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4MaterialPropertiesTable* MPT4 = new G4MaterialPropertiesTable();
  MPT4->AddProperty("REFLECTIVITY", mPhotonEnergyD, mReflMCP, nBins);
  MPT4->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbsMCP, nBins);
  MPT4->AddProperty("EFFICIENCY", mPhotonEnergyD, mEffMCP, nBins);
  
  G4cout << "GaOx G4MaterialPropertiesTable:" << G4endl;
  MPT4->DumpTable();

  GaOx->SetMaterialPropertiesTable(MPT4);
  
  
  G4ThreeVector pos7 = G4ThreeVector(0.0*mm, 0.0*mm, 21.01*mm);    //correct:21.03  test:21.01
           
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
 /* 
  G4OpticalSurface* opMCPSurface = new G4OpticalSurface("MCPSurface");
  opMCPSurface->SetType(dielectric_dielectric);       //dielectric_dielectric
  opMCPSurface->SetFinish(groundbackpainted);		//groundbackpainted
  opMCPSurface->SetModel(unified);			//unified
  													
  
   G4LogicalBorderSurface* mcpSurface = new G4LogicalBorderSurface(
    "MCPSurface", logicPhys5, worldPhysical, opMCPSurface);

  G4OpticalSurface* opticalSurface5 = dynamic_cast<G4OpticalSurface*>(
    mcpSurface->GetSurface(logicPhys5, worldPhysical)
      ->GetSurfaceProperty());
  if(opticalSurface5) {
    opticalSurface5->DumpInfo();
  }
  */
  
                                                                                           
 
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;  //worldPhysical
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  auto hadCalorimeter = new B5HadCalorimeterSD(SDname="/HadCalorimeter");
  sdManager->AddNewDetector(hadCalorimeter);
  fWirePlane1Logical->SetSensitiveDetector(hadCalorimeter);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B5DetectorConstruction::DefineCommands()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
