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

#include "G4Cerenkov.hh"
#include "G4OpticalPhoton.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "B5DetectorConstruction.hh"
#include "B5HadCalorimeterSD.hh"
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
  // Construct materials
  ConstructMaterials();
  auto nist = G4NistManager::Instance();

  // Air 
  auto air = nist->FindOrBuildMaterial("G4_AIR");
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*cm, env_sizeZ = 10*cm;
  
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
               
  G4VPhysicalVolume* logicEnvPhys = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    worldLogical,            //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
 G4double mPhotonEnergyD[354];
 G4double mEfficMet[354];
 G4double mReflMet[354];
 G4double mAbs[354];
 G4double mRefractiveIndex2[354];
 int nBins = sizeof(mPhotonEnergyD)/sizeof(mPhotonEnergyD[0]);
 
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
 mReflMet[num] = ref;
 mEfficMet[num] = eff;  
 std::cout << mPhotonEnergyD[num] << " " << mAbs[num] << " " << mReflMet[num] << " " << mEfficMet[num] << "\n";
 ++num;
      } 
    
  
  for (auto u=0; u< 354; u++) {
  std::cout << mPhotonEnergyD[u] << " " << mAbs[u] << " " << mReflMet[u] << " " << mEfficMet[u] << "\n";
  }


 for (auto i = 0; i < nBins; i++) {
    mRefractiveIndex2[i] = 1.0;
  }


  //     
  // Quartz radiator 1
  //  
  G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.0855*g/mole);
  G4Element* elO = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Material* SiO2 = new G4Material("Silicon_dioxide", 2.533*g/cm3, 2);
  SiO2->AddElement(elO, 2);
  SiO2->AddElement(elSi, 1);
  
  
  G4MaterialPropertiesTable *MPT = new G4MaterialPropertiesTable();
  MPT->AddProperty("RINDEX", mPhotonEnergyD, mReflMet, nBins)->SetSpline(true);
  MPT->AddProperty("ABSLENGTH", mPhotonEnergyD, mAbs, nBins)->SetSpline(true);
  MPT->AddProperty("EFFICIENCY", mPhotonEnergyD, mEfficMet, nBins)->SetSpline(true);
  
  
  G4cout << "Quartz G4MaterialPropertiesTable:" << G4endl;
  MPT->DumpTable();  
  SiO2->SetMaterialPropertiesTable(MPT);

  //Air properties
  
  G4MaterialPropertiesTable* MPT2 = new G4MaterialPropertiesTable();
  MPT2->AddProperty("RINDEX", mPhotonEnergyD, mRefractiveIndex2, nBins);

  G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
  MPT2->DumpTable();

  air->SetMaterialPropertiesTable(MPT2);
 

  G4ThreeVector pos1 = G4ThreeVector(-13.255*mm, 13.255*mm, 0*mm);
             

  G4Box* solidShape1 =    
    new G4Box("Shape1", 
    0.5*26.5*mm, 0.5*26.5*mm, 0.5*20.0*mm);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        SiO2,          //its material
                        "Shape1");           //its name                  

 G4VPhysicalVolume* logicPhys1;   //interpreted as last or all 4 volumes? 
   G4int p = 0;
  for (auto i=0;i<2;i++) {
     G4double x1 = -13.255*mm+i*2*13.255*mm;
     for (auto j=0;j<2; j++){
     p++;
     G4double y1 = -13.255*mm+j*2*13.255*mm;
    logicPhys1 =  new G4PVPlacement(0,G4ThreeVector(x1,y1,0.0*mm),logicShape1,
                        "Shape1",logicEnv,
                        false,p,checkOverlaps);
  	}
  }
  
  
  G4OpticalSurface* opQuartzSurface = new G4OpticalSurface("QuartzSurface");
  opQuartzSurface->SetType(dielectric_LUTDAVIS);
  opQuartzSurface->SetFinish(Rough_LUT);
  opQuartzSurface->SetModel(DAVIS);
  
  
   G4LogicalBorderSurface* quartzSurface = new G4LogicalBorderSurface(
    "QuartzSurface", logicPhys1, logicEnvPhys, opQuartzSurface);

  G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
    quartzSurface->GetSurface(logicPhys1, logicEnvPhys)
      ->GetSurfaceProperty());
  if(opticalSurface)
    opticalSurface->DumpInfo();
   
  
  logicShape1->SetVisAttributes(blueVis);              
 
  G4ThreeVector pos5 = G4ThreeVector(0*mm, 0*mm, 11*mm);
      
        
  // Quartz window (monolithic)      

  G4Box* solidShape5 =    
    new G4Box("Shape5", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*2.0*mm);
                      
  G4LogicalVolume* logicShape5 =                         
    new G4LogicalVolume(solidShape5,         //its solid
                        SiO2,          //its material
                        "Shape5");           //its name
               
  logicShape5->SetVisAttributes(purpleVis);           
               
  
  G4VPhysicalVolume* logicPhys2 = new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logicShape5,             //its logical volume
                    "Shape5",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
                
  G4OpticalSurface* opQuartzSurface2 = new G4OpticalSurface("QuartzSurface2");
  opQuartzSurface2->SetType(dielectric_LUTDAVIS);
  opQuartzSurface2->SetFinish(Rough_LUT);
  opQuartzSurface2->SetModel(DAVIS);
  
  
   G4LogicalBorderSurface* quartzSurface2 = new G4LogicalBorderSurface(
    "QuartzSurface2", logicPhys2, logicEnvPhys, opQuartzSurface2);

  G4OpticalSurface* opticalSurface2 = dynamic_cast<G4OpticalSurface*>(
    quartzSurface2->GetSurface(logicPhys2, logicEnvPhys)
      ->GetSurfaceProperty());
  if(opticalSurface2)
    opticalSurface2->DumpInfo();         
          
                                       
 //Photocathode                    
  
  G4Material* GaAr = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");

        
 // first arm
  
  auto firstArmSolid 
    = new G4Box("firstArmBox",1*27.5*mm, 1*27.5*mm, 1*0.01*mm);
  auto firstArmLogical
    = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,12.015*mm),firstArmLogical,
                    "firstArmPhysical",logicEnv,
                    false,0,checkOverlaps);
  
  // hodoscopes in first arm
  G4int o = 0;
  auto chamber1Solid
    = new G4Box("chamber1Box",0.5*26.5*mm, 0.5*26.5*mm, 0.5*0.01*mm);
  auto chamber1Logical
    = new G4LogicalVolume(chamber1Solid,GaAr,"chamber1Logical");

  for (auto i=0;i<2;i++) {
      G4double x1 = -13.255*mm+i*2*13.255*mm;
      for (auto j=0;j<2; j++){
      o++;
      G4double y1 = -13.255*mm+j*2*13.255*mm;
      new G4PVPlacement(0,G4ThreeVector(x1,y1,0.0*mm),chamber1Logical,
                        "chamber1Physical",firstArmLogical,
                        false,o,checkOverlaps);
  	}
  }
    chamber1Logical->SetVisAttributes(yellowVis);  
    

  //"virtual" wire plane
  auto wirePlane1Solid 
    = new G4Box("wirePlane1Box", 0.5*26.5*mm, 0.5*26.5*mm, 0.25*0.01*mm);
  fWirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,air,"wirePlane1Logical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0),fWirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,
                    false,0,checkOverlaps);
                    
 
 
 
   
  //MCP-PMT (ceramic)
  G4Material* shape7_mat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4ThreeVector pos7 = G4ThreeVector(0.0*mm, 0.0*mm, 21.03*mm);
           
  G4Box* solidShape7 =    
    new G4Box("Shape7", 
    0.5*59.0*mm, 0.5*59.0*mm, 0.5*18.0*mm);
                      
  G4LogicalVolume* logicShape7 =                         
    new G4LogicalVolume(solidShape7,         //its solid
                        shape7_mat,          //its material
                        "Shape7");           //its name
                        
  logicShape7->SetVisAttributes(greyVis);                        
               
  new G4PVPlacement(0,                       //no rotation
                    pos7,                    //at position
                    logicShape7,             //its logical volume
                    "Shape7",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking                 
                                                                                           
 
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;
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
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
