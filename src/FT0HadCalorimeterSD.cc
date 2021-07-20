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
/// \file FT0HadCalorimeterSD.cc
/// \brief Implementation of the FT0HadCalorimeterSD class
#include <iostream>
#include <fstream>

#include "FT0HadCalorimeterSD.hh"
#include "FT0HadCalorimeterHit.hh"
#include "FT0EventAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4VHitsCollection.hh"
#include "G4SystemOfUnits.hh"
#include "g4analysis.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterSD::FT0HadCalorimeterSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1), fCerenkovCounter(0)
{
  collectionName.insert("HadCalorimeterColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterSD::~FT0HadCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0HadCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new FT0HadCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
  
  // fill calorimeter hits with zero energy deposition
  for (auto column=0;column<2;column++) {
    for (auto row=0;row<2;row++) {
      fHitsCollection->insert(new FT0HadCalorimeterHit());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FT0HadCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  
  auto analysisManager = G4AnalysisManager::Instance();

  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto touchable = step->GetPreStepPoint()->GetTouchable(); //step->GetPreStepPoint()->GetTouchable(); 
  auto hit = (*fHitsCollection)[0];
  auto encoding = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  auto track = step->GetTrack();
  G4int runManager = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  std::ofstream myfile;
  myfile.open("FT0output.txt", std::ofstream::app);
  
    if (track->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    hit->AddCerenkov(1);

     //Store data to .txt	
      myfile << runManager << "_";
      myfile << touchable->GetVolume(0)->GetCopyNo() << "_";   
      myfile << encoding << "_";
      myfile << step->GetTrack()->GetMomentum()(0) << "_";
      myfile << step->GetTrack()->GetMomentum()(1) << "_";
      myfile << step->GetTrack()->GetMomentum()(2) << "_";
      myfile << step->GetTrack()->GetPosition()(0) << "_";
      myfile << step->GetTrack()->GetPosition()(1) << "_";
      myfile << step->GetTrack()->GetPosition()(2) << "\n"; 
        
      
     //Process hit data to root file
      hit->SetEvent(runManager);
      hit->SetDetectorID(touchable->GetVolume(0)->GetCopyNo()); 
      hit->SetPDG(encoding);  
      hit->SetPX(step->GetTrack()->GetMomentum()(0));                       
      hit->SetPY(step->GetTrack()->GetMomentum()(1));
      hit->SetPZ(step->GetTrack()->GetMomentum()(2));      
      hit->SetX(step->GetTrack()->GetPosition()(0));                       
      hit->SetY(step->GetTrack()->GetPosition()(1));
      hit->SetZ(step->GetTrack()->GetPosition()(2));
          
     //Write to root file
      analysisManager->FillNtupleDColumn(0, hit->GetEvent());
      analysisManager->FillNtupleDColumn(1, hit->GetDetectorID());
      analysisManager->FillNtupleDColumn(2, hit->GetPDG());
      analysisManager->FillNtupleDColumn(3, hit->GetPX());
      analysisManager->FillNtupleDColumn(4, hit->GetPY());
      analysisManager->FillNtupleDColumn(5, hit->GetPZ());
      analysisManager->FillNtupleDColumn(6, hit->GetX());
      analysisManager->FillNtupleDColumn(7, hit->GetY());
      analysisManager->FillNtupleDColumn(8, hit->GetZ());
      analysisManager->AddNtupleRow();    
    }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
