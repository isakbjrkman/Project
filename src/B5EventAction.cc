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
/// \file B5EventAction.cc
/// \brief Implementation of the B5EventAction class
/// \brief Implementation of the B5EventAction class
#include <iostream>
#include <fstream>

#include "B5Run.hh"

#include "B5EventAction.hh"

#include "B5HadCalorimeterHit.hh"
#include "B5Constants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "g4analysis.hh"

using std::array;
using std::vector;


namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found 


G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()",
                  "B5Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl; 
    G4Exception("B5EventAction::EndOfEventAction()",
                "B5Code001", JustWarning, msg);
  }
  return hc;  
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::B5EventAction()
: G4UserEventAction(), 


  fCalHCID  {{ -1, -1 }},
  fDriftHistoID{{ {{ -1, -1 }}, {{ -1, -1 }} }},
  fCalEdep{{ vector<G4double>(kNofEmCells, 0.), vector<G4double>(kNofHadCells, 0.) }}
      // std::array<T, N> is an aggregate that contains a C array. 
    // To initialize it, we need outer braces for the class itself 
      // and inner braces for the C array
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EventAction::~B5EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EventAction::BeginOfEventAction(const G4Event*)
{
    auto sdManager = G4SDManager::GetSDMpointer();
    // hits collections names    
    array<G4String, kDim> cHCName 
      = {{ "EMcalorimeter/EMcalorimeterColl", "HadCalorimeter/HadCalorimeterColl" }};

    for (G4int iDet = 0; iDet < kDim; ++iDet) {
      // hit collections IDs
      fCalHCID[iDet]   = sdManager->GetCollectionID(cHCName[iDet]);
    }

}     


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B5EventAction::EndOfEventAction(const G4Event* event)
{

 // B5Run* run = static_cast<B5Run*>(
   // G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //
  // Fill histograms & ntuple
  // 
  std::ofstream myfile;
  myfile.open("filename.txt", std::ofstream::app);
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

      
  // Em/Had Calorimeters hits
  array<G4int, kDim> totalCalHit = {{ 0, 0 }}; 
  array<G4double, kDim> totalCalEdep = {{ 0., 0. }}; 

  for (G4int iDet = 1; iDet < kDim; ++iDet) {
    auto hc = GetHC(event, fCalHCID[iDet]);
    if ( ! hc ) return;

    totalCalHit[iDet] = 0;
    totalCalEdep[iDet] = 0.;
    for (unsigned long i = 0; i < hc->GetSize(); ++i) {
      G4double edep = 0.;
      // The EM and Had calorimeter hits are of different types
     
        auto hit = static_cast<B5HadCalorimeterHit*>(hc->GetHit(i));   //check i & det to be correct!
        edep = hit->GetEdep(); //GetEdep
      
      if ( edep > 0. ) {
        totalCalHit[iDet]++;
        totalCalEdep[iDet] += edep;
      }
      fCalEdep[iDet][i] = edep;
    }
    // columns 0, 1
    analysisManager->FillNtupleDColumn(0, totalCalEdep[iDet]);	/*iDet + 0*/ 
  }

 
    auto hc = GetHC(event, fCalHCID[1]);
    if ( ! hc ) {
    return;
   } else {   
   
      auto hit = static_cast<B5HadCalorimeterHit*>(hc->GetHit(0));
      // columns 6->10
      //auto hceID = event->GetHCofThisEvent()->GetHC(collId);
      hit->SetEvent(event->GetEventID());
      analysisManager->FillNtupleDColumn(1 + 5, hit->GetPX());
      G4cout << hit->GetPX() << G4endl;
      analysisManager->FillNtupleDColumn(1 + 6, hit->GetPY());
      G4cout << hit->GetPY() << G4endl;
      analysisManager->FillNtupleDColumn(1 + 7, hit->GetPZ());
      G4cout << hit->GetPZ() << G4endl;
      analysisManager->FillNtupleDColumn(1 + 8, hit->GetDetectorID());
      G4cout << hit->GetDetectorID() << G4endl;
      myfile << hit->GetDetectorID() << "_";
      myfile << "Cerenkov" << hit->GetCerenkov() << "Cerenkov";
        // HadCalorimeter hits
  for (G4int iDet = 1; iDet < kDim; ++iDet) {
    auto hc2 = GetHC(event, fCalHCID[iDet]);
    if ( ! hc2 ) return;

    for (unsigned int i = 0; i<1; ++i) {
      auto hit2 = static_cast<B5HadCalorimeterHit*>(hc2->GetHit(i));
      // columns 2
      analysisManager->FillNtupleDColumn(iDet + 1, hit2->GetPDG());
      G4cout << hit2->GetPDG() << G4endl;
      myfile << hit2->GetPDG() << "_";
    }
  }      
      myfile << hit->GetPX() << "_";
      myfile << hit->GetPY() << "_";
      myfile << hit->GetPZ() << "_";       
  }
  
    auto hc3 = GetHC(event, fCalHCID[1]);
    if ( ! hc3 ) return;

    for (unsigned int i = 0; i<1; ++i) {
      auto hit = static_cast<B5HadCalorimeterHit*>(hc3->GetHit(i));
      // columns 3,4,5
   
      analysisManager->FillNtupleDColumn(1 + 2, hit->GetX());
      G4cout << hit->GetX() << G4endl;
      analysisManager->FillNtupleDColumn(1 + 3, hit->GetY());
      G4cout << hit->GetY() << G4endl;
      analysisManager->FillNtupleDColumn(1 + 4, hit->GetZ());
      G4cout << hit->GetZ() << G4endl;
      myfile << hit->GetX() << "_";
      myfile << hit->GetY() << "_";
      myfile << hit->GetZ() << "\n";
    }
    
   myfile.close();
  


  analysisManager->AddNtupleRow();


  //
  // Print diagnostics
  // 
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( printModulo == 0 || event->GetEventID() % printModulo != 0) return;
  auto primary = event->GetPrimaryVertex(0)->GetPrimary(0);
  G4cout 
    << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " " << primary->GetMomentum() << G4endl;


  // Calorimeters

    G4cout << "Hadron Calorimeter has " << totalCalHit[1] << " hits." 			
           << " Total Edep is " << totalCalEdep[1]/MeV << " (MeV)" << G4endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
