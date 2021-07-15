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


  fCalHCID  {{ -1 }},
  fCalEdep{{ vector<G4double>(kNofHadCells, 0.) }}
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
    array<G4String, 1> cHCName 
      = {{ "HadCalorimeter/HadCalorimeterColl" }};

 
      // hit collections IDs
      fCalHCID[0]   = sdManager->GetCollectionID(cHCName[0]);
   

}     


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B5EventAction::EndOfEventAction(const G4Event* event)
{

  // Fill ntuple

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

      
  // Had Calorimeters hits
  array<G4int, 1> totalCalHit = {{ 0 }}; 
  array<G4double, 1> totalCalEdep = {{ 0. }}; 

  
    auto hc = GetHC(event, fCalHCID[0]);
    if ( ! hc ) return;

    totalCalHit[0] = 0;
    totalCalEdep[0] = 0.;
    for (unsigned long i = 0; i < hc->GetSize(); ++i) {
      G4double edep = 0.;
     
        auto hit = static_cast<B5HadCalorimeterHit*>(hc->GetHit(i));   //check i & det to be correct!
        edep = hit->GetEdep(); //GetEdep
      
      if ( edep > 0. ) {
        totalCalHit[0]++;
        totalCalEdep[0] += edep;
      }
      fCalEdep[0][i] = edep;
    }
    // columns 0, 1
    analysisManager->FillNtupleDColumn(0, totalCalEdep[0]);	/*iDet + 0*/ 
  
    if ( ! hc ) {
    return;
   } else {   
   
      auto hit = static_cast<B5HadCalorimeterHit*>(hc->GetHit(0));
      // columns 6->10
      //auto hceID = event->GetHCofThisEvent()->GetHC(collId);
      hit->SetEvent(event->GetEventID());
      analysisManager->FillNtupleDColumn(5, hit->GetPX());
      G4cout << hit->GetPX() << G4endl;
      analysisManager->FillNtupleDColumn(6, hit->GetPY());
      G4cout << hit->GetPY() << G4endl;
      analysisManager->FillNtupleDColumn(7, hit->GetPZ());
      G4cout << hit->GetPZ() << G4endl;
      analysisManager->FillNtupleDColumn(8, hit->GetDetectorID());
      G4cout << hit->GetDetectorID() << G4endl; 
        // HadCalorimeter hits

    auto hc2 = GetHC(event, fCalHCID[0]);
    if ( ! hc2 ) return;

    
      auto hit2 = static_cast<B5HadCalorimeterHit*>(hc2->GetHit(0));
      // columns 2
      analysisManager->FillNtupleDColumn(1, hit2->GetPDG());
      G4cout << hit2->GetPDG() << G4endl;
      
  }
  
    auto hc3 = GetHC(event, fCalHCID[0]);
    if ( ! hc3 ) return;


      auto hit = static_cast<B5HadCalorimeterHit*>(hc3->GetHit(0));
      // columns 3,4,5
   
      analysisManager->FillNtupleDColumn(2, hit->GetX());
      G4cout << hit->GetX() << G4endl;
      analysisManager->FillNtupleDColumn(3, hit->GetY());
      G4cout << hit->GetY() << G4endl;
      analysisManager->FillNtupleDColumn(4, hit->GetZ());
      G4cout << hit->GetZ() << G4endl;
 
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

    G4cout << "Hadron Calorimeter has " << totalCalHit[0] << " hits." 			
           << " Total Edep is " << totalCalEdep[0]/MeV << " (MeV)" << G4endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
