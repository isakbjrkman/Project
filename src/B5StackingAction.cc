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
/// \file OpNovice/src/OpNoviceStackingAction.cc
/// \brief Implementation of the OpNoviceStackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B5HadCalorimeterSD.hh"
#include "B5HadCalorimeterHit.hh"
#include "B5StackingAction.hh"
#include "B5Run.hh"   
#include "B5EventAction.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4VHitsCollection.hh"
#include "G4SystemOfUnits.hh"
#include "g4analysis.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"   
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5StackingAction::B5StackingAction()
  : G4UserStackingAction()
  , fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5StackingAction::~B5StackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack B5StackingAction::ClassifyNewTrack(
  const G4Track* aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  {  // particle is optical photon      
	if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        ++fCerenkovCounter;        
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5StackingAction::NewStage()
{
  // G4cout << "Number of Cerenkov photons produced in this event : "
  //       << fCerenkovCounter << G4endl;

  B5Run* run = static_cast<B5Run*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddCerenkov((G4double) fCerenkovCounter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5StackingAction::PrepareNewEvent()
{
  fCerenkovCounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
