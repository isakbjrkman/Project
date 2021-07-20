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
/// \file FT0RunAction.cc
/// \brief Implementation of the FT0RunAction class


#include "FT0RunAction.hh"
#include "FT0EventAction.hh"
#include "FT0Run.hh"
#include "G4ParticleDefinition.hh"
#include "FT0PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericAnalysisManager.hh"

using G4AnalysisManager = G4GenericAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0RunAction::FT0RunAction(FT0PrimaryGeneratorAction* prim)    
 : G4UserRunAction(),
   fRun(nullptr),
   fPrimary(prim)
{ 

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);

  // Default settings
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("FT0");
  // Book ntuple

  if ( fPrimary ) {
    analysisManager->CreateNtuple("FT0", "Hits");
    analysisManager->CreateNtupleDColumn("EventID");  // column Id = 0
    analysisManager->CreateNtupleDColumn("DetectorID");  // column Id = 1
    analysisManager->CreateNtupleDColumn("PDG");      // column Id = 2
    analysisManager->CreateNtupleDColumn("PX");      // column Id = 3
    analysisManager->CreateNtupleDColumn("PY");      // column Id = 4
    analysisManager->CreateNtupleDColumn("PZ");      // column Id = 5
    analysisManager->CreateNtupleDColumn("X");      // column Id = 6
    analysisManager->CreateNtupleDColumn("Y");      // column Id = 7
    analysisManager->CreateNtupleDColumn("Z");      // column Id = 8
    analysisManager->FinishNtuple();
  }

  // Set ntuple output file
  analysisManager->SetNtupleFileName(0, "FT0ntuple");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0RunAction::~FT0RunAction()
{
  delete G4AnalysisManager::Instance();  
}

G4Run* FT0RunAction::GenerateRun()
{
  fRun = new FT0Run();
  return fRun;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 

  if(fPrimary)
  {
    G4ParticleDefinition* particle =
    fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Open an output file 
   analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
 auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  if(isMaster)
    fRun->EndOfRun();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
