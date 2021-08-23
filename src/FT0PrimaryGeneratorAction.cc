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
/// \file FT0PrimaryGeneratorAction.cc
/// \brief Implementation of the FT0PrimaryGeneratorAction class

#include <cstdlib> 
#include <cmath>

#include "FT0Constants.hh"
#include "FT0Run.hh"
#include "FT0EventAction.hh"
#include "FT0PrimaryGeneratorAction.hh"
#include "FT0HadCalorimeterHit.hh"

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "g4analysis.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0PrimaryGeneratorAction::FT0PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
 , fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun  = new G4ParticleGun(nofParticles);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0PrimaryGeneratorAction::~FT0PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  auto particleTable = G4ParticleTable::GetParticleTable();
fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
  
  //  fParticleGun->SetParticlePosition(G4ThreeVector(0,74.30*mm,3800*mm));
  
  fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu+"));
  fParticleGun->SetParticleEnergy(1.*GeV);
  
  G4double x = (-26.50-dist/2)*mm + ((double) rand()/RAND_MAX)*(53.00+dist)*mm;
  G4double y = (74.30-26.50-dist/2)*mm + ((double) rand()/RAND_MAX)*(53.00+dist)*mm;
  G4double z = 3332.945*mm;
  
 fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.01,0.01,-1));
  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0PrimaryGeneratorAction::DefineCommands()
{}     

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
