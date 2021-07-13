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
/// \file B5PrimaryGeneratorAction.cc
/// \brief Implementation of the B5PrimaryGeneratorAction class

#include "B5Run.hh"

#include "B5EventAction.hh"

#include "B5HadCalorimeterHit.hh"
#include "B5Constants.hh"

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "g4analysis.hh"

#include "B5PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5PrimaryGeneratorAction::B5PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
 , fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun  = new G4ParticleGun(nofParticles);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5PrimaryGeneratorAction::~B5PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  auto particleTable = G4ParticleTable::GetParticleTable();
  fParticleGun->SetParticlePosition(G4ThreeVector(-1.*cm,-1.*cm,-12.*cm));
  fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu+"));
  fParticleGun->SetParticleEnergy(1.*GeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  fParticleGun->GeneratePrimaryVertex(event);
 
  //auto angle = (G4UniformRand()-0.5)*fSigmaAngle;
  //fParticleGun->SetParticleMomentumDirection(
  //              G4ThreeVector(std::sin(angle),0.,std::cos(angle)));
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5PrimaryGeneratorAction::DefineCommands()
{}     

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
