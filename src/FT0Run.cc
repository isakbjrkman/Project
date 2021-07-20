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
/// \file OpNovice/src/OpNoviceRun.cc
/// \brief Implementation of the OpNoviceRun class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include <iostream>
#include <fstream>
#include "FT0Run.hh"

#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0Run::FT0Run()
  : G4Run()
{
  fParticle             = nullptr;
  fEnergy               = -1.;
  fCerenkovCounter      = 0.;
  fCerenkov2            = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0Run::~FT0Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEnergy   = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0Run::Merge(const G4Run* run)
{
  const FT0Run* localRun = static_cast<const FT0Run*>(run);

  fParticle = localRun->fParticle;
  fEnergy   = localRun->fEnergy;
  fCerenkovCounter += localRun->fCerenkovCounter;
  fCerenkov2 += localRun->fCerenkov2;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FT0Run::EndOfRun()
{
  if(numberOfEvent == 0)
    return;
  G4double TotNbofEvents = G4double(numberOfEvent);
  fCerenkovCounter /= TotNbofEvents;
  fCerenkov2 /= TotNbofEvents;
  G4double rmsCerenkov = fCerenkov2 - fCerenkovCounter * fCerenkovCounter;
  if(rmsCerenkov > 0.)
    rmsCerenkov = std::sqrt(rmsCerenkov);
  else
    rmsCerenkov = 0.;


  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";

  G4cout << "Primary particle was: " << fParticle->GetParticleName()
         << " with energy " << G4BestUnit(fEnergy, "Energy") << "." << G4endl;
  G4cout << "Number of events: " << numberOfEvent << G4endl;

  G4cout << "Average number of Cerenkov photons created per event: "
         << fCerenkovCounter << " +- " << rmsCerenkov << G4endl;


  G4cout << G4endl;
  G4cout.precision(prec);
}
