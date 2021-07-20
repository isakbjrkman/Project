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
/// \file FT0HadCalorimeterHit.cc
/// \brief Implementation of the FT0HadCalorimeterHit class
#include "G4ParticleDefinition.hh"
#include "FT0HadCalorimeterHit.hh"
#include "FT0DetectorConstruction.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<FT0HadCalorimeterHit>* FT0HadCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterHit::FT0HadCalorimeterHit()
: G4VHit(), 
  fEdep(0.), fCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
FT0HadCalorimeterHit::FT0HadCalorimeterHit(G4int columnID,G4int rowID)
: G4VHit(), 
  fEdep(0.), fCerenkovCounter(0)
{}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterHit::~FT0HadCalorimeterHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterHit::FT0HadCalorimeterHit(const FT0HadCalorimeterHit &right)
: G4VHit(),
  fEdep(right.fEdep),
  fCerenkovCounter(right.fCerenkovCounter)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const FT0HadCalorimeterHit& FT0HadCalorimeterHit::operator=(
        const FT0HadCalorimeterHit &right)
{
  fEdep = right.fEdep;
  fCerenkovCounter = right.fCerenkovCounter;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
G4bool FT0HadCalorimeterHit::operator==(const FT0HadCalorimeterHit &right) const
{
  //return (fColumnID==right.fColumnID&&fRowID==right.fRowID);
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* FT0HadCalorimeterHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("FT0HadCalorimeterHit",isNew);

  if (isNew) {
    (*store)["HitType"] 
      = G4AttDef("HitType","Hit Type","Physics","","G4String");

    (*store)["Energy"] 
      = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit",
                 "G4double");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* FT0HadCalorimeterHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  
  values
    ->push_back(G4AttValue("HitType","HadCalorimeterHit",""));
  values
    ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),"")); 
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0HadCalorimeterHit::Print()
{
  G4cout << "Cell hit \n" << G4endl;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
