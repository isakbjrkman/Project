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
/// \file FT0PrimaryGeneratorAction.hh
/// \brief Definition of the FT0PrimaryGeneratorAction class

#ifndef FT0PrimaryGeneratorAction_h
#define FT0PrimaryGeneratorAction_h 1

#include "G4UserEventAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"

#include <vector>
#include <array>

// named constants

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

/// Primary generator
///


class FT0PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    FT0PrimaryGeneratorAction();
    virtual ~FT0PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event*);
    
    G4ParticleGun* GetParticleGun() { return fParticleGun; }
    std::vector<G4double>& GetHadCalEdep() { return fCalEdep[1]; }
    
  private:
    void DefineCommands();

    G4ParticleGun* fParticleGun;
    G4GenericMessenger* fMessenger;
    
    
    std::array<G4int, 2> fCalHCID;						//delete if not filled x3, found in eventaction class
    // energy deposit in calorimeters cells
    std::array<std::vector<G4double>, 2> fCalEdep;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
