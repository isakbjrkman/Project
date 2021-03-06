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
/// \file FT0RunAction.hh
/// \brief Definition of the FT0RunAction class

#ifndef FT0RunAction_h
#define FT0RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class FT0EventAction;
class FT0Run;
class FT0PrimaryGeneratorAction;

class G4Run;

/// Run action class

class FT0RunAction : public G4UserRunAction
{
  public:
    FT0RunAction(FT0PrimaryGeneratorAction* = nullptr);
    virtual ~FT0RunAction();

    G4Run* GenerateRun() override;
    virtual void BeginOfRunAction(const G4Run*) override;
    virtual void   EndOfRunAction(const G4Run*) override;

  private:
    FT0Run* fRun;
    FT0PrimaryGeneratorAction* fPrimary; 		
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
