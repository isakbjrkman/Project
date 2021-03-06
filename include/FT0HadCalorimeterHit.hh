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
/// \file FT0HadCalorimeterHit.hh
/// \brief Definition of the FT0HadCalorimeterHit class

#ifndef FT0HadCalorimeterHit_h
#define FT0HadCalorimeterHit_h 1
#include "G4ParticleDefinition.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hadron Calorimeter hit
///

class FT0HadCalorimeterHit : public G4VHit
{
  public:
    FT0HadCalorimeterHit();
    FT0HadCalorimeterHit(G4int iCol,G4int iRow);
    FT0HadCalorimeterHit(const FT0HadCalorimeterHit &right);
    virtual ~FT0HadCalorimeterHit();

    const FT0HadCalorimeterHit& operator=(const FT0HadCalorimeterHit &right);
    G4bool operator==(const FT0HadCalorimeterHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetPDG(G4int p) { fpdg = p; }
    G4int GetPDG() const { return fpdg; }
    
    void SetX(G4double x) { fx = x; }
    G4double GetX() const { return fx; }
    
    void SetY(G4double y) { fy = y; }
    G4double GetY() const { return fy; }
    
    void SetZ(G4double z) { fz = z; }
    G4double GetZ() const { return fz; }
    
    void SetPX(G4double px) { fpx = px; }
    G4double GetPX() const { return fpx; }
    
    void SetPY(G4double py) { fpy = py; }
    G4double GetPY() const { return fpy; }
    
    void SetPZ(G4double pz) { fpz = pz; }
    G4double GetPZ() const { return fpz; }
    
    void SetEvent(G4int e) { fEvent = e; }
    G4int GetEvent() const { return fEvent; }
    
    void SetDetectorID(G4int id) { fDetectorID = id; }
    G4int GetDetectorID() const { return fDetectorID; }
    
    void SetCerenkov(G4int cer) { fCerenkovCounter = cer; }
    void AddCerenkov(G4int cer) { fCerenkovCounter += cer; }
    G4int GetCerenkov() const { return fCerenkovCounter; }
    
  private:
    G4double fEdep;
    G4int fpdg;
    G4double fx;
    G4double fy;
    G4double fz;
    G4double fpx;
    G4double fpy;
    G4double fpz;
    G4int fEvent;
    G4int fDetectorID;
    G4int fCerenkovCounter;
};

using FT0HadCalorimeterHitsCollection = G4THitsCollection<FT0HadCalorimeterHit>;

extern G4ThreadLocal G4Allocator<FT0HadCalorimeterHit>* FT0HadCalorimeterHitAllocator;

inline void* FT0HadCalorimeterHit::operator new(size_t)
{
  if (!FT0HadCalorimeterHitAllocator) {
       FT0HadCalorimeterHitAllocator = new G4Allocator<FT0HadCalorimeterHit>;
  }
  return (void*)FT0HadCalorimeterHitAllocator->MallocSingle();
}

inline void FT0HadCalorimeterHit::operator delete(void* aHit)
{
  FT0HadCalorimeterHitAllocator->FreeSingle((FT0HadCalorimeterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
