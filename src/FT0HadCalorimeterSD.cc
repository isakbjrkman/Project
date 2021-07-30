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
/// \file FT0HadCalorimeterSD.cc
/// \brief Implementation of the FT0HadCalorimeterSD class
#include <iostream>
#include <fstream>

#include <cstdlib> 
#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <stdexcept>
#include <typeinfo>
#include <sstream>
#include <limits>

#include "FT0HadCalorimeterSD.hh"
#include "FT0HadCalorimeterHit.hh"
#include "FT0EventAction.hh"


#include <G4LinInterpolation.hh>
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4VHitsCollection.hh"
#include "G4SystemOfUnits.hh"
#include "g4analysis.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterSD::FT0HadCalorimeterSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1), fCerenkovCounter(0)
{
  collectionName.insert("HadCalorimeterColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FT0HadCalorimeterSD::~FT0HadCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FT0HadCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new FT0HadCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
  
  // fill calorimeter hits with zero energy deposition
  for (auto column=0;column<2;column++) {
    for (auto row=0;row<2;row++) {
      fHitsCollection->insert(new FT0HadCalorimeterHit());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FT0HadCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
int nbins = 354; 
  
 double photonEnergy[354] = { 2.38882e-06*MeV,2.39344e-06*MeV,2.39807e-06*MeV,2.40271e-06*MeV,2.40738e-06*MeV,2.41206e-06*MeV,2.41676e-06*MeV,2.42148e-06*MeV,2.42622e-06*MeV,2.43098e-06*MeV,2.43576e-06*MeV,2.44055e-06*MeV,2.44536e-06*MeV,2.4502e-06*MeV,2.45505e-06*MeV,2.45992e-06*MeV,2.46481e-06*MeV,2.46972e-06*MeV,2.47465e-06*MeV,2.4796e-06*MeV,2.48457e-06*MeV,2.48956e-06*MeV,2.49457e-06*MeV,2.4996e-06*MeV,2.50465e-06*MeV,2.50972e-06*MeV,2.51481e-06*MeV,2.51992e-06*MeV,2.52505e-06*MeV,2.5302e-06*MeV,2.53538e-06*MeV,2.54057e-06*MeV,2.54579e-06*MeV,2.55103e-06*MeV,2.55629e-06*MeV,2.56157e-06*MeV,2.56687e-06*MeV,2.5722e-06*MeV,2.57755e-06*MeV,2.58292e-06*MeV,2.58831e-06*MeV,2.59372e-06*MeV,2.59916e-06*MeV,2.60462e-06*MeV,2.61011e-06*MeV,2.61561e-06*MeV,2.62114e-06*MeV,2.62669e-06*MeV,2.63227e-06*MeV,2.63787e-06*MeV,2.6435e-06*MeV,2.64915e-06*MeV,2.65482e-06*MeV,2.66052e-06*MeV,2.66624e-06*MeV,2.67198e-06*MeV,2.67775e-06*MeV,2.68355e-06*MeV,2.68937e-06*MeV,2.69522e-06*MeV,2.70109e-06*MeV,2.70699e-06*MeV,2.71291e-06*MeV,2.71886e-06*MeV,2.72484e-06*MeV,2.73084e-06*MeV,2.73687e-06*MeV,2.74292e-06*MeV,2.749e-06*MeV,2.75511e-06*MeV,2.76125e-06*MeV,2.76741e-06*MeV,2.7736e-06*MeV,2.77982e-06*MeV,2.78607e-06*MeV,2.79234e-06*MeV,2.79865e-06*MeV,2.80498e-06*MeV,2.81134e-06*MeV,2.81773e-06*MeV,2.82415e-06*MeV,2.83059e-06*MeV,2.83707e-06*MeV,2.84358e-06*MeV,2.85011e-06*MeV,2.85668e-06*MeV,2.86328e-06*MeV,2.86991e-06*MeV,2.87657e-06*MeV,2.88326e-06*MeV,2.88998e-06*MeV,2.89673e-06*MeV,2.90351e-06*MeV,2.91033e-06*MeV,2.91718e-06*MeV,2.92406e-06*MeV,2.93097e-06*MeV,2.93791e-06*MeV,2.94489e-06*MeV,2.9519e-06*MeV,2.95895e-06*MeV,2.96603e-06*MeV,2.97314e-06*MeV,2.98029e-06*MeV,2.98747e-06*MeV,2.99469e-06*MeV,3.00194e-06*MeV,3.00922e-06*MeV,3.01655e-06*MeV,3.0239e-06*MeV,3.0313e-06*MeV,3.03873e-06*MeV,3.04619e-06*MeV,3.05369e-06*MeV,3.06123e-06*MeV,3.06881e-06*MeV,3.07643e-06*MeV,3.08408e-06*MeV,3.09177e-06*MeV,3.0995e-06*MeV,3.10727e-06*MeV,3.11508e-06*MeV,3.12292e-06*MeV,3.13081e-06*MeV,3.13873e-06*MeV,3.1467e-06*MeV,3.15471e-06*MeV,3.16276e-06*MeV,3.17084e-06*MeV,3.17897e-06*MeV,3.18715e-06*MeV,3.19536e-06*MeV,3.20362e-06*MeV,3.21192e-06*MeV,3.22026e-06*MeV,3.22865e-06*MeV,3.23708e-06*MeV,3.24555e-06*MeV,3.25407e-06*MeV,3.26263e-06*MeV,3.27124e-06*MeV,3.27989e-06*MeV,3.28859e-06*MeV,3.29734e-06*MeV,3.30613e-06*MeV,3.31497e-06*MeV,3.32386e-06*MeV,3.3328e-06*MeV,3.34178e-06*MeV,3.35081e-06*MeV,3.35989e-06*MeV,3.36902e-06*MeV,3.3782e-06*MeV,3.38743e-06*MeV,3.39671e-06*MeV,3.40604e-06*MeV,3.41543e-06*MeV,3.42486e-06*MeV,3.43435e-06*MeV,3.44389e-06*MeV,3.45348e-06*MeV,3.46313e-06*MeV,3.47283e-06*MeV,3.48258e-06*MeV,3.49239e-06*MeV,3.50226e-06*MeV,3.51218e-06*MeV,3.52216e-06*MeV,3.53219e-06*MeV,3.54229e-06*MeV,3.55244e-06*MeV,3.56264e-06*MeV,3.57291e-06*MeV,3.58324e-06*MeV,3.59362e-06*MeV,3.60407e-06*MeV,3.61458e-06*MeV,3.62515e-06*MeV,3.63578e-06*MeV,3.64647e-06*MeV,3.65723e-06*MeV,3.66805e-06*MeV,3.67893e-06*MeV,3.68988e-06*MeV,3.7009e-06*MeV,3.71198e-06*MeV,3.72312e-06*MeV,3.73434e-06*MeV,3.74562e-06*MeV,3.75697e-06*MeV,3.76839e-06*MeV,3.77988e-06*MeV,3.79144e-06*MeV,3.80307e-06*MeV,3.81477e-06*MeV,3.82654e-06*MeV,3.83839e-06*MeV,3.85031e-06*MeV,3.86231e-06*MeV,3.87437e-06*MeV,3.88652e-06*MeV,3.89874e-06*MeV,3.91104e-06*MeV,3.92342e-06*MeV,3.93587e-06*MeV,3.94841e-06*MeV,3.96102e-06*MeV,3.97372e-06*MeV,3.9865e-06*MeV,3.99935e-06*MeV,4.0123e-06*MeV,4.02532e-06*MeV,4.03844e-06*MeV,4.05163e-06*MeV,4.06492e-06*MeV,4.07829e-06*MeV,4.09175e-06*MeV,4.1053e-06*MeV,4.11894e-06*MeV,4.13267e-06*MeV,4.14649e-06*MeV,4.1604e-06*MeV,4.17441e-06*MeV,4.18851e-06*MeV,4.20271e-06*MeV,4.21701e-06*MeV,4.2314e-06*MeV,4.24589e-06*MeV,4.26048e-06*MeV,4.27517e-06*MeV,4.28997e-06*MeV,4.30486e-06*MeV,4.31986e-06*MeV,4.33497e-06*MeV,4.35018e-06*MeV,4.36549e-06*MeV,4.38092e-06*MeV,4.39645e-06*MeV,4.4121e-06*MeV,4.42786e-06*MeV,4.44373e-06*MeV,4.45971e-06*MeV,4.47581e-06*MeV,4.49203e-06*MeV,4.50836e-06*MeV,4.52482e-06*MeV,4.54139e-06*MeV,4.55809e-06*MeV,4.57491e-06*MeV,4.59185e-06*MeV,4.60892e-06*MeV,4.62612e-06*MeV,4.64345e-06*MeV,4.6609e-06*MeV,4.67849e-06*MeV,4.69621e-06*MeV,4.71407e-06*MeV,4.73206e-06*MeV,4.75019e-06*MeV,4.76846e-06*MeV,4.78687e-06*MeV,4.80543e-06*MeV,4.82412e-06*MeV,4.84297e-06*MeV,4.86196e-06*MeV,4.8811e-06*MeV,4.9004e-06*MeV,4.91984e-06*MeV,4.93944e-06*MeV,4.9592e-06*MeV,4.97912e-06*MeV,4.99919e-06*MeV,5.01943e-06*MeV,5.03984e-06*MeV,5.06041e-06*MeV,5.08115e-06*MeV,5.10206e-06*MeV,5.12314e-06*MeV,5.1444e-06*MeV,5.16583e-06*MeV,5.18745e-06*MeV,5.20924e-06*MeV,5.23122e-06*MeV,5.25339e-06*MeV,5.27574e-06*MeV,5.29829e-06*MeV,5.32103e-06*MeV,5.34397e-06*MeV,5.3671e-06*MeV,5.39043e-06*MeV,5.41397e-06*MeV,5.43772e-06*MeV,5.46167e-06*MeV,5.48584e-06*MeV,5.51022e-06*MeV,5.53482e-06*MeV,5.55964e-06*MeV,5.58468e-06*MeV,5.60995e-06*MeV,5.63545e-06*MeV,5.66119e-06*MeV,5.68716e-06*MeV,5.71336e-06*MeV,5.73981e-06*MeV,5.76651e-06*MeV,5.79346e-06*MeV,5.82066e-06*MeV,5.84811e-06*MeV,5.87583e-06*MeV,5.90381e-06*MeV,5.93206e-06*MeV,5.96058e-06*MeV,5.98937e-06*MeV,6.01845e-06*MeV,6.0478e-06*MeV,6.07745e-06*MeV,6.10739e-06*MeV,6.13762e-06*MeV,6.16816e-06*MeV,6.199e-06*MeV,6.23015e-06*MeV,6.26162e-06*MeV,6.2934e-06*MeV,6.32551e-06*MeV,6.35795e-06*MeV,6.39072e-06*MeV,6.42383e-06*MeV,6.45729e-06*MeV,6.4911e-06*MeV,6.52526e-06*MeV,6.55979e-06*MeV,6.59468e-06*MeV,6.62995e-06*MeV,6.66559e-06*MeV,6.70162e-06*MeV,6.73804e-06*MeV,6.77486e-06*MeV,6.81209e-06*MeV,6.84972e-06*MeV,6.88778e-06*MeV,6.92626e-06*MeV,6.96517e-06*MeV,7.00452e-06*MeV,7.04432e-06*MeV,7.08457e-06*MeV,7.12529e-06*MeV,7.16647e-06*MeV,7.20814e-06*MeV,7.25029e-06*MeV,7.29294e-06*MeV,7.33609e-06*MeV,7.37976e-06*MeV,7.42395e-06*MeV,7.46867e-06*MeV};
  
double eff[354] = {0.0966295,0.0998724,0.103128,0.106396,0.109677,0.11218,0.114692,0.117215,0.119747,0.12229,0.123772,
0.12526,0.126754,0.128254,0.12976,0.130429,0.1311,0.131774,0.132451,0.13313,0.133825,0.134521,0.135221,0.135924,0.136629,0.13674,0.136851,0.136963,0.137075,0.137188,
0.137946,0.138708,0.139472,0.14024,0.141011,0.143063,0.145124,0.147193,0.149271,0.151357,0.152411,0.153469,0.154532,0.155599,0.15667,0.158349,0.160034,0.161727,0.163427,
0.165135,0.165215,0.165296,0.165377,0.165458,0.16554,0.168134,0.170739,0.173355,0.175983,0.178622,0.180839,0.183065,0.1853,0.187546,0.189801,0.190908,0.19202,0.193137,
0.194259,0.195385,0.197107,0.198836,0.200572,0.202317,0.204069,0.205326,0.206588,0.207856,0.20913,0.21041,0.210022,0.209633,0.209241,0.208848,0.208453,0.209131,0.209811,
0.210495,0.211182,0.211873,0.213663,0.215462,0.217269,0.219085,0.220909,0.221313,0.221719,0.222127,0.222537,0.222949,0.224175,0.225406,0.226644,0.227887,0.229137,0.229239,
0.229342,0.229446,0.22955,0.229654,0.230149,0.230646,0.231146,0.231649,0.232153,0.233597,0.235049,0.236507,0.237973,0.239446,0.239794,0.240144,0.240496,0.240849,0.241205,
0.241242,0.241279,0.241317,0.241355,0.241393,0.24078,0.240164,0.239545,0.238923,0.238297,0.238017,0.237736,0.237453,0.237168,0.236882,0.236787,0.236691,0.236595,0.236498,
0.236401,0.236576,0.236752,0.236929,0.237106,0.237285,0.237076,0.236867,0.236656,0.236444,0.236231,0.236046,0.23586,0.235673,0.235485,0.235296,0.234754,0.234209,0.233661,
0.23311,0.232556,0.232441,0.232326,0.23221,0.232093,0.231976,0.232028,0.232081,0.232134,0.232188,0.232241,0.232159,0.232076,0.231993,0.231909,0.231825,0.231847,0.23187,
0.231893,0.231916,0.231939,0.232056,0.232174,0.232292,0.232412,0.232532,0.232449,0.232366,0.232282,0.232198,0.232114,0.231788,0.23146,0.23113,0.230798,0.230463,0.229834,
0.229201,0.228564,0.227922,0.227277,0.226647,0.226014,0.225376,0.224734,0.224088,0.223465,0.222838,0.222206,0.221571,0.220931,0.219714,0.218488,0.217255,0.216013,0.214764,
0.213563,0.212355,0.211139,0.209914,0.208681,0.207758,0.206829,0.205894,0.204952,0.204003,0.202881,0.201751,0.200613,0.199468,0.198314,0.197222,0.196122,0.195014,0.193898,
0.192774,0.19185,0.19092,0.189982,0.189038,0.188087,0.186929,0.185763,0.184588,0.183405,0.182212,0.181477,0.180735,0.179989,0.179236,0.178478,0.177884,0.177285,0.176682,
0.176074,0.175461,0.174629,0.173789,0.172944,0.172091,0.171232,0.170733,0.17023,0.169723,0.169212,0.168697,0.168375,0.168051,0.167725,0.167395,0.167063,0.166787,0.166509,
0.166228,0.165945,0.16566,0.165735,0.16581,0.165886,0.165962,0.166039,0.16609,0.166142,0.166194,0.166247,0.1663,0.166164,0.166027,0.165889,0.16575,0.165609,0.165417,0.165223,
0.165027,0.16483,0.164631,0.164954,0.165281,0.16561,0.165943,0.166278,0.167092,0.167914,0.168743,0.169581,0.170426,0.171972,0.173534,0.17511,0.176702,0.178309,0.180543,0.182799,
0.185077,0.187379,0.189703,0.192684,0.195696,0.198738,0.201811,0.204916,0.208053,0.211703,0.215597,0.219533,0.22351,0.227529,0.231591,0.233567,0.235325,0.237103,0.2389,0.240717,
0.241684,0.240623,0.239549,0.238464,0.237367,0.236257,0.235135,0.234,0.225402,0.216705,0.207907,0.199006,0.19,0.165996,0.141707,0.117127,0.0922502};  

  
  auto analysisManager = G4AnalysisManager::Instance();

  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto touchable = step->GetPreStepPoint()->GetTouchable(); //step->GetPreStepPoint()->GetTouchable(); 
  auto hit = (*fHitsCollection)[0];
  auto encoding = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  auto track = step->GetTrack();
  G4int runManager = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  std::ofstream myfile;
  myfile.open("FT0output.txt", std::ofstream::app);
  
    if (track->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    hit->AddCerenkov(1);

	//Quantum efficiency property added
	double momentX = step->GetTrack()->GetMomentum()(0);
	double momentY = step->GetTrack()->GetMomentum()(1);
	double momentZ = step->GetTrack()->GetMomentum()(2);
        double energy = abs(momentX) + abs(momentY) + abs(momentZ);
        int i = 0;
        /*
    while ( i<(nbins-1) && energy > photonEnergy[i] ) {
     i++;          
    }
    
     if ((float) rand()/RAND_MAX > eff[i]) {
        return true;     //kill photon e.g. don't register photon        
     } 
     */
     //Store data to .txt	
      myfile << runManager << "_";
      myfile << touchable->GetVolume(0)->GetCopyNo() << "_";   
      myfile << encoding << "_";
      myfile << momentX << "_";
      myfile << momentY << "_";
      myfile << momentZ << "_";
      myfile << step->GetTrack()->GetPosition()(0) << "_";
      myfile << step->GetTrack()->GetPosition()(1) << "_";
      myfile << step->GetTrack()->GetPosition()(2) << "\n"; 
        
      
     //Process hit data to root file
      hit->SetEvent(runManager);
      hit->SetDetectorID(touchable->GetVolume(0)->GetCopyNo()); 
      hit->SetPDG(encoding);  
      hit->SetPX(step->GetTrack()->GetMomentum()(0));                       
      hit->SetPY(step->GetTrack()->GetMomentum()(1));
      hit->SetPZ(step->GetTrack()->GetMomentum()(2));      
      hit->SetX(step->GetTrack()->GetPosition()(0));                       
      hit->SetY(step->GetTrack()->GetPosition()(1));
      hit->SetZ(step->GetTrack()->GetPosition()(2));
          
     //Write to root file
      analysisManager->FillNtupleDColumn(0, hit->GetEvent());
      analysisManager->FillNtupleDColumn(1, hit->GetDetectorID());
      analysisManager->FillNtupleDColumn(2, hit->GetPDG());
      analysisManager->FillNtupleDColumn(3, hit->GetPX());
      analysisManager->FillNtupleDColumn(4, hit->GetPY());
      analysisManager->FillNtupleDColumn(5, hit->GetPZ());
      analysisManager->FillNtupleDColumn(6, hit->GetX());
      analysisManager->FillNtupleDColumn(7, hit->GetY());
      analysisManager->FillNtupleDColumn(8, hit->GetZ());
      analysisManager->AddNtupleRow();               
    
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
