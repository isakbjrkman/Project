#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <array>
#include <cstdlib> 
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>


int rootMacro() {
//enter name command
int number;
cout << "Type root file number to be read. \nIf file name is FT0ntupleX.root, type X.\nEnter number:" << endl;
cin >> number;
string numberStr = to_string(number);

//read file

string init = "FT0ntuple";
string last = ".root";
string inp = init + numberStr + last;

cout << "Reading filename: " << inp << endl;

TFile* f1 = new TFile( inp.c_str(), "READ");
if (!f1) {
cout << "File read failed" << endl;
return -1;
}
TTree* t = (TTree*)f1->Get("FT0");

//Output filename according to input number, i.e. dataProcessingX belongs to FT0ntupleX. 
string outputName = "dataProcessing";
string outp = outputName + numberStr + last;

TFile* f = new TFile( outp.c_str(), "RECREATE");
if (!f) {
cout << "File read failed" << endl;
return -1;
}

double eventID,detectorID,PDG;
double x,y,z,px,py,pz;

t->SetBranchAddress("EventID",&eventID);
t->SetBranchAddress("DetectorID",&detectorID);
t->SetBranchAddress("PDG",&PDG);
t->SetBranchAddress("PX",&px);
t->SetBranchAddress("PY",&py);
t->SetBranchAddress("PZ",&pz);
t->SetBranchAddress("X",&x);
t->SetBranchAddress("Y",&y);
t->SetBranchAddress("Z",&z);

int nentries = t->GetEntries(); 

double h = 6.62607015*pow(10.,-34);
double c = 299792458;
double E;

TH1D* h1ev = new TH1D("h1ev", "EventID; Event ID; frequency", 10000, 0, nentries/2000);
TH1D* h1de = new TH1D("h1de", "DetectorID; Detector ID; frequency", 50, 0.5, 4.5);
TH1D* h1pdg = new TH1D("h1pdg", "PDG; pdg; frequency", 100, -50, 50);
TH1D* h1px = new TH1D("h1px", "PX; px (MeV); frequency", 100, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1py = new TH1D("h1py", "PY; py (MeV); frequency", 100, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1pz = new TH1D("h1pz", "PZ; pz (MeV); frequency", 100, 0, 8*pow(10., -6));
TH1D* h1x = new TH1D("h1x", "X; x (mm); frequency", 200, -30, 30);
TH1D* h1y = new TH1D("h1y", "Y; y (mm); frequency", 200, -30, 30);
TH1D* h1z = new TH1D("h1z", "Z; z (mm); frequency", 200, 11.99, 12.015);
TH1D* h1w = new TH1D("h1w", "Wavelength distribution; wavelength (nm); frequency", 100, 0.06*pow(10.,3), 0.53*pow(10., 3));

double convert = 1.60217662*pow(10.,-13);   //MeV to Joules

for (int i=0; i<nentries; i++) {
  t->GetEntry(i);
  E = (abs(px)+abs(py)+abs(pz))*convert; 
  h1w->Fill(pow(10,9)*h*c/E);
  
  h1ev->Fill(eventID);
  h1de->Fill(detectorID);
  h1pdg->Fill(PDG);
  h1px->Fill(px);
  h1py->Fill(py);
  h1pz->Fill(pz);
  h1x->Fill(x);
  h1y->Fill(y);
  h1z->Fill(z);
}  


h1w->Write();
h1ev->Write();
h1de->Write();
h1pdg->Write();
h1px->Write();
h1py->Write();
h1pz->Write();
h1x->Write();
h1y->Write();
h1z->Write();

f->Close();
f1->Close();

cout << "Dataprocessing completed successfully. \nOutput file " << outp << " created." << endl;
return 1;

}
