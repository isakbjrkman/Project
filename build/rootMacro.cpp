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

int help(int number) {
//enter name command
//int number;
//cout << "Type root file number to be read. \nIf file name is FT0ntupleX.root, type X.\nEnter number:" << endl;
//cin >> number;
string numberStr = to_string(number);

//read file

string init = "FT0ntuple";
string last = ".root";
string inp = "ntupleSave/" + init + numberStr + last;

cout << "Reading filename: " << inp << endl;

TFile* f1 = new TFile( inp.c_str(), "READ");
if (!f1) {
cout << "File read failed" << endl;
return -1;
}
TTree* t = (TTree*)f1->Get("FT0");

//Output filename according to input number, i.e. dataProcessingX belongs to FT0ntupleX. 
string outputName = "dataProcessing";
string outp = outputName + "/" + outputName + numberStr + last;
cout << outp << endl;

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

TH1D* h1ev = new TH1D("h1ev", "EventID; Event ID; Frequency", 10000, 0, nentries/320);
TH1D* h1de = new TH1D("h1de", "DetectorID; Detector ID; Frequency", 50, 0.5, 4.5);
TH1D* h1pdg = new TH1D("h1pdg", "PDG; pdg; Frequency", 100, -50, 50);
TH1D* h1px = new TH1D("h1px", "PX; px (MeV); Frequency", 200, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1py = new TH1D("h1py", "PY; py (MeV); Frequency", 200, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1pz = new TH1D("h1pz", "PZ; pz (MeV); Frequency", 200, 0, 8*pow(10., -6));
TH1D* h1x = new TH1D("h1x", "X; x (mm); Frequency", 200, -30, 30);
TH1D* h1y = new TH1D("h1y", "Y; y (mm); Frequency", 200, 47, 102);
TH1D* h1z = new TH1D("h1z", "Z; z (mm); Frequency", 200, 3354.995, 3355.015);
TH1D* h1w = new TH1D("h1w", "Wavelength distribution; Wavelength (nm); frequency", 200, 0.06*pow(10.,3), 0.53*pow(10., 3));

TH2F* h2xy = new TH2F("h2xy", "X-Y hits",  100, -31, 31,1000, 44, 105);
TH2F* h2xz = new TH2F("h2xz", "X-Z hits",  100, -30, 30,1000, 3354.993, 3355.010);
h2xy->SetStats(0);
h2xz->SetStats(0);

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
  h2xy->Fill(x,y);
  h2xz->Fill(x,z);  
}  

h2xy->SetContour(1000);
h2xz->SetContour(1000);
h2xy->SetOption("colz");
h2xz->SetOption("colz");

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

h2xy->Write();
h2xz->Write();

f->Close();
f1->Close();

cout << "Dataprocessing completed successfully. \nOutput file " << outp << " created." << endl;
return 1;

}



void rootMacro() {

int number;
char p;
cout << "Type root file number to be read. \nIf file name is FT0ntupleX.root, type X.\nEnter number:" << endl;
cin >> number;
cout << "Do you want to include preceding ntuple files starting from FT0ntuple0.root?\n [y] for yes, [n] for no." << endl;
cin >> p;

if (p == 'y') {

for (int i = 0; i <= number; i++) {
   help(i);
 }

} else {
   help(number);
}

return;
}






