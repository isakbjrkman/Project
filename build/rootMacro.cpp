#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <array>

void rootMacro() {
TFile* f = new TFile("data.root", "RECREATE");

TFile* f1 = new TFile("Results.root", "READ");
TTree* t = (TTree*)f1->Get("B5");


double px,py,pz;

t->SetBranchAddress("PX",&px);
t->SetBranchAddress("PY",&py);
t->SetBranchAddress("PZ",&pz);

int nentries = t->GetEntries(); 

double h = 6.62607015*pow(10.,-34);
double c = 299792458;
double E;

//double lambda[nentries];

TH1D* h1 = new TH1D("h1", "Wavelength distribution; wavelength (m); frequency", 100, pow(10., -9), 2*pow(10., -6));

double convert = 1.60217662*pow(10.,-13);   //MeV to Joules

for (int i=0; i<nentries; i++) {
  t->GetEntry(i);
  E = (px+py+pz)*convert; 
  h1->Fill(h*c/E);
}  


h1->Draw();
h1->Write();

}
