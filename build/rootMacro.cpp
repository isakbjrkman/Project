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
double lambda;
double qe;
double temp;


double waveE[354] = {2.38882,2.39344,2.39807,2.40271,2.40738,2.41206,2.41676,2.42148,2.42622,2.43098,2.43576,2.44055,2.44536,2.4502,2.45505,2.45992,2.46481,2.46972,2.47465,2.4796,2.48457,2.48956,
2.49457,2.4996,2.50465,2.50972,2.51481,2.51992,2.52505,2.5302,2.53538,2.54057,2.54579,2.55103,2.55629,2.56157,2.56687,2.5722,2.57755,2.58292,2.58831,2.59372,2.59916,2.60462,2.61011,2.61561,
2.62114,2.62669,2.63227,2.63787,2.6435,2.64915,2.65482,2.66052,2.66624,2.67198,2.67775,2.68355,2.68937,2.69522,2.70109,2.70699,2.71291,2.71886,2.72484,2.73084,2.73687,2.74292,2.749,2.75511,
2.76125,2.76741,2.7736,2.77982,2.78607,2.79234,2.79865,2.80498,2.81134,2.81773,2.82415,2.83059,2.83707,2.84358,2.85011,2.85668,2.86328,2.86991,2.87657,2.88326,2.88998,2.89673,2.90351,2.91033,
2.91718,2.92406,2.93097,2.93791,2.94489,2.9519,2.95895,2.96603,2.97314,2.98029,2.98747,2.99469,3.00194,3.00922,3.01655,3.0239,3.0313,3.03873,3.04619,3.05369,3.06123,3.06881,3.07643,3.08408,
3.09177,3.0995,3.10727,3.11508,3.12292,3.13081,3.13873,3.1467,3.15471,3.16276,3.17084,3.17897,3.18715,3.19536,3.20362,3.21192,3.22026,3.22865,3.23708,3.24555,3.25407,3.26263,3.27124,3.27989,
3.28859,3.29734,3.30613,3.31497,3.32386,3.3328,3.34178,3.35081,3.35989,3.36902,3.3782,3.38743,3.39671,3.40604,3.41543,3.42486,3.43435,3.44389,3.45348,3.46313,3.47283,3.48258,3.49239,3.50226,
3.51218,3.52216,3.53219,3.54229,3.55244,3.56264,3.57291,3.58324,3.59362,3.60407,3.61458,3.62515,3.63578,3.64647,3.65723,3.66805,3.67893,3.68988,3.7009,3.71198,3.72312,3.73434,3.74562,3.75697,
3.76839,3.77988,3.79144,3.80307,3.81477,3.82654,3.83839,3.85031,3.86231,3.87437,3.88652,3.89874,3.91104,3.92342,3.93587,3.94841,3.96102,3.97372,3.9865,3.99935,4.0123,4.02532,4.03844,4.05163,
4.06492,4.07829,4.09175,4.1053,4.11894,4.13267,4.14649,4.1604,4.17441,4.18851,4.20271,4.21701,4.2314,4.24589,4.26048,4.27517,4.28997,4.30486,4.31986,4.33497,4.35018,4.36549,4.38092,4.39645,
4.4121,4.42786,4.44373,4.45971,4.47581,4.49203,4.50836,4.52482,4.54139,4.55809,4.57491,4.59185,4.60892,4.62612,4.64345,4.6609,4.67849,4.69621,4.71407,4.73206,4.75019,4.76846,4.78687,4.80543,
4.82412,4.84297,4.86196,4.8811,4.9004,4.91984,4.93944,4.9592,4.97912,4.99919,5.01943,5.03984,5.06041,5.08115,5.10206,5.12314,5.1444,5.16583,5.18745,5.20924,5.23122,5.25339,5.27574,5.29829,
5.32103,5.34397,5.3671,5.39043,5.41397,5.43772,5.46167,5.48584,5.51022,5.53482,5.55964,5.58468,5.60995,5.63545,5.66119,5.68716,5.71336,5.73981,5.76651,5.79346,5.82066,5.84811,5.87583,
5.90381,5.93206,5.96058,5.98937,6.01845,6.0478,6.07745,6.10739,6.13762,6.16816,6.199,6.23015,6.26162,6.2934,6.32551,6.35795,6.39072,6.42383,6.45729,6.4911,6.52526,6.55979,6.59468,6.62995,
6.66559,6.70162,6.73804,6.77486,6.81209,6.84972,6.88778,6.92626,6.96517,7.00452,7.04432,7.08457,7.12529,7.16647,7.20814,7.25029,7.29294,7.33609,7.37976,7.42395,7.46867};

double conv2 = 1.60217662*pow(10.,-19);   //eV to Joules

double quantumEff[354] = {0.0966295,0.0998724,0.103128,0.106396,0.109677,0.11218,0.114692,0.117215,0.119747,0.12229,0.123772,0.12526,0.126754,0.128254,0.12976,0.130429,0.1311,0.131774,0.132451,0.13313,0.133825,
0.134521,0.135221,0.135924,0.136629,0.13674,0.136851,0.136963,0.137075,0.137188,0.137946,0.138708,0.139472,0.14024,0.141011,0.143063,0.145124,0.147193,0.149271,0.151357,0.152411,0.153469,
0.154532,0.155599,0.15667,0.158349,0.160034,0.161727,0.163427,0.165135,0.165215,0.165296,0.165377,0.165458,0.16554,0.168134,0.170739,0.173355,0.175983,0.178622,0.180839,0.183065,0.1853,
0.187546,0.189801,0.190908,0.19202,0.193137,0.194259,0.195385,0.197107,0.198836,0.200572,0.202317,0.204069,0.205326,0.206588,0.207856,0.20913,0.21041,0.210022,0.209633,0.209241,0.208848,
0.208453,0.209131,0.209811,0.210495,0.211182,0.211873,0.213663,0.215462,0.217269,0.219085,0.220909,0.221313,0.221719,0.222127,0.222537,0.222949,0.224175,0.225406,0.226644,0.227887,0.229137,
0.229239,0.229342,0.229446,0.22955,0.229654,0.230149,0.230646,0.231146,0.231649,0.232153,0.233597,0.235049,0.236507,0.237973,0.239446,0.239794,0.240144,0.240496,0.240849,0.241205,0.241242,
0.241279,0.241317,0.241355,0.241393,0.24078,0.240164,0.239545,0.238923,0.238297,0.238017,0.237736,0.237453,0.237168,0.236882,0.236787,0.236691,0.236595,0.236498,0.236401,0.236576,0.236752,
0.236929,0.237106,0.237285,0.237076,0.236867,0.236656,0.236444,0.236231,0.236046,0.23586,0.235673,0.235485,0.235296,0.234754,0.234209,0.233661,0.23311,0.232556,0.232441,0.232326,0.23221,
0.232093,0.231976,0.232028,0.232081,0.232134,0.232188,0.232241,0.232159,0.232076,0.231993,0.231909,0.231825,0.231847,0.23187,0.231893,0.231916,0.231939,0.232056,0.232174,0.232292,0.232412,
0.232532,0.232449,0.232366,0.232282,0.232198,0.232114,0.231788,0.23146,0.23113,0.230798,0.230463,0.229834,0.229201,0.228564,0.227922,0.227277,0.226647,0.226014,0.225376,0.224734,0.224088,
0.223465,0.222838,0.222206,0.221571,0.220931,0.219714,0.218488,0.217255,0.216013,0.214764,0.213563,0.212355,0.211139,0.209914,0.208681,0.207758,0.206829,0.205894,0.204952,0.204003,0.202881,
0.201751,0.200613,0.199468,0.198314,0.197222,0.196122,0.195014,0.193898,0.192774,0.19185,0.19092,0.189982,0.189038,0.188087,0.186929,0.185763,0.184588,0.183405,0.182212,0.181477,0.180735,
0.179989,0.179236,0.178478,0.177884,0.177285,0.176682,0.176074,0.175461,0.174629,0.173789,0.172944,0.172091,0.171232,0.170733,0.17023,0.169723,0.169212,0.168697,0.168375,0.168051,0.167725,
0.167395,0.167063,0.166787,0.166509,0.166228,0.165945,0.16566,0.165735,0.16581,0.165886,0.165962,0.166039,0.16609,0.166142,0.166194,0.166247,0.1663,0.166164,0.166027,0.165889,0.16575,0.165609,
0.165417,0.165223,0.165027,0.16483,0.164631,0.164954,0.165281,0.16561,0.165943,0.166278,0.167092,0.167914,0.168743,0.169581,0.170426,0.171972,0.173534,0.17511,0.176702,0.178309,0.180543,
0.182799,0.185077,0.187379,0.189703,0.192684,0.195696,0.198738,0.201811,0.204916,0.208053,0.211703,0.215597,0.219533,0.22351,0.227529,0.231591,0.233567,0.235325,0.237103,0.2389,0.240717,
0.241684,0.240623,0.239549,0.238464,0.237367,0.236257,0.235135,0.234,0.225402,0.216705,0.207907,0.199006,0.19,0.165996,0.141707,0.117127,0.0922502};


TH1D* h1ev = new TH1D("h1ev", "EventID; Event ID; Count", nentries/320, 0, nentries/320);
TH1D* h1de = new TH1D("h1de", "DetectorID; Detector ID; Count", 50, 0.5, 4.5);
TH1D* h1pdg = new TH1D("h1pdg", "PDG; pdg; Count", 100, -50, 50);
TH1D* h1px = new TH1D("h1px", "PX; px (MeV); Count", 200, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1py = new TH1D("h1py", "PY; py (MeV); Count", 200, 8*pow(10., -6), 8*pow(10., -6));
TH1D* h1pz = new TH1D("h1pz", "PZ; pz (MeV); Count", 200, 0, 8*pow(10., -6));
TH1D* h1x = new TH1D("h1x", "X; x (mm); Count", 200, -35, 35);
TH1D* h1y = new TH1D("h1y", "Y; y (mm); Count", 200, 42, 107);
TH1D* h1z = new TH1D("h1z", "Z; z (mm); Count", 500, 3354.993, 3355.008);
TH1D* h1w = new TH1D("h1w", "Wavelength distribution; Wavelength (nm); Count", 500, 0.06*pow(10.,3), 0.53*pow(10., 3));

TH2D* h2qe = new TH2D("h2qe", "Quantum efficiency as a function of wavelength; Wavelength (nm); Efficiency", 500, 0.06*pow(10.,3), 0.53*pow(10., 3), 500, 0, 0.5);

TH2D* h2xy = new TH2D("h2xy", "X-Y hits; x (mm); y (mm)",  200, -36, 36,200, 39, 110);
TH2D* h2xz = new TH2D("h2xz", "X-Z hits; x (mm); z (mm)",  200, -35, 35,200, 3354.993, 3355.008); 
TH2D* h2xpx = new TH2D("h2xpx", "X-PX hits; x (mm); px (MeV)",  200, -35, 35,200, 8*pow(10., -6), 8*pow(10., -6)); 
 
h2qe->SetStats(0);
h2xy->SetStats(0);
h2xz->SetStats(0);
h2xpx->SetStats(0);

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
  h2xpx->Fill(x,px);  
}  

for (int j = 0; j<354; j++) {
  temp = waveE[j];
  lambda = pow(10,9)*h*c/(temp*conv2);
  qe = quantumEff[j];
  h2qe->Fill(lambda, qe);
 } 


h1ev->SetFillColor(kBlue+2);

h2xy->SetContour(10000);
h2xz->SetContour(10000);
h2xpx->SetContour(10000);
h2xy->SetOption("colz");
h2xz->SetOption("colz");
h2xpx->SetOption("colz");

TCanvas canvas("Canvas","",950,600);

//gStyle->SetStatX(1.0);

h1w->Draw();
canvas.SaveAs("h1w.pdf","pdf");
h1w->Write();
canvas.Clear();

h1ev->Draw();
canvas.SaveAs("h1ev.pdf","pdf");
h1ev->Write();
canvas.Clear();

h1de->Draw();
gPad->Update();
auto statDE = dynamic_cast<TPaveStats*>(h1de->FindObject("stats"));
if (statDE) {
statDE->SetX1NDC(0.96);
statDE->SetX2NDC(0.81);
statDE->SetY1NDC(1.00);
statDE->SetY2NDC(0.85);
statDE->Draw();
} else {
cerr << "Error, no stats box found. \n";
}
canvas.SaveAs("h1de.pdf","pdf");
h1de->Write();
canvas.Clear();

h1pdg->Draw();
canvas.SaveAs("h1pdg.pdf","pdf");
h1pdg->Write();
canvas.Clear();

h1px->Draw();
canvas.SaveAs("h1px.pdf","pdf");
h1px->Write();
canvas.Clear();

h1py->Draw();
canvas.SaveAs("h1py.pdf","pdf");
h1py->Write();
canvas.Clear();

h1pz->Draw();
canvas.SaveAs("h1pz.pdf","pdf");
h1pz->Write();
canvas.Clear();

h1x->Draw();
gPad->Update();
auto statX = dynamic_cast<TPaveStats*>(h1x->FindObject("stats"));
if (statX) {
statX->SetX1NDC(0.76);
statX->SetX2NDC(0.61);
statX->SetY1NDC(1.00);
statX->SetY2NDC(0.86);
statX->Draw();
} else {
cerr << "Error, no stats box found. \n";
}
canvas.SaveAs("h1x.pdf","pdf");
h1x->Write();
canvas.Clear();

h1y->Draw();
gPad->Update();
auto statY = dynamic_cast<TPaveStats*>(h1y->FindObject("stats"));
if (statY) {
statY->SetX1NDC(0.76);
statY->SetX2NDC(0.61);
statY->SetY1NDC(1.00);
statY->SetY2NDC(0.86);
statY->Draw();
} else {
cerr << "Error, no stats box found. \n";
}
canvas.SaveAs("h1y.pdf","pdf");
h1y->Write();
canvas.Clear();

h1z->Draw();
gPad->Update();
gPad->SetLogy();   //test
canvas.SaveAs("h1z.pdf","pdf");
h1z->Write();
canvas.Clear();

h2xy->Draw();
gPad->Update();
gPad->SetLogy(0);
canvas.SaveAs("h2xy.pdf","pdf");
h2xy->Write();
canvas.Clear();

h2qe->Draw();
canvas.SaveAs("h2qe.pdf","pdf");
h2qe->Write();
canvas.Clear();

h2xz->Draw();
canvas.SaveAs("h2xz.pdf","pdf");
h2xz->Write();
canvas.Clear();

h2xpx->Draw();
canvas.SaveAs("h2xpx.pdf","pdf");
h2xpx->Write();
canvas.Clear();


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






