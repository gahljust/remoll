// Usage:
// ./build/reroot -l -b -q 'processRate.C("Moller.list",1)' --> for beam generator
// ./build/reroot -l -b -q 'processRate.C("Moller.list",0)' --> for physics generators
//
// Arguments:
// finName:  Path to a single remoll ROOT file OR a text file listing ROOT files (one path per line).
// beamGenerator: 1 --> BEAM generator mode (weight per event, rate is 1).
// other physics generators: 0 --> use event weight from TTree branch 'rate'
//
// Detectors: histograms are defined for a small set of main-detector tiles (see detH/dtM).
// Adjust detector lists and histogram binning if you use a different detector plane.
// Particle Species (nSp=5):
//     0: e, pi (all energies)
//     1: e, pi with E > 1 MeV
//     2: gamma
//     3: neutron
//     4: primary e with E > 1 MeV  (mtrid==0 and trid==1 or 2 for Moller and beam; for ep-elastic/inelastic use mtrid==0 and trid==1)
//
// Output:
// A ROOT file named "<finNameStem>_det_proc.root" with directories per detector:
//     hRateLog_* : log10(rate) distributions (x-axis in log10)
//     hRate_*    : rate distributions (log-binned edges via niceLogBins; change if you prefer linear)
//     hrRate_*   : rate-weighted radial distributions, rate vs r [mm]
//     *XY_*      : rate-weighted 2D x–y hit distributions [mm]
//
// Normalization:
// - BEAM generator (rate=1): histograms scaled by 1/nTotEv.
// - Physics mode (tree weights, 'rate'): histograms scaled by 1/nFiles.


#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

TFile *fout = nullptr;

const int nSp = 5;// [e/pi,e/pi E>1,gamma,neutron]
const std::string spTit[nSp] = {"e/#pi","e/#pi E>1","#gamma","neutron","primary e E>1"};
const std::string spH[nSp]   = {"e","e1","g","n","eP1"};
std::map<int,int> spM {{11,1},{211,1},{22,3},{2112,4}};
//std::map<int,int> spM {{11,1},{-211,1},{22,3},{2112,4},{-11,6},{211,6}};

const int nDet = 9;
const std::string detH[nDet] = {"det150851", "det150852", "det150853","det150451", "det150452", "det150453","det110451", "det110452", "det110453"};
std::map<int,int> dtM = {{150851,1}, {150852,2}, {150853,3},{150451,4}, {150452,5}, {150453,6},{110451,7}, {110452,8}, {110453,9}};

TH1D *hRateLog[nSp][nDet], *hRate[nSp][nDet], *hrRate[nSp][nDet];
TH2D *hXYrate[nSp][nDet];

std::string fin;
int   nFiles = 0;
int beamGen = 1;
long nTotEv = 0;
long  currentEvNr = 0;

void niceLogBins(TH1*);

void initHisto();
long processOne(const std::string &fnm);
void process();
void writeOutput();

void processRate(const std::string& finName = "./remollout.root", int beamGenerator = 1){
  fin = finName;
  beamGen = beamGenerator;
  initHisto();
  process();
  writeOutput();
}

void process(){


  if(fin==""){
    std::cout << "\t did not find input file. Quitting!" << std::endl;
    return;
  }

  if( fin.find(".root") < fin.size() ){
    std::cout << "Processing single file:\n\t" << fin << std::endl;
    nTotEv += processOne(fin);
    nFiles = 1;
  }else{
    std::cout << "Attempting to process list of output from\n\t" << fin << std::endl;
    std::ifstream ifile(fin.c_str());
    std::string data;
    while(ifile >> data){
      std::cout << "processing: " << data << std::endl;
      nTotEv += processOne(data);
      nFiles++;
    }
  }


  std::cout << "\nFinished processing a total of " << nTotEv
            << " events across " << nFiles << " file(s).\n";
}

void initHisto(){


  std::string foutNm = Form("%s_det_proc.root", fin.substr(0, fin.find(".")).c_str());
  fout = new TFile(foutNm.c_str(), "RECREATE");

  for(int j=0;j<nDet;j++){
    fout->mkdir(detH[j].c_str(),Form("%s plane",detH[j].c_str()));
    fout->cd(detH[j].c_str());
    for(int i=0;i<nSp;i++){

      hRateLog[i][j] = new TH1D(Form("hRateLog_%s_%s",detH[j].c_str(),spH[i].c_str()),Form("Rate distribution (log10) %s;log10(rate);counts", spTit[i].c_str()),100,0,20);
      hRateLog[i][j]->Sumw2();

      hRate[i][j] = new TH1D(Form("hRate_%s_%s",detH[j].c_str(),spH[i].c_str()),Form("Rate distribution %s;rate;counts", spTit[i].c_str()),120,-1,20);
      hRate[i][j]->Sumw2();
      niceLogBins(hRate[i][j]);

      hrRate[i][j] = new TH1D(Form("hrRate_%s_%s",detH[j].c_str(),spH[i].c_str()),Form("Rate-weighted radial distribution %s;r [mm];rate", spTit[i].c_str()),175,600,1300);
      hrRate[i][j]->Sumw2();
        
      hXYrate[i][j] = new TH2D(Form("%sXY_%s",detH[j].c_str(),spH[i].c_str()), Form("Rate-weighted hits for %s;x [mm];y [mm]",spTit[i].c_str()), 1000,-2000,2000, 1000,-2000,2000);
      hXYrate[i][j]->Sumw2();


    }
  }
}

long processOne(const std::string &fnm){
  TFile *fin = TFile::Open(fnm.c_str(),"READ");
  if(!fin || !fin->IsOpen() || fin->IsZombie()){
    std::cout << "Problem: can't find file: " << fnm << std::endl;
    if (fin){ fin->Close(); delete fin; }
    return 0;
  }
  TTree *t = (TTree*)fin->Get("T");
  if (!t) { fin->Close(); delete fin; return 0; }
  Double_t rate = 0;
  std::vector<remollGenericDetectorHit_t> *hit = 0;

  t->SetBranchAddress("rate", &rate);
  t->SetBranchAddress("hit", &hit);

  Long64_t nEntries = t->GetEntries();
  std::cout << "\tTotal events: " << nEntries << std::endl;
  float currentProc = 1,procStep = 60;
  std::vector<int> procID;

  for (Long64_t event = 0; event < nEntries; t->GetEntry(event++)) {
    currentEvNr++;
    if( float(event+1)/nEntries*100 > currentProc){
      std::cout << "at tree entry\t" << event << "\t" << float(event+1)/nEntries*100 << std::endl;
      currentProc+=procStep;
    }

    if (beamGen) rate = 1; //this is for beam simulations only (right now)
    if (std::isnan(rate) || std::isinf(rate) || rate <= 0.0) continue;
    double lgRate = log10(rate);
      
    for(int j=0;j<hit->size();j++){

      //Particle species map
      int sp = spM[int(abs(hit->at(j).pid))] - 1;
      if(sp==-1) continue;

      //Detector map
      int dt = dtM[int(hit->at(j).det)]-1;
      if(dt==-1) continue;

      // e or pi with E>1 and without E>1 MeV
      if(sp == 0 ){
	    
	// all e or pi
	hRateLog[0][dt]->Fill(lgRate);
	hRate[0][dt]->Fill(rate);
	hrRate[0][dt]->Fill(hit->at(j).r,rate);
	hXYrate[0][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
          
	// all e or pi with E>1 MeV
	if((hit->at(j).k>1)){
	  hRateLog[1][dt]->Fill(lgRate);
	  hRate[1][dt]->Fill(rate);
	  hrRate[1][dt]->Fill(hit->at(j).r,rate);
	  hXYrate[1][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);

	  // primary electrons only (trid==1 or 2) with E>1 MeV
	  // for epelastic,epinelastic use trid==1 only
	  if(hit->at(j).mtrid == 0 && (hit->at(j).trid == 1 || hit->at(j).trid == 2 )){
	    hRateLog[4][dt]->Fill(lgRate);
	    hRate[4][dt]->Fill(rate);
	    hrRate[4][dt]->Fill(hit->at(j).r,rate);
	    hXYrate[4][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);
	  }
	}
      } else{
	//gamma and neutrons
	hRateLog[sp][dt]->Fill(lgRate);
	hRate[sp][dt]->Fill(rate);
	hrRate[sp][dt]->Fill(hit->at(j).r,rate);
	hXYrate[sp][dt]->Fill(hit->at(j).x,hit->at(j).y,rate);

      }
    }
  }
  fin->Close();
  delete fin;
  return nEntries;
};

void writeOutput(){
  for(int j=0;j<nDet;j++){
    fout->cd(detH[j].c_str());
    for(int i=0;i<nSp;i++){

      double scaleFactor = 1./nFiles;
      if(beamGen) scaleFactor = 1./nTotEv;
    
      hRateLog[i][j]->Scale(scaleFactor);
      hRateLog[i][j]->Write();
     
      hRate[i][j]->Scale(scaleFactor);
      hRate[i][j]->Write();
        
      hrRate[i][j]->Scale(scaleFactor);
      hrRate[i][j]->Write();
        
      hXYrate[i][j]->Scale(scaleFactor);
      hXYrate[i][j]->Write();
    }
  }
  fout->Close();
}

void niceLogBins(TH1*h)
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  double from = axis->GetXmin();
  double to = axis->GetXmax();
  double width = (to - from) / bins;
  double *new_bins = new double[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = pow(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}
