#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TH1D.h>                     // 1D histogram class
#include <TLorentzVector.h>           // 4-vector class
#include <TVector2.h>                 // 2-vector class
#include <vector>                     // STL vector class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <cmath>                      // C++ math library
#include <cassert>

#include "CSample.hh"                 // helper class to manage samples
#endif

using namespace std;


void DefineListHistos(vector< pair <string, vector<float> > > *list, vector <bool> *settingslogscale, vector < pair< float, float> > *settingsrange, vector<string> *settingslabel, bool SaveGenLevel, bool phosel, bool MVARESTOPTAGGER){
  //list of histo names and their binning (nbins,xmin,xmax) 
  string hname; vector<float> hbinning;
  hname = "hNPU"; hbinning.push_back(100); hbinning.push_back(-0.5); hbinning.push_back(95.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of PU");
  hname = "hNPUlinear"; hbinning.push_back(100); hbinning.push_back(-0.5); hbinning.push_back(95.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of PU");
  hname = "hNPV"; hbinning.push_back(100); hbinning.push_back(-0.5); hbinning.push_back(95.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of PV");
  hname = "hNPVlinear"; hbinning.push_back(100); hbinning.push_back(-0.5); hbinning.push_back(95.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of PV");
  hname = "hMET"; hbinning.push_back(112); hbinning.push_back(0); hbinning.push_back(560); //name, nbins, xmin, xmax
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} [GeV]");
  hname = "hMETlinear"; hbinning.push_back(112); hbinning.push_back(0); hbinning.push_back(560);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} [GeV]");
  hname = "hMETphi"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("#MET #Phi");
  hname = "hMETphilinear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("#MET #Phi");
  hname = "hNJets"; hbinning.push_back(16); hbinning.push_back(-0.5); hbinning.push_back(15.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);// ,make_pair(2e-4,2)));
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of jets");
  hname = "hNBJets"; hbinning.push_back(10); hbinning.push_back(-0.5); hbinning.push_back(9.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of b-tagged jets");
  hname = "hNBLJets"; hbinning.push_back(10); hbinning.push_back(-0.5); hbinning.push_back(9.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of loose b-tagged jets");
  hname = "res_hJet1csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 1 CSV");
  hname = "res_hJet2csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 2 CSV");
  hname = "res_hJet3csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 3 CSV");
  hname = "res_hJet4csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 4 CSV");
  hname = "res_hJet5csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 5 CSV");
  hname = "res_hJet6csv"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 6 CSV");
  hname = "res_hJet1qgid"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 1 QGID");
  hname = "res_hJet2qgid"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 2 QGID");
  hname = "res_hJet3qgid"; hbinning.push_back(48); hbinning.push_back(-0.2); hbinning.push_back(1.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 3 QGID");
  hname = "res_hJet1Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 1 #eta");
  hname = "res_hJet2Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 2 #eta");
  hname = "res_hJet3Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 3 #eta");
  hname = "res_hJet4Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 4 #eta");
  hname = "res_hJet5Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 5 #eta");
  hname = "res_hJet6Eta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 6 #eta");
  hname = "res_hJet1Pt"; hbinning.push_back(100); hbinning.push_back(0.0); hbinning.push_back(600.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 1 P_{T} [GeV/c]");
  hname = "res_hJet2Pt"; hbinning.push_back(70); hbinning.push_back(0.0); hbinning.push_back(420.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 2 P_{T} [GeV/c]");
  hname = "res_hJet3Pt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 3 P_{T} [GeV/c]");
  hname = "res_hJet4Pt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 4 P_{T} [GeV/c]");
  hname = "res_hJet5Pt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 5 P_{T} [GeV/c]");
  hname = "res_hJet6Pt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Jet 6 P_{T} [GeV/c]");
  hname = "hDeltaPhiMETmultijets"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
settingsrange->push_back(make_pair(2e-4,2));
 settingslabel->push_back("#Delta #phi (#slash{E}_{T}, multijets)");
  hname = "hMinDeltaPhiMETjets"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet");
  hname = "hMinDeltaPhiMETjetslinear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet");
  hname = "hMinDeltaPhiMETjets12"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet12");
  hname = "hMinDeltaPhiMETjets12linear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet12");
  hname = "hMinDeltaPhiMETjets123"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet123");
  hname = "hMinDeltaPhiMETjets123linear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet123");
  hname = "hMinDeltaPhiMETjets1234"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet1234");
  hname = "hMinDeltaPhiMETjets1234linear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet1234");
  hname = "hMinDeltaPhiMETjets12345"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet12345");
  hname = "hMinDeltaPhiMETjets12345linear"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("mindphimetjet12345");
  //  hname = "hTopMass"; hbinning.push_back(500); hbinning.push_back(0.0); hbinning.push_back(2000);
  hname = "hTopMass"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(2000);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("TopMass [GeV/c]");
  // hname = "hTopMasslinear"; hbinning.push_back(500); hbinning.push_back(0.0); hbinning.push_back(2000);
  hname = "hTopMasslinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(2000);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("TopMass [GeV/c]");
  hname = "hWMass"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(400);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("WMass [GeV/c]");
  hname = "hWMasslinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(400);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("WMass [GeV/c]");

  if(MVARESTOPTAGGER){
    hname = "hResTopMVA"; hbinning.push_back(50); hbinning.push_back(-1.0); hbinning.push_back(1.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("Res Top TMVA");
    hname = "hResTopMVAlinear"; hbinning.push_back(50); hbinning.push_back(-1.0); hbinning.push_back(1.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("Res Top TMVA");
    hname = "hResProb"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(1.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("Res Prob");
    hname = "hResChisq"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(100.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("Res Chisq");
    hname = "hResCost"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(1000.);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("Res Cost");
    hname = "hResFit_Mass"; hbinning.push_back(70); hbinning.push_back(0.0); hbinning.push_back(350);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("TMVA Fit TopMass [GeV/c]");
    hname = "hResFit_MassW"; hbinning.push_back(70); hbinning.push_back(0.0); hbinning.push_back(350);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("TMVA Fit WMass [GeV/c]");
    hname = "hResBDTDPhij1b"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#delta #phi j1,b");
    hname = "hResBDTDPhij2b"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#delta #phi j2,b");
    hname = "hResBDTDRj1b"; hbinning.push_back(80); hbinning.push_back(0.0); hbinning.push_back(8.);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#delta R j1,b");
    hname = "hResBDTDRj2b"; hbinning.push_back(80); hbinning.push_back(0.0); hbinning.push_back(8.);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#delta R j2,b");
  }//if(MVARESTOPTAGGER){


  hname = "hNEle"; hbinning.push_back(5); hbinning.push_back(-0.5); hbinning.push_back(4.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);// ,make_pair(2e-4,2)));
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of electron");
  hname = "hNMuo"; hbinning.push_back(5); hbinning.push_back(-0.5); hbinning.push_back(4.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);// ,make_pair(2e-4,2)));
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of muon");
  hname = "hNTau"; hbinning.push_back(5); hbinning.push_back(-0.5); hbinning.push_back(4.5);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);// ,make_pair(2e-4,2)));
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Number of tau");
  hname = "LowestPtEle"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Ele P_{T} [GeV/c]");
  hname = "LowestPtElelinear"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Ele P_{T} [GeV/c]");
  hname = "HighestPtEle"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Ele P_{T} [GeV/c]");
  hname = "HighestPtElelinear"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Ele P_{T} [GeV/c]");
  hname = "LowestPtMuo"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Muo P_{T} [GeV/c]");
  hname = "LowestPtMuolinear"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Muo P_{T} [GeV/c]");
  hname = "HighestPtMuo"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Muo P_{T} [GeV/c]");
  hname = "HighestPtMuolinear"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Muo P_{T} [GeV/c]");
  hname = "LowestPtTau"; hbinning.push_back(20); hbinning.push_back(0.0); hbinning.push_back(200.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("Tau P_{T} [GeV/c]");
  hname = "hMT"; hbinning.push_back(100); hbinning.push_back(0); hbinning.push_back(500);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("M_{T} [GeV]");
  hname = "hMTlinear"; hbinning.push_back(100); hbinning.push_back(0); hbinning.push_back(500);
  //hbinning.push_back(40); hbinning.push_back(0); hbinning.push_back(400); //name, nbins, xmin, xmax
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("M_{T} [GeV]");
  hname = "hMT2W"; hbinning.push_back(100); hbinning.push_back(0); hbinning.push_back(700);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("MT2W [GeV]");
  hname = "hMT2Wlinear"; hbinning.push_back(100); hbinning.push_back(0); hbinning.push_back(700);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("MT2W [GeV]");
  hname = "hDeltaPhiLepMET"; hbinning.push_back(64); hbinning.push_back(0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("#Delta #Phi Lep-MET");
  hname = "hMTnolepinmet"; hbinning.push_back(40); hbinning.push_back(0); hbinning.push_back(400); //name, nbins, xmin, xmax
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("M_{T} - nolep in met [GeV]");
  hname = "hDeltaPhiLepMETnolepinmet"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("#Delta #Phi Lep-MET - no lep in met");
  hname = "hMETnoLep"; hbinning.push_back(50); hbinning.push_back(0); hbinning.push_back(100); //name, nbins, xmin, xmax
  //hbinning.push_back(51); hbinning.push_back(50); hbinning.push_back(560); //name, nbins, xmin, xmax
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} - no lep [GeV]");
  hname = "hMETnoLeplinear"; hbinning.push_back(50); hbinning.push_back(0); hbinning.push_back(100); 
  //hbinning.push_back(51); hbinning.push_back(50); hbinning.push_back(560); //name, nbins, xmin, xmax
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} - no lep [GeV]");
  hname = "hMETphinoLep"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("#MET #Phi - no lep");
  hname = "hMinDeltaRLepJet"; hbinning.push_back(64); hbinning.push_back(0.0); hbinning.push_back(3.2);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("MinDeltaREleorMuoJet");
  hname = "hDilepMass"; hbinning.push_back(400); hbinning.push_back(0.0); hbinning.push_back(200);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("DilepMass");
  hname = "hDilepMasslinear"; hbinning.push_back(400); hbinning.push_back(0.0); hbinning.push_back(200);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,3));
  settingslabel->push_back("DilepMass");



  if(SaveGenLevel && !phosel){
  //   hname = "hNLepG"; hbinning.push_back(5); hbinning.push_back(-0.5); hbinning.push_back(4.5);
  //   list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  //   settingslogscale->push_back(false);
  //   settingsrange->push_back(make_pair(2e-4,2));
  //   settingslabel->push_back("Number of generator lep (status 3)");
    hname = "hEleorMuGEta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("EleorMuG #eta");
    hname = "hEleorMuGEtalinear"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("EleorMuG #eta");
    hname = "hNuGEta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#{nu}G #eta");
    hname = "hNuGEtalinear"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#{nu}G #eta");
    hname = "hWGEta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("WG #eta");
    hname = "hWGEtalinear"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("WG #eta");
  //   hname = "hTauGEta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  //   list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  //   settingslogscale->push_back(true);
  //   settingsrange->push_back(make_pair(2e-4,2));
  //   settingslabel->push_back("TauG #eta");
    hname = "hEleorMuoPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("EleorMuoG P_{T}");
    hname = "hEleorMuoPtlinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("EleorMuoG P_{T}");
    hname = "hNuPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#{nu}G P_{T}");
    hname = "hNuPtlinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#{nu}G P_{T}");
    hname = "hWPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("WG P_{T}");
    hname = "hWPtlinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("WG P_{T}");
    hname = "hGLepPtmLepPt"; hbinning.push_back(100); hbinning.push_back(-10.); hbinning.push_back(10.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("GLepPtmLepPt");
    hname = "hGLepPtmLepPtlinear"; hbinning.push_back(100); hbinning.push_back(-10.); hbinning.push_back(10.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("GLepPtmLepPt");




    //   hname = "hTauPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(300.0);
    //   list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    //   settingslogscale->push_back(true);
    //   settingsrange->push_back(make_pair(2e-4,2));
    //   settingslabel->push_back("TauG P_{T}");
    //   hname = "hLepgPdgId"; hbinning.push_back(10); hbinning.push_back(9.5); hbinning.push_back(19.5);
    //   list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    //   settingslogscale->push_back(false);
    //   settingsrange->push_back(make_pair(2e-4,2));
    //   settingslabel->push_back("LepG PDGID");
  }//if(SaveGenLevel && !phosel)

  if(SaveGenLevel && phosel){
    hname = "GPhotonPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("G#gamma P_{T}");
    hname = "GPhotonPtlinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(false);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("G#gamma P_{T}");
    hname = "GPhotonPtmRecoPhotonPt"; hbinning.push_back(200); hbinning.push_back(-10.0); hbinning.push_back(10.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("G#gamma P_{T} - Reco#gamma P_{T}");
    hname = "GPhotonRecoPhotondR"; hbinning.push_back(200); hbinning.push_back(-5.0); hbinning.push_back(5.0);
    list->push_back(make_pair(hname, hbinning)); hbinning.clear();
    settingslogscale->push_back(true);
    settingsrange->push_back(make_pair(2e-4,2));
    settingslabel->push_back("#Delta R (Gen - Reco #gamma)");
  }

  hname = "hphotonPt"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("photon P_{T} [GeV/c]");
  hname = "hphotonPtlinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(500.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("photon P_{T} [GeV/c]");
  hname = "hphotoneta"; hbinning.push_back(60); hbinning.push_back(-6.0); hbinning.push_back(6.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("photon #eta");
  hname = "hPFMETnopho"; hbinning.push_back(51); hbinning.push_back(50.0); hbinning.push_back(560.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} - no photon [GeV]");
  hname = "hPFMETnopholinear"; hbinning.push_back(50); hbinning.push_back(0.0); hbinning.push_back(1000.0);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(false);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("#slash{E}_{T} - no photon [GeV]");

  hname = "hisHadrTopHighestMVA"; hbinning.push_back(5); hbinning.push_back(-1.25); hbinning.push_back(1.25);
  list->push_back(make_pair(hname, hbinning)); hbinning.clear();
  settingslogscale->push_back(true);
  settingsrange->push_back(make_pair(2e-4,2));
  settingslabel->push_back("isHadrTopHighestMVA");


   
}

void DefineListHistos2D(vector< pair < string, vector< pair <float, float> > > > *list, vector<string> *settingslabelX, vector<string> *settingslabelY){
  //list of histo names and their binning (nbins,xmin,xmax) 
  string hname; vector< pair<float, float> >hbinningXY;
  hname = "hMETVsWpt"; 
  hbinningXY.push_back(make_pair(51,51)); hbinningXY.push_back(make_pair(50,50)); hbinningXY.push_back(make_pair(560,500)); //name, nbinsX, xmin, xmax
  settingslabelX->push_back("W gen P_T [GeV]");
  settingslabelY->push_back("#slash{E}_{T} [GeV]");
  list->push_back(make_pair(hname, hbinningXY)); hbinningXY.clear();
  hname = "hMETnolepVsWpt"; 
  hbinningXY.push_back(make_pair(51,51)); hbinningXY.push_back(make_pair(50,50)); hbinningXY.push_back(make_pair(560,500)); //name, nbinsX, xmin, xmax
  settingslabelX->push_back("W gen P_T [GeV]");
  settingslabelY->push_back("#slash{E}_{T}-nolep [GeV]");
  list->push_back(make_pair(hname, hbinningXY)); hbinningXY.clear();
}

void DefineInputSamples(vector<CSample*> *samplev, string inputdir, string inputdirQCD, bool splitttbar, int singlelepdata, bool blind, int mask){
  //do not touch the default list of samples. If you wanna change something go in the main macro and change them from there

  cout<<"mask "<<mask<<"(mask & 0x1) "<<(mask & 0x1)<<endl;


  //mask 
  //0x0000001 : diboson
  //0x0000010 : Zjets
  //0x0000100 : single top
  //0x0001000 : Wjets
  //0x0010000 : ttbar
  //0x0100000 : QCD
  //0x1000000 : ttDM
  //

   // data 

  if(!blind){
    samplev->push_back(new CSample("data",0,0));

    if(singlelepdata==1)
      samplev->back()->fnamev.push_back("../"+inputdir+"/SingleMuon_Run2015D_PrReco_MINIAOD_allhadrbits.root");
    else if(singlelepdata==2)
      samplev->back()->fnamev.push_back("../"+inputdir+"/SingleElectron_2012-22Jan2013_allhadrbits.root");
    else
      samplev->back()->fnamev.push_back("../"+inputdir+"/MET_2012-22Jan2013_allhadrbits.root");
  }          

  // dibosons 
  if((mask & 0000001) == 0000001){
    samplev->push_back(new CSample("diboson",kRed-6,kRed+2));
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WW_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WZ_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_ZZ_MINIAOD_allhadrbits.root");
  }
  // Z+jets
  if((mask & 0x0000010) == 0x0000010){
    samplev->push_back(new CSample("Z+jets",kOrange-2,kOrange-3));
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_DYJetsToLL_M-50_MINIAOD_allhadrbits.root");
  }            
  //single top
  if((mask & 0x0000100) == 0x0000100){
    samplev->push_back(new CSample("single t",kAzure-8,kBlue-6));
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_ST_tW_antitop_5f_inclusiveDecays_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_ST_tW_top_5f_inclusiveDecays_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_ST_t-channel_antitop_4f_leptonDecays_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_ST_t-channel_top_4f_leptonDecays_MINIAOD_allhadrbits.root");

  }
  //W+jets 
  if((mask & 0x0001000) == 0x0001000){
    samplev->push_back(new CSample("W+jets",kOrange+7,kOrange+10));
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-100To200_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-200To400_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-400To600_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-600To800_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-800To1200_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-1200To2500_MINIAOD_allhadrbits.root");
    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_WJetsToLNu_HT-2500ToInf_MINIAOD_allhadrbits.root");
  }   
 
  // ttbar 
  if((mask & 0x0010000) == 0x0010000){
    if(splitttbar){
      samplev->push_back(new CSample("t#bar{t} Dil",kGray+2,kGray+3)); 
      samplev->back()->fnamev.push_back("../"+inputdir+"/baconbits_1tightmuon_3ormorejets_ttbarsplit_27Oct2015_2/2L/Spring15_a25ns_TTJets_amcatnlo_MINIAOD_allhadrbits.root");
      samplev->push_back(new CSample("t#bar{t} Semi",kGreen+2,kGreen+3));
      samplev->back()->fnamev.push_back("../"+inputdir+"/baconbits_1tightmuon_3ormorejets_ttbarsplit_27Oct2015_2/1L/Spring15_a25ns_TTJets_amcatnlo_MINIAOD_allhadrbits.root");
      samplev->push_back(new CSample("t#bar{t} AllHadr,t#bar{t}XJets",kCyan+2,kCyan+3));
      samplev->back()->fnamev.push_back("../"+inputdir+"/baconbits_1tightmuon_3ormorejets_ttbarsplit_27Oct2015_2/HAD/Spring15_a25ns_TTJets_amcatnlo_MINIAOD_allhadrbits.root");    
    }
    else{
      samplev->push_back(new CSample("t#bar{t}",kGreen+2,kGreen+3));

    //    samplev->back()->fnamev.push_back("../"+inputdir+"/Spring15_a25ns_TTJets_amcatnlo_MINIAOD_allhadrbits.root");
    }
    

    cout<<"splitttbar "<<splitttbar<<endl;
  }
  
 // //QCD 
 if((mask & 0x0100000) == 0x0100000){
   samplev->push_back(new CSample("QCD ",kMagenta+2,kMagenta+3));
   samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Spring15_a25ns_QCD_HT700to1000_MINIAOD_allhadrbits.root");
   samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Spring15_a25ns_QCD_HT1000to1500_MINIAOD_allhadrbits.root");
   samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Spring15_a25ns_QCD_HT1500to2000_MINIAOD_allhadrbits.root");
   samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Spring15_a25ns_QCD_HT2000toInf_MINIAOD_allhadrbits.root");

 }
 
  
  
  
 
 //old
 // // // DM 
 // if((mask & 0x1000000) == 0x1000000){
 //   samplev->push_back(new CSample("DM 1 GeV",kBlue+2,kBlue+3));
 //   samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_TTChiChiJets_M-1_TuneZ2star_allhadrbits.root");
 // }


}


void DefineInputSamplesGammaJets(vector<CSample*> *samplev, string inputdir, string inputdirQCD, bool blind){
  //do not touch the default list of samples. If you wanna change something go in the main macro and change them from there

   // data

  samplev->push_back(new CSample("data",0,0));
  /*  if(!blind){
    samplev->back()->fnamev.push_back("../"+inputdir+"/SinglePhoton_2012-22Jan2013_bacon_allhadrbits.root");
    } */         

  
 // //QCD
 samplev->push_back(new CSample("QCD",kMagenta+2,kMagenta+3));
 /* samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Summer12_QCD_HT-100To250_TuneZ2star_allhadrbits.root"); //only 14 entries
 samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Summer12_QCD_HT-250To500_TuneZ2star_allhadrbits.root");
 samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Summer12_QCD_HT-500To1000_TuneZ2star_allhadrbits.root");
 samplev->back()->fnamev.push_back("../"+inputdirQCD+"/Summer12_QCD_HT-1000ToInf_TuneZ2star_allhadrbits.root");*/
 
 //Gamma+jets
 samplev->push_back(new CSample("#gamma+jet",kGreen+2,kGreen+3));
 samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_GJets_HT-40To100_allhadrbits.root");
 samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_GJets_HT-100To200_allhadrbits.root");
 samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_GJets_HT-200To400_allhadrbits.root");
 samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_GJets_HT-400ToInf_allhadrbits.root");

  
  
  
 

 // // DM
 samplev->push_back(new CSample("DM 1 GeV",kBlue+2,kBlue+3));
 //samplev->back()->fnamev.push_back("../"+inputdir+"/Summer12_TTChiChiJets_M-1_TuneZ2star_allhadrbits.root");
  


}
