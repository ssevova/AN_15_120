//FOUND BUG WITH GetEntriesFast() for TClonesArray. Always count entries of the tclonearray by yourself until I find where the bug is

//================================================================================================
//
// Plot various distributions in all-hadronic ttbar events
//
//________________________________________________________________________________________________

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
#include <sstream>
#include <fstream>      // std::ofstream
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <cmath>                      // C++ math library
#include <cassert>

#include "CPlot.hh"                   // helper class for plots
#include "KStyle.hh"                  // style settings for drawing
#include "CSample.hh"                 // helper class to manage samples
#include "DefineInputSamples.C"
#include "Utils.C"
#include <TClonesArray.h>
#endif

using namespace std;



//
void plotAllHadrPresel()//bool blind=false, int masksamples = 0x0111111, float res_topmva_cut=-999., int isHadrTopHighestMVA_cut=-999, string outputDir = "tmpoutputdir")
{
  //res_topmva>res_topmva_cut; isHadrTopHighestMVA = isHadrTopHighestMVA_cut (if isHadrTopHighestMVA_cut used masksamples=0x0010000)

  //  cout<<"blind, masksamples , res_topmva_cut, isHadrTopHighestMVA_cut, outputDir "<<blind<<", "<< masksamples <<", "<< res_topmva_cut<<", "<< isHadrTopHighestMVA_cut<<", "<< outputDir<<endl;


  int useturnonpar=0;
  //forMC
  //-1 no turnon no emulation for MC
  //0 no turnon - relying on MC emulation of the trigger
  //1 multiply the weights for turnon pfmet150 for muo
  //2 multiply the weights for turnon pfmet150 for ele

  if(useturnonpar!=-1 && useturnonpar!=0 && useturnonpar!=1 && useturnonpar!=2){
    cout<<"useturnonpar "<<useturnonpar<<" not allowed"<<endl;
    exit(1);
  }

  unsigned int verbose = 0; //0,1,2 depending on the level of verbosity
  bool SaveGenLevel=false;
  bool SplitTtbar=true; //separate in the legend and yields different modes in the ttbar
  bool blind=false; int masksamples = 0x0111111; //0x0000010;
  //bool blind=true; int masksamples = 0x0010000; //Default 0x1111111
  string extvarweighting="";//NLOCorrPT_GtoZ";//"WGPT", "NUGPT_ELE", "NUGPT_MUO"; //if you wanna sculpt a given distribution and reweight events accordingly.. Suppoerted only wpt 
  bool phosel=false; //true photon+jets, false=SR, CRttbar,Zjets,Z, etc..
  int singlelepdata=1;  //1 for singlemu, 2 for singleele, 0 for  met. 1,2 for  trigger turnon studies or lepton energy studies (Z->ll)
  bool MVARESTOPTAGGER=true;
  bool SaveGenLevelTMVA=false; //if true will be used only for ttbar semileptonic

  //
  // Cuts
  //
  /*  const unsigned int NLEP_CUT = 0;
  const unsigned int NJETS_CUT  = 4;
  const unsigned int NBJETS_CUT = 2;
  const double       MET_CUT    = 320;
  const double       MET_CUT_PS    = 120;*/

  //CR cuts
  //  const unsigned int NBJETS_CUT_CR = 0;
  const double BTAG_WPL      = 0.244;  // CSV loose working point 

  //CR cuts (W+jets)
  //const unsigned int NLEP_CUT = 1;





  //charactherization of tmva tagger
  //string inputdir = "baconbits_tightmuo_3ormorejets_27Oct2015"; 
  //  string inputdir = "baconbits_2tightmuon_27Oct2015_2"; 
  //  string inputdir = "baconbits_2tightmuon_27Oct2015_wnegtweightdealt";


  //  string inputdir = "baconbits_1tightmuon_3ormorejets_MVARESTOPTAGGERfalse_27Oct2015";
  //  string inputdir = "baconbits_1tightmuon_3ormorejets_MVARESTOPTAGGERfalse_27Oct2015_2"; //with ggsmear
  
  string inputdir = "baconbits_1tightmuon_3ormorejets_27Oct2015"; //with ggsmear
  
  
  string inputdirQCD = inputdir;


  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_Oct282015");
  //string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_noPUreweigthing_Oct282015");
    //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_PUreweigthingcorrect_Oct282015");

  //  string outputDir("semilep_2tightmuo_trigrequestonMC_Oct282015");
  //  string outputDir("semilep_2tightmuo_trigrequestonMC_wnegtweightdealt_Oct282015");
  //string outputDir("semilep_2tightmuo_trigrequestonMC_wnegtweightdealt_first058fb_Oct282015");



  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_Oct292015");
  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_qgsmearedwith8TeVrecipe_Oct292015");
  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_qgsmearedwith8TeVrecipe_andtmvatoplot_Oct292015");
  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_Oct292015_muoptmore100");

  //  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_qgsmearedwith8TeVrecipe_nostop_Oct302015");
  string outputDir("semilep_1tightmuo_resolved_3ormorejets_2ormorejetWPm_pfmetmore100_pfmtmore40_trigrequestonMC_qgsmearedwith8TeVrecipe_Oct302015");




  cout<<"qui"<<endl;

  // Create output directory
 {
   string tmp = outputDir+"_tobedeleted";
   gSystem->Rename(outputDir.c_str(),tmp.c_str());
 }
 gSystem->mkdir(outputDir.c_str(), true);
 {
   string tmp = outputDir+"/plotAllHadrPresel.C";
   gSystem->CopyFile("./plotAllHadrPresel.C",tmp.c_str());
   tmp = outputDir+"/DefineInputSamples.C";
   gSystem->CopyFile("./DefineInputSamples.C",tmp.c_str());
   tmp = outputDir+"/Utils.C";
   gSystem->CopyFile("./Utils.C",tmp.c_str());


 }
 CPlot::sOutDir = outputDir;

 //
 // samples
 // Note: macro assumes samplev[0] is data
 //
 vector<CSample*> samplev;

 if(SaveGenLevel){//erase samples where no gen particles where saved
   blind = true; //trick not to plot the data without screweing up things
 }

 if(!phosel){
   cout<<"DefineInputSamples"<<endl;
   DefineInputSamples(&samplev, inputdir, inputdirQCD, SplitTtbar, singlelepdata, blind, masksamples); //default list of samples

 }
 else{
   DefineInputSamplesGammaJets(&samplev, inputdir, inputdirQCD, blind);
 }
 if(verbose)
   cout<<"after DefineInputSamplesGammaJets"<<endl;

 if(SaveGenLevel && !phosel){//erase samples where no gen particles where saved
     for(int iss=0; iss<samplev.size(); iss++){
     if((samplev[iss]->label).compare("diboson")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("Z+jets")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("single t")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("W+jets")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("t#bar{t} Dil")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("t#bar{t} AllHadr,t#bar{t}XJets")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("QCD")==0)
       samplev.erase(samplev.begin()+iss);
     if((samplev[iss]->label).compare("DM 1 GeV")==0)
       samplev.erase(samplev.begin()+iss);
   }
 }

   
 

 //   }
 // }
 for(int iss=0; iss<samplev.size(); iss++)
   cout<<"After process "<<samplev[iss]->label<<endl;
 

 // integrated luminosity to scale MC
 const double LUMI = 1.267;
  
 // histograms for various corrections
 const string puWeightFilename("/tthome/ksung/cms/Analysis/05/CMSSW_7_4_12/src/DMSAna/Utils/data/PVWeights.root");
 const string muonEffFilename("");
 const string eleEffFilename("");
  

 //--------------------------------------------------------------------------------------------------------------
 // Main analysis code
 //==============================================================================================================

 //
 // Declare histograms
 //  

 vector<double> neventsv, nerrv; 
 vector< vector<double> > neventssamplev, nentriessamplev;
 vector< vector<TH1D*> > hv;
 vector< vector<TH2D*> > hv2D;
 double neventsMC, nerrMC;
 vector<TH1D*> hMC;
 vector<TH2D*> hMC2D;
 vector<TH1D*> hRatio;


 //listhisto[ilist][isam] //for each var one entry per each sample
 //listvar float to fill, string for xaxis label
 vector< pair <string, vector<float> > > listhisto; 
 vector< pair < string, vector< pair <float, float> > > > listhisto2D; //For 2D var 
 vector <bool> settingslogscale; //(bool for log scale, 
 vector < pair< float, float > > settingsrange; //the two floats for ymin., ymax)
 vector<string> settingslabel, settingslabelX, settingslabelY; //xlabel
 vector<float> listvar; //those two lists must be kept in synch since the latter is the list of variable used the fill the former
 vector< pair<float, float> > listvar2D; //those two lists must be kept in synch since the latter is the list of variable used the fill the former
 //all var assumed to be float: does it matter?
if(verbose==2)
   cout<<"Before DefineListHistos"<<endl;
 DefineListHistos(&listhisto,&settingslogscale,&settingsrange,&settingslabel, SaveGenLevel, phosel, MVARESTOPTAGGER);  

 if(listhisto.size() != settingslogscale.size()){
   cerr<<"listhisto and settinglogscale should have the same size. Instead they have "<<listhisto.size()<<", "<<settingslogscale.size()<<endl;
   exit(1);
 }
 if(listhisto.size() != settingsrange.size()){
   cerr<<"listhisto and settingrange should have the same size. Instead they have "<<listhisto.size()<<", "<<settingsrange.size()<<endl;
   exit(1);
 }
 if(listhisto.size() != settingslabel.size()){
   cerr<<"listhisto and settinglabel should have the same size. Instead they have "<<listhisto.size()<<", "<<settingslabel.size()<<endl;
   exit(1);
 }
 if(verbose==2)
   cout<<"Before DefineListHistos"<<endl;
 string hname; vector<float> hbinning; //should be an array
 vector< pair<float, float> > hbinningXY;


 DefineListHistos2D(&listhisto2D, &settingslabelX, &settingslabelY);  
 if(listhisto2D.size() != settingslabelX.size()){
   cerr<<"listhisto2D and settinglabelX should have the same size. Instead they have "<<listhisto2D.size()<<", "<<settingslabelX.size()<<endl;
   exit(1);
 }
 if(listhisto2D.size() != settingslabelY.size()){
   cerr<<"listhisto2D and settinglabelY should have the same size. Instead they have "<<listhisto2D.size()<<", "<<settingslabelY.size()<<endl;
   exit(1);
 }


 for(unsigned int isam=0; isam<samplev.size(); isam++) {   
   neventsv.push_back(0);
 }
 
 for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {
   vector<TH1D*> htmp;
   for(unsigned int isam=0; isam<samplev.size(); isam++) {   
     ostringstream ostr_isam;
      hname = listhisto[ilist].first; hbinning = listhisto[ilist].second;
      if(verbose==2){
	cout<<"ilist "<<ilist<<" isam "<<isam<<endl;     
	cout<<"hname "<<hname<<endl;
	cout<<"nbins, xmin, xmax "<<hbinning[0]<<", "<<hbinning[1]<<", "<<hbinning[2]<<endl;
      }
      ostr_isam << isam; hname = hname + "v_"+ostr_isam.str();
      //     htmp.push_back(new TH1D(hname.c_str(),"",hbinning[0],hbinning[1],hbinning[2])); 
      htmp.push_back(new TH1D(hname.c_str(),"",hbinning[0],hbinning[1],hbinning[2])); 
      htmp[isam]->Sumw2();
      hbinning.clear();
   }//for(unsigned int isam=0; isam<samplev.size(); isam++)
   hv.push_back(htmp);
 }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++)

 for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++) {
   vector<TH2D*> htmp;
   for(unsigned int isam=0; isam<samplev.size(); isam++) {   
     ostringstream ostr_isam;
     hname = listhisto2D[ilist2D].first; hbinningXY = listhisto2D[ilist2D].second; 
      if(verbose==2){
	cout<<"ilist2D "<<ilist2D<<" isam "<<isam<<endl;     
	cout<<"hname "<<hname<<endl;
	cout<<"nbinsX, xmin, xmax "<<hbinningXY[0].first<<", "<<hbinningXY[1].first<<", "<<hbinningXY[2].first<<endl;
	cout<<"nbinsY, ymin, ymax "<<hbinningXY[0].second<<", "<<hbinningXY[1].second<<", "<<hbinningXY[2].second<<endl;
      }
      ostr_isam << isam; hname = hname + "v_"+ostr_isam.str();
      htmp.push_back(new TH2D(hname.c_str(),"",hbinningXY[0].first,hbinningXY[1].first,hbinningXY[2].first,hbinningXY[0].second,hbinningXY[1].second,hbinningXY[2].second));    
      hbinningXY.clear();
   }//for(unsigned int isam=0; isam<samplev.size(); isam++)
   hv2D.push_back(htmp);
 }//for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++)
   

 neventsMC=0;
 for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {
   hname = listhisto[ilist].first; //hbinning = listhisto[ilist].second;
   hname += "MC";
   TH1D *htmp = (TH1D*)hv[ilist][0]->Clone(hname.c_str());
   hMC.push_back(htmp);
 }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {

 for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++) {
   hname = listhisto2D[ilist2D].first; 
   hname += "MC";
   TH2D *htmp = (TH2D*)hv2D[ilist2D][0]->Clone(hname.c_str());
   hMC2D.push_back(htmp);
 }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {

 //
 // variables to read in bacon bits
 //
  //TO BE CHANGED ACCORDING TO runALLHadronicPreselection.cc
  unsigned int runNum, lumiSec, evtNum;          // event ID
  //  unsigned int selType;                          // selection category bits
  unsigned int metfilter;                        // MET filter bits
  unsigned int npv, npu;                         // number of PV / PU
  unsigned int njets, nfjets, /* njets08,*/ nbjets, nbljets;         // jet multiplicity
  float rho;                                     // event energy density
  float scale1fb;                                // cross section scale factor per 1/fb
  float evtWeight;                               // event weight (NOT from cross section normalization)
  float pfmet, pfmetphi, pfmt; 
  float pfmetraw, pfmetphiraw;                   // PF RAW MET
  float pfmetnolep, pfmetphinolep, pfmtnolepinmet;                   // PF MET removing lep (should be equal to pfmet when no lep)
  float pfmetnolep_in, pfmetphinolep_in; //computed at .cc must be equal to above
  float pfmetnopho, pfmetphinopho; //PFMET no pho, considering leading tight photon
  float mvamet, mvametphi;                // MVA MET


  float scEt, scEta, sieie;
  float chHadIso, neuHadIso, phoIso, sthovere;
  bool passElectronVeto;
  unsigned int photypeBits;
  unsigned int nloosepho, ntightpho;
  TLorentzVector *phovec=0;



  float deltaphilepmet;
  float deltaphilepmetnolepinmet;
  float mindeltarlepjet;
  bool passtrigger;
  float dilepmass; //nele==2 or nmuo==2 compite dilepmass
  //  int lepId;                                     // lepton PDG ID
  //  TLorentzVector *electron=0; TLorentzVector *muon=0; TLorentzVector *tau=0;                     // lepton 4-vector
  TClonesArray electron("TLorentzVector"); TClonesArray muon("TLorentzVector");  TClonesArray tau("TLorentzVector");  
  TClonesArray *electron_arr  = 0; TBranch *electron_br = 0; //don't remove the 0 or it wil crash
  TClonesArray *muon_arr = 0; TBranch *muon_br = 0;
  TClonesArray *tau_arr = 0; TBranch *tau_br = 0;
  unsigned int nele, nmuo, ntau, nlep;

  //gen level stuff for tmva
  int             genId1, genId2, genId3;
  TLorentzVector *genpar1=0, *genpar2=0, *genpar3=0;




  //info of the generated lep from W with status 3 
  unsigned int nlepg, nqg, ndmg, ntopg, nbg; 
  int lepgpdgid[4];
  int qgpdgid[4], dmgpdgid[2];
  int topgpdgid[2], bgpdgid[2];
  int lepgdaugtlep[4]; 
  TClonesArray lepg("TLorentzVector");
  TClonesArray qg("TLorentzVector");
  TClonesArray dmg("TLorentzVector");
  TClonesArray topg("TLorentzVector");
  TClonesArray bg("TLorentzVector");
  TClonesArray *lepg_arr = 0; TBranch *lepg_br = 0;
  TClonesArray *qg_arr = 0; TBranch *qg_br = 0;
  TClonesArray *dmg_arr = 0; TBranch *dmg_br = 0;
  TClonesArray *topg_arr = 0; TBranch *topg_br = 0;
  TClonesArray *bg_arr = 0; TBranch *bg_br = 0;
  TClonesArray *phog_arr = 0; TBranch *phog_br = 0;

  TClonesArray phog("TLorentzVector");

  unsigned int nphog;




  ///// boosted /////
  float bst_tau21, bst_tau32;
  // float  bst_prunedm, bst_subjetcsv;
  // TLorentzVector *bst_jet=0, *bst_topjet=0;
  // TLorentzVector *bst_subjet1=0, *bst_subjet2=0, *bst_subjet3=0;

  float bst_jetcsv;
  float bst_jetmprun, bst_jetmtrim, bst_jetmsd;
  TLorentzVector *bst_jet=0;

  float           bst_sub1csv,    bst_sub2csv,    bst_sub3csv,    bst_sub4csv;
  TLorentzVector *bst_subjet1=0, *bst_subjet2=0, *bst_subjet3=0, *bst_subjet4=0;


  //from the least to the most significant bit: MET, PV, Trigger, lepton veto; njet selection (either boosted or resolved)
  unsigned int cutflow; //keep track of which cut have been passed



  ///// resolved /////
  float           res_wmva, res_topmva;
  float           res_chisq, res_prob, res_cost, res_fitmass, res_fitmassW;
  float res_bdt_dphij1b_out, res_bdt_dphij2b_out, res_bdt_drj1b_out, res_bdt_drj2b_out;
  float           res_jet1csv,     res_jet2csv,     res_jet3csv;
  float           res_jet4csv,     res_jet5csv,     res_jet6csv;
  float           res_jet1qgid,    res_jet2qgid,    res_jet3qgid;
  float           res_jet4qgid,    res_jet5qgid,    res_jet6qgid;
  TVector2       *res_jet1pull=0, *res_jet2pull=0, *res_jet3pull=0;
  TLorentzVector *res_jet1=0,     *res_jet2=0,     *res_jet3=0;
  TLorentzVector *res_jet4=0,     *res_jet5=0,     *res_jet6=0;
  //

  float           inc_mT2W=-999.;



  //variable to be computed in this macro
  float isHadrTopHighestMVA; //0,1 is semileptonic ttbar event, -1 not semileptonic ttbar event; 0: highest MVA score is not top hadronic, 1 is top hadronic


 float dphimetmultijets; //DeltaPhi between met and the 4 jets system
 float mindphimetjet; //dPhi between met and the closest jet (out of max 6)
 float mindphimetjet12; //dPhi between met and the closest jet (out of the two leading jets)
 float mindphimetjet123; //dPhi between met and the closest jet (out of the three leading jets)
 float mindphimetjet1234; //dPhi between met and the closest jet (out of the four leading jets)
 float mindphimetjet12345; //dPhi between met and the closest jet (out of the five leading jets)
 unsigned int nbljets_mt;   
 vector<TLorentzVector> ele;
 vector<TLorentzVector> muo;
 vector<TLorentzVector> ta;
 vector<TLorentzVector> lepgg;

 vector<TLorentzVector> phogg;

 float topmass, wmass; //take all the jets in the event above threshold and build the mass from the combination of jets (3 and 2 respectively) which give the mass closer to the nominal top and W mass respectively


float eleormuogpt=-999.,  eleormuogeta=-999.; //filled only if SaveGenLevel = true
float nugpt,  nugeta; //filled only if SaveGenLevel = true
float wgpt,  wgeta; //filled only if SaveGenLevel = true

 float phogpt=-999.;  float phogdR=999.;//dR wrt to the matching reco pho

 TFile puWeightFile(puWeightFilename.c_str());
 TH1D *hPUWeights = (TH1D*)puWeightFile.Get("PVWeights");//pileup");
 hPUWeights->SetDirectory(0);
 puWeightFile.Close();

 TFile *infile=0;
 TTree *intree=0;
 vector<double> neventssamplev_tmp, nentriessamplev_tmp;

   for(int iss=0; iss<samplev.size(); iss++)
     cout<<"Before processing process "<<samplev[iss]->label<<endl;


 for(unsigned int isam=0; isam<samplev.size(); isam++) {
   cout<<"isam "<<isam<<endl;
   nentriessamplev_tmp.clear();
   neventssamplev_tmp.clear();
   if(verbose>=1){
     cout<<"nentriessamplev_tmp.size() "<<nentriessamplev_tmp.size()<<endl;
     cout<<"neventssamplev_tmp.size() "<<neventssamplev_tmp.size()<<endl;
   }
   CSample *sample = samplev[isam];
   if(verbose>=1)
     cout << "Sample: " << sample->label << endl;
   bool isData   = (isam==0);
   bool isSignal = (isam==samplev.size()-1);
    
   cout<<"qui "<<sample->fnamev.size()<<endl;


   for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++) {
     nentriessamplev_tmp.push_back(0);
     neventssamplev_tmp.push_back(0);      
     string infilename = sample->fnamev[ifile];
     cout << " ==> Processing " << infilename << "..." << endl;
     /*     if(infilename.find("Summer12_TTJets_FullLeptMGDecays_allhadrbits") != string::npos){
       cout<<"setting the verbose to 2"<<endl;
       verbose = 2;
       }*/
     infile = new TFile(infilename.c_str()); assert(infile);
     intree = (TTree*)infile->Get("Events"); assert(intree);


     intree->SetBranchAddress("runNum",    &runNum);
     intree->SetBranchAddress("lumiSec",   &lumiSec);
     intree->SetBranchAddress("evtNum",    &evtNum);
     intree->SetBranchAddress("metfilter", &metfilter);
     intree->SetBranchAddress("npv",       &npv);
     intree->SetBranchAddress("npu",       &npu);
     intree->SetBranchAddress("njets",   &njets);
     //    intree->SetBranchAddress("njets08",   &njets08);
     intree->SetBranchAddress("nfjets",        &nfjets);
     intree->SetBranchAddress("nbjets",    &nbjets);
     intree->SetBranchAddress("nbljets",    &nbljets);
     intree->SetBranchAddress("rho",       &rho);
     intree->SetBranchAddress("scale1fb",  &scale1fb);
     intree->SetBranchAddress("evtWeight", &evtWeight);
     intree->SetBranchAddress("pfmet",     &pfmet);
     intree->SetBranchAddress("pfmetphi",  &pfmetphi);
     intree->SetBranchAddress("pfmetraw",     &pfmetraw);
     intree->SetBranchAddress("pfmetphiraw",  &pfmetphiraw);
     intree->SetBranchAddress("pfmetnolep",     &pfmetnolep_in);
     intree->SetBranchAddress("pfmetphinolep",  &pfmetphinolep_in);
     if(phosel){
       intree->SetBranchAddress("pfmetnopho",     &pfmetnopho);
       intree->SetBranchAddress("pfmetphinopho",  &pfmetphinopho);
     }
     intree->SetBranchAddress("mvamet",      &mvamet);
     intree->SetBranchAddress("mvametphi",   &mvametphi);
     if(phosel){
       intree->SetBranchAddress("scEt", &scEt);
       intree->SetBranchAddress("nloosepho", &nloosepho);
       intree->SetBranchAddress("ntightpho", &ntightpho);
       intree->SetBranchAddress("photon", &phovec); //highest ET tight photon
       intree->SetBranchAddress("scEta", &scEta);
       intree->SetBranchAddress("sieie", &sieie);
       intree->SetBranchAddress("chHadIso", &chHadIso);
       intree->SetBranchAddress("neuHadIso", &neuHadIso);
       intree->SetBranchAddress("phoIso", &phoIso);
       intree->SetBranchAddress("sthovere", &sthovere);
       intree->SetBranchAddress("passElectronVeto", &passElectronVeto);
       intree->SetBranchAddress("photypeBits", &photypeBits);
     }

     // lepton variables
     //     cout<<"qui3"<<endl;
     intree->SetBranchAddress("electron",    &electron_arr,    &electron_br);
     intree->SetBranchAddress("muon",    &muon_arr,    &muon_br);
     intree->SetBranchAddress("tau",    &tau_arr,    &tau_br);
     intree->SetBranchAddress("nlep",   &nlep);
     intree->SetBranchAddress("nele",   &nele);
     intree->SetBranchAddress("nmuo",   &nmuo);
     intree->SetBranchAddress("ntau",   &ntau);
     cout<<"qui4"<<endl;
     
     if(SaveGenLevelTMVA && infilename.find("Summer12_TTJets_SemiLeptMGDecays_allhadrbits") != string::npos){
      intree->SetBranchAddress("genId1", &genId1);
      intree->SetBranchAddress("genId2", &genId2);
      intree->SetBranchAddress("genId3", &genId3);
      intree->SetBranchAddress("genpar1", &genpar1);
      intree->SetBranchAddress("genpar2", &genpar2);
      intree->SetBranchAddress("genpar3", &genpar3);
    }



     if(SaveGenLevel && !phosel){
       intree->SetBranchAddress("lepg",    &lepg_arr,    &lepg_br);
       intree->SetBranchAddress("qg",    &qg_arr,    &qg_br);
       intree->SetBranchAddress("dmg",    &dmg_arr,    &dmg_br);
       intree->SetBranchAddress("topg",    &topg_arr,    &topg_br);
       intree->SetBranchAddress("bg",    &bg_arr,    &bg_br);
       intree->SetBranchAddress("nlepg",   &nlepg);
       intree->SetBranchAddress("nqg",   &nqg);
       intree->SetBranchAddress("ndmg",   &ndmg);
       intree->SetBranchAddress("ntopg",   &ntopg);
       intree->SetBranchAddress("nbg",   &nbg);
       intree->SetBranchAddress("lepgpdgid",   lepgpdgid);
       intree->SetBranchAddress("qgpdgid",   &qgpdgid);  
       intree->SetBranchAddress("dmgpdgid",   &dmgpdgid);
       intree->SetBranchAddress("topgpdgid",   &topgpdgid);
       intree->SetBranchAddress("bgpdgid",   &bgpdgid);
       intree->SetBranchAddress("lepgdaugtlep",&lepgdaugtlep); //keeps track of whether is a tau leptonic(12) or tau hadronic (1)
     }
    if(SaveGenLevel && phosel){
      intree->SetBranchAddress("phog", &phog_arr, &phog_br);
      intree->SetBranchAddress("nphog", &nphog);
    }

     //cout<<"qui5"<<endl;
     // boosted
    intree->SetBranchAddress("bst_tau21",    &bst_tau21);
    intree->SetBranchAddress("bst_tau32",     &bst_tau32);
    
    intree->SetBranchAddress("bst_jetcsv",   &bst_jetcsv);
    intree->SetBranchAddress("bst_jetmprun", &bst_jetmprun);
    intree->SetBranchAddress("bst_jetmtrim", &bst_jetmtrim);
    intree->SetBranchAddress("bst_jetmsd",   &bst_jetmsd);
    intree->SetBranchAddress("bst_jet", &bst_jet);

    intree->SetBranchAddress("bst_sub1csv", &bst_sub1csv);
    intree->SetBranchAddress("bst_sub2csv", &bst_sub2csv);
    intree->SetBranchAddress("bst_sub3csv", &bst_sub3csv);
    intree->SetBranchAddress("bst_sub4csv", &bst_sub4csv);
    intree->SetBranchAddress("bst_subjet1", &bst_subjet1);
    intree->SetBranchAddress("bst_subjet2", &bst_subjet2);
    intree->SetBranchAddress("bst_subjet3", &bst_subjet3);
    intree->SetBranchAddress("bst_subjet4", &bst_subjet4);


     // resolved
     intree->SetBranchAddress("res_wmva",     &res_wmva);
     intree->SetBranchAddress("res_topmva",   &res_topmva);
     intree->SetBranchAddress("res_jet1csv",  &res_jet1csv);
     intree->SetBranchAddress("res_jet2csv",  &res_jet2csv);
     intree->SetBranchAddress("res_jet3csv",  &res_jet3csv);
     if(!MVARESTOPTAGGER){
       intree->SetBranchAddress("res_jet4csv",  &res_jet4csv);
       intree->SetBranchAddress("res_jet5csv",  &res_jet5csv);
       intree->SetBranchAddress("res_jet6csv",  &res_jet6csv);
     }
     intree->SetBranchAddress("res_jet1qgid", &res_jet1qgid);
     intree->SetBranchAddress("res_jet2qgid", &res_jet2qgid);
     intree->SetBranchAddress("res_jet3qgid", &res_jet3qgid);
     if(!MVARESTOPTAGGER){
       intree->SetBranchAddress("res_jet4qgid", &res_jet4qgid);
       intree->SetBranchAddress("res_jet5qgid", &res_jet5qgid);
       intree->SetBranchAddress("res_jet6qgid", &res_jet6qgid);
     }
     intree->SetBranchAddress("res_jet1pull", &res_jet1pull);
     intree->SetBranchAddress("res_jet2pull", &res_jet2pull);
     intree->SetBranchAddress("res_jet3pull", &res_jet3pull);
     intree->SetBranchAddress("res_jet1", &res_jet1);
     intree->SetBranchAddress("res_jet2", &res_jet2);
     intree->SetBranchAddress("res_jet3", &res_jet3);
     if(!MVARESTOPTAGGER){
       intree->SetBranchAddress("res_jet4", &res_jet4);
       intree->SetBranchAddress("res_jet5", &res_jet5);
       intree->SetBranchAddress("res_jet6", &res_jet6);
     }
     //     cout<<"qui6"<<endl;

    if(MVARESTOPTAGGER){
      intree->SetBranchAddress("res_prob",     &res_prob);
      intree->SetBranchAddress("res_chisq",    &res_chisq);
      intree->SetBranchAddress("res_cost",     &res_cost);
      intree->SetBranchAddress("res_fitmass",  &res_fitmass);
      intree->SetBranchAddress("res_fitmassW", &res_fitmassW);
      intree->SetBranchAddress("res_bdt_dphij1b", &res_bdt_dphij1b_out);
      intree->SetBranchAddress("res_bdt_dphij2b", &res_bdt_dphij2b_out);
      intree->SetBranchAddress("res_bdt_drj1b", &res_bdt_drj1b_out);
      intree->SetBranchAddress("res_bdt_drj2b", &res_bdt_drj2b_out);
    }

    //    intree->SetBranchAddress("inc_mT2W",    &inc_mT2W);

    
    intree->SetBranchAddress("passtrigger",       &passtrigger);
    intree->SetBranchAddress("cutflow",    &cutflow);
     //cout<<"qui7"<<endl;




     if(verbose>=1)
       cout<<"intree->GetEntries() "<<intree->GetEntries()<<endl;
     vector<float> elept, muopt, taupt, lepggpt;
     for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
       //       if(ientry>10000) continue;
       //     cout<<"****ientry "<<ientry<<endl;
       intree->GetEntry(ientry);

       if(verbose)
	 cout<<"runNum "<<runNum<<" evtNum "<<evtNum<<endl;

       /////////////////////////////////////////////
       ///////////INITIALIZE POINTERS//////////
       /////////////////////////////////////////////
       //Don't remove this: re-initialize pointers to 0 if no jet is found -> without this it will bias dphiminjet
       if(njets<6){
	 res_jet6 = 0;
	 if(njets<5){
	   res_jet5 = 0;
	   if(njets<4){
	     res_jet4 = 0;
	     if(njets<3){
	       res_jet3 = 0;
	       res_jet3pull=0;
	       if(njets<2){
		 res_jet2 = 0;
		   res_jet2pull=0;
		 if(njets<1){
		   res_jet1 = 0;
		   res_jet1pull=0;
		 }//if(njets<1)		 
	       }//if(njets<2)
	     }//if(njets<3)
	   }//if(njets<4)
	 }//if(njets<5)
	 }//if(njets<6)
       /////////////////////////////////////////////
       /////////////////////////////////////////////

	 

       //       cout<<"BBB njets "<<njets<<" res_jet "<<res_jet1<<", "<<res_jet2<<", "<<res_jet3<<", "<<res_jet4<<", "<<res_jet5<<", "<<res_jet6<<", "<<endl;



       //     cout<<"qui8"<<endl;


       //       if(evtNum != 2074643) continue;
       //       if(ientry>3) continue;
       //cout<<"****ientry "<<ientry<<" evtNum "<<evtNum<<endl;

     if(verbose>=1){
      cout<<"electron_br "<<electron_br<<endl;
      cout<<"muon_br "<<muon_br<<endl;
      cout<<"tau_br "<<tau_br<<endl;
     }
     if(electron_arr){
       electron_arr->Clear(); 
       electron_br->GetEntry(ientry);
     }
     else{
       cout<<"*******WARNING****** electron branch not set"<<endl;
     }
     //     cout<<"AAA muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
     if(muon_arr){
       muon_arr->Clear(); 
       //cout<<"BBB muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
       muon_br->GetEntry(ientry);
       //cout<<"CCC muon_arr->GetEntriesFast() "<<muon_arr->GetEntriesFast()<<endl;
     }
     else{
       cout<<"*******WARNING****** muon branch not set"<<endl;
     }

     if(tau_arr){
       tau_arr->Clear();
       tau_br->GetEntry(ientry); 
     }
     else{
       cout<<"*******WARNING****** tau branch not set"<<endl;
     }

     //cout<<"qui9"<<endl;
     
     elept.clear(); muopt.clear(); taupt.clear(); lepggpt.clear();
     //cout<<"BBBBBBB "<<muo.size()<<endl; 
     ele.clear(); /*muo.clear();*/ ta.clear(); lepgg.clear(); phogg.clear();
     //cout<<"CCCCCC "<<muo.size()<<endl; 
     muo.erase(muo.begin(),muo.begin()+muo.size());
     //cout<<"DDDD "<<muo.size()<<endl; 
     muo.clear();
     //cout<<"EEEE "<<muo.size()<<endl; 
     if(electron_arr){ 
       //       for(int ie=0; ie<electron_arr->GetEntriesFast(); ie++) { //FOUND BUG WITH GetEntriesFast(). in runAllHadronic.cc I have a check nele != vEle.sizee() exit(1)
       for(int ie=0; ie<nele; ie++) {
	 TLorentzVector *electron_tmp = (TLorentzVector*)electron_arr->At(ie);
	 ele.push_back(*electron_tmp);
	 elept.push_back(electron_tmp->Pt());
	 if(verbose)
	   cout<<"electron Pt "<<electron_tmp->Pt()<<endl;
       }      
     }//if(electron_arr)
     std::sort (ele.begin(), ele.begin()+ele.size(), lorentzvecptorder);
     std::sort (elept.begin(), elept.begin()+elept.size(), floatorder);

     if(nele != ele.size()){
       cout<<"nele "<<nele<<" ele.size() "<<ele.size()<<" are different..gonna cause problems"<<endl;
       cout<<"runNum, evtNum "<<runNum<<" "<<evtNum<<endl;
       for (std::vector<TLorentzVector>::iterator it = ele.begin() ; it < ele.end(); ++it){
	 cout<<"Pt, eta "<<it->Pt()<<" "<<it->Eta()<<endl;
       }
       exit(1);
     }

     if(ele.size()>=2){
       if(ele[0].Pt()<ele[1].Pt()){ 
	 cerr<<"AAAAAAA Wrong ordering "<<ele.size()<<" ele 0, 1 "<<ele[0].Pt()<<" "<<ele[1].Pt()<<endl;
	 exit(1);
       }
     }

     if(muon_arr){
       //       for(int im=0; im<muon_arr->GetEntriesFast(); im++) { //FOUND BUG WITH GetEntriesFast(). nmuo != vMuo.size() exit(1): for event with nmuo=1 returned two muons, where the second muon ws the second muon of the previous stored event
       for(int im=0; im<nmuo; im++) {
	 TLorentzVector *muon_tmp = (TLorentzVector*)muon_arr->At(im);
	 //	 cout<<"muon Pt "<<muon_tmp->Pt()<<endl;
	 muo.push_back(*muon_tmp);
	 muopt.push_back(muon_tmp->Pt());
	 if(verbose)
	   cout<<"muon Pt "<<muon_tmp->Pt()<<endl;


	 	 
	 //cout<<"FFFF "<<muo.size()<<" "<<muo[0].Pt()<<" "<< muo[1].Pt()<<endl;
       }
     }//if(muon_arr)
     std::sort (muo.begin(), muo.begin()+muo.size(), lorentzvecptorder);
     std::sort (muopt.begin(), muopt.begin()+muopt.size(), floatorder);

     if(nmuo != muo.size()){
       cout<<"nmuo "<<nmuo<<" muo.size() "<<muo.size()<<" are different..gonna cause problems"<<endl;
       cout<<"runNum, evtNum "<<runNum<<" "<<evtNum<<endl;
       for (std::vector<TLorentzVector>::iterator it = muo.begin() ; it < muo.end(); ++it){
	 cout<<"Pt, eta "<<it->Pt()<<" "<<it->Eta()<<endl;
       }
       exit(1);
     }

     if(muo.size()>=2){
       if(muo[0].Pt()<muo[1].Pt()){ 
	 cerr<<"AAAAAAA Wrong ordering "<<muo.size()<<" muo 0, 1 "<<muo[0].Pt()<<" "<<muo[1].Pt()<<endl;
	 exit(1);
       }
     }
     //cout<<"qui2"<<endl;
     if(tau_arr){
       //for(int it=0; it<tau_arr->GetEntriesFast(); it++) { //FOUND BUG WITH GetEntriesFast(). nmuo != vMuo.size() exit(1): for event with nmuo=1 returned two muons, where the second muon ws the second muon of the previous store d event                            
       for(int it=0; it<ntau; it++) { 
	 TLorentzVector *tau_tmp = (TLorentzVector*)tau_arr->At(it);
	 ta.push_back(*tau_tmp);
	 taupt.push_back(tau_tmp->Pt());
	 if(verbose)
	  cout<<"tau Pt "<<tau_tmp->Pt()<<endl;
       }
     }//if(tau_arr)
     std::sort (ta.begin(), ta.begin()+ta.size(), lorentzvecptorder);
     std::sort (taupt.begin(), taupt.begin()+taupt.size(), floatorder);
     //     cout<<"after muo0 muo1"<<muo[0].Pt()<<", "<<muo[1].Pt()<<endl;

     if(ntau != ta.size()){
       cout<<"ntau "<<nmuo<<" ta.size() "<<ta.size()<<" are different..gonna cause problems"<<endl;
       cout<<"runNum, evtNum "<<runNum<<" "<<evtNum<<endl;
       for (std::vector<TLorentzVector>::iterator it = ta.begin() ; it < ta.end(); ++it){
	 cout<<"Pt, eta "<<it->Pt()<<" "<<it->Eta()<<endl;
       }
       exit(1);
     }

     if(ta.size()>=2){
       if(ta[0].Pt()<ta[1].Pt()){ 
	 cerr<<"AAAAAAA Wrong ordering "<<ta.size()<<" ta 0, 1 "<<ta[0].Pt()<<" "<<ta[1].Pt()<<endl;
	 exit(1);
       }
     }

     //TMVA top tagger: assuming that res_jet1, 2, 3 are the triplet associated to highest MVA score
     //assuming the genpar1, genpar2, genpar3 are from the hadronic top
     //count the number of match -> if three found hadronic top
     if(SaveGenLevelTMVA && infilename.find("Summer12_TTJets_SemiLeptMGDecays_allhadrbits") != string::npos)
       isHadrTopHighestMVA = isHadronicTop(res_jet1, res_jet2, res_jet3, genpar1, genpar2, genpar3, genId1, genId2, genId3, verbose);
     else
       isHadrTopHighestMVA = -1.;

     //     cout<<"isHadrTopHighestMVA "<<isHadrTopHighestMVA<<endl;
     //



     if(SaveGenLevel && !phosel){
       if(lepg_arr){
	 //	 for(int il=0; il<lepg_arr->GetEntriesFast(); il++) { //FOUND BUG WITH GetEntriesFast()
	 for(int il=0; il<nlepg; il++) {
	   TLorentzVector *lepgg_tmp = (TLorentzVector*)lepg_arr->At(il);
	   lepgg.push_back(*lepgg_tmp);
	   lepggpt.push_back(lepgg_tmp->Pt());
	   if(verbose)
	     cout<<"lepgg Pt "<<lepgg_tmp->Pt()<<endl;
	 }
       }//if(lepg_arr)
     }

     if(SaveGenLevel && phosel){
       if(phog_arr){
	 //	 cout<<"nphog "<<nphog<<endl;
	 for(int ipho=0; ipho<nphog; ipho++) {
	   TLorentzVector *phog_tmp = (TLorentzVector*)phog_arr->At(ipho);
	   phogg.push_back(*phog_tmp);
	   if(verbose)
	     cout<<"phogg Pt "<<phog_tmp->Pt()<<endl;
	 }
       }//if(phog_arr)



       for(int ipho=0; ipho<phog_arr->GetEntriesFast(); ipho++) {
	 TLorentzVector *phog_tmp = (TLorentzVector*)phog_arr->At(ipho);
	 new(phog[ipho]) TLorentzVector(*phog_tmp);
	 if(verbose)
	   cout<<" phog E "<<phog_tmp->E()<<endl;
       }      
       
     }//if(SaveGenLevel && !phosel)

     if(verbose)
       cout<<"Starting to calculate my own var"<<endl;
        
     //CALC NEW VARIABLE WHICH SHOULD BE PORTED IN .CC 
     //compute the amount of bjets passing the medium work point - to be transferred in the .cc
     nbljets_mt=0;
     if(njets>0){
       if(res_jet1 != 0){
	 if(res_jet1csv>BTAG_WPL && fabs(res_jet1->Eta())<2.4)
	   nbljets_mt++;
       }
     }
     if(njets>1){
       if(res_jet2 != 0){
	 if(res_jet2csv>BTAG_WPL  && fabs(res_jet2->Eta())<2.4)
	   nbljets_mt++;
       }
     }
     if(njets>2){
       if(res_jet3 != 0){
	 if(res_jet3csv>BTAG_WPL  && fabs(res_jet3->Eta())<2.4)
	   nbljets_mt++;
       }
     }
       //	 cout<<"qui5 njets "<<njets<<" "<<res_jet4<<" "<<res_jet5<<" "<<res_jet6<<endl;
     if(njets>3){
       if(res_jet4 != 0){
	 if(res_jet4csv>BTAG_WPL  && fabs(res_jet4->Eta())<2.4)
	   nbljets_mt++;
       }
     }
     if(njets>4){
       if(res_jet5 != 0){
	 if(res_jet5csv>BTAG_WPL  && fabs(res_jet5->Eta())<2.4)
	   nbljets_mt++;
       }
     }
       //       cout<<"qui6"<<endl;
     if(njets>5){
       if(res_jet6 != 0){
	 if(res_jet6csv>BTAG_WPL  && fabs(res_jet6->Eta())<2.4)
	   nbljets_mt++;
       }
     }
     //cal deltaphi between met and the closest jet
     if(res_jet1 != 0 && res_jet2 !=0 && njets>=2)
       mindphimetjet12 = CalcMinDPhiJet(pfmetphi,res_jet1,res_jet2);
     else
       mindphimetjet12 = 999.;
     
     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 != 0 && njets>=3)
	 mindphimetjet123 = CalcMinDPhiJet(pfmetphi,res_jet1,res_jet2,res_jet3);
     else
       mindphimetjet123 = 999.;

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 != 0 && res_jet4 != 0 && njets>=4)
       mindphimetjet1234 = CalcMinDPhiJet(pfmetphi,res_jet1,res_jet2,res_jet3,res_jet4);
     else
       mindphimetjet1234 = 999.;

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 != 0 && res_jet4 != 0 && res_jet5 != 0 && njets>=5)
       mindphimetjet12345 = CalcMinDPhiJet(pfmetphi,res_jet1,res_jet2,res_jet3,res_jet4,res_jet5); 
     else
       mindphimetjet12345 = 999.;

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 != 0 && res_jet4 != 0 && res_jet5 != 0 && res_jet6 != 0 && njets>=6)
       mindphimetjet = CalcMinDPhiJet(pfmetphi,res_jet1,res_jet2,res_jet3,res_jet4,res_jet5,res_jet6);
     else
       mindphimetjet = 999.;



     

     //     cout<<"mindphimetjet "<<mindphimetjet<<endl;

	      
     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 !=0 && njets>=3)
       topmass = calcmass(res_jet1,res_jet2,res_jet3); //Calctopmass(njets, res_jet1,res_jet2,res_jet3,res_jet4,res_jet5,res_jet6);
     else
       topmass =  -999.;  

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 !=0 && njets>=3) //assuming I have a top with at least 3 jets and want to compute the W mass
       wmass = calcwmass(res_jet1,res_jet2,res_jet3,res_jet1csv,res_jet2csv,res_jet3csv); 
     else
       wmass =  -999.;  
     

     
       //pfmt meant only for nlep=1
     TLorentzVector vLep, vLep2; //highest energy lep
     vLep.SetPtEtaPhiM(0.,0,0,0);
     vLep2.SetPtEtaPhiM(0.,0,0,0);
     //     TLorentzVector vEle[nele], vMuo[nmuo];
     TLorentzVector vEle[10], vMuo[10];
     //       cout<<"---------"<<endl;
     //cout<<"nele "<<nele<<" nmuo "<<nmuo<<endl;
     for(int ie=0; ie<nele; ie++){
       //cout<<"ele "<<ie<<" Pt, Eta "<<ele[ie].Pt()<<", "<<ele[ie].Eta()<<endl;
       vEle[ie] = ele[ie];
       if(vEle[ie].Pt()>vLep.Pt() || vLep.Pt()==0){
	 vLep2 = vLep;
	   vLep = vEle[ie];	   
       }
       else if (vEle[ie].Pt()>vLep2.Pt() || vLep2.Pt()==0)	   
	 vLep2 = vEle[ie];
     }
     for(int im=0; im<nmuo; im++){
       //cout<<"muo "<<im<<" Pt, Eta "<<muo[im].Pt()<<", "<<muo[im].Eta()<<endl;
       vMuo[im] = muo[im];
       if(vMuo[im].Pt()>vLep.Pt() || vLep.Pt()==0){
	 vLep2 = vLep;
	 vLep = vMuo[im];
       }
       else if (vMuo[im].Pt()>vLep2.Pt() || vLep2.Pt()==0)
	 vLep2 = vMuo[im]; 
     }
     //if(vLep.Pt() != 0)
     // cout<<"highest Et lep Pt, Eta "<<vLep.Pt()<<", "<<vLep.Eta()<<endl;
     //if(vLep2.Pt() != 0)
     //	 cout<<"2nd highest Et lep Pt, Eta "<<vLep2.Pt()<<", "<<vLep2.Eta()<<endl;
     //cout<<"---------"<<endl;
     
     
     /*       if(ta.size()>0){
	      if(ele.size()>0 && ta[0].Pt() > ele[0].Pt())
	      vLep = ta[0];
	      if(muo.size()>0 && vLep.Pt() < muo[0].Pt())
	      vLep = muo[0];
	      }*/
     //       cout<<"ele.size() "<<ele.size()<<"muo.size() "<<muo.size()<<"ta.size() "<<ta.size()<<endl;
     
     if(verbose)
       cout<<"nmuo "<<nmuo<<endl;
     
     

     if(nele==1 || nmuo==1){//should be nlep but be careful about ntau
       pfmt = computeMt(vLep.Pt(), vLep.Phi(), pfmet, pfmetphi);
       deltaphilepmet = deltaPhi(vLep.Phi(),pfmetphi);
     }
     else{
       pfmt = -999.;
       deltaphilepmet = -999.;
     }

     //     cout<<"quiiiii"<<endl;

     TVector2 lepvect2(-999.,-999.);
     TVector2 pfmetvect2(pfmet*TMath::Cos(pfmetphi),pfmet*TMath::Sin(pfmetphi));
     for(int ie=0; ie<nele; ie++){
       lepvect2.Set(vEle[ie].Px(),vEle[ie].Py());
       pfmetvect2 += lepvect2;
     }
     for(int im=0; im<nmuo; im++){
       lepvect2.Set(vMuo[im].Px(),vMuo[im].Py());
       pfmetvect2 += lepvect2;
     }
     
     //cout<<"quiiiii2"<<endl;     

     
     
     pfmetnolep = pfmetvect2.Mod(); //should be equal to pfmet when no lep - check
     pfmetphinolep = pfmetvect2.Phi();

     //       cout<<"nele "<<nele<<" nmuo "<<nmuo<<endl;
     /*     if(TMath::Abs(pfmetnolep-pfmetnolep_in)>0.001){ 
       cout<<"ERROR"<<endl;
       cout<<"pfmetnolep "<<pfmetnolep<<" pfmetphinolep "<<pfmetphinolep<<endl;
       cout<<"pfmetnolep_in "<<pfmetnolep_in<<" pfmetphinolep_in "<<pfmetphinolep_in<<endl;

       }*/

     //cout<<"quiiiii3"<<endl;     

       if(nele>=1 || nmuo>=1)//should be nlep but be careful about ntau
	 pfmtnolepinmet = computeMt(vLep.Pt(), vLep.Phi(), pfmetnolep, pfmetphinolep);
       else
	 pfmtnolepinmet = -999.;
       deltaphilepmetnolepinmet = deltaPhi(vLep.Phi(),pfmetphinolep);


       //       cout<<"qui2 "<<nmuo<<endl; 
       //       cout<<"quiii "<<njets<<endl;

       //creates problem with the new run for TMVA tagger 
       /*       if(nele==1 || nmuo==1){
	 float mindeltarlepjettmp = 999.;

	 if(njets>0 && res_jet1 != 0){

	   mindeltarlepjettmp = vLep.DeltaR(*res_jet1); 
	   if(njets>1 && res_jet2 != 0){

	     if(vLep.DeltaR(*res_jet2)<mindeltarlepjettmp)
	       mindeltarlepjettmp = vLep.DeltaR(*res_jet2);
	     if(njets>2 && res_jet3 != 0){

	       if(vLep.DeltaR(*res_jet3)<mindeltarlepjettmp)
		 mindeltarlepjettmp = vLep.DeltaR(*res_jet3);
	       cout<<"qui3 "<<res_jet4<<endl; 
	       if(njets>3 && res_jet4 != 0){
	       cout<<"res_jet4->Pt()"<<res_jet4->Pt()<<endl;
	       cout<<"---"<<endl;

		        cout<<"quiii4 "<<endl;
		 if(vLep.DeltaR(*res_jet4)<mindeltarlepjettmp)
		   mindeltarlepjettmp = vLep.DeltaR(*res_jet4);
		 		 cout<<"qui4 "<<endl; 
		 if(njets>4 && res_jet5 != 0){
		        cout<<"quiii5 "<<endl;
		   if(vLep.DeltaR(*res_jet5)<mindeltarlepjettmp)
		     mindeltarlepjettmp = vLep.DeltaR(*res_jet5);
		   if(njets>5 && res_jet6 != 0){
		            cout<<"quiii6 "<<endl;
		     if(vLep.DeltaR(*res_jet6)<mindeltarlepjettmp)
		       mindeltarlepjettmp = vLep.DeltaR(*res_jet6);
		   }	   
		 }//if(njets>4)
	       }//if(njets>3)
	     }//f(njets>2)
	   }//if(njets>1)
	 }//if(njets>0)
	 mindeltarlepjet = mindeltarlepjettmp;
	 }//if(nele==1 || nmuo==1)*/

       //     cout<<"quiiiii4"<<endl;     

       if(nele==2 || nmuo==2){ //to check
	 dilepmass = (vLep+vLep2).M();
	 //	 cout<<"dilepmass "<<dilepmass<<endl;
       }
       else
	 dilepmass = -999.;
       //END CALC NEW VARIABLE WHICH SHOULD BE PORTED IN .CC
       //       cout<<"quiiifine"<<endl;


       //       cout<<"sample->label "<<sample->label<<endl;
       //cout<<"sample->label.compare QCD "<<(sample->label).compare("QCD")<<endl;

       if(verbose)
	 cout<<"Selecting events"<<endl;


       //if(metfilter!=0)    continue;
       
   
       if(useturnonpar==-1); //do nothing for MC
       if(useturnonpar==0){
       	 if(!passtrigger) continue;
       }
       else{//assuming I will be using my own trigger efficiency curve for MC
       	 if(isData /*|| (sample->label).compare("QCD")==0*/){
       	   if(!passtrigger) continue;
       	 }
       }

   
//       if((cutflow & 0x0000000000000100) != 0x0000000000000100) continue;
       //if(pfmet < 120) continue;
	
//       if( (cutflow & 0x0000000000050000) == 0x0000000000050000); //resolved
//       else continue;


     // if(evtNum != 571980070) continue;

     // cout<<"Lep Pt, Phi, pfmet, pfmetphi "<<vLep.Pt()<<", "<<vLep.Phi()<<", "<<pfmet<<", "<< pfmetphi<<endl;
     // cout<<"passtrigger, njets, nbjets, pfmet, nmuo, pfmt "<<passtrigger<<", "<<njets<<", "<<nbjets<<", "<<pfmet<<", "<<nmuo<<", "<<pfmt<<endl;



       //if(njets  < 4)  continue;	
       if(njets  < 3)  continue;	
       

       //SR
       //if(nlep != 0) continue;
       if(nbjets < 2) continue;
       //if(nbjets !=  1) continue;
       //if(TMath::Abs(lepg->Eta())>2.5) continue; 
       //if(lepg->Pt()<30) continue;
       //       if(mindphimetjet<1.0) continue;
       if(pfmet<100) continue;

       //CR - QCD no btag WlP	
       //if(nbljets_mt>0) continue;
       //       if(mindphimetjet>0.1) continue;
       //	if(nbjets!=0) continue; 
       //if(mindphimetjet>0.2) continue;

       //CR tt 1 lep & >=2btag 
       //if(pfmetnolep  < 180)    continue;
       if(nmuo != 1) continue;
       //       if(muo[0].Pt()<100) continue;
       //if(nele != 1) continue;
       //if(nlep != 1) continue;
       //       if(ntau != 0) continue;
       //if(pfmtnolepinmet>100) continue;
       if(pfmt<40) continue;
       //       if(isHadrTopHighestMVA_cut == 0 || isHadrTopHighestMVA_cut == 1)
       //if(isHadrTopHighestMVA!=isHadrTopHighestMVA_cut) continue;
       //if(res_topmva<res_topmva_cut) continue;	

       //W+jets 1 lep & 0btag  
       //if(nlep != 1) continue;
       //	if(nbljets_mt!=0) continue; 
       //if(nbjets!=0) continue; 
       //if(nbljets_mt!=1) continue; 
       //if(mindphimetjet<1) continue;
       //       if(nbjets <  1) continue;
       //if(topmass>160 && topmass<200) continue;

       //Z(ll)+jets
       //       if(nmuo != 2) continue;
       //if(nele != 2) continue;
       //if(nbljets_mt>0) continue;
       //if(pfmetnolep  < 120)    continue;

       //dilepton for Z peak
       //if(nmuo != 2) continue; 
       //       if(nele != 2) continue; 
       //if(pfmet>50) continue;

       //gamma+jets
       // if(ntightpho != 1) continue;
       // if((scEta <= 1.479 && sieie <= 0.011)||(scEta > 1.479 && sieie <= 0.033));
       // else continue;
       // if(phovec->Pt()<170) continue;
       // if(TMath::Abs(phovec->Eta())>2.5 || (TMath::Abs(phovec->Eta())>1.442 && TMath::Abs(phovec->Eta())<1.56)) continue;
       // if(pfmetnopho<180) continue;
       //

       if(verbose)
	 cout<<"Event Passed selection evtNum "<<evtNum<<endl;


      
       TLorentzVector multijet; 
       if(res_jet1 != 0 ) {
	 multijet = (*res_jet1);
	 if(res_jet2 != 0)
	   multijet += (*res_jet2);
	 if(res_jet3 != 0)
	   multijet += (*res_jet3);
	 if(res_jet4 != 0)
	   multijet += (*res_jet4);
	 if(res_jet5 != 0)
	   multijet += (*res_jet5);
	 if(res_jet6 != 0)
	   multijet += (*res_jet6);
	 dphimetmultijets = deltaPhi(pfmetphi,multijet.Phi());
       }
       else
	 dphimetmultijets = -999.;
	   
       //       cout<<"qui5"<<endl;

       if(SaveGenLevel && !phosel){
	 TLorentzVector lep_tmp, nu_tmp; 
	 int pdglep_tmp=-999.;
	 int pdgnu_tmp=-999.; //assuming only one leptonically-decaying W, and looking at only e,mu
	 for(int il=0; il<4; il++){
	   //cout<<"il "<<il<<" lepgpdgid "<<lepgpdgid[il]<<endl;
	   if(TMath::Abs(lepgpdgid[il]) == 11 || TMath::Abs(lepgpdgid[il]) == 13 || TMath::Abs(lepgpdgid[il]) == 15){
	     lep_tmp = lepgg[il];
	     //	     cout<<"A lep_tmp il "<<il<<" "<<lep_tmp.Pt()<<endl;
	     pdglep_tmp = lepgpdgid[il];
	   }
	   else if(TMath::Abs(lepgpdgid[il]) == 12 || TMath::Abs(lepgpdgid[il]) == 14 || TMath::Abs(lepgpdgid[il]) == 16){
	     nu_tmp = lepgg[il];
	     //cout<<"A nu_tmp il "<<il<<" "<<nu_tmp.Pt()<<endl;
	     pdgnu_tmp = lepgpdgid[il];
	   }	   
	 }//for(int il=0; il<4; il++)	 
	 //cout<<"pdglep_tmp "<<pdglep_tmp<<endl;
	 //cout<<"pdgnu_tmp "<<pdgnu_tmp<<endl;
	 if((pdglep_tmp == 11 && pdgnu_tmp == -12) || (pdglep_tmp == 13 && pdgnu_tmp == -14) || (pdglep_tmp == -11 && pdgnu_tmp == 12) ||  (pdglep_tmp == -13 && pdgnu_tmp == 14) ||  (pdglep_tmp == -15 && pdgnu_tmp == 16) ||  (pdglep_tmp == 15 && pdgnu_tmp == -16)){
	   //cout<<"lep_tmp "<<lep_tmp.Pt()<<endl;
	   eleormuogpt = lep_tmp.Pt();
	   eleormuogeta = lep_tmp.Eta();
	   //cout<<"nu_tmp "<<nu_tmp.Pt()<<endl;
	   nugpt = nu_tmp.Pt();
	   nugeta = nu_tmp.Eta();
	   wgpt = (lep_tmp+nu_tmp).Pt();
	   wgeta = (lep_tmp+nu_tmp).Eta();
	 }
	 else{
	 eleormuogpt = -999.;
	 eleormuogeta = -999.; 
	 nugpt = -999.;
	 nugeta = -999.;
	 wgpt = -999.;
	 wgeta = -999.;

	 //taugpt = -999.;
	 //	 taugeta = -999.;
	 }
       }//if(SaveGenLevel && !phosel)
	 
       if(SaveGenLevel && phosel){
	 TLorentzVector pho_tmp;
	 pho_tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	 //	 cout<<"nphog "<<nphog<<endl;
	 for(int ipho=0; ipho<nphog; ipho++){
	   phogpt=-999.; phogdR=999.;
	   if(phovec){ 
	     float dR = phovec->DeltaR(phogg[ipho]);
	     //  cout<<"dR gen pho reco pho "<<dR<<endl;
	     //cout<<"gen, reco pho pt "<<phogg[ipho].Pt()<<", "<<phovec->Pt()<<endl;
	     if(dR<0.3 && phogg[ipho].Pt()>phogpt){//saving the highest pt gen photon which matches with reco pho
	       //cout<<"here"<<endl;
	       phogpt = phogg[ipho].Pt();
	       phogdR = dR;
	     }
	   }//if(phovec)
	   //cout<<"BBBB gen, reco pho pt "<<phogpt<<", "<<phovec->Pt()<<endl;
	 }//for(int ipho=0; ipho<nphog; ipho++)
       }//if(SaveGenLevel && phosel)

       if(verbose)
	 cout<<"before the weights"<<endl;


       double weight = 1;
       //       if(verbose==2)
       //	 cout<<"isData "<<isData<<endl;
       if(!isData) {     
	 weight *= LUMI*scale1fb;//*evtWeight;
	 weight *= hPUWeights->GetBinContent(hPUWeights->FindBin(npv));
	 //	 weight *= hPUWeights->GetBinContent(hPUWeights->FindBin(npu));
	 //	 cout<<"turnon("<<pfmetnolep<<","<<useturnonpar<<") "<<turnon(pfmetnolep,useturnonpar)<<endl;
	 //if((sample->label).compare("QCD")!=0){ //no trigger weighting for QCD
	 if(useturnonpar==1)
	   weight *= turnon(pfmetnolep,useturnonpar);
	 else if(useturnonpar==2)
	   weight *= turnon(pfmet,useturnonpar);
	   //}

	 //re-weighting to artifically sculpt a given distribution (Wpt in this case)
	 //	 cout<<"extvarweighting "<<extvarweighting<<" "<<extvarweighting.find("WGPT")<<endl;
	 if(extvarweighting.find("WGPT") == 0){
	   weight *= Wptweight(wgpt);
	 }
	 else if(extvarweighting.find("NUGPT_ELE") == 0){
	   weight *= Nuptweight_ele(nugpt);
	 }
	 else if(extvarweighting.find("NUGPT_MUO") == 0){
	   weight *= Nuptweight_muo(nugpt);
	 }
	 else if(extvarweighting.find("NLOCorrPT_GtoZ") == 0){
	   if(phogpt>=0)
	     weight *= NLOcorrphotoZ(phogpt);
	 }
	 
       }
       //       weight = 1;



       if(verbose)
	 cout<<"Define variables to be plotted"<<endl;

       //*******************
       //*******************variables to be plotted
       //*******************
       //make sure the order is the same as DefineListHistos
       listvar.clear();
       listvar.push_back(npu);        listvar.push_back(npu);
       listvar.push_back(npv);        listvar.push_back(npv);
       listvar.push_back(pfmet); listvar.push_back(pfmet); //one log one linear
       listvar.push_back(pfmetphi);       listvar.push_back(pfmetphi);
       listvar.push_back(njets); listvar.push_back(nbjets); listvar.push_back(nbljets_mt);
       //             cout<<"res_jet1csv "<<res_jet1csv<<endl;
       listvar.push_back(res_jet1csv); listvar.push_back(res_jet2csv); listvar.push_back(res_jet3csv); listvar.push_back(res_jet4csv); listvar.push_back(res_jet5csv); listvar.push_back(res_jet6csv);
       listvar.push_back(res_jet1qgid); listvar.push_back(res_jet2qgid); listvar.push_back(res_jet3qgid);
       if(njets>0 && res_jet1 != 0  && res_jet1->Pt() != 0)
       	 listvar.push_back(res_jet1->Eta());
       else
       	 listvar.push_back(-999.);
       if(njets>1 && res_jet2 != 0 && res_jet2->Pt() != 0)
       	 listvar.push_back(res_jet2->Eta()); 
       else
       	 listvar.push_back(-999.);       
       if(njets>2 && res_jet3 != 0 && res_jet3->Pt() != 0)
       	 listvar.push_back(res_jet3->Eta());
       else
       	 listvar.push_back(-999.);
       if(njets>3 && res_jet4 != 0 && res_jet4->Pt() != 0) 
       	 listvar.push_back(res_jet4->Eta());
       else
       	 listvar.push_back(-999.);
       if(njets>4 && res_jet5 != 0 && res_jet5->Pt() != 0) 
       	 listvar.push_back(res_jet5->Eta());
       else
       	 listvar.push_back(-999.);
       if(njets>5 && res_jet6 != 0 && res_jet6->Pt() != 0) 
       	 listvar.push_back(res_jet6->Eta());
       else
       	 listvar.push_back(-999);
       if(njets>0 && res_jet1 != 0)
       	 listvar.push_back(res_jet1->Pt()); 
       else
       	 listvar.push_back(-999.);
       if(njets>1 && res_jet2 != 0)
       	 listvar.push_back(res_jet2->Pt()); 
       else
       	 listvar.push_back(-999.);
       if(njets>2 && res_jet3 != 0)
       	 listvar.push_back(res_jet3->Pt()); 
       else
       	 listvar.push_back(-999.);


       if(njets>3 && res_jet4 != 0) 
       	 listvar.push_back(res_jet4->Pt());
       else
       	 listvar.push_back(-999.);
       if(njets>4 && res_jet5 != 0) 
       	 listvar.push_back(res_jet5->Pt());
       else
       	 listvar.push_back(-999.);
       if(njets>5 && res_jet6 != 0) 
       	 listvar.push_back(res_jet6->Pt());
       else
       	 listvar.push_back(-999.);
       listvar.push_back(TMath::Abs(dphimetmultijets));
       listvar.push_back(TMath::Abs(mindphimetjet));        listvar.push_back(TMath::Abs(mindphimetjet)); 
       listvar.push_back(TMath::Abs(mindphimetjet12));        listvar.push_back(TMath::Abs(mindphimetjet12)); 
       listvar.push_back(TMath::Abs(mindphimetjet123));        listvar.push_back(TMath::Abs(mindphimetjet123)); 
       listvar.push_back(TMath::Abs(mindphimetjet1234));        listvar.push_back(TMath::Abs(mindphimetjet1234)); 
       listvar.push_back(TMath::Abs(mindphimetjet12345));       listvar.push_back(TMath::Abs(mindphimetjet12345));
       listvar.push_back(topmass);        listvar.push_back(topmass);
       listvar.push_back(wmass);        listvar.push_back(wmass);
       // //cout<<"qui9"<<endl;
       if(MVARESTOPTAGGER){
       	 listvar.push_back(res_topmva);	 listvar.push_back(res_topmva);
       	 listvar.push_back(res_prob);
       	 listvar.push_back(res_chisq);
       	 listvar.push_back(res_cost);
       	 listvar.push_back(res_fitmass);
       	 listvar.push_back(res_fitmassW);
       	 listvar.push_back(res_bdt_dphij1b_out);
       	 listvar.push_back(res_bdt_dphij2b_out);
       	 listvar.push_back(res_bdt_drj1b_out);
       	 listvar.push_back(res_bdt_drj2b_out);
       }

       
       //cout<<"qui10"<<endl;
       // if(ele.size() != nele){
       // 	 cerr<<"***************** elesize "<<ele.size()<<" nele "<<nele<<" are different...check"<<endl;
       // 	 exit(1);
       // }
       // if(muo.size() != nmuo){
       // 	 cerr<<"***************** muosize "<<muo.size()<<" nmuo "<<nmuo<<" are different...check"<<endl;
       // 	 exit(1);
       // }
       // if(ta.size() != ntau){
       // 	 cerr<<"***************** tausize "<<ta.size()<<" ntau "<<ntau<<" are different...check"<<endl;
       // 	 exit(1);
       // }
       listvar.push_back(nele);
       listvar.push_back(nmuo);
       listvar.push_back(ntau);
       //         listvar.push_back(nlep);
       if(ele.size()>0){
       	 listvar.push_back(ele[ele.size()-1].Pt());//plotting the lowest pt electron
       	 listvar.push_back(ele[ele.size()-1].Pt());//plotting the lowest pt electron
       	 listvar.push_back(ele[0].Pt());
       	 listvar.push_back(ele[0].Pt());
       }
       else{
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       }
       //       cout<<"lowes muo "<<muo[muo.size()-1].Pt()<<" muo[0].Pt() " << muo[0].Pt()<<endl;
       if(muo.size()>0){
       	 listvar.push_back(muo[muo.size()-1].Pt()); //need to sort the muons (sort it when you build the muon)
       	 listvar.push_back(muo[muo.size()-1].Pt()); //need to sort the muons (sort it when you build the muon)
       	 listvar.push_back(muo[0].Pt()); //need to sort the muons
       	 listvar.push_back(muo[0].Pt()); //need to sort the muons
       }
       else{
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       	 listvar.push_back(-999);
       }
       if(ta.size()>0)
       	 listvar.push_back(ta[ta.size()-1].Pt());
       else
       	 listvar.push_back(-999);

       if(nele==1 || nmuo==1 /*|| ntau==1*/){ //should be nlep but be careful about ntaup
       	 listvar.push_back(pfmt);	 listvar.push_back(pfmt);
       	 listvar.push_back(inc_mT2W);	 listvar.push_back(inc_mT2W);
       	 listvar.push_back(TMath::Abs(deltaphilepmet));
       	 listvar.push_back(pfmtnolepinmet);
       	 listvar.push_back(TMath::Abs(deltaphilepmetnolepinmet));
       }
       else{
       	 listvar.push_back(-999.);        	 listvar.push_back(-999.);
       	 listvar.push_back(-999.);       	 listvar.push_back(-999.);
       	 listvar.push_back(-999.);
       	 listvar.push_back(-999.);
       	 listvar.push_back(-999.);
       }

       listvar.push_back(pfmetnolep);       listvar.push_back(pfmetnolep);
       listvar.push_back(pfmetphinolep);
       listvar.push_back(mindeltarlepjet);
       listvar.push_back(dilepmass);        listvar.push_back(dilepmass);
	 
       if(SaveGenLevel && !phosel){
       // 	 listvar.push_back(nlepg);
       	 listvar.push_back(eleormuogeta);	 listvar.push_back(eleormuogeta);
       	 listvar.push_back(nugeta);	 listvar.push_back(nugeta);
       	 listvar.push_back(wgeta);	 listvar.push_back(wgeta);
       // 	 listvar.push_back(taugeta);
       	 listvar.push_back(eleormuogpt);	 listvar.push_back(eleormuogpt);
       	 listvar.push_back(nugpt);	 listvar.push_back(nugpt);
       	 listvar.push_back(wgpt);	 listvar.push_back(wgpt);
       	 if((nele==1 && nmuo==0) || (nmuo==1 && nele==0)){

       	   //cout<<"runNum "<<runNum<<" evtNum "<<evtNum<<" nele "<<nele<<" nmuo "<<nmuo<<endl;
       	   //	   cout<<"gen eleormuogpt "<<eleormuogpt<<" vLep pt"<<vLep.Pt()<<" "<<eleormuogpt-vLep.Pt()<<endl;
       	   //cout<<"*************************"<<endl;
       	   listvar.push_back(eleormuogpt-vLep.Pt());	 
       	   listvar.push_back(eleormuogpt-vLep.Pt());	 
       	 }
       	 else{
       	   listvar.push_back(-9.);	 
       	   listvar.push_back(-9.);	 
       	 }

       // 	 listvar.push_back(taugpt);
       // 	 listvar.push_back(lepgpdgid);
       // 	 //add here other gen lvel variables that you wanna plot
	 
       }

       if(SaveGenLevel && phosel){
       	 //	 cout<<"AAA phogpt "<<phogpt<<endl; 
       	 listvar.push_back(phogpt);	   listvar.push_back(phogpt);
       	 if(phovec){	     
       	   listvar.push_back(phogpt-phovec->Pt());	   
       	   listvar.push_back(phogdR);
       	 }
       	 else{
       	   listvar.push_back(-999.);	   
       	   listvar.push_back(-999.);
       	 }
       }//if(SaveGenLevel && phosel) 
	 

       if(phovec){
       	 listvar.push_back(phovec->Pt());	 listvar.push_back(phovec->Pt());
       	 listvar.push_back(phovec->Eta());
       }
       else{
       	 listvar.push_back(-999.);	 listvar.push_back(-999.);
       	 listvar.push_back(-999.);
       }
       listvar.push_back(pfmetnopho);       listvar.push_back(pfmetnopho);


       listvar.push_back(isHadrTopHighestMVA);


       if(listvar.size() != listhisto.size()){
	 cerr<<"listvar and listhisto should have the same size. Instead they have "<<listvar.size()<<", "<<listhisto.size()<<endl;
	 cout<<"listhisto"<<endl;
	 for(int il=0; il<listhisto.size(); il++)
	   cout<<""<<listhisto[il].first<<endl;
	 exit(1);
       }
       //*******************

       //*******VARIABLE 2D to be plotted
       listvar2D.clear();
       listvar2D.push_back(make_pair(wgpt,pfmet));
       listvar2D.push_back(make_pair(wgpt,pfmetnolep));

       //*****

       if(listvar2D.size() != listhisto2D.size()){
	 cerr<<"listvar2D and listhisto2D should have the same size. Instead they have "<<listvar2D.size()<<", "<<listhisto2D.size()<<endl;
	 for(int il2D=0; il2D<listhisto2D.size(); il2D++)
	   cout<<""<<listhisto2D[il2D].first<<endl;
	 exit(1);
       }




       //*******************
       //*******************





       if(verbose==2)
	 cout<<"before making hv"<<endl;


       neventsv[isam]+=weight;
       nentriessamplev_tmp[ifile]++;
       neventssamplev_tmp[ifile]+=weight;
       for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) { 
	 if(verbose>=1)
	   cout<<"ilist "<<ilist<<" hv "<<hv[ilist][isam]->GetName()<<endl;	 
	 //	 cout<<"ilist "<<ilist<<"listvar "<<listvar[ilist].first<<" "<<listvar[ilist].second<<endl;
	 //	   hv[ilist][isam] ->Fill(TMath::Min((double)listvar[ilist].first,hv[ilist][isam]->GetXaxis()->GetXmax()-1.0), weight); //need to have the variable as a list as wll
	   hv[ilist][isam] ->Fill((double)listvar[ilist], weight); 
       }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++)
       for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++) { 
	 if(verbose>=1)
	   cout<<"ilist2D "<<ilist2D<<" hv "<<hv2D[ilist2D][isam]->GetName()<<endl;	 
	 hv2D[ilist2D][isam] ->Fill((double)listvar2D[ilist2D].first,(double)listvar2D[ilist2D].second, weight); 
       }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++)



       if(verbose==2)
	 cout<<"before making hMC"<<endl;

       if(!isData && !isSignal) {
	 //cout<<"qui13"<<endl;
	 neventsMC+=weight;
	 for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {
	   if(verbose>=1)
	     cout<<"ilist "<<ilist<<" hMC "<<hMC[ilist]->GetName()<<endl;
	   //	     hMC[ilist]   ->Fill(TMath::Min((double)listvar[ilist].first,hMC[ilist]->GetXaxis()->GetXmax()-1.0), weight);
	   hMC[ilist]   ->Fill((double)listvar[ilist], weight);
	 }//       for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {
	 for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++) {
	   if(verbose>=1)
	     cout<<"ilist2D "<<ilist2D<<" hMC "<<hMC2D[ilist2D]->GetName()<<endl;
	   hMC2D[ilist2D]   ->Fill((double)listvar2D[ilist2D].first,(double)listvar2D[ilist2D].second, weight);
	 }//	 for(unsigned int ilist2D=0; ilist<listhisto2D.size(); ilist2D++) {
       }
       //	  cout<<"qui21"<<endl;
     }//for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++)

      //cout<<"qui15"<<endl;
      
     delete infile;
     infile=0;
     intree=0;
   }//for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++)
   nentriessamplev.push_back(nentriessamplev_tmp);        
   neventssamplev.push_back(neventssamplev_tmp);        
 }//for(unsigned int isam=0; isam<samplev.size(); isam++)

 if(verbose>=1)
   cout<<"end filling"<<endl;

 for(unsigned int isam=0; isam<samplev.size(); isam++){
   CSample *sample = samplev[isam];
   if(verbose>=1)
     cout << "Sample: " << sample->label << endl;
   for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++){
     if(verbose>=1)
       cout<<"isam, ifile, nentries, nevents "<<isam<<", "<<ifile<<", "<<nentriessamplev[isam][ifile]<<", "<<neventssamplev[isam][ifile]<<endl;
   }
   if(verbose>=1)
     cout<<"nevents "<<neventsv[isam]<<endl;
 }


 //
 // Make ratio histograms
 //
 for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) { 
   hname = listhisto[ilist].first; hbinning = listhisto[ilist].second;
   hname += "Ratio";
   if(verbose>=1)
     cout<<"hname "<<hname<<endl;
   TH1D *htmp = makeRatioHist(hv[ilist][0], hMC[ilist], hname, verbose);
   hRatio.push_back(htmp);   

 }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++)





 //--------------------------------------------------------------------------------------------------------------
 // Make plots
 //==============================================================================================================

 TCanvas *c = MakeCanvas("c","c",800,800);
 c->Divide(1,2,0,0);
 c->cd(1)->SetPad(0,0.3,1.0,1.0);
 c->cd(1)->SetTopMargin(0.1);
 c->cd(1)->SetBottomMargin(0.01);
 c->cd(1)->SetLeftMargin(0.15);
 c->cd(1)->SetRightMargin(0.07);
 c->cd(1)->SetTickx(1);
 c->cd(1)->SetTicky(1);
 c->cd(2)->SetPad(0,0,1.0,0.3);
 c->cd(2)->SetTopMargin(0.05);
 c->cd(2)->SetBottomMargin(0.45);
 c->cd(2)->SetLeftMargin(0.15);
 c->cd(2)->SetRightMargin(0.07);
 c->cd(2)->SetTickx(1);
 c->cd(2)->SetTicky(1);

 char pname[100];
 char ylabel[100];
 string suffix;
 
 if(verbose==2) 
   cout<<"makeplot"<<endl;

 for(unsigned int ilist=0; ilist<listhisto.size(); ilist++) {

   sprintf(pname,(listhisto[ilist].first).c_str());
   sprintf(ylabel,"Events / %i GeV",int(hv[ilist][0]->GetBinWidth(1))); //cambiare string ylabel = "Events /"+hv[ilist][0]->GetBinWidth(1) + unity[ilist] (use ostringstream, and unity in DefineListSample.C)
   makePlot(c, pname, settingslabel[ilist], ylabel, hv[ilist], samplev, hMC[ilist], hRatio[ilist], outputDir+"/out.root", LUMI, settingslogscale[ilist], 0.05, -0.03,
	    (settingsrange[ilist].first)*(hMC[ilist]->GetBinContent(hMC[ilist]->GetMaximumBin())), (settingsrange[ilist].second)*(hMC[ilist]->GetBinContent(hMC[ilist]->GetMaximumBin())), verbose);
 }//for(unsigned int ilist=0; ilist<listhisto.size(); ilist++)


 for(unsigned int ilist2D=0; ilist2D<listhisto2D.size(); ilist2D++){
   for(unsigned int isam=0; isam<samplev.size(); isam++){
    sprintf(pname,(listhisto2D[ilist2D].first).c_str());
    makePlot2D(c, pname, settingslabelX[ilist2D], settingslabelY[ilist2D], hv2D[ilist2D][isam], samplev[isam], hMC2D[ilist2D], outputDir+"/out2D.root", LUMI, verbose);   
  }
 }//	 for(unsigned int ilist2d=0; ilist2D<listhisto2D.size(); ilist2D++) 
 


 //--------------------------------------------------------------------------------------------------------------
 // Output
 //==============================================================================================================
 std::ofstream yieldfile;
 string yieldstring = outputDir+"/"+"Yields.txt";
 yieldfile.open (yieldstring.c_str(), std::ofstream::out | std::ofstream::app);

 //SCRIVERE QUA I TAGLI NELLE YIELDS.txt
 /*  cout<<"******Cuts (on top of preselection applied in runAllHadronicPreselection****"<<endl;
     cout<<"metfilter"<<endl;
     cout<<""
     cout<<"trigger "<<cutflow_cut<<endl;
     cout<<"MET_CUT "<<MET_CUT<<endl;
     cout<<"njtes05 "<<njets_cut<<endl;*/
   
 cout << "*" << endl;
 cout << "* SUMMARY" << endl;
 cout << "*--------------------------------------------------" << endl;
 cout << endl;
 yieldfile << "Yields" << endl;
 yieldfile << "--------------------------------------------------" << endl;
 yieldfile << endl;


 //compute errors
 cout<<"Computing errors "<<nentriessamplev.size()<<endl;
 if(nentriessamplev.size() > 0)
   if(nentriessamplev[0].size() > 0)
     nerrv.push_back(TMath::Sqrt(nentriessamplev[0][0])); //data has one only file
 nerrMC = 0;
 for(unsigned int isam=1; isam<samplev.size(); isam++) {
   nerrv.push_back(0);
   CSample *sample = samplev[isam];
   if(verbose==1)
     cout << "Sample: " << sample->label << endl;
   for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++){
     if(verbose==1)
       cout<<"isam, ifile, nentries, nevents "<<isam<<", "<<ifile<<", "<<nentriessamplev[isam][ifile]<<", "<<neventssamplev[isam][ifile]<<endl;
     if(nentriessamplev[isam][ifile]>0)
       nerrv[isam] += neventssamplev[isam][ifile]*neventssamplev[isam][ifile]/nentriessamplev[isam][ifile];
   }
   nerrMC += nerrv[isam];
   nerrv[isam] = TMath::Sqrt(nerrv[isam]);
   if(verbose==1)
     cout<<"isam, err "<<isam<<", "<<nerrv[isam]<<endl;
 }

 nerrMC = TMath::Sqrt(nerrMC);

 for(unsigned int isam=1; isam<samplev.size()-1; isam++) {
   cout << "  * " << setw(10) << samplev[isam]->label;
   cout << setprecision(2) << fixed << setw(10) << neventsv[isam]<<" +/- "<<nerrv[isam];    
   cout << endl;
   yieldfile << "  * " << setw(10) << samplev[isam]->label;
   yieldfile << setprecision(2) << fixed << setw(10) << neventsv[isam]<<" +/- "<<nerrv[isam];    
   yieldfile << endl;
 }
 cout << "  ============" << endl;
 cout << "    " << setw(10) << "total bkg";
 cout << setprecision(2) << fixed << setw(10) << neventsMC<<" +/- "<<nerrMC;

 cout << endl;
 cout << "    " << setw(10) << "observed";
 cout << setprecision(2) << fixed << setw(10) << neventsv[0];

 cout << endl;
 cout << "    " << setw(10) << samplev.back()->label;
 cout << setprecision(2) << fixed << setw(10) << neventsv[samplev.size()-1]<<" +/- "<<nerrv[samplev.size()-1];
  
 cout << endl;

 cout << endl;
 cout << " <> Output saved in " << outputDir << "/" << endl;
 cout << endl; 

 yieldfile << "  ============" << endl;
 yieldfile << "    " << setw(10) << "total bkg";
 yieldfile << setprecision(2) << fixed << setw(10) << neventsMC<<" +/- "<<nerrMC;

 yieldfile << endl;
 yieldfile << "    " << setw(10) << "observed";
 yieldfile << setprecision(2) << fixed << setw(10) << neventsv[0];

 yieldfile << endl;
 yieldfile << "    " << setw(10) << samplev.back()->label;
 yieldfile << setprecision(2) << fixed << setw(10) << neventsv[samplev.size()-1]<<" +/- "<<nerrv[samplev.size()-1];
  
 yieldfile << endl;

 // yieldfile << endl;
 // yieldfile << " <> Output saved in " << outputDir << "/" << endl;
 // yieldfile << endl; 


 yieldfile.close();

 for(unsigned int isam=1; isam<samplev.size(); isam++){
   nentriessamplev[isam].clear();
   neventssamplev[isam].clear();
 }
 nentriessamplev.clear();
 neventssamplev.clear();
}

