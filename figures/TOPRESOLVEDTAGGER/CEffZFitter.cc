#include "CEffZFitter.hh"
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <cassert>
#include <sstream>
#include <iomanip>

#include "CPlot.hh"
#include "KStyle.hh"
#include "CEffUser1D.hh"


// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"


#include "RooAbsPdf.h"
#include "RooHistPdf.h"




// bin size constants
#define BIN_SIZE_PASS 20
#define BIN_SIZE_FAIL 20
#define LUMI 19.7 //MT



//--------------------------------------------------------------------------------------------------
CEffZFitter::CEffZFitter():
fIsInitialized(false),
fSigPass      (0),
fBkgPass      (0),
fSigFail      (0),
fBkgFail      (0),
fMassLo       (60),
fMassHi       (120),
fFitMassLo    (60),
fFitMassHi    (120),
fOutputDir    ("."),
topmvacut     (-999.),
maketemplates     (true)

{}

//--------------------------------------------------------------------------------------------------
CEffZFitter::~CEffZFitter()
{
  delete fPassTree;     fPassTree=0;  
  
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::initialize(const std::string infname, const std::string outdir, const std::string temfname1, const std::string temfname2, const std::string temfname3,
                             const double massLo, const double massHi, const double fitMassLo, const double fitMassHi, 
		             const unsigned int runNumLo, const unsigned int runNumHi,
			     const double TOPMVACUT, const bool MAKETEMPLATES)
{
  std::cout << "   [CEffZFitter] Initializing... " << std::endl;



  fMassLo    = massLo;
  fMassHi    = massHi;
  fFitMassLo = fitMassLo;
  fFitMassHi = fitMassHi;
  topmvacut = TOPMVACUT;
  maketemplates = MAKETEMPLATES;
  
  // set up output directory
  fOutputDir = outdir;
  gSystem->mkdir(fOutputDir.c_str(),true);
  //CPlot::sOutDir = TString(outdir.c_str()) + TString("/plots");
      
  /*TFile *pufile=0;
    TH1D *puWeights=0;
  if(pufname.compare("none")!=0) {
    pufile = new TFile(pufname.c_str());      assert(pufile);
    puWeights = (TH1D*)pufile->Get("pileup"); assert(puWeights); 
    }*/

  if(maketemplates){
    makeBinnedTemplates(temfname1, "signal");
    makeBinnedTemplates(temfname2, "bkg1");
    makeBinnedTemplates(temfname3, "bkg2");
  }

  
  //------------------------------------------------------------------------------------------------
  // Read in probes data
  //===============================================================================================

  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv, npu;                       // number of primary vertices  
  
  float        scale1fb, evtWeight;                  // event weight per 1/fb and other non xsec-lumi weights
  TLorentzVector *res_jet1=0,     *res_jet2=0,     *res_jet3=0;
  unsigned int nele, nmuo, nbjets, njets05;
  unsigned int  metfilter;
  float pfmet;
  bool passtrigger;

  float res_topmva;



  float topmass; //will be computed on the fly


  
  TFile *infile = new TFile(infname.c_str());    assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);

  //from plotAllHadronic.C: keeping the minimum amount of info for selection and fit
  intree->SetBranchAddress("runNum",    &runNum);
  intree->SetBranchAddress("lumiSec",   &lumiSec);
  intree->SetBranchAddress("evtNum",    &evtNum);
  intree->SetBranchAddress("metfilter", &metfilter);
  intree->SetBranchAddress("npv",       &npv);
  intree->SetBranchAddress("npu",       &npu);
  intree->SetBranchAddress("njets05",   &njets05);
  intree->SetBranchAddress("nbjets",    &nbjets);
  intree->SetBranchAddress("scale1fb",  &scale1fb);
  intree->SetBranchAddress("evtWeight", &evtWeight);
  intree->SetBranchAddress("pfmet",     &pfmet);
  intree->SetBranchAddress("nele",   &nele);
  intree->SetBranchAddress("nmuo",   &nmuo);
  intree->SetBranchAddress("passtrigger",       &passtrigger);
  intree->SetBranchAddress("res_topmva",   &res_topmva);
  intree->SetBranchAddress("res_jet1", &res_jet1);
  intree->SetBranchAddress("res_jet2", &res_jet2);
  intree->SetBranchAddress("res_jet3", &res_jet3);



  unsigned int pass;
  
 

  
  char tname[50];
  float wgt;

  {
    sprintf(tname,"pass");
    fPassTree = new TTree(tname,"");
    fPassTree->Branch("m",&topmass,"m/F");
    fPassTree->Branch("w",&wgt, "w/F");
    fPassTree->SetDirectory(0);
    sprintf(tname,"fail");
    fFailTree = new TTree(tname,"");
    fFailTree->Branch("m",&topmass,"m/F");
    fFailTree->Branch("w",&wgt, "w/F");
    fFailTree->SetDirectory(0);
  }


  //
  // loop over probes
  //
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
    
    //    if(qprobe*charge < 0) continue; //MT

     if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 !=0 && njets05>=3)
       topmass = calcmass(res_jet1,res_jet2,res_jet3); //Calctopmass(njets05, res_jet1,res_jet2,res_jet3,res_jet4,res_jet5,res_jet6);
     else
       topmass =  -999.;  


    //PUTHERE SELECTION
    if(topmass < fFitMassLo) continue;
    if(topmass > fFitMassHi) continue;
    if(runNum < runNumLo) continue;
    if(runNum > runNumHi) continue;

    if(res_topmva > topmvacut) pass = 1;
    else pass = 0;



    wgt = scale1fb;
    //if(puWeights) {
    // wgt *= LUMI*evtWeight; //MT why Kevin doesn't include LUMI
      //wgt *= puWeights->GetBinContent(puWeights->FindBin(npu));
      // }
    

    //
    // Fill trees
    //
    if(pass) {
      fPassTree->Fill();
    
    } else {
      fFailTree->Fill();
      
    }
  }
  delete infile;
  infile=0, intree=0;
  
  //  delete pufile;
  //pufile=0, puWeights=0;
  
  fIsInitialized = true;
}

//--------------------------------------------------------------------------------------------------
void CEffZFitter::computeEff()
{
  assert(fIsInitialized);
  
  std::cout << "   [CEffZFitter] Computing efficiencies..." << std::endl;


  //------------------------------------------------------------------------------------------------
  // Efficiency calculation
  //================================================================================================

  
  double eff, errl, errh;  
  performFit(eff, errl, errh,
	     fPassTree, fFailTree); 


  std::cout<<"eff, errl, errh"<<errh<<" "<<errl<<" "<<errh<<std::endl;  
  
}


//--------------------------------------------------------------------------------------------------
void CEffZFitter::makeBinnedTemplates(const std::string temfname, string name)
{
  std::cout << "   [CEffZFitter] Creating binned templates... "; std::cout.flush();

  char hname[50];
  
  
  TH1D* h_pass;
  TH1D* h_fail;

  sprintf(hname,"h_pass_%s",name.c_str());
  h_pass = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi);
  h_pass->SetDirectory(0);
  sprintf(hname,"h_fail_%s",name.c_str());
  h_fail = new TH1D(hname,"",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
  h_fail->SetDirectory(0);

  
    
  
  unsigned int runNum, lumiSec, evtNum;   // event ID
  unsigned int npv, npu;                       // number of primary vertices  

  float        scale1fb, evtWeight;                  // event weight per 1/fb and other non xsec-lumi weights
  TLorentzVector *res_jet1=0,     *res_jet2=0,     *res_jet3=0;
  unsigned int nele, nmuo, nbjets, njets05;
  unsigned int  metfilter;
  float pfmet;
  bool passtrigger;

  float res_topmva;
  
  float topmass; //will be computed on the fly



  TFile *infile = new TFile(temfname.c_str());   assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);


  //from plotAllHadronic.C: keeping the minimum amount of info for selection and fit
  intree->SetBranchAddress("runNum",    &runNum);
  intree->SetBranchAddress("lumiSec",   &lumiSec);
  intree->SetBranchAddress("evtNum",    &evtNum);
  intree->SetBranchAddress("metfilter", &metfilter);
  intree->SetBranchAddress("npv",       &npv);
  intree->SetBranchAddress("npu",       &npu);
  intree->SetBranchAddress("njets05",   &njets05);
  intree->SetBranchAddress("nbjets",    &nbjets);
  intree->SetBranchAddress("scale1fb",  &scale1fb);
  intree->SetBranchAddress("evtWeight", &evtWeight);
  intree->SetBranchAddress("pfmet",     &pfmet);
  intree->SetBranchAddress("nele",   &nele);
  intree->SetBranchAddress("nmuo",   &nmuo);
  intree->SetBranchAddress("passtrigger",       &passtrigger);
  intree->SetBranchAddress("res_topmva",   &res_topmva);
  intree->SetBranchAddress("res_jet1", &res_jet1);
  intree->SetBranchAddress("res_jet2", &res_jet2);
  intree->SetBranchAddress("res_jet3", &res_jet3);

  unsigned int pass;
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);


    if(res_jet1 != 0 && res_jet2 !=0 && res_jet3 !=0 && njets05>=3)
      topmass = calcmass(res_jet1,res_jet2,res_jet3); //Calctopmass(njets05, res_jet1,res_jet2,res_jet3,res_jet4,res_jet5,res_jet6);
     else
       topmass =  -999.;  
    
    
    double weight = scale1fb;
    weight *= LUMI*evtWeight; //MT why Kevin doesn't include LUMI
    /*if(puWeights)
      weight *= puWeights->GetBinContent(puWeights->FindBin(npu));*/
    
    
    //    if(qprobe*charge < 0) continue; //MT
    //same selection as in initialize
    if(topmass < fFitMassLo) continue;
    if(topmass > fFitMassHi) continue;

    if(res_topmva > topmvacut) pass = 1;
    else pass = 0;
    


    if(pass) {
      h_pass->Fill(topmass,weight);
    } else {
      h_fail->Fill(topmass,weight);
    }    
  }
  infile->Close();
 
  string str_tmp = fOutputDir+"/binnedTemplates_"+name+".root";
  TFile outfile(str_tmp.c_str(), "RECREATE");

  h_pass->Write();
  h_fail->Write();
  delete h_pass;
  delete h_fail;

  outfile.Write();
  outfile.Close(); 

  cout << "Done!" << endl;
}




//--------------------------------------------------------------------------------------------------
void CEffZFitter::performFit(double &resEff, double &resErrl, double &resErrh,
                             TTree *passTree, TTree *failTree)

{

  std::cout << " ...performFit..." << std::endl; 

  RooRealVar m("m","Top Mass [GeV/c^2]",fFitMassLo,fFitMassHi);
  m.setBins(10000);
  
  
  TFile *histfile_signal = 0;
 TFile *histfile_bkg1 = 0;  TFile *histfile_bkg2 = 0;
  {
    string str_tmp = fOutputDir+"/binnedTemplates_signal.root";
    histfile_signal = new TFile(str_tmp.c_str());
    assert(histfile_signal);

   
    str_tmp = fOutputDir+"/binnedTemplates_bkg1.root";
    histfile_bkg1 = new TFile(str_tmp.c_str());  
    str_tmp = fOutputDir+"/binnedTemplates_bkg2.root";
    histfile_bkg2 = new TFile(str_tmp.c_str());
    assert(histfile_bkg1);  assert(histfile_bkg2);
  }
  
  // Define categories
  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);
  
  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  TH1D histPass("histPass","",int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS,fFitMassLo,fFitMassHi); 
  TH1D histFail("histFail","",int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL,fFitMassLo,fFitMassHi);
  RooAbsData *dataCombined=0;
  
  std::cout<<"passtree entries "<<passTree->Draw("m>>histPass","w")<<std::endl;
  std::cout<<"failtree entries "<<failTree->Draw("m>>histFail","w")<<std::endl;
  dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),&histPass);
  dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),&histFail);
  //m.setBins(100);  

  dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                   RooFit::Index(sample),
                                   RooFit::Import("Pass",*((RooDataHist*)dataPass)),
                                   RooFit::Import("Fail",*((RooDataHist*)dataFail)));



  std::cout<<"qui"<<std::endl;
  RooHistPdf *sigModPass;
  RooHistPdf *bkg1ModPass;
  RooHistPdf *bkg2ModPass;
  RooHistPdf *sigModFail;
  RooHistPdf *bkg1ModFail;
  RooHistPdf *bkg2ModFail;
  RooDataHist *sigModPass_rdh; 
  RooDataHist *bkg1ModPass_rdh;
  RooDataHist *bkg2ModPass_rdh;  
  RooDataHist *sigModFail_rdh;
  RooDataHist *bkg1ModFail_rdh;
  RooDataHist *bkg2ModFail_rdh;   

  TH1D *h_pass_signal, *h_pass_bkg1, *h_pass_bkg2; 
  TH1D *h_fail_signal, *h_fail_bkg1, *h_fail_bkg2; 
  {
    char hname_signal[50];
    sprintf(hname_signal,"h_pass_%s","signal");
    h_pass_signal = (TH1D*)histfile_signal->Get(hname_signal);
    cout<<"looking for histo "<<hname_signal<<" in file "<<histfile_signal->GetName()<<std::endl;
    assert(h_pass_signal);
    char hname_bkg1[50];
    sprintf(hname_bkg1,"h_pass_%s","bkg1");
    h_pass_bkg1 = (TH1D*)histfile_bkg1->Get(hname_bkg1);
    assert(h_pass_bkg1);
    char hname_bkg2[50];
    sprintf(hname_bkg2,"h_pass_%s","bkg2");
    h_pass_bkg2 = (TH1D*)histfile_bkg2->Get(hname_bkg2);
    assert(h_pass_bkg2);
    





    sigModPass_rdh= new RooDataHist("sigModPass_rdh","sigModPass_rdh",RooArgSet(m),h_pass_signal);
    bkg1ModPass_rdh = new RooDataHist("bkg1ModPass_rdh","bkg1ModPass_rdh",RooArgSet(m),h_pass_bkg1);
    bkg2ModPass_rdh = new RooDataHist("bkg2ModPass_rdh","bkg2ModPass_rdh",RooArgSet(m),h_pass_bkg2);

    sigModPass = new RooHistPdf("sigModPass","sigModPass",RooArgSet(m),*sigModPass_rdh);
    bkg1ModPass = new RooHistPdf("bkg1ModPass","bkg1ModPass",RooArgSet(m),*bkg1ModPass_rdh);
    bkg2ModPass = new RooHistPdf("bkg2ModPass","bkg2ModPass",RooArgSet(m),*bkg2ModPass_rdh);



                                                                                                                    
    // Make PDF from MC histograms                                                                                                                                                                                                                
  }  

  std::cout<<"qui3"<<std::endl;

  {
    char hname_signal[50];
    sprintf(hname_signal,"h_fail_%s","signal");
    h_fail_signal = (TH1D*)histfile_signal->Get(hname_signal);
    assert(h_fail_signal);
    char hname_bkg1[50];
    sprintf(hname_bkg1,"h_fail_%s","bkg1");
    h_fail_bkg1 = (TH1D*)histfile_bkg1->Get(hname_bkg1);
    assert(h_fail_bkg1);
    char hname_bkg2[50];
    sprintf(hname_bkg2,"h_fail_%s","bkg2");
    h_fail_bkg2 = (TH1D*)histfile_bkg2->Get(hname_bkg2);
    assert(h_fail_bkg2);





                                                                                                                    
    // Make PDF from MC histograms                                                                                                                                                                                                          

      
    sigModFail_rdh = new RooDataHist("sigModFail_rdh","sigModFail_rdh",RooArgSet(m),h_fail_signal);
    bkg1ModFail_rdh = new RooDataHist("bkg1ModFail_rdh","bkg1ModFail_rdh",RooArgSet(m),h_fail_bkg1);
    bkg2ModFail_rdh = new RooDataHist("bkg2ModFail_rdh","bkg2ModFail_rdh",RooArgSet(m),h_fail_bkg2);

    sigModFail = new RooHistPdf("sigModFail","sigModFail",RooArgSet(m),*sigModFail_rdh);
    bkg1ModFail = new RooHistPdf("bkg1ModFail","bkg1ModFail",RooArgSet(m),*bkg1ModFail_rdh);
    bkg2ModFail = new RooHistPdf("bkg2ModFail","bkg2ModFail",RooArgSet(m),*bkg2ModFail_rdh);





  }


    std::cout<<"**********Fail**********"<<std::endl;
    std::cout<<"h_signal, h_bkg1, h_bkg2 "<<h_fail_signal->Integral()<<" "<<h_fail_bkg1->Integral()<<" "<<h_fail_bkg2->Integral()<<std::endl;
    std::cout<<"**********Pass********"<<std::endl;
    std::cout<<"h_signal, h_bkg1, h_bkg2 "<<h_pass_signal->Integral()<<" "<<h_pass_bkg1->Integral()<<" "<<h_pass_bkg2->Integral()<<std::endl;
    std::cout<<"*********TOT**********"<<std::endl;
    std::cout<<"h_signal, h_bkg1, h_bkg2 "<<h_fail_signal->Integral()+h_pass_signal->Integral()<<" "<<h_fail_bkg1->Integral()+h_pass_bkg1->Integral()<<" "<<h_fail_bkg2->Integral()+h_pass_bkg2->Integral()<<std::endl;
    std::cout<<"*********TOTMC**********"<<std::endl;
    std::cout<<"Pass, Fail "<<h_pass_signal->Integral()+h_pass_bkg1->Integral()+h_pass_bkg2->Integral()<<" "<<h_fail_signal->Integral()+h_fail_bkg1->Integral()+h_fail_bkg2->Integral()<<std::endl;
    std::cout<<"*******Data*******"<<std::endl;
    std::cout<<"Pass, Fail, Tot "<<histPass.Integral()<<" "<<histFail.Integral()<<" "<<histPass.Integral()+histFail.Integral()<<std::endl;
    std::cout<<"**************"<<std::endl;



  

  // Define free parameters
  double NsigMax     = histPass.Integral()+histFail.Integral();
  //  double Nbkg1FailMax = histFail.Integral();
  double Nbkg1Max = histPass.Integral()+histFail.Integral();
  double Nbkg2Max = histPass.Integral()+histFail.Integral();
  //  double Nbkg1PassMax = histPass.Integral();
  //double Nbkg2FailMax = histFail.Integral();
  //double Nbkg2PassMax = histPass.Integral();

  std::cout<<"NsigMax, Nbkg1Max, Nbkg2Max "<<NsigMax<<" "<<Nbkg1Max<<" "<<" "<<Nbkg2Max<<std::endl;


  RooRealVar Nsig("Nsig","Signal Yield",0.0,0,NsigMax);
  RooRealVar eff("eff","Efficiency",0.0,0,1.0);
  RooRealVar Nbkg1("Nbkg1","Bkg1 Yield",0.0,0,Nbkg1Max);
  RooRealVar FR1("FR1","FR Bkg1",0.0,0,1.0);
  RooRealVar Nbkg2("Nbkg2","Bkg2 Yield",0.0,0,Nbkg2Max);
  RooRealVar FR2("FR2","FR Bkg2",0.0,0,1.0);

  RooFormulaVar Nbkg1Pass("Nbkg1Pass","FR1*Nbkg1",RooArgList(FR1,Nbkg1));
  RooFormulaVar Nbkg1Fail("Nbkg1Fail","(1.0-FR1)*Nbkg1",RooArgList(FR1,Nbkg1));

  RooFormulaVar Nbkg2Pass("Nbkg2Pass","FR2*Nbkg2",RooArgList(FR2,Nbkg2));
  RooFormulaVar Nbkg2Fail("Nbkg2Fail","(1.0-FR2)*Nbkg2",RooArgList(FR2,Nbkg2));
    
  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));   
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));
 
  RooAddPdf *modelPass=0, *modelFail=0;


  modelPass = new RooAddPdf("modelPass","Model for PASS sample",RooArgList(*sigModPass,*bkg1ModPass,*bkg2ModPass),RooArgList(NsigPass,Nbkg1Pass,Nbkg2Pass));
  
  modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*sigModFail,*bkg1ModFail,*bkg2ModFail),RooArgList(NsigFail,Nbkg1Fail,Nbkg2Fail)); 
  
  
  
  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass"); 
  totalPdf.addPdf(*modelFail,"Fail");




  // Plot the imported histogram(s)                                                                                   
  TH1D *h_pass_tot = (TH1D *)h_pass_signal->Clone("h_pass_tot");
  h_pass_tot->Add(h_pass_signal);
  h_pass_tot->Add(h_pass_bkg1);
  h_pass_tot->Add(h_pass_bkg2);
  h_pass_tot->SetLineColor(kBlue);
  h_pass_signal->SetLineColor(kGreen);
  h_pass_bkg1->SetLineColor(kRed);
  h_pass_bkg2->SetLineColor(kGray);
  h_pass_signal->SetLineStyle(kDashed);
  h_pass_bkg1->SetLineStyle(kDashed);
  h_pass_bkg2->SetLineStyle(kDashed);

  TH1D *h_fail_tot = (TH1D *)h_fail_signal->Clone("h_fail_tot");
  h_fail_tot->Add(h_fail_signal);
  h_fail_tot->Add(h_fail_bkg1);
  h_fail_tot->Add(h_fail_bkg2);
  h_fail_tot->SetLineColor(kBlue);
  h_fail_signal->SetLineColor(kGreen);
  h_fail_bkg1->SetLineColor(kRed);
  h_fail_bkg2->SetLineColor(kGray);
  h_fail_signal->SetLineStyle(kDashed);
  h_fail_bkg1->SetLineStyle(kDashed);
  h_fail_bkg2->SetLineStyle(kDashed);


  TLegend *legend = new TLegend(0.65,0.54,0.99,0.99);
  legend->SetFillColor(0);
  legend->AddEntry(&histPass,"Data" , "lp");
  legend->AddEntry(h_pass_tot, "TOT MC" , "f");
  legend->AddEntry(h_pass_signal, "Signal" , "f");
  legend->AddEntry(h_pass_bkg1, "tt1l Combinatorial BG" , "f");
  legend->AddEntry(h_pass_bkg2, "Non-tt1l BG" , "f");
                                                                                                                        {
    std::cout<<"plotting data and model prefit via TH1"<<std::endl;

    TCanvas* cpass_bf_mt = new TCanvas("cpass_bf_mt","cpass_bf_mt",800,1200);
    cpass_bf_mt->SetWindowPosition(cpass_bf_mt->GetWindowTopX()+cpass_bf_mt->GetBorderSize()+800,0);
    cpass_bf_mt->SetTickx(1);
    cpass_bf_mt->SetTicky(1);
    if(histPass.GetMaximum() < h_pass_tot->GetMaximum())
      histPass.SetMaximum(h_pass_tot->GetMaximum()*1.1);
    histPass.Draw("E");
    h_pass_tot->Draw("histsame");
    h_pass_signal->Draw("samehist");
    h_pass_bkg1->Draw("samehist");
    h_pass_bkg2->Draw("samehist");

    legend->Draw();

    string str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_mt.png"; 
    cpass_bf_mt->SaveAs(str_tmp.c_str());
    str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit_mt.C"; 
    cpass_bf_mt->SaveAs(str_tmp.c_str());

    TCanvas* cfail_bf_mt = new TCanvas("cfail_bf_mt","cfail_bf_mt",800,1200);
    cfail_bf_mt->SetWindowPosition(cfail_bf_mt->GetWindowTopX()+cfail_bf_mt->GetBorderSize()+800,0);
    cfail_bf_mt->SetTickx(1);
    cfail_bf_mt->SetTicky(1);
    if(histFail.GetMaximum() < h_fail_tot->GetMaximum())
      histFail.SetMaximum(h_fail_tot->GetMaximum()*1.1);

    histFail.Draw("E");
    h_fail_tot->Draw("histsame");
    h_fail_signal->Draw("samehist");
    h_fail_bkg1->Draw("samehist");
    h_fail_bkg2->Draw("samehist");

    legend->Draw();


    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_mt.png"; 
    cfail_bf_mt->SaveAs(str_tmp.c_str());
    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit_mt.C"; 
    cfail_bf_mt->SaveAs(str_tmp.c_str());


  }
 


  {
    std::cout<<"plotting data and model prefit"<<std::endl;
    RooPlot* mframePass_bf = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS));
    mframePass_bf->SetTitle("");
    dataPass->plotOn(mframePass_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
    modelPass->plotOn(mframePass_bf);
    modelPass->plotOn(mframePass_bf,Components("sigModPass"),LineStyle(kDashed),LineColor(kGreen));
    modelPass->plotOn(mframePass_bf,Components("bkg1ModPass"),LineStyle(kDashed),LineColor(kRed));
    modelPass->plotOn(mframePass_bf,Components("bkg2ModPass"),LineStyle(kDashed),LineColor(kGray));
    TCanvas* cpass_bf = new TCanvas("cpass_bf","cpass_bf",800,1200);
    cpass_bf->SetWindowPosition(cpass_bf->GetWindowTopX()+cpass_bf->GetBorderSize()+800,0);
    cpass_bf->SetTickx(1);
    cpass_bf->SetTicky(1);
    mframePass_bf->Draw();

    legend->Draw();

    string str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit.png"; 
    cpass_bf->SaveAs(str_tmp.c_str());
    str_tmp = fOutputDir+"/TemplatesandData_Passsample_beforefit.C"; 
    cpass_bf->SaveAs(str_tmp.c_str());



    RooPlot* mframeFail_bf = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS));
    mframeFail_bf->SetTitle("");
    dataFail->plotOn(mframeFail_bf,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
    modelFail->plotOn(mframeFail_bf);
    modelFail->plotOn(mframeFail_bf,Components("sigModFail"),LineStyle(kDashed),LineColor(kGreen));
    modelFail->plotOn(mframeFail_bf,Components("bkg1ModFail"),LineStyle(kDashed),LineColor(kRed));
    modelFail->plotOn(mframeFail_bf,Components("bkg2ModFail"),LineStyle(kDashed),LineColor(kGray));
    TCanvas* cfail_bf = new TCanvas("cfail_bf","cfail_bf",800,1200);
    cfail_bf->SetWindowPosition(cfail_bf->GetWindowTopX()+cfail_bf->GetBorderSize()+800,0);
    cfail_bf->SetTickx(1);
    cfail_bf->SetTicky(1);
    mframeFail_bf->Draw();

    legend->Draw();

    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit.png"; 
    cfail_bf->SaveAs(str_tmp.c_str());
    str_tmp = fOutputDir+"/TemplatesandData_Failsample_beforefit.C"; 
    cfail_bf->SaveAs(str_tmp.c_str());
  }






  RooFitResult *fitResult=0;
  std::cout<<"before fit"<<std::endl;
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::Extended(),
                             //RooFit::Strategy(2),
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::NumCPU(4),
                             RooFit::Save());
 
  std::cout<<"after fit"<<std::endl;

  // Refit w/o MINOS if MINOS errors are strange...
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  resEff  = eff.getVal(); 
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();






  //plot

  RooPlot *mframePass = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_PASS));
  mframePass->SetTitle("");
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));    
  modelPass->plotOn(mframePass);
  modelPass->plotOn(mframePass,Components("sigModPass"),LineStyle(kDashed),LineColor(kGreen));
  modelPass->plotOn(mframePass,Components("bkg1ModPass"),LineStyle(kDashed),LineColor(kRed));
  modelPass->plotOn(mframePass,Components("bkg2ModPass"),LineStyle(kDashed),LineColor(kGray));

  
  RooPlot *mframeFail = m.frame(Bins(int(fFitMassHi-fFitMassLo)/BIN_SIZE_FAIL));
  mframeFail->SetTitle("");
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail);
  modelFail->plotOn(mframeFail,Components("sigModFail"),LineStyle(kDashed),LineColor(kGreen));
  modelFail->plotOn(mframeFail,Components("bkg1ModFail"),LineStyle(kDashed),LineColor(kRed));
  modelFail->plotOn(mframeFail,Components("bkg2ModFail"),LineStyle(kDashed),LineColor(kGray));


  TCanvas *cpass = MakeCanvas("cpass","cpass",800,1200);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  cpass->SetTickx(1);
  cpass->SetTicky(1);
  mframePass->Draw();

  legend->Draw();

  string str_tmp = fOutputDir+"/TemplatesandData_Passsample_afterfit.png"; 
  cpass->SaveAs(str_tmp.c_str());
  str_tmp = fOutputDir+"/TemplatesandData_Passsample_afterfit.C"; 
  cpass->SaveAs(str_tmp.c_str());

  TCanvas *cfail = MakeCanvas("cfail","cfail",800,1200); 
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,0);
  cfail->SetTickx(1);
  cfail->SetTicky(1);  
  mframeFail->Draw();

  legend->Draw();
  str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit.png"; 
  cfail->SaveAs(str_tmp.c_str());
  str_tmp = fOutputDir+"/TemplatesandData_Failsample_afterfit.C"; 
  cfail->SaveAs(str_tmp.c_str());



  std::cout<<"++++++++++Fit++++++++++"<<std::endl;
  std::cout<<"Eff "<<resEff<<" - "<<resErrl<<" + "<<resErrh<<std::endl;
  std::cout<<"FR1 "<<FR1.getVal()<<" - "<<fabs(FR1.getErrorLo())<<" + "<<FR1.getErrorHi()<<std::endl;
  std::cout<<"FR2 "<<FR2.getVal()<<" - "<<fabs(FR2.getErrorLo())<<" + "<<FR2.getErrorHi()<<std::endl;

  std::cout<<"++++++++++Expect++++++++++"<<std::endl; //compute errors
  std::cout<<"Eff "<<h_pass_signal->Integral()/(h_pass_signal->Integral()+h_fail_signal->Integral())<<std::endl;
  std::cout<<"FR1 "<<h_pass_bkg1->Integral()/(h_pass_bkg1->Integral()+h_fail_bkg1->Integral())<<std::endl;
  std::cout<<"FR2 "<<h_pass_bkg2->Integral()/(h_pass_bkg2->Integral()+h_fail_bkg2->Integral())<<std::endl;

  std::cout<<"++++++++++Fit (extras)++++++++++"<<std::endl;
  std::cout<<"#chi^{2}/dof (Pass)"<<mframePass->chiSquare()<<std::endl;
  std::cout<<"#chi^{2}/dof (Fail)"<<mframeFail->chiSquare()<<std::endl;


    
  std::cout<<"*********Fitted Pass**********"<<std::endl;
  std::cout<<"Yields (data) "<<(int)passTree->GetEntries()<<std::endl;
  std::cout<<"Nsig "<<NsigPass.getVal()<<" +/- "<<NsigPass.getPropagatedError(*fitResult)<<std::endl;
  std::cout<<"Nbkg1 "<<Nbkg1Pass.getVal()<<" +/- "<<Nbkg1Pass.getPropagatedError(*fitResult)<<std::endl;
  std::cout<<"Nbkg2 "<<Nbkg2Pass.getVal()<<" +/- "<<Nbkg2Pass.getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
  std::cout<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigPass.getVal()/histPass.Integral()<<" "<<Nbkg1Pass.getVal()/histPass.Integral()<<" "<<Nbkg2Pass.getVal()/histPass.Integral()<<std::endl;
  std::cout<<"*********EXPECTED FRACMC Pass**********"<<std::endl;
  std::cout<<"Signal, Bkg1, Bkg2 "<<h_pass_signal->Integral()/histPass.Integral()<<" "<<h_pass_bkg1->Integral()/histPass.Integral()<<" "<<h_pass_bkg2->Integral()/histPass.Integral()<<std::endl;
  std::cout<<"*********Fitted Fail**********"<<std::endl;
  std::cout<<"Yields (data) "<<(int)failTree->GetEntries()<<std::endl;
  std::cout<<"Nsig "<<NsigFail.getVal()<<" +/- "<<NsigFail.getPropagatedError(*fitResult)<<std::endl;
  std::cout<<"Nbkg1 "<<Nbkg1Fail.getVal()<<" +/- "<<Nbkg1Fail.getPropagatedError(*fitResult)<<std::endl;
  std::cout<<"Nbkg2 "<<Nbkg2Fail.getVal()<<" +/- "<<Nbkg2Fail.getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
  std::cout<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigFail.getVal()/histFail.Integral()<<" "<<Nbkg1Fail.getVal()/histFail.Integral()<<" "<<Nbkg2Fail.getVal()/histFail.Integral()<<std::endl;

  std::cout<<"*********EXPECTED FRACMC Fail**********"<<std::endl;
  std::cout<<"Signal, Bkg1, Bkg2 "<<h_fail_signal->Integral()/histFail.Integral()<<" "<<h_fail_bkg1->Integral()/histFail.Integral()<<" "<<h_fail_bkg2->Integral()/histFail.Integral()<<std::endl;
 

  //
  // Write fit results
  //
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/fitres.txt",fOutputDir.c_str());
  txtfile.open(txtfname);
  assert(txtfile.is_open());



  txtfile<<"++++++++++Fit++++++++++"<<std::endl;
  txtfile<<"Eff "<<resEff<<" - "<<resErrl<<" + "<<resErrh<<std::endl;
  txtfile<<"FR1 "<<FR1.getVal()<<" - "<<fabs(FR1.getErrorLo())<<" + "<<FR1.getErrorHi()<<std::endl;
  txtfile<<"FR2 "<<FR2.getVal()<<" - "<<fabs(FR2.getErrorLo())<<" + "<<FR2.getErrorHi()<<std::endl;

  txtfile<<"++++++++++Expect++++++++++"<<std::endl; //compute errors
  txtfile<<"Eff "<<h_pass_signal->Integral()/(h_pass_signal->Integral()+h_fail_signal->Integral())<<std::endl;
  txtfile<<"FR1 "<<h_pass_bkg1->Integral()/(h_pass_bkg1->Integral()+h_fail_bkg1->Integral())<<std::endl;
  txtfile<<"FR2 "<<h_pass_bkg2->Integral()/(h_pass_bkg2->Integral()+h_fail_bkg2->Integral())<<std::endl;

  txtfile<<"++++++++++Fit (extras)++++++++++"<<std::endl;
  txtfile<<"#chi^{2}/dof (Pass)"<<mframePass->chiSquare()<<std::endl;
  txtfile<<"#chi^{2}/dof (Fail)"<<mframeFail->chiSquare()<<std::endl;

    
  txtfile<<"*********Fitted Pass**********"<<std::endl;
  txtfile<<"Yields (data) "<<(int)passTree->GetEntries()<<std::endl;
  txtfile<<"Nsig "<<NsigPass.getVal()<<" +/- "<<NsigPass.getPropagatedError(*fitResult)<<std::endl;
  txtfile<<"Nbkg1 "<<Nbkg1Pass.getVal()<<" +/- "<<Nbkg1Pass.getPropagatedError(*fitResult)<<std::endl;
  txtfile<<"Nbkg2 "<<Nbkg2Pass.getVal()<<" +/- "<<Nbkg2Pass.getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
  txtfile<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigPass.getVal()/histPass.Integral()<<" "<<Nbkg1Pass.getVal()/histPass.Integral()<<" "<<Nbkg2Pass.getVal()/histPass.Integral()<<std::endl;
  txtfile<<"*********EXPECTED FRACMC Pass**********"<<std::endl;
  txtfile<<"Signal, Bkg1, Bkg2 "<<h_pass_signal->Integral()/histPass.Integral()<<" "<<h_pass_bkg1->Integral()/histPass.Integral()<<" "<<h_pass_bkg2->Integral()/histPass.Integral()<<std::endl;
  txtfile<<"*********Fitted Fail**********"<<std::endl;
  txtfile<<"Yields (data) "<<(int)failTree->GetEntries()<<std::endl;
  txtfile<<"Nsig "<<NsigFail.getVal()<<" +/- "<<NsigFail.getPropagatedError(*fitResult)<<std::endl;
  txtfile<<"Nbkg1 "<<Nbkg1Fail.getVal()<<" +/- "<<Nbkg1Fail.getPropagatedError(*fitResult)<<std::endl;
  txtfile<<"Nbkg2 "<<Nbkg2Fail.getVal()<<" +/- "<<Nbkg2Fail.getPropagatedError(*fitResult)<<std::endl;
  //propagate error fraction  
  txtfile<<"FracNsig, FracNbkg1, FracNbkg2 "<<NsigFail.getVal()/histFail.Integral()<<" "<<Nbkg1Fail.getVal()/histFail.Integral()<<" "<<Nbkg2Fail.getVal()/histFail.Integral()<<std::endl;

  txtfile<<"*********EXPECTED FRACMC Fail**********"<<std::endl;
  txtfile<<"Signal, Bkg1, Bkg2 "<<h_fail_signal->Integral()/histFail.Integral()<<" "<<h_fail_bkg1->Integral()/histFail.Integral()<<" "<<h_fail_bkg2->Integral()/histFail.Integral()<<std::endl;

  txtfile <<"DONE WITH PRINTING MY RESULTS"<< endl;
  txtfile << endl;  txtfile << endl;
  txtfile <<"NOW DUMPING fit results from fitResult->printStream"<< endl;
  fitResult->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  txtfile.close();



  //
  // Clean up
  //

  
  delete modelPass;
  delete modelFail;  
  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigModPass_rdh;
  delete bkg1ModPass_rdh;
  delete bkg2ModPass_rdh;
  delete sigModFail_rdh;
  delete bkg1ModFail_rdh;
  delete bkg2ModFail_rdh;
  delete sigModPass;
  delete bkg1ModPass;  delete bkg2ModPass;
  delete sigModFail;
  delete bkg1ModFail;          delete bkg2ModFail;        
  delete histfile_signal;  delete histfile_bkg1;  delete histfile_bkg2;
  delete legend;
}

float CEffZFitter::calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3){
  float invmass = -999.;
  TLorentzVector jetSumLV;

  jetSumLV = (*obj1)+(*obj2)+(*obj3);

  invmass = jetSumLV.Mag();

  return invmass;

}

 
