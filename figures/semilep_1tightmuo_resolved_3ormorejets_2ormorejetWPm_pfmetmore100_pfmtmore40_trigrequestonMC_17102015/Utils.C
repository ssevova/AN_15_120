#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TH1D.h>                     // 1D histogram class
#include <TLorentzVector.h>           // 4-vector class
#include <TVector2.h>                 // 2-vector class
#include <TF1.h>
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

bool lorentzvecptorder (TLorentzVector i,TLorentzVector j) { return (i.Pt()>j.Pt()); }
bool floatorder (float i,float j) { return (i>j); }

float turnon(float var, int muorele){
  float output=1;

  string func = "([2])/(1+TMath::Exp(([0]-x)/[1]))";
  float xminfunc=0; float xmaxfunc=500;
  TF1 *f1 = new TF1("f1",func.c_str(),xminfunc,xmaxfunc);

  if(muorele==1){
    //HLT_MET120_HBHENoiseCleaned Vs pfmetnolep - mu
    f1->SetParameter(0,179.); //+/- 0.1762
    f1->SetParameter(1,18.41); //+/- 0.084
    f1->SetParameter(2,0.9972); //+/- 0.0008
  }
  else if(muorele==2){
    //HLT_MET120_HBHENoiseCleaned Vs pfmet -ele
    f1->SetParameter(0,181.2); //+/- 0.4978
    f1->SetParameter(1,20.94); //+/- 0.1467
    f1->SetParameter(2,0.9742); //+/- 0.005
  }

  if(muorele==1 || muorele==2)
    output = f1->Eval(var);

  return output;
}

double deltaPhi(const double phi1, const double phi2) {
  double result = phi1 - phi2;
  if     (result >  TMath::Pi()) { result = result - 2*TMath::Pi(); }
  else if(result < -TMath::Pi()) { result = result + 2*TMath::Pi(); }
  return result;
}

double CalcMinDPhiJet(double pfmetphi, TLorentzVector *res_jet1, TLorentzVector *res_jet2=0, TLorentzVector *res_jet3=0, TLorentzVector *res_jet4=0, TLorentzVector *res_jet5=0, TLorentzVector *res_jet6=0, bool verbose=false){
  if(verbose){
    cout<<"---CalcMinDPhiJet---"<<endl;
    cout<<"res_jet "<<res_jet1<<", "<<res_jet2<<", "<<res_jet3<<", "<<res_jet4<<", "<<res_jet5<<", "<<res_jet6<<", "<<endl;
  }
  double dphimetjet[7] = {999.,999.,999.,999.,999.,999.,999.}; //dphimetjet[0] is the minimum of the others
  if(res_jet1)
    dphimetjet[1] = TMath::Abs(deltaPhi(pfmetphi,res_jet1->Phi())); 
  if(res_jet2)
    dphimetjet[2] = TMath::Abs(deltaPhi(pfmetphi,res_jet2->Phi())); 
  if(res_jet3){
    //    cout<<"qui"<<endl;
    dphimetjet[3] = TMath::Abs(deltaPhi(pfmetphi,res_jet3->Phi())); 
  }
  if(res_jet4)
    dphimetjet[4] = TMath::Abs(deltaPhi(pfmetphi,res_jet4->Phi()));
  if(res_jet5)
    dphimetjet[5] = TMath::Abs(deltaPhi(pfmetphi,res_jet5->Phi()));
  if(res_jet6)
    dphimetjet[6] = TMath::Abs(deltaPhi(pfmetphi,res_jet6->Phi()));
  for(int ij=1; ij<=6; ij++){
    //cout<<"dphimetjet["<<ij<<"] "<<dphimetjet[ij]<<endl;
    if(dphimetjet[ij] < dphimetjet[0])
	     dphimetjet[0] = dphimetjet[ij];
  }

  if(verbose){
    for(int i=1; i<7; i++)
      cout<<"dphimetjet["<<i<<"] "<<dphimetjet[i]<<endl;
    cout<<"dphimetjet "<<dphimetjet[0]<<endl;
  }
  return dphimetjet[0];      
}

float calcmass(TLorentzVector *obj1, TLorentzVector *obj2, TLorentzVector *obj3){
  float invmass = -999.;
  TLorentzVector jetSumLV;

  jetSumLV = (*obj1)+(*obj2)+(*obj3);

  invmass = jetSumLV.Mag();

  return invmass;

}


float Calctopmass(unsigned int njets05, TLorentzVector *res_jet1, TLorentzVector *res_jet2=0, TLorentzVector *res_jet3=0, TLorentzVector *res_jet4=0, TLorentzVector *res_jet5=0, TLorentzVector *res_jet6=0){
  //next: can make use of of btagging info

  //assuming res_jet ordered in decreasing pt
  //  cout<<"Calctopmass "<<endl;

  float topmassbest = -999.;
  float topmass[20];
  for(int i=0; i<20; i++)
    topmass[i] = -999.;
  unsigned int countavailablejet = 0;


  if(res_jet1)
    countavailablejet++;
  if(res_jet2)
    countavailablejet++;
  if(res_jet3)
    countavailablejet++;
  if(res_jet4)
    countavailablejet++;
  if(res_jet5)
    countavailablejet++;
  if(res_jet6)
    countavailablejet++;

  if(countavailablejet > njets05){ //I save up to 6 Tlorentzvector jets, but can count more
    cout<<"countavailablejet, njets05 "<<countavailablejet<<", "<<njets05<<" are different...exit"<<endl;
    exit(1);
    }


  //  cout<<" res_jet1 "<<res_jet1<<" res_jet2 "<<res_jet2<<" res_jet3 "<<res_jet3<<" res_jet4 "<<res_jet4<<" res_jet5 "<<res_jet5<<" res_jet6 "<<res_jet6<<endl;

  //cout<<" res_jet6 "<<res_jet6->Pt()<<endl;

  if(countavailablejet >= 3){
    topmass[0] = calcmass(res_jet1,res_jet2,res_jet3);
    if(countavailablejet >= 4){
      topmass[1] = calcmass(res_jet4,res_jet2,res_jet3); //1 sub wrt 1,2,3
      topmass[2] = calcmass(res_jet1,res_jet4,res_jet3);
      topmass[3] = calcmass(res_jet1,res_jet2,res_jet4);
   
      if(countavailablejet >= 5){
	topmass[4] = calcmass(res_jet5,res_jet2,res_jet3); //1sub
	topmass[5] = calcmass(res_jet1,res_jet5,res_jet3);
	topmass[6] = calcmass(res_jet1,res_jet2,res_jet5);
	
	topmass[7] = calcmass(res_jet4,res_jet5,res_jet3); //2sub
	topmass[8] = calcmass(res_jet1,res_jet4,res_jet5);
	topmass[9] = calcmass(res_jet5,res_jet2,res_jet4);
    
	if(countavailablejet == 6){
	  topmass[10] = calcmass(res_jet6,res_jet2,res_jet3);//1sub
	  topmass[11] = calcmass(res_jet1,res_jet6,res_jet3);
	  topmass[12] = calcmass(res_jet1,res_jet2,res_jet6);

	  topmass[13] = calcmass(res_jet5,res_jet6,res_jet3); //2sub
	  topmass[14] = calcmass(res_jet1,res_jet5,res_jet6);
	  topmass[15] = calcmass(res_jet6,res_jet2,res_jet5);
	  topmass[16] = calcmass(res_jet4,res_jet6,res_jet3); //2sub
	  topmass[17] = calcmass(res_jet1,res_jet4,res_jet6);
	  topmass[18] = calcmass(res_jet6,res_jet2,res_jet4);

	  topmass[19] = calcmass(res_jet4,res_jet5,res_jet6); //3 sub
	}//if(countavailablejet == 6)
      }//if(countavailablejet >= 5)
    }//if(countavailablejet >= 4)
  }//if(countavailablejet >= 3)

  //select mass which is closer to nominal top mass
  for(int i=0; i<20; i++){
    if(TMath::Abs(topmass[i]-172.5)<TMath::Abs(topmassbest-172.5))
      topmassbest = topmass[i];
  }

  return topmassbest;      
}



//--------------------------------------------------------------------------------------------------

void makePlot(TCanvas *c, const string outname, const string xlabel, const string ylabel,
              const vector<TH1D*>& histv, const vector<CSample*>& samplev, TH1D* hExp, TH1D* hRatio, string str_rootfile, const double lumi, const bool doLogy, const double legdx, const double legdy,
              const double ymin, const double ymax, unsigned int verbose)
{
  if(verbose==1){
    cout<<"-----makePlot----"<<endl;
    cout<<"outname "<<outname<<endl;
  }
  TFile *outfile = NULL;
  if(str_rootfile != "")
    outfile = new TFile(str_rootfile.c_str(),"UPDATE");

  const int uncColor = kGray+3;
  
  //add overflow/underflow in the last bin/first bin
  for(unsigned int i=0; i<histv.size(); i++){
    int nbin = histv[i]->GetNbinsX();
    float overflow = histv[i]->GetBinContent(nbin+1); float underflow = histv[i]->GetBinContent(0);
    float overflowerr = histv[i]->GetBinError(nbin+1);    float underflowerr = histv[i]->GetBinError(0);
    float lastbin = histv[i]->GetBinContent(nbin);    float firstbin = histv[i]->GetBinContent(1);
    float lastbinerr = histv[i]->GetBinError(nbin);    float firstbinerr = histv[i]->GetBinError(1);
    float lastbinmod = lastbin + overflow;    float firstbinmod = firstbin + underflow;
    float lastbinerrmod = TMath::Sqrt(lastbinerr*lastbinerr + overflowerr*overflowerr);    float firstbinerrmod = TMath::Sqrt(firstbinerr*firstbinerr + underflowerr*underflowerr);
    histv[i]->SetBinContent(nbin,lastbinmod);    histv[i]->SetBinContent(1,firstbinmod);
    histv[i]->SetBinError(nbin,lastbinerrmod);    histv[i]->SetBinError(1,firstbinerrmod);
    
  }

  //cout<<"qui"<<endl;

  histv[0]->SetMarkerSize(0.9);
  
  CPlot plot(outname.c_str(),"",xlabel.c_str(),ylabel.c_str());  
  plot.AddHist1D(hExp,"E2",uncColor,1,3004);
  plot.AddHist1D(histv[0],samplev[0]->label,"E"); //blinded
  for(unsigned int i=1; i<histv.size()-1; i++) {
    if(verbose==1){
      cout<<"Sample name "<<samplev[i]->label<<endl;
      cout<<"i "<<i<<" histv[i]->Integral() "<<histv[i]->Integral()<<endl;
      cout<<"i "<<i<<" histv[i]->GetEntries() "<<histv[i]->GetEntries()<<endl;
    }
    plot.AddToStack(histv[i],samplev[i]->label,samplev[i]->fillcolor,samplev[i]->linecolor);
  }
  plot.AddHist1D(histv.back(),samplev.back()->label,"hist",samplev.back()->linecolor);
  
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1} (8 TeV)",lumi);
  plot.AddTextBox(lumitext,0.66,0.99,0.95,0.925,0,kBlack);
  //  plot.AddTextBox("CMS",0.18,0.88,0.30,0.82,0,kBlack,62);
  //  plot.AddTextBox("Preliminary",0.18,0.82,0.37,0.77,0,kBlack,52);
  
  plot.TransLegend(legdx, legdy);
  //  cout<<"qui2"<<endl;
  if(doLogy) {
    plot.SetLogy();
  }
  
  if(ymin!=ymax) {
    //cout<<"ymin "<<ymin<<" ymax "<<ymax<<endl;
    plot.SetYRange(ymin,ymax);
  }
  //cout<<"qui3"<<endl;

  if(verbose==1){
    for(int i=1; i<hRatio->GetNbinsX(); i++)
      cout<<"bin, hratio bc "<<i<<" "<<hRatio->GetBinContent(i)<<endl;
  }
  hRatio->SetMarkerSize(0.8);  
  const double xmin = histv[0]->GetXaxis()->GetBinLowEdge(1);
  const double xmax = histv[0]->GetXaxis()->GetBinUpEdge(histv[0]->GetNbinsX());
  CPlot plotRatio(outname.c_str(),"",xlabel.c_str(),"Data/Bkg");
  plotRatio.AddHist1D(hRatio,"EX0",kBlack);
  plotRatio.SetYRange(0,2);
  plotRatio.AddLine(xmin,1.5,xmax,1.5,kBlack,3);
  plotRatio.AddLine(xmin,1.0,xmax,1.0,kBlack,3);
  plotRatio.AddLine(xmin,0.5,xmax,0.5,kBlack,3);
  //cout<<"qui4"<<endl;
  plot.Draw(c,false,"png",1);
  //cout<<"qui5"<<endl;
  plotRatio.Draw(c,true,"png",2);
  //cout<<"qui6"<<endl;
  //saving into the rootfile
  if(outfile){
    //    cout<<"Writing into the root file"<<endl;
    for(unsigned int i=0; i<histv.size(); i++)
      histv[i]->Write();
    //    cout<<"Done writing"<<endl;
    outfile->Close();
  }
  //  cout<<"end makePlot"<<endl;

  // cout<<"qui4 "<<c<<endl;
  // plot.Draw(c,true,"png");
  //cout<<"qui5"<<endl;
  //  plotRatio.Draw(c,true,"pdf");
}

void makePlot_nostack(TCanvas *c, const string outname, const string xlabel, const string ylabel,
		      const vector<TH1D*>& histv, bool norm, const vector<CSample*>& samplev, TH1D* hExp,
		      const double lumi, const bool doLogy, const double legdx, const double legdy,
		      const double ymin, const double ymax) //[4] give it dynamically
{

  cout<<"---makePlot_nostack---"<<endl;
  const int uncColor = kGray+3;
  TString outname2 = TString(histv[0]->GetName())+"_"+outname;
  cout<<"histv[0]->GetName() "<<histv[0]->GetName()<<endl;
  cout<<"outname "<<outname<<endl;
  CPlot plot(outname2,"",xlabel.c_str(),ylabel.c_str());   
  //  TString rootfile = plot.GetsOutDir()+"/"+outname2+".root";
  //  TFile *fout = new TFile(rootfile,"RECREATE","");
  //cout<<"plot.GetsOutDir() "<<plot.GetsOutDir()<<endl;
  for(unsigned int i=0; i<histv.size(); i++) {
    histv[i]->SetMarkerSize(0.9);
    cout<<"name histo "<<histv[i]->GetName()<<endl;  
    plot.AddHist1D(hExp,"E2",uncColor,1,3004);
    if(norm){ //to be implemented better
      TH1D *htemp = (TH1D*)histv[i]->Clone(histv[i]->GetName());
      htemp->Scale(1./htemp->Integral());
      plot.AddHist1D(htemp,samplev[i]->label,"E",samplev[i]->linecolor,1,0);
    }
    else
      plot.AddHist1D(histv[i],samplev[i]->label,"E",samplev[i]->linecolor,1,0);
    //    fout->cd();
    //histv[i]->Write();

  }//for(unsigned int i=0; i<histv.size(); i++)
  
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1} (8 TeV)",lumi);
  plot.AddTextBox(lumitext,0.66,0.99,0.95,0.925,0,kBlack);
  //  plot.AddTextBox("CMS",0.18,0.88,0.30,0.82,0,kBlack,62);
  //  plot.AddTextBox("Preliminary",0.18,0.82,0.37,0.77,0,kBlack,52);
  
  plot.TransLegend(legdx, legdy);
  
  if(doLogy) {
    plot.SetLogy();
  }
  
  if(ymin!=ymax) {
    plot.SetYRange(ymin,ymax);
  }
  
  

  //  plotRatio.Draw(c,true,"pdf");
  

  plot.Draw(c,true,"png");
  //fout->Close();
}


TH1D* makeRatioHist(TH1D* hData, TH1D* hMC, const string name, unsigned int verbose)
{
  if(verbose==2)
    cout<<"---makeRatioHist--- "<<hData->GetName()<<endl;
  TH1D *hRatio = new TH1D(name.c_str(),"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    double numer = hData->GetBinContent(ibin);
    double denom = hMC->GetBinContent(ibin);
    double ratio = (denom>0) ? numer/denom : 0;
    if(verbose==2)
      cout<<"ibin "<<ibin<<" numer "<<numer<<" denom "<<denom<<" ratio "<<ratio<<endl;

    double errn = hData->GetBinError(ibin);
    double errd = hMC->GetBinError(ibin);        
    double err = (denom>0 && numer>0) ? ratio*sqrt(errn*errn/numer/numer + errd*errd/denom/denom) : 0;

    hRatio->SetBinContent(ibin,ratio);
    hRatio->SetBinError(ibin,err);
  }
  hRatio->GetYaxis()->SetTitleOffset(0.42);
  hRatio->GetYaxis()->SetTitleSize(0.13);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetNdivisions(104);
  hRatio->GetYaxis()->CenterTitle();
  hRatio->GetXaxis()->SetTitleOffset(1.2);
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetLabelSize(0.12);
  //hRatio->GetXaxis()->CenterTitle();

  return hRatio;
}


void makePlot2D(TCanvas *c, const string outname, const string xlabel, const string ylabel,
		TH2D* histv, CSample* samplev, TH2D* hExp, string str_rootfile, const double lumi, unsigned int verbose)
{
  //  if(verbose==1){
    cout<<"-----makePlot2D----"<<endl;
    cout<<"outname "<<outname<<" histv->GetName() "<<histv->GetName()<<endl;
    // }
  TFile *outfile = NULL;
  if(str_rootfile != "")
    outfile = new TFile(str_rootfile.c_str(),"UPDATE");


  ////weird plot to be fixed
  /*  const int uncColor = kGray+3;
  

  histv->SetMarkerSize(0.9);
  
  cout<<"histv->Integral() "<<histv->Integral()<<endl;

  CPlot plot(outname.c_str(),"",xlabel.c_str(),ylabel.c_str());  
  plot.AddHist2D(hExp,"E2",uncColor,1);
  plot.AddHist2D(histv,"E");
  
  char lumitext[100];
  sprintf(lumitext,"%.1f fb^{-1} (8 TeV)",lumi);
  plot.AddTextBox(lumitext,0.66,0.99,0.95,0.925,0,kBlack);
  //  plot.AddTextBox("CMS",0.18,0.88,0.30,0.82,0,kBlack,62);
  //  plot.AddTextBox("Preliminary",0.18,0.82,0.37,0.77,0,kBlack,52);
  

  cout<<"qui4"<<endl;
  plot.Draw(c,true,"png"); */

  if(outfile){
    cout<<"Writing into the root file"<<endl;
    histv->Write();
    cout<<"Done writing"<<endl;
    outfile->Close();
  }
  cout<<"end makePlot2D"<<endl;

  // cout<<"qui4 "<<c<<endl;
  // plot.Draw(c,true,"png");
  //cout<<"qui5"<<endl;
  //  plotRatio.Draw(c,true,"pdf");
}



float computeMt(const float pt, const float phi, const float met, const float metphi)
{
  TLorentzVector vLep; vLep.SetPtEtaPhiM(pt,0,phi,0);
  TLorentzVector vMet; vMet.SetPtEtaPhiM(met,0,metphi,0);
  double dphi = vLep.DeltaPhi(vMet);
  return sqrt( 2.0* pt * met *(1.0-cos(dphi)) );
}

float Wptweight(float wgpt){

  //derived as follows
  //compareshapes("ttbarsemileCR/muo/allhadronic_1tightmu_resolved_4ormorejets_2ormorejetWPm_pfmet180_pfmtless100_trigrequestonMC_genlevel_21092015/","SR/allhadronic_0tightlep_resolved_4ormorejets_2ormorejetWPm_pfmet180_mindphimetjetmore1_notrigrequestonMC_parmu_genlevel_21092015/","hWPtlinearv_1")
  float weight;

  if (wgpt < 0)
    weight = 0;
  if (wgpt >= 0 && wgpt < 10)
    weight = 0;
  if (wgpt >= 10 && wgpt < 20)
    weight = 0;
  if (wgpt >= 20 && wgpt < 30)
    weight = 0;
  if (wgpt >= 30 && wgpt < 40)
    weight = 0;
  if (wgpt >= 40 && wgpt < 50)
    weight = 0;
  if (wgpt >= 50 && wgpt < 60)
    weight = 1.01657;
  if (wgpt >= 60 && wgpt < 70)
    weight = 0.0603251;
  if (wgpt >= 70 && wgpt < 80)
    weight = 0.13954;
  if (wgpt >= 80 && wgpt < 90)
    weight = 0.147986;
  if (wgpt >= 90 && wgpt < 100)
    weight = 0.106036;
  if (wgpt >= 100 && wgpt < 110)
    weight = 0.610301;
  if (wgpt >= 110 && wgpt < 120)
    weight = 1.35605;
  if (wgpt >= 120 && wgpt < 130)
    weight = 2.6694;
  if (wgpt >= 130 && wgpt < 140)
    weight = 3.87264;
  if (wgpt >= 140 && wgpt < 150)
    weight = 5.71802;
  if (wgpt >= 150 && wgpt < 160)
    weight = 4.77864;
  if (wgpt >= 160 && wgpt < 170)
    weight = 4.65914;
  if (wgpt >= 170 && wgpt < 180)
    weight = 3.99841;
  if (wgpt >= 180 && wgpt < 190)
    weight = 3.3691;
  if (wgpt >= 190 && wgpt < 200)
    weight = 2.68324;
  if (wgpt >= 200 && wgpt < 210)
    weight = 2.07273;
  if (wgpt >= 210 && wgpt < 220)
    weight = 1.56676;
  if (wgpt >= 220 && wgpt < 230)
    weight = 1.2055;
  if (wgpt >= 230 && wgpt < 240)
    weight = 1.01517;
  if (wgpt >= 240 && wgpt < 250)
    weight = 0.803723;
  if (wgpt >= 250 && wgpt < 260)
    weight = 0.70951;
  if (wgpt >= 260 && wgpt < 270)
    weight = 0.587388;
  if (wgpt >= 270 && wgpt < 280)
    weight = 0.56441;
  if (wgpt >= 280 && wgpt < 290)
    weight = 0.531869;
  if (wgpt >= 290 && wgpt < 300)
    weight = 0.438924;
  if (wgpt >= 300 && wgpt < 310)
    weight = 0.420812;
  if (wgpt >= 310 && wgpt < 320)
    weight = 0.32246;
  if (wgpt >= 320 && wgpt < 330)
    weight = 0.379156;
  if (wgpt >= 330 && wgpt < 340)
    weight = 0.278408;
  if (wgpt >= 340 && wgpt < 350)
    weight = 0.231482;
  if (wgpt >= 350 && wgpt < 360)
    weight = 0.251567;
  if (wgpt >= 360 && wgpt < 370)
    weight = 0.189346;
  if (wgpt >= 370 && wgpt < 380)
    weight = 0.145539;
  if (wgpt >= 380 && wgpt < 390)
    weight = 0.207189;
  if (wgpt >= 390 && wgpt < 400)
    weight = 0.105004;
  if (wgpt >= 400 && wgpt < 410)
    weight = 0.177698;
  if (wgpt >= 410 && wgpt < 420)
    weight = 0.182866;
  if (wgpt >= 420 && wgpt < 430)
    weight = 0.18189;
  if (wgpt >= 430 && wgpt < 440)
    weight = 0.147903;
  if (wgpt >= 440 && wgpt < 450)
    weight = 0.104074;
  if (wgpt >= 450 && wgpt < 460)
    weight = 0.0790584;
  if (wgpt >= 460 && wgpt < 470)
    weight = 0.1482;
  if (wgpt >= 470 && wgpt < 480)
    weight = 0.0981846;
  if (wgpt >= 480 && wgpt < 490)
    weight = 0.0528554;
  if (wgpt >= 490)
    weight = 0.0640049;

  //  cout<<"wgpt "<<wgpt<<" weight "<<weight<<endl;

  return weight;
}



float Nuptweight_muo(float var){                                                                                      //    root [1] compareshapes("SR/allhadronic_0tightlep_resolved_4ormorejets_2ormorejetWPm_pfmet180_mindphimetjetmore1_notrigrequestonMC_parmu_genlevel_21092015/","ttbarsemileCR/muo/allhadronic_1tightmu_resolved_4ormorejets_2ormorejetWPm_pfmet180_pfmtless100_trigrequestonMC_genlevel_21092015/","hNuPtlinearv_1")
  
  float weight;  
  
  if (var < 0)
    weight = 0;
  if (var >= 0 && var < 10)
    weight = 0.602232;
  if (var >= 10 && var < 20)
    weight = 1.02456;
  if (var >= 20 && var < 30)
    weight = 0.722626;
  if (var >= 30 && var < 40)
    weight = 1.07706;
  if (var >= 40 && var < 50)
    weight = 1.23765;
  if (var >= 50 && var < 60)
    weight = 1.35448;
  if (var >= 60 && var < 70)
    weight = 1.24221;
  if (var >= 70 && var < 80)
    weight = 1.29993;
  if (var >= 80 && var < 90)
    weight = 1.15848;
  if (var >= 90 && var < 100)
    weight = 1.74289;
  if (var >= 100 && var < 110)
    weight = 1.16585;
  if (var >= 110 && var < 120)
    weight = 1.17082;
  if (var >= 120 && var < 130)
    weight = 1.12403;
  if (var >= 130 && var < 140)
    weight = 0.978057;
  if (var >= 140 && var < 150)
    weight = 1.02212;
  if (var >= 150 && var < 160)
    weight = 1.08675;
  if (var >= 160 && var < 170)
    weight = 1.13336;
  if (var >= 170 && var < 180)
    weight = 1.06129;
  if (var >= 180 && var < 190)
    weight = 1.04891;
  if (var >= 190 && var < 200)
    weight = 1.0859;
  if (var >= 200 && var < 210)
    weight = 1.03307;
  if (var >= 210 && var < 220)
    weight = 0.925138;
  if (var >= 220 && var < 230)
    weight = 1.05391;
  if (var >= 230 && var < 240)
    weight = 0.821367;
  if (var >= 240 && var < 250)
    weight = 0.923456;
  if (var >= 250 && var < 260)
    weight = 0.988058;
  if (var >= 260 && var < 270)
    weight = 0.823316;
  if (var >= 270 && var < 280)
    weight = 0.805477;
  if (var >= 280 && var < 290)
    weight = 0.688109;
  if (var >= 290 && var < 300)
    weight = 0.726746;
  if (var >= 300 && var < 310)
    weight = 0.678557;
  if (var >= 310 && var < 320)
    weight = 0.53124;
  if (var >= 320 && var < 330)
    weight = 0.486138;
  if (var >= 330 && var < 340)
    weight = 0.463498;
  if (var >= 340 && var < 350)
    weight = 0.450141;
  if (var >= 350 && var < 360)
    weight = 0.522469;
  if (var >= 360 && var < 370)
    weight = 0.314685;
  if (var >= 370 && var < 380)
    weight = 0.581282;
  if (var >= 380 && var < 390)
    weight = 0.579348;
  if (var >= 390 && var < 400)
    weight = 0.450015;
  if (var >= 400 && var < 410)
    weight = 0.321663;
  if (var >= 410 && var < 420)
    weight = 0.726501;
  if (var >= 420 && var < 430)
    weight = 0.645924;
  if (var >= 430 && var < 440)
    weight = 0.200121;
  if (var >= 440 && var < 450)
    weight = 0;
  if (var >= 450 && var < 460)
    weight = 0.177381;
  if (var >= 460 && var < 470)
    weight = 0.5279;
  if (var >= 470 && var < 480)
    weight = 0.590355;
  if (var >= 480 && var < 490)
    weight = 0;
  if (var >= 490)
    weight = 0.103886;

  return weight;
}

float Nuptweight_ele(float var){

  //compareshapes("SR/allhadronic_0tightlep_resolved_4ormorejets_2ormorejetWPm_pfmet180_mindphimetjetmore1_notrigrequestonMC_parmu_genlevel_21092015/","ttbarsemileCR/ele/allhadronic_1tightele_resolved_4ormorejets_2ormorejetWPm_pfmet180_pfmtless100_notrigrequestonMC_parele_genlevel_21092015","hNuPtlinearv_1")
  float weight;

  if (var < 0)
      weight = 0;
  if (var >= 0 && var < 10)
    weight = 0.827279;
  if (var >= 10 && var < 20)
    weight = 0.83067;
  if (var >= 20 && var < 30)
    weight = 1.07761;
  if (var >= 30 && var < 40)
    weight = 1.07246;
  if (var >= 40 && var < 50)
    weight = 1.13964;
  if (var >= 50 && var < 60)
    weight = 1.32284;
  if (var >= 60 && var < 70)
    weight = 1.14489;
  if (var >= 70 && var < 80)
    weight = 1.28684;
  if (var >= 80 && var < 90)
    weight = 1.13264;
  if (var >= 90 && var < 100)
    weight = 1.8438;
  if (var >= 100 && var < 110)
    weight = 1.50657;
  if (var >= 110 && var < 120)
    weight = 1.42022;
  if (var >= 120 && var < 130)
    weight = 1.38404;
  if (var >= 130 && var < 140)
    weight = 1.19959;
  if (var >= 140 && var < 150)
    weight = 1.31523;
  if (var >= 150 && var < 160)
    weight = 1.38668;
  if (var >= 160 && var < 170)
    weight = 1.35524;
  if (var >= 170 && var < 180)
    weight = 1.17955;
  if (var >= 180 && var < 190)
    weight = 1.15303;
  if (var >= 190 && var < 200)
    weight = 1.09293;
  if (var >= 200 && var < 210)
    weight = 1.06926;
  if (var >= 210 && var < 220)
    weight = 0.890696;
  if (var >= 220 && var < 230)
    weight = 0.9383;
  if (var >= 230 && var < 240)
    weight = 0.743145;
  if (var >= 240 && var < 250)
    weight = 0.774953;
  if (var >= 250 && var < 260)
    weight = 0.71939;
  if (var >= 260 && var < 270)
    weight = 0.706429;
  if (var >= 270 && var < 280)
    weight = 0.610216;
  if (var >= 280 && var < 290)
    weight = 0.588262;
  if (var >= 290 && var < 300)
    weight = 0.501927;
  if (var >= 300 && var < 310)
    weight = 0.498859;
  if (var >= 310 && var < 320)
    weight = 0.434174;
  if (var >= 320 && var < 330)
    weight = 0.312839;
  if (var >= 330 && var < 340)
    weight = 0.350397;
  if (var >= 340 && var < 350)
    weight = 0.248306;
  if (var >= 350 && var < 360)
    weight = 0.344745;
  if (var >= 360 && var < 370)
    weight = 0.258001;
  if (var >= 370 && var < 380)
    weight = 0.414676;
  if (var >= 380 && var < 390)
    weight = 0.382955;
  if (var >= 390 && var < 400)
    weight = 0.250412;
  if (var >= 400 && var < 410)
    weight = 0.2029;
  if (var >= 410 && var < 420)
    weight = 0.245635;
  if (var >= 420 && var < 430)
    weight = 0.205321;
  if (var >= 430 && var < 440)
    weight = 0.236067;
  if (var >= 440 && var < 450)
    weight = 0;
  if (var >= 450 && var < 460)
    weight = 0.137666;
  if (var >= 460 && var < 470)
    weight = 0.644619;
  if (var >= 470 && var < 480)
    weight = 0.787524;
  if (var >= 480 && var < 490)
    weight = 0;
  if (var >= 490)
    weight = 0.0598389;

  return weight;
}

float NLOcorrphotoZ(float var){

  //from Photon_Z_NLO_kfactors.root Z_pho_NLO
  float weight;
  if (var < 110)
    weight = 0.124126;
  if (var >= 150 && var< 175)
    weight = 0.171022;
  if (var>= 175 && var< 200)
    weight = 0.186607;
  if (var>= 200 && var< 225)
    weight = 0.197397;
  if ( var>= 225 && var< 250 )
    weight = 0.206985;
  if (var>= 250 && var< 275 )
    weight = 0.212915;
  if (var>= 275 && var< 300)
    weight = 0.218689;
  if (var>= 300 && var< 350)
    weight = 0.225314;
  if (var>= 350 && var< 400 )
    weight = 0.23158;
  if ( var>= 400 && var< 550 )
    weight = 0.234917;		
  if (var>= 550 ) 
    weight = 0.24054;
					    
  return weight;


}

float isHadronicTop(TLorentzVector *res_jet1, TLorentzVector *res_jet2, TLorentzVector *res_jet3, TLorentzVector *genpar1, TLorentzVector *genpar2, TLorentzVector *genpar3, int genId1, int genId2, int genId3, bool verbose){

     //TMVA top tagger: assuming that res_jet1, 2, 3 are the triplet associated to highest MVA score
     //assuming the genpar1, genpar2, genpar3 are from the hadronic top
     //count the number of match -> if three found hadronic top
  if(verbose)
    cout<<"------isHadronicTop------"<<endl;


  float isHadrTopHighestMVA=0.;


  TLorentzVector resjet_tmp[3] = {*res_jet1, *res_jet2, *res_jet3};
  TLorentzVector genpar_tmp[3] = {*genpar1, *genpar2, *genpar3};
  int genpid_tmp[3] = {genId1, genId2, genId3};
  int match_jtop[3] = {0,0,0};
  int counttop = 0;

  for(int ij = 0; ij<3; ij++ ){
    for(int ip = 0; ip<3; ip++ ){
      if(resjet_tmp[ij].DeltaR(genpar_tmp[ip])<0.3){
	if(verbose)
	  cout<<"ij "<<ij<<" matches part "<<ip<<", whose pdg is "<<genpid_tmp[ip]<<endl;
	match_jtop[ij] += 1;	       
      }
	 }//for(int ip = 0; ip<3; ip++ ){
  }//for(int ij = 0; ij<3; ij++ ){
  
  //a jet matching two or more particles is allowd (counting as one though)
  for(int ij = 0; ij<3; ij++ ){
    if(match_jtop[ij] > 0)
      counttop += 1; 
    /*    else{ //if a jet matdouble match is a failure
      cout<<"********double match found for jet***** "<<ij<<endl;
      cout<<"deltaR w/ 3 part "<<resjet_tmp[ij].DeltaR(genpar_tmp[0])<<", "<<resjet_tmp[ij].DeltaR(genpar_tmp[1])<<", "<<resjet_tmp[ij].DeltaR(genpar_tmp[2])<<endl;
      //      exit(1);
      }*/
  }//for(int ij = 0; ij<3; ij++ )
  
  if(counttop == 3){
    isHadrTopHighestMVA=1.;
    if(verbose)
      cout<<"hadronic top is associated to highest MVA score"<<endl;
  }
  else{
    isHadrTopHighestMVA=0.;
    if(verbose)
      cout<<"hadronic top is not associated to highest MVA score "<<counttop<<endl;
  }


  if(verbose)
    cout<<"------end of isHadronicTop------"<<endl;


  return isHadrTopHighestMVA;	
 
}
