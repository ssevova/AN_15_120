void parseFitResults(std::ifstream& ifs, vector<double> *fiteff, vector<double> *expeff, bool verbose)
{

  //assumin txtfile is
// ++++++++++Fit++++++++++
// Eff 1 - 0.0097367 + 0.0097367
// FR1 0.560602 - 0.0662585 + 0.0662585
// FR2 0.387178 - 0.0763985 + 0.0763985
// ++++++++++Expect++++++++++
// Eff 0.947539
// FR1 0.584033
// FR2 0.372671
// ...

  if(verbose)
    cout<<"--------parseFitResults------"<<endl;

  std::string line;
  int fit_or_expect = -1; //0=fit, 1=expect
  int counterfit=0;   int counterexp=0;
  while(getline(ifs,line)) {
    //    cout<<"line "<<line<<endl;
    size_t found = line.find("+Fit+");
    if(verbose)
      cout<<"found "<<found<<endl;
    if(found!=string::npos){ 
      fit_or_expect=0;
      continue;
    }
    found = line.find("+Expect+");
    if(found!=string::npos) {
      fit_or_expect=1;
      continue;
    }
    if(verbose)
      cout<<"fit_or_expect "<<fit_or_expect<<endl;
    std::string varname, value, signlow, errlow, signhigh, errhigh;
    std::stringstream ss(line);
    ss >> varname >> value >> signlow >> errlow >> signhigh >> errhigh;    

    if(verbose){
      cout<<"varname "<<varname<<" value "<<value<<" signlow  "<<signlow<<" errlow  "<<errlow<<" signhigh  "<<signhigh<<" errhigh "<<errhigh<<endl;
 
      cout<<"fit_or_expect "<<fit_or_expect<<endl;
      cout<<"counterfit "<<counterfit<<" counterexp "<<counterexp<<endl;
    }
    if(fit_or_expect==0 && counterfit<3){
      counterfit++;
      fiteff->push_back(atof(value.c_str())); fiteff->push_back(atof(errlow.c_str())); fiteff->push_back(atof(errhigh.c_str()));
    }
    else if(fit_or_expect==1 && counterexp<3){
      counterexp++;
      expeff->push_back(atof(value.c_str())); expeff->push_back(atof(errlow.c_str())); expeff->push_back(atof(errhigh.c_str()));
    }
    else
      continue;
  }//if(found!=string::npos)


  if(verbose)
    cout<<"--------Done parseFitResults------"<<endl;

}


MakeGraph(){


  const int n = 9; //#files

  double xval[n], xerr[n];
  double yval[n], yerrl[n], yerrh[n];

  vector<double> fiteff,  expeff;

    string str_tmvacut[n] = {"-0.8","-0.6","-0.4","-0.2","0.0","0.2","0.4","0.6","0.8"}; //make dir to second digit
  //string str_tmvacut[n] = {"-0.6"}; //make dir to second digit
  double tmvacut[n], tmvacuterr[n];
  for(int i=0; i<n; i++){
    tmvacut[i] = atof(str_tmvacut[i].c_str());
    tmvacuterr[i] = 0.1;
  }
  double fit_eff[n], fit_eff_errlow[n], fit_eff_errhigh[n]; 
  double fit_fr1[n], fit_fr1_errlow[n], fit_fr1_errhigh[n]; 
  double fit_fr2[n], fit_fr2_errlow[n], fit_fr2_errhigh[n]; 
  double exp_eff[n], exp_eff_errlow[n], exp_eff_errhigh[n]; 
  double exp_fr1[n], exp_fr1_errlow[n], exp_fr1_errhigh[n]; 
  double exp_fr2[n], exp_fr2_errlow[n], exp_fr2_errhigh[n]; 

  //read txt files
  for(int i=0; i<n; i++){
    fiteff.clear();    expeff.clear();
    cout<<"i "<<i<<" str_tmvacut "<<str_tmvacut[i]<<endl;
    ifstream rfile;
    char rname[100];
    string fOutputDir = "TOPRESOLVEDTAGGER_"+str_tmvacut[i];
    sprintf(rname,"%s/fitres.txt",fOutputDir.c_str());
    rfile.open(rname);

    if(rfile.is_open()) {
      parseFitResults(rfile, &fiteff, &expeff,false);
      rfile.close();
    }
    else{
      cout<<"Error file "<<rname<<" not found "<<endl;
      exit(1);
    }
    

    fit_eff[i]  = fiteff[0];
    fit_eff_errlow[i]  = TMath::Abs(fiteff[1]);
    fit_eff_errhigh[i]  = fiteff[2];
    fit_fr1[i]  = fiteff[3];
    fit_fr1_errlow[i]  = TMath::Abs(fiteff[4]);
    fit_fr1_errhigh[i]  = fiteff[5];
    fit_fr2[i]  = fiteff[6];
    fit_fr2_errlow[i]  = TMath::Abs(fiteff[7]);
    fit_fr2_errhigh[i]  = fiteff[8];

    exp_eff[i]  = expeff[0];
    exp_eff_errlow[i]  = TMath::Abs(expeff[1]);
    exp_eff_errhigh[i]  = expeff[2];
    exp_fr1[i]  = expeff[3];
    exp_fr1_errlow[i]  = TMath::Abs(expeff[4]);
    exp_fr1_errhigh[i]  = expeff[5];
    exp_fr2[i]  = expeff[6];
    exp_fr2_errlow[i]  = TMath::Abs(expeff[7]);
    exp_fr2_errhigh[i]  = expeff[8];


  }
  
  for(int i=0; i<n; i++){
    cout<<"tmvacut["<<i<<"] "<<tmvacut[i]<<" err "<<tmvacuterr[i]<<endl;
    cout<<"fit_eff["<<i<<"] "<<fit_eff[i]<<" - "<<fit_eff_errlow[i]<<" + "<<fit_eff_errhigh[i]<<endl;
    cout<<"fit_fr1["<<i<<"] "<<fit_fr1[i]<<" - "<<fit_fr1_errlow[i]<<" + "<<fit_fr1_errhigh[i]<<endl;
    cout<<"fit_fr2["<<i<<"] "<<fit_fr2[i]<<" - "<<fit_fr2_errlow[i]<<" + "<<fit_fr2_errhigh[i]<<endl;
    cout<<"exp_eff["<<i<<"] "<<exp_eff[i]<<" - "<<exp_eff_errlow[i]<<" + "<<exp_eff_errhigh[i]<<endl;
    cout<<"exp_fr1["<<i<<"] "<<exp_fr1[i]<<" - "<<exp_fr1_errlow[i]<<" + "<<exp_fr1_errhigh[i]<<endl;
    cout<<"exp_fr2["<<i<<"] "<<exp_fr2[i]<<" - "<<exp_fr2_errlow[i]<<" + "<<exp_fr2_errhigh[i]<<endl;

  }

  TGraphAsymmErrors *fit_eff_gr = new TGraphAsymmErrors(n,tmvacut,fit_eff,tmvacuterr,tmvacuterr,fit_eff_errlow,fit_eff_errhigh);
  TGraphAsymmErrors *fit_fr1_gr = new TGraphAsymmErrors(n,tmvacut,fit_fr1,tmvacuterr,tmvacuterr,fit_fr1_errlow,fit_fr1_errhigh);
  TGraphAsymmErrors *fit_fr2_gr = new TGraphAsymmErrors(n,tmvacut,fit_fr2,tmvacuterr,tmvacuterr,fit_fr2_errlow,fit_fr2_errhigh);

  TGraphAsymmErrors *exp_eff_gr = new TGraphAsymmErrors(n,tmvacut,exp_eff,tmvacuterr,tmvacuterr,exp_eff_errlow,exp_eff_errhigh);
  TGraphAsymmErrors *exp_fr1_gr = new TGraphAsymmErrors(n,tmvacut,exp_fr1,tmvacuterr,tmvacuterr,exp_fr1_errlow,exp_fr1_errhigh);
  TGraphAsymmErrors *exp_fr2_gr = new TGraphAsymmErrors(n,tmvacut,exp_fr2,tmvacuterr,tmvacuterr,exp_fr2_errlow,exp_fr2_errhigh);


  fit_eff_gr->SetTitle("");
  fit_eff_gr->GetXaxis()->SetTitle("TMVA cut");
  fit_eff_gr->GetYaxis()->SetTitle("Obs. Efficiency [%]");
  fit_eff_gr->GetXaxis()->SetTitleOffset(1.4);
  fit_fr1_gr->GetYaxis()->SetTitleOffset(1.2);
  fit_eff_gr->SetLineColor(kBlue);
  fit_eff_gr->SetLineWidth(4);

  fit_fr1_gr->SetTitle("");
  fit_fr1_gr->GetXaxis()->SetTitle("TMVA cut");
  fit_fr1_gr->GetYaxis()->SetTitle("Obs. Fake Rate tt1l Combinatorial Bkg [%]");
  fit_fr1_gr->GetXaxis()->SetTitleOffset(1.4);
  fit_fr1_gr->GetYaxis()->SetTitleOffset(1.2);
  fit_fr1_gr->SetLineColor(kBlue);
  fit_fr1_gr->SetLineWidth(4);

  fit_fr2_gr->SetTitle("");
  fit_fr2_gr->GetXaxis()->SetTitle("TMVA cut");
  fit_fr2_gr->GetYaxis()->SetTitle("Obs. Fake Rate non-tt1l Bkg [%]");
  fit_fr2_gr->GetXaxis()->SetTitleOffset(1.4);
  fit_fr2_gr->GetYaxis()->SetTitleOffset(1.2);
  fit_fr2_gr->SetLineColor(kBlue);
  fit_fr2_gr->SetLineWidth(4);

  exp_eff_gr->SetTitle("");
  exp_eff_gr->GetXaxis()->SetTitle("TMVA cut");
  exp_eff_gr->GetYaxis()->SetTitle("Obs. Efficiency [%]");
  exp_eff_gr->GetXaxis()->SetTitleOffset(1.4);
  exp_fr1_gr->GetYaxis()->SetTitleOffset(1.2);
  exp_eff_gr->SetLineColor(kRed);
  exp_eff_gr->SetLineWidth(4);

  exp_fr1_gr->SetTitle("");
  exp_fr1_gr->GetXaxis()->SetTitle("TMVA cut");
  exp_fr1_gr->GetYaxis()->SetTitle("Obs. Fake Rate tt1l Combinatorial Bkg [%]");
  exp_fr1_gr->GetXaxis()->SetTitleOffset(1.4);
  exp_fr1_gr->GetYaxis()->SetTitleOffset(1.2);
  exp_fr1_gr->SetLineColor(kRed);
  exp_fr1_gr->SetLineWidth(4);

  exp_fr2_gr->SetTitle("");
  exp_fr2_gr->GetXaxis()->SetTitle("TMVA cut");
  exp_fr2_gr->GetYaxis()->SetTitle("Obs. Fake Rate non-tt1l Bkg [%]");
  exp_fr2_gr->GetXaxis()->SetTitleOffset(1.4);
  exp_fr2_gr->GetYaxis()->SetTitleOffset(1.2);
  exp_fr2_gr->SetLineColor(kRed);
  exp_fr2_gr->SetLineWidth(4);

  TLegend *legend = new TLegend(0.65,0.24,0.75,0.54);
  legend->SetFillColor(0);
  legend->AddEntry(fit_eff_gr,"Obs." , "lp");
  legend->AddEntry(exp_eff_gr,"Exp." , "lp");

  TLegend *legend2 = new TLegend(0.65,0.54,0.75,0.84);
  legend2->SetFillColor(0);
  legend2->AddEntry(fit_eff_gr,"Obs." , "lp");
  legend2->AddEntry(exp_eff_gr,"Exp." , "lp");




  TMultiGraph *eff_gr = new TMultiGraph();
  eff_gr->Add(fit_eff_gr);
  eff_gr->Add(exp_eff_gr);
  TCanvas *c_eff = new TCanvas("c_eff","c_eff",800,1200);
  c_eff->SetWindowPosition(c_eff->GetWindowTopX()+c_eff->GetBorderSize()+800,0);
  c_eff->SetTickx(1);
  c_eff->SetTicky(1);

  eff_gr->Draw("AP");
  eff_gr->GetXaxis()->SetTitle("TMVA cut");
  eff_gr->GetYaxis()->SetTitle("Efficiency [%]");
  legend->Draw();
  c_eff->SaveAs("c_eff.png");
  
  TMultiGraph *fr1_gr = new TMultiGraph();
  fit_fr1_gr->SetLineColor(kBlue);
  exp_fr1_gr->SetLineColor(kRed);
  fr1_gr->Add(fit_fr1_gr);
  fr1_gr->Add(exp_fr1_gr);
  TCanvas *c_fr1 = new TCanvas("c_fr1");
  fr1_gr->Draw("AP");
  fr1_gr->GetXaxis()->SetTitle("TMVA cut");
  fr1_gr->GetYaxis()->SetTitle("Fake Rate tt1l Combinatorial Bkg [%]");
  legend2->Draw(); 
  c_fr1->SaveAs("c_fr1.png");

  TMultiGraph *fr2_gr = new TMultiGraph();
  fit_fr2_gr->SetLineColor(kBlue);
  exp_fr2_gr->SetLineColor(kRed);
  fr2_gr->Add(fit_fr2_gr);
  fr2_gr->Add(exp_fr2_gr);
  TCanvas *c_fr2 = new TCanvas("c_fr2");
  fr2_gr->Draw("AP");
  fr2_gr->GetXaxis()->SetTitle("TMVA cut");
  fr2_gr->GetYaxis()->SetTitle("Fake rate non-tt1l Bkg [%]");
  legend2->Draw(); 
  c_fr2->SaveAs("c_fr2.png");
  

}
