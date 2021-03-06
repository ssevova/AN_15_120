{
//=========Macro generated from canvas: cfail_bf_mt/cfail_bf_mt
//=========  (Sat Oct 31 01:57:58 2015) by ROOT version5.32/00
   TCanvas *cfail_bf_mt = new TCanvas("cfail_bf_mt", "cfail_bf_mt",1,1,800,1176);
   gStyle->SetOptStat(0);
   cfail_bf_mt->SetHighLightColor(2);
   cfail_bf_mt->Range(-150,-26.23429,850,175.5679);
   cfail_bf_mt->SetFillColor(0);
   cfail_bf_mt->SetBorderMode(0);
   cfail_bf_mt->SetBorderSize(10);
   cfail_bf_mt->SetTickx(1);
   cfail_bf_mt->SetTicky(1);
   cfail_bf_mt->SetLeftMargin(0.15);
   cfail_bf_mt->SetRightMargin(0.05);
   cfail_bf_mt->SetTopMargin(0.08);
   cfail_bf_mt->SetBottomMargin(0.13);
   cfail_bf_mt->SetFrameFillStyle(0);
   cfail_bf_mt->SetFrameLineStyle(0);
   cfail_bf_mt->SetFrameBorderMode(0);
   cfail_bf_mt->SetFrameBorderSize(10);
   cfail_bf_mt->SetFrameFillStyle(0);
   cfail_bf_mt->SetFrameLineStyle(0);
   cfail_bf_mt->SetFrameBorderMode(0);
   cfail_bf_mt->SetFrameBorderSize(10);
   
   TH1D *histFail = new TH1D("histFail","",40,0,800);
   histFail->SetBinContent(4,1);
   histFail->SetBinContent(5,5);
   histFail->SetBinContent(6,4);
   histFail->SetBinContent(7,23);
   histFail->SetBinContent(8,44);
   histFail->SetBinContent(9,51);
   histFail->SetBinContent(10,73);
   histFail->SetBinContent(11,91);
   histFail->SetBinContent(12,129);
   histFail->SetBinContent(13,140);
   histFail->SetBinContent(14,104);
   histFail->SetBinContent(15,107);
   histFail->SetBinContent(16,93);
   histFail->SetBinContent(17,81);
   histFail->SetBinContent(18,100);
   histFail->SetBinContent(19,70);
   histFail->SetBinContent(20,65);
   histFail->SetBinContent(21,65);
   histFail->SetBinContent(22,49);
   histFail->SetBinContent(23,51);
   histFail->SetBinContent(24,47);
   histFail->SetBinContent(25,36);
   histFail->SetBinContent(26,35);
   histFail->SetBinContent(27,27);
   histFail->SetBinContent(28,30);
   histFail->SetBinContent(29,25);
   histFail->SetBinContent(30,21);
   histFail->SetBinContent(31,19);
   histFail->SetBinContent(32,23);
   histFail->SetBinContent(33,11);
   histFail->SetBinContent(34,14);
   histFail->SetBinContent(35,13);
   histFail->SetBinContent(36,14);
   histFail->SetBinContent(37,12);
   histFail->SetBinContent(38,9);
   histFail->SetBinContent(39,14);
   histFail->SetBinContent(40,6);
   histFail->SetEntries(1702);
   histFail->SetStats(0);
   histFail->SetFillColor(2);
   histFail->SetFillStyle(0);
   histFail->SetLineStyle(0);
   histFail->SetLineWidth(2);
   histFail->SetMarkerStyle(20);
   histFail->SetMarkerSize(1.2);
   histFail->GetXaxis()->SetNdivisions(505);
   histFail->GetXaxis()->SetLabelFont(42);
   histFail->GetXaxis()->SetLabelSize(0.05);
   histFail->GetXaxis()->SetTitleSize(0.055);
   histFail->GetXaxis()->SetTitleOffset(1.2);
   histFail->GetXaxis()->SetTitleFont(42);
   histFail->GetYaxis()->SetLabelFont(42);
   histFail->GetYaxis()->SetLabelOffset(0.01);
   histFail->GetYaxis()->SetLabelSize(0.05);
   histFail->GetYaxis()->SetTitleSize(0.055);
   histFail->GetYaxis()->SetTitleOffset(1.4);
   histFail->GetYaxis()->SetTitleFont(42);
   histFail->GetZaxis()->SetLabelFont(42);
   histFail->GetZaxis()->SetLabelSize(0.035);
   histFail->GetZaxis()->SetTitleSize(0.035);
   histFail->GetZaxis()->SetTitleFont(42);
   histFail->Draw("E");
   
   TH1D *h_fail_tot = new TH1D("h_fail_tot","",40,0,800);
   h_fail_tot->SetBinContent(4,0.8597657);
   h_fail_tot->SetBinContent(5,4.375461);
   h_fail_tot->SetBinContent(6,13.52181);
   h_fail_tot->SetBinContent(7,33.34119);
   h_fail_tot->SetBinContent(8,51.64651);
   h_fail_tot->SetBinContent(9,57.74548);
   h_fail_tot->SetBinContent(10,86.48542);
   h_fail_tot->SetBinContent(11,108.5158);
   h_fail_tot->SetBinContent(12,126.6499);
   h_fail_tot->SetBinContent(13,114.3985);
   h_fail_tot->SetBinContent(14,119.7653);
   h_fail_tot->SetBinContent(15,102.3827);
   h_fail_tot->SetBinContent(16,94.17219);
   h_fail_tot->SetBinContent(17,95.92078);
   h_fail_tot->SetBinContent(18,71.49468);
   h_fail_tot->SetBinContent(19,70.04584);
   h_fail_tot->SetBinContent(20,61.7069);
   h_fail_tot->SetBinContent(21,52.71534);
   h_fail_tot->SetBinContent(22,45.9019);
   h_fail_tot->SetBinContent(23,40.36973);
   h_fail_tot->SetBinContent(24,41.47618);
   h_fail_tot->SetBinContent(25,40.14997);
   h_fail_tot->SetBinContent(26,29.66679);
   h_fail_tot->SetBinContent(27,24.06112);
   h_fail_tot->SetBinContent(28,27.69686);
   h_fail_tot->SetBinContent(29,27.97494);
   h_fail_tot->SetBinContent(30,25.84192);
   h_fail_tot->SetBinContent(31,13.0971);
   h_fail_tot->SetBinContent(32,18.23663);
   h_fail_tot->SetBinContent(33,13.88696);
   h_fail_tot->SetBinContent(34,15.27068);
   h_fail_tot->SetBinContent(35,12.50362);
   h_fail_tot->SetBinContent(36,7.226198);
   h_fail_tot->SetBinContent(37,11.91174);
   h_fail_tot->SetBinContent(38,7.940075);
   h_fail_tot->SetBinContent(39,8.844361);
   h_fail_tot->SetBinContent(40,6.553836);
   h_fail_tot->SetEntries(81447);
   h_fail_tot->SetStats(0);
   h_fail_tot->SetFillColor(2);
   h_fail_tot->SetFillStyle(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   h_fail_tot->SetLineColor(ci);
   h_fail_tot->SetLineStyle(0);
   h_fail_tot->SetLineWidth(2);
   h_fail_tot->SetMarkerStyle(20);
   h_fail_tot->SetMarkerSize(1.2);
   h_fail_tot->GetXaxis()->SetNdivisions(505);
   h_fail_tot->GetXaxis()->SetLabelFont(42);
   h_fail_tot->GetXaxis()->SetLabelSize(0.05);
   h_fail_tot->GetXaxis()->SetTitleSize(0.055);
   h_fail_tot->GetXaxis()->SetTitleOffset(1.2);
   h_fail_tot->GetXaxis()->SetTitleFont(42);
   h_fail_tot->GetYaxis()->SetLabelFont(42);
   h_fail_tot->GetYaxis()->SetLabelOffset(0.01);
   h_fail_tot->GetYaxis()->SetLabelSize(0.05);
   h_fail_tot->GetYaxis()->SetTitleSize(0.055);
   h_fail_tot->GetYaxis()->SetTitleOffset(1.4);
   h_fail_tot->GetYaxis()->SetTitleFont(42);
   h_fail_tot->GetZaxis()->SetLabelFont(42);
   h_fail_tot->GetZaxis()->SetLabelSize(0.035);
   h_fail_tot->GetZaxis()->SetTitleSize(0.035);
   h_fail_tot->GetZaxis()->SetTitleFont(42);
   h_fail_tot->Draw("histsame");
   
   TH1D *h_fail_signal = new TH1D("h_fail_signal","",40,0,800);
   h_fail_signal->SetBinContent(5,0.07164714);
   h_fail_signal->SetBinContent(6,0.2865886);
   h_fail_signal->SetBinContent(7,0.8597657);
   h_fail_signal->SetBinContent(8,3.00918);
   h_fail_signal->SetBinContent(9,2.794239);
   h_fail_signal->SetBinContent(10,2.865886);
   h_fail_signal->SetBinContent(11,1.50459);
   h_fail_signal->SetBinContent(12,1.791179);
   h_fail_signal->SetBinContent(13,1.50459);
   h_fail_signal->SetBinContent(14,1.00306);
   h_fail_signal->SetBinContent(15,0.7164714);
   h_fail_signal->SetBinContent(16,0.2865886);
   h_fail_signal->SetBinContent(17,0.50153);
   h_fail_signal->SetBinContent(18,0.50153);
   h_fail_signal->SetBinContent(20,0.2149414);
   h_fail_signal->SetBinContent(21,0.07164714);
   h_fail_signal->SetBinContent(23,0.07164714);
   h_fail_signal->SetBinContent(24,0.2149414);
   h_fail_signal->SetBinContent(25,0.07164714);
   h_fail_signal->SetBinContent(26,-0.1432943);
   h_fail_signal->SetBinContent(27,0.1432943);
   h_fail_signal->SetBinContent(32,0.07164714);
   h_fail_signal->SetBinContent(37,-0.07164714);
   h_fail_signal->SetEntries(762);
   h_fail_signal->SetStats(0);
   h_fail_signal->SetFillColor(2);
   h_fail_signal->SetFillStyle(0);

   ci = TColor::GetColor("#00ff00");
   h_fail_signal->SetLineColor(ci);
   h_fail_signal->SetLineStyle(2);
   h_fail_signal->SetLineWidth(2);
   h_fail_signal->SetMarkerStyle(20);
   h_fail_signal->SetMarkerSize(1.2);
   h_fail_signal->GetXaxis()->SetNdivisions(505);
   h_fail_signal->GetXaxis()->SetLabelFont(42);
   h_fail_signal->GetXaxis()->SetLabelSize(0.05);
   h_fail_signal->GetXaxis()->SetTitleSize(0.055);
   h_fail_signal->GetXaxis()->SetTitleOffset(1.2);
   h_fail_signal->GetXaxis()->SetTitleFont(42);
   h_fail_signal->GetYaxis()->SetLabelFont(42);
   h_fail_signal->GetYaxis()->SetLabelOffset(0.01);
   h_fail_signal->GetYaxis()->SetLabelSize(0.05);
   h_fail_signal->GetYaxis()->SetTitleSize(0.055);
   h_fail_signal->GetYaxis()->SetTitleOffset(1.4);
   h_fail_signal->GetYaxis()->SetTitleFont(42);
   h_fail_signal->GetZaxis()->SetLabelFont(42);
   h_fail_signal->GetZaxis()->SetLabelSize(0.035);
   h_fail_signal->GetZaxis()->SetTitleSize(0.035);
   h_fail_signal->GetZaxis()->SetTitleFont(42);
   h_fail_signal->Draw("samehist");
   
   TH1D *h_fail_bkg1 = new TH1D("h_fail_bkg1","",40,0,800);
   h_fail_bkg1->SetBinContent(5,0.6448243);
   h_fail_bkg1->SetBinContent(6,4.51377);
   h_fail_bkg1->SetBinContent(7,11.75013);
   h_fail_bkg1->SetBinContent(8,20.41944);
   h_fail_bkg1->SetBinContent(9,24.43168);
   h_fail_bkg1->SetBinContent(10,35.96686);
   h_fail_bkg1->SetBinContent(11,53.73536);
   h_fail_bkg1->SetBinContent(12,60.54183);
   h_fail_bkg1->SetBinContent(13,56.31465);
   h_fail_bkg1->SetBinContent(14,53.30547);
   h_fail_bkg1->SetBinContent(15,46.78558);
   h_fail_bkg1->SetBinContent(16,36.25345);
   h_fail_bkg1->SetBinContent(17,41.19711);
   h_fail_bkg1->SetBinContent(18,33.10098);
   h_fail_bkg1->SetBinContent(19,32.5278);
   h_fail_bkg1->SetBinContent(20,26.93933);
   h_fail_bkg1->SetBinContent(21,23.07038);
   h_fail_bkg1->SetBinContent(22,20.34779);
   h_fail_bkg1->SetBinContent(23,15.90567);
   h_fail_bkg1->SetBinContent(24,14.90261);
   h_fail_bkg1->SetBinContent(25,15.83402);
   h_fail_bkg1->SetBinContent(26,13.18307);
   h_fail_bkg1->SetBinContent(27,9.099187);
   h_fail_bkg1->SetBinContent(28,13.97119);
   h_fail_bkg1->SetBinContent(29,11.82178);
   h_fail_bkg1->SetBinContent(30,10.24554);
   h_fail_bkg1->SetBinContent(31,4.083887);
   h_fail_bkg1->SetBinContent(32,4.585417);
   h_fail_bkg1->SetBinContent(33,6.161654);
   h_fail_bkg1->SetBinContent(34,7.236361);
   h_fail_bkg1->SetBinContent(35,4.083887);
   h_fail_bkg1->SetBinContent(36,3.654004);
   h_fail_bkg1->SetBinContent(37,5.445183);
   h_fail_bkg1->SetBinContent(38,2.077767);
   h_fail_bkg1->SetBinContent(39,1.647884);
   h_fail_bkg1->SetBinContent(40,3.080827);
   h_fail_bkg1->SetEntries(33381);
   h_fail_bkg1->SetStats(0);
   h_fail_bkg1->SetFillColor(2);
   h_fail_bkg1->SetFillStyle(0);

   ci = TColor::GetColor("#ff0000");
   h_fail_bkg1->SetLineColor(ci);
   h_fail_bkg1->SetLineStyle(2);
   h_fail_bkg1->SetLineWidth(2);
   h_fail_bkg1->SetMarkerStyle(20);
   h_fail_bkg1->SetMarkerSize(1.2);
   h_fail_bkg1->GetXaxis()->SetNdivisions(505);
   h_fail_bkg1->GetXaxis()->SetLabelFont(42);
   h_fail_bkg1->GetXaxis()->SetLabelSize(0.05);
   h_fail_bkg1->GetXaxis()->SetTitleSize(0.055);
   h_fail_bkg1->GetXaxis()->SetTitleOffset(1.2);
   h_fail_bkg1->GetXaxis()->SetTitleFont(42);
   h_fail_bkg1->GetYaxis()->SetLabelFont(42);
   h_fail_bkg1->GetYaxis()->SetLabelOffset(0.01);
   h_fail_bkg1->GetYaxis()->SetLabelSize(0.05);
   h_fail_bkg1->GetYaxis()->SetTitleSize(0.055);
   h_fail_bkg1->GetYaxis()->SetTitleOffset(1.4);
   h_fail_bkg1->GetYaxis()->SetTitleFont(42);
   h_fail_bkg1->GetZaxis()->SetLabelFont(42);
   h_fail_bkg1->GetZaxis()->SetLabelSize(0.035);
   h_fail_bkg1->GetZaxis()->SetTitleSize(0.035);
   h_fail_bkg1->GetZaxis()->SetTitleFont(42);
   h_fail_bkg1->Draw("samehist");
   
   TH1D *h_fail_bkg2 = new TH1D("h_fail_bkg2","",40,0,800);
   h_fail_bkg2->SetBinContent(4,0.8597657);
   h_fail_bkg2->SetBinContent(5,3.587343);
   h_fail_bkg2->SetBinContent(6,8.434864);
   h_fail_bkg2->SetBinContent(7,19.87153);
   h_fail_bkg2->SetBinContent(8,25.20872);
   h_fail_bkg2->SetBinContent(9,27.72533);
   h_fail_bkg2->SetBinContent(10,44.78678);
   h_fail_bkg2->SetBinContent(11,51.77129);
   h_fail_bkg2->SetBinContent(12,62.52572);
   h_fail_bkg2->SetBinContent(13,55.07471);
   h_fail_bkg2->SetBinContent(14,64.4537);
   h_fail_bkg2->SetBinContent(15,54.16418);
   h_fail_bkg2->SetBinContent(16,57.34556);
   h_fail_bkg2->SetBinContent(17,53.72062);
   h_fail_bkg2->SetBinContent(18,37.39065);
   h_fail_bkg2->SetBinContent(19,37.51804);
   h_fail_bkg2->SetBinContent(20,34.3377);
   h_fail_bkg2->SetBinContent(21,29.50166);
   h_fail_bkg2->SetBinContent(22,25.55411);
   h_fail_bkg2->SetBinContent(23,24.32077);
   h_fail_bkg2->SetBinContent(24,26.14369);
   h_fail_bkg2->SetBinContent(25,24.17265);
   h_fail_bkg2->SetBinContent(26,16.7703);
   h_fail_bkg2->SetBinContent(27,14.67535);
   h_fail_bkg2->SetBinContent(28,13.72567);
   h_fail_bkg2->SetBinContent(29,16.15316);
   h_fail_bkg2->SetBinContent(30,15.59637);
   h_fail_bkg2->SetBinContent(31,9.013215);
   h_fail_bkg2->SetBinContent(32,13.50792);
   h_fail_bkg2->SetBinContent(33,7.725305);
   h_fail_bkg2->SetBinContent(34,8.034324);
   h_fail_bkg2->SetBinContent(35,8.419729);
   h_fail_bkg2->SetBinContent(36,3.572194);
   h_fail_bkg2->SetBinContent(37,6.609855);
   h_fail_bkg2->SetBinContent(38,5.862308);
   h_fail_bkg2->SetBinContent(39,7.196476);
   h_fail_bkg2->SetBinContent(40,3.473008);
   h_fail_bkg2->SetEntries(46542);
   h_fail_bkg2->SetStats(0);
   h_fail_bkg2->SetFillColor(2);
   h_fail_bkg2->SetFillStyle(0);

   ci = TColor::GetColor("#cccccc");
   h_fail_bkg2->SetLineColor(ci);
   h_fail_bkg2->SetLineStyle(2);
   h_fail_bkg2->SetLineWidth(2);
   h_fail_bkg2->SetMarkerStyle(20);
   h_fail_bkg2->SetMarkerSize(1.2);
   h_fail_bkg2->GetXaxis()->SetNdivisions(505);
   h_fail_bkg2->GetXaxis()->SetLabelFont(42);
   h_fail_bkg2->GetXaxis()->SetLabelSize(0.05);
   h_fail_bkg2->GetXaxis()->SetTitleSize(0.055);
   h_fail_bkg2->GetXaxis()->SetTitleOffset(1.2);
   h_fail_bkg2->GetXaxis()->SetTitleFont(42);
   h_fail_bkg2->GetYaxis()->SetLabelFont(42);
   h_fail_bkg2->GetYaxis()->SetLabelOffset(0.01);
   h_fail_bkg2->GetYaxis()->SetLabelSize(0.05);
   h_fail_bkg2->GetYaxis()->SetTitleSize(0.055);
   h_fail_bkg2->GetYaxis()->SetTitleOffset(1.4);
   h_fail_bkg2->GetYaxis()->SetTitleFont(42);
   h_fail_bkg2->GetZaxis()->SetLabelFont(42);
   h_fail_bkg2->GetZaxis()->SetLabelSize(0.035);
   h_fail_bkg2->GetZaxis()->SetTitleSize(0.035);
   h_fail_bkg2->GetZaxis()->SetTitleFont(42);
   h_fail_bkg2->Draw("samehist");
   
   TLegend *leg = new TLegend(0.65,0.54,0.99,0.99,NULL,"brNDC");
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("histPass","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("h_pass_tot","TOT MC","f");
   entry->SetFillColor(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_pass_signal","Signal","f");
   entry->SetFillColor(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_pass_bkg1","tt1l Combinatorial BG","f");
   entry->SetFillColor(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h_pass_bkg2","Non-tt1l BG","f");
   entry->SetFillColor(2);

   ci = TColor::GetColor("#cccccc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   cfail_bf_mt->Modified();
   cfail_bf_mt->cd();
   cfail_bf_mt->SetSelected(cfail_bf_mt);
}
