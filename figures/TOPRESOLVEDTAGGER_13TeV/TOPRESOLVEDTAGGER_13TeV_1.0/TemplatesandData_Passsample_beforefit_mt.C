{
//=========Macro generated from canvas: cpass_bf_mt/cpass_bf_mt
//=========  (Sat Oct 31 02:03:27 2015) by ROOT version5.32/00
   TCanvas *cpass_bf_mt = new TCanvas("cpass_bf_mt", "cpass_bf_mt",1,1,800,1176);
   gStyle->SetOptStat(0);
   cpass_bf_mt->SetHighLightColor(2);
   cpass_bf_mt->Range(-150,-0.1727848,850,1.156329);
   cpass_bf_mt->SetFillColor(0);
   cpass_bf_mt->SetBorderMode(0);
   cpass_bf_mt->SetBorderSize(10);
   cpass_bf_mt->SetTickx(1);
   cpass_bf_mt->SetTicky(1);
   cpass_bf_mt->SetLeftMargin(0.15);
   cpass_bf_mt->SetRightMargin(0.05);
   cpass_bf_mt->SetTopMargin(0.08);
   cpass_bf_mt->SetBottomMargin(0.13);
   cpass_bf_mt->SetFrameFillStyle(0);
   cpass_bf_mt->SetFrameLineStyle(0);
   cpass_bf_mt->SetFrameBorderMode(0);
   cpass_bf_mt->SetFrameBorderSize(10);
   cpass_bf_mt->SetFrameFillStyle(0);
   cpass_bf_mt->SetFrameLineStyle(0);
   cpass_bf_mt->SetFrameBorderMode(0);
   cpass_bf_mt->SetFrameBorderSize(10);
   
   TH1D *histPass = new TH1D("histPass","",40,0,800);
   histPass->SetStats(0);
   histPass->SetFillColor(2);
   histPass->SetFillStyle(0);
   histPass->SetLineStyle(0);
   histPass->SetLineWidth(2);
   histPass->SetMarkerStyle(20);
   histPass->SetMarkerSize(1.2);
   histPass->GetXaxis()->SetNdivisions(505);
   histPass->GetXaxis()->SetLabelFont(42);
   histPass->GetXaxis()->SetLabelSize(0.05);
   histPass->GetXaxis()->SetTitleSize(0.055);
   histPass->GetXaxis()->SetTitleOffset(1.2);
   histPass->GetXaxis()->SetTitleFont(42);
   histPass->GetYaxis()->SetLabelFont(42);
   histPass->GetYaxis()->SetLabelOffset(0.01);
   histPass->GetYaxis()->SetLabelSize(0.05);
   histPass->GetYaxis()->SetTitleSize(0.055);
   histPass->GetYaxis()->SetTitleOffset(1.4);
   histPass->GetYaxis()->SetTitleFont(42);
   histPass->GetZaxis()->SetLabelFont(42);
   histPass->GetZaxis()->SetLabelSize(0.035);
   histPass->GetZaxis()->SetTitleSize(0.035);
   histPass->GetZaxis()->SetTitleFont(42);
   histPass->Draw("E");
   
   TH1D *h_pass_tot = new TH1D("h_pass_tot","",40,0,800);
   h_pass_tot->SetStats(0);
   h_pass_tot->SetFillColor(2);
   h_pass_tot->SetFillStyle(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   h_pass_tot->SetLineColor(ci);
   h_pass_tot->SetLineStyle(0);
   h_pass_tot->SetLineWidth(2);
   h_pass_tot->SetMarkerStyle(20);
   h_pass_tot->SetMarkerSize(1.2);
   h_pass_tot->GetXaxis()->SetNdivisions(505);
   h_pass_tot->GetXaxis()->SetLabelFont(42);
   h_pass_tot->GetXaxis()->SetLabelSize(0.05);
   h_pass_tot->GetXaxis()->SetTitleSize(0.055);
   h_pass_tot->GetXaxis()->SetTitleOffset(1.2);
   h_pass_tot->GetXaxis()->SetTitleFont(42);
   h_pass_tot->GetYaxis()->SetLabelFont(42);
   h_pass_tot->GetYaxis()->SetLabelOffset(0.01);
   h_pass_tot->GetYaxis()->SetLabelSize(0.05);
   h_pass_tot->GetYaxis()->SetTitleSize(0.055);
   h_pass_tot->GetYaxis()->SetTitleOffset(1.4);
   h_pass_tot->GetYaxis()->SetTitleFont(42);
   h_pass_tot->GetZaxis()->SetLabelFont(42);
   h_pass_tot->GetZaxis()->SetLabelSize(0.035);
   h_pass_tot->GetZaxis()->SetTitleSize(0.035);
   h_pass_tot->GetZaxis()->SetTitleFont(42);
   h_pass_tot->Draw("histsame");
   
   TH1D *h_pass_signal = new TH1D("h_pass_signal","",40,0,800);
   h_pass_signal->SetStats(0);
   h_pass_signal->SetFillColor(2);
   h_pass_signal->SetFillStyle(0);

   ci = TColor::GetColor("#00ff00");
   h_pass_signal->SetLineColor(ci);
   h_pass_signal->SetLineStyle(2);
   h_pass_signal->SetLineWidth(2);
   h_pass_signal->SetMarkerStyle(20);
   h_pass_signal->SetMarkerSize(1.2);
   h_pass_signal->GetXaxis()->SetNdivisions(505);
   h_pass_signal->GetXaxis()->SetLabelFont(42);
   h_pass_signal->GetXaxis()->SetLabelSize(0.05);
   h_pass_signal->GetXaxis()->SetTitleSize(0.055);
   h_pass_signal->GetXaxis()->SetTitleOffset(1.2);
   h_pass_signal->GetXaxis()->SetTitleFont(42);
   h_pass_signal->GetYaxis()->SetLabelFont(42);
   h_pass_signal->GetYaxis()->SetLabelOffset(0.01);
   h_pass_signal->GetYaxis()->SetLabelSize(0.05);
   h_pass_signal->GetYaxis()->SetTitleSize(0.055);
   h_pass_signal->GetYaxis()->SetTitleOffset(1.4);
   h_pass_signal->GetYaxis()->SetTitleFont(42);
   h_pass_signal->GetZaxis()->SetLabelFont(42);
   h_pass_signal->GetZaxis()->SetLabelSize(0.035);
   h_pass_signal->GetZaxis()->SetTitleSize(0.035);
   h_pass_signal->GetZaxis()->SetTitleFont(42);
   h_pass_signal->Draw("samehist");
   
   TH1D *h_pass_bkg1 = new TH1D("h_pass_bkg1","",40,0,800);
   h_pass_bkg1->SetStats(0);
   h_pass_bkg1->SetFillColor(2);
   h_pass_bkg1->SetFillStyle(0);

   ci = TColor::GetColor("#ff0000");
   h_pass_bkg1->SetLineColor(ci);
   h_pass_bkg1->SetLineStyle(2);
   h_pass_bkg1->SetLineWidth(2);
   h_pass_bkg1->SetMarkerStyle(20);
   h_pass_bkg1->SetMarkerSize(1.2);
   h_pass_bkg1->GetXaxis()->SetNdivisions(505);
   h_pass_bkg1->GetXaxis()->SetLabelFont(42);
   h_pass_bkg1->GetXaxis()->SetLabelSize(0.05);
   h_pass_bkg1->GetXaxis()->SetTitleSize(0.055);
   h_pass_bkg1->GetXaxis()->SetTitleOffset(1.2);
   h_pass_bkg1->GetXaxis()->SetTitleFont(42);
   h_pass_bkg1->GetYaxis()->SetLabelFont(42);
   h_pass_bkg1->GetYaxis()->SetLabelOffset(0.01);
   h_pass_bkg1->GetYaxis()->SetLabelSize(0.05);
   h_pass_bkg1->GetYaxis()->SetTitleSize(0.055);
   h_pass_bkg1->GetYaxis()->SetTitleOffset(1.4);
   h_pass_bkg1->GetYaxis()->SetTitleFont(42);
   h_pass_bkg1->GetZaxis()->SetLabelFont(42);
   h_pass_bkg1->GetZaxis()->SetLabelSize(0.035);
   h_pass_bkg1->GetZaxis()->SetTitleSize(0.035);
   h_pass_bkg1->GetZaxis()->SetTitleFont(42);
   h_pass_bkg1->Draw("samehist");
   
   TH1D *h_pass_bkg2 = new TH1D("h_pass_bkg2","",40,0,800);
   h_pass_bkg2->SetStats(0);
   h_pass_bkg2->SetFillColor(2);
   h_pass_bkg2->SetFillStyle(0);

   ci = TColor::GetColor("#cccccc");
   h_pass_bkg2->SetLineColor(ci);
   h_pass_bkg2->SetLineStyle(2);
   h_pass_bkg2->SetLineWidth(2);
   h_pass_bkg2->SetMarkerStyle(20);
   h_pass_bkg2->SetMarkerSize(1.2);
   h_pass_bkg2->GetXaxis()->SetNdivisions(505);
   h_pass_bkg2->GetXaxis()->SetLabelFont(42);
   h_pass_bkg2->GetXaxis()->SetLabelSize(0.05);
   h_pass_bkg2->GetXaxis()->SetTitleSize(0.055);
   h_pass_bkg2->GetXaxis()->SetTitleOffset(1.2);
   h_pass_bkg2->GetXaxis()->SetTitleFont(42);
   h_pass_bkg2->GetYaxis()->SetLabelFont(42);
   h_pass_bkg2->GetYaxis()->SetLabelOffset(0.01);
   h_pass_bkg2->GetYaxis()->SetLabelSize(0.05);
   h_pass_bkg2->GetYaxis()->SetTitleSize(0.055);
   h_pass_bkg2->GetYaxis()->SetTitleOffset(1.4);
   h_pass_bkg2->GetYaxis()->SetTitleFont(42);
   h_pass_bkg2->GetZaxis()->SetLabelFont(42);
   h_pass_bkg2->GetZaxis()->SetLabelSize(0.035);
   h_pass_bkg2->GetZaxis()->SetTitleSize(0.035);
   h_pass_bkg2->GetZaxis()->SetTitleFont(42);
   h_pass_bkg2->Draw("samehist");
   
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
   cpass_bf_mt->Modified();
   cpass_bf_mt->cd();
   cpass_bf_mt->SetSelected(cpass_bf_mt);
}
