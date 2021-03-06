{
//=========Macro generated from canvas: cfail_bf/cfail_bf
//=========  (Tue Oct 27 00:00:46 2015) by ROOT version5.32/00
   TCanvas *cfail_bf = new TCanvas("cfail_bf", "cfail_bf",1,1,800,1176);
   gStyle->SetOptStat(0);
   cfail_bf->SetHighLightColor(2);
   cfail_bf->Range(-150,-228.6596,850,1530.261);
   cfail_bf->SetFillColor(0);
   cfail_bf->SetBorderMode(0);
   cfail_bf->SetBorderSize(10);
   cfail_bf->SetTickx(1);
   cfail_bf->SetTicky(1);
   cfail_bf->SetLeftMargin(0.15);
   cfail_bf->SetRightMargin(0.05);
   cfail_bf->SetTopMargin(0.08);
   cfail_bf->SetBottomMargin(0.13);
   cfail_bf->SetFrameFillStyle(0);
   cfail_bf->SetFrameLineStyle(0);
   cfail_bf->SetFrameBorderMode(0);
   cfail_bf->SetFrameBorderSize(10);
   cfail_bf->SetFrameFillStyle(0);
   cfail_bf->SetFrameLineStyle(0);
   cfail_bf->SetFrameBorderMode(0);
   cfail_bf->SetFrameBorderSize(10);
   
   TH1D *frame_179c990__3 = new TH1D("frame_179c990__3","",40,0,800);
   frame_179c990__3->SetMaximum(1389.547);
   frame_179c990__3->SetDirectory(0);
   frame_179c990__3->SetStats(0);
   frame_179c990__3->SetFillColor(2);
   frame_179c990__3->SetFillStyle(0);
   frame_179c990__3->SetLineStyle(0);
   frame_179c990__3->SetLineWidth(2);
   frame_179c990__3->SetMarkerStyle(20);
   frame_179c990__3->SetMarkerSize(1.2);
   frame_179c990__3->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_179c990__3->GetXaxis()->SetNdivisions(505);
   frame_179c990__3->GetXaxis()->SetLabelFont(42);
   frame_179c990__3->GetXaxis()->SetLabelSize(0.05);
   frame_179c990__3->GetXaxis()->SetTitleSize(0.055);
   frame_179c990__3->GetXaxis()->SetTitleOffset(1.2);
   frame_179c990__3->GetXaxis()->SetTitleFont(42);
   frame_179c990__3->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_179c990__3->GetYaxis()->SetLabelFont(42);
   frame_179c990__3->GetYaxis()->SetLabelOffset(0.01);
   frame_179c990__3->GetYaxis()->SetLabelSize(0.05);
   frame_179c990__3->GetYaxis()->SetTitleSize(0.055);
   frame_179c990__3->GetYaxis()->SetTitleOffset(1.4);
   frame_179c990__3->GetYaxis()->SetTitleFont(42);
   frame_179c990__3->GetZaxis()->SetLabelFont(42);
   frame_179c990__3->GetZaxis()->SetLabelSize(0.035);
   frame_179c990__3->GetZaxis()->SetTitleSize(0.035);
   frame_179c990__3->GetZaxis()->SetTitleFont(42);
   frame_179c990__3->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(40);
   grae->SetName("h_dataFail");
   grae->SetTitle("Histogram of dataFail_plot__m");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   grae->SetPoint(0,10,0);
   grae->SetPointError(0,10,10,0,1.147874);
   grae->SetPoint(1,30,0);
   grae->SetPointError(1,10,10,0,1.147874);
   grae->SetPoint(2,50,1);
   grae->SetPointError(2,10,10,0.8272462,2.299527);
   grae->SetPoint(3,70,12);
   grae->SetPointError(3,10,10,3.415266,4.559819);
   grae->SetPoint(4,90,80);
   grae->SetPointError(4,10,10,8.925539,9.981567);
   grae->SetPoint(5,110,210);
   grae->SetPointError(5,10,10,14,15);
   grae->SetPoint(6,130,412);
   grae->SetPointError(6,10,10,19.80394,20.80394);
   grae->SetPoint(7,150,703);
   grae->SetPointError(7,10,10,26.01886,27.01886);
   grae->SetPoint(8,170,1166);
   grae->SetPointError(8,10,10,33.6504,34.6504);
   grae->SetPoint(9,190,1287);
   grae->SetPointError(9,10,10,35.37827,36.37827);
   grae->SetPoint(10,210,985);
   grae->SetPointError(10,10,10,30.88869,31.88869);
   grae->SetPoint(11,230,794);
   grae->SetPointError(11,10,10,27.68244,28.68244);
   grae->SetPoint(12,250,608);
   grae->SetPointError(12,10,10,24.16272,25.16272);
   grae->SetPoint(13,270,540);
   grae->SetPointError(13,10,10,22.74328,23.74328);
   grae->SetPoint(14,290,463);
   grae->SetPointError(14,10,10,21.02324,22.02324);
   grae->SetPoint(15,310,414);
   grae->SetPointError(15,10,10,19.85313,20.85313);
   grae->SetPoint(16,330,343);
   grae->SetPointError(16,10,10,18.02701,19.02701);
   grae->SetPoint(17,350,307);
   grae->SetPointError(17,10,10,17.02855,18.02855);
   grae->SetPoint(18,370,272);
   grae->SetPointError(18,10,10,16,17);
   grae->SetPoint(19,390,235);
   grae->SetPointError(19,10,10,14.83786,15.83786);
   grae->SetPoint(20,410,213);
   grae->SetPointError(20,10,10,14.10308,15.10308);
   grae->SetPoint(21,430,191);
   grae->SetPointError(21,10,10,13.32932,14.32932);
   grae->SetPoint(22,450,162);
   grae->SetPointError(22,10,10,12.23774,13.23774);
   grae->SetPoint(23,470,159);
   grae->SetPointError(23,10,10,12.11943,13.11943);
   grae->SetPoint(24,490,135);
   grae->SetPointError(24,10,10,11.1297,12.1297);
   grae->SetPoint(25,510,108);
   grae->SetPointError(25,10,10,9.904326,10.90433);
   grae->SetPoint(26,530,124);
   grae->SetPointError(26,10,10,10.64675,11.64675);
   grae->SetPoint(27,550,93);
   grae->SetPointError(27,10,10,9.626284,10.67824);
   grae->SetPoint(28,570,70);
   grae->SetPointError(28,10,10,8.346566,9.406468);
   grae->SetPoint(29,590,71);
   grae->SetPointError(29,10,10,8.406258,9.465736);
   grae->SetPoint(30,610,75);
   grae->SetPointError(30,10,10,8.640903,9.698771);
   grae->SetPoint(31,630,48);
   grae->SetPointError(31,10,10,6.903979,7.97633);
   grae->SetPoint(32,650,69);
   grae->SetPointError(32,10,10,8.286444,9.346779);
   grae->SetPoint(33,670,43);
   grae->SetPointError(33,10,10,6.531834,7.608278);
   grae->SetPoint(34,690,36);
   grae->SetPointError(34,10,10,5.971996,7.055545);
   grae->SetPoint(35,710,28);
   grae->SetPointError(35,10,10,5.259711,6.354446);
   grae->SetPoint(36,730,35);
   grae->SetPointError(36,10,10,5.887675,6.97241);
   grae->SetPoint(37,750,35);
   grae->SetPointError(37,10,10,5.887675,6.97241);
   grae->SetPoint(38,770,22);
   grae->SetPointError(38,10,10,4.654502,5.761366);
   grae->SetPoint(39,790,25);
   grae->SetPointError(39,10,10,4.966335,6.066589);
   
   TH1F *Graph_h_dataFail2 = new TH1F("Graph_h_dataFail2","Histogram of dataFail_plot__m",100,0,880);
   Graph_h_dataFail2->SetMinimum(0);
   Graph_h_dataFail2->SetMaximum(1455.716);
   Graph_h_dataFail2->SetDirectory(0);
   Graph_h_dataFail2->SetStats(0);
   Graph_h_dataFail2->SetFillColor(2);
   Graph_h_dataFail2->SetFillStyle(0);
   Graph_h_dataFail2->SetLineStyle(0);
   Graph_h_dataFail2->SetLineWidth(2);
   Graph_h_dataFail2->SetMarkerStyle(20);
   Graph_h_dataFail2->SetMarkerSize(1.2);
   Graph_h_dataFail2->GetXaxis()->SetNdivisions(505);
   Graph_h_dataFail2->GetXaxis()->SetLabelFont(42);
   Graph_h_dataFail2->GetXaxis()->SetLabelSize(0.05);
   Graph_h_dataFail2->GetXaxis()->SetTitleSize(0.055);
   Graph_h_dataFail2->GetXaxis()->SetTitleOffset(1.2);
   Graph_h_dataFail2->GetXaxis()->SetTitleFont(42);
   Graph_h_dataFail2->GetYaxis()->SetLabelFont(42);
   Graph_h_dataFail2->GetYaxis()->SetLabelOffset(0.01);
   Graph_h_dataFail2->GetYaxis()->SetLabelSize(0.05);
   Graph_h_dataFail2->GetYaxis()->SetTitleSize(0.055);
   Graph_h_dataFail2->GetYaxis()->SetTitleOffset(1.4);
   Graph_h_dataFail2->GetYaxis()->SetTitleFont(42);
   Graph_h_dataFail2->GetZaxis()->SetLabelFont(42);
   Graph_h_dataFail2->GetZaxis()->SetLabelSize(0.035);
   Graph_h_dataFail2->GetZaxis()->SetTitleSize(0.035);
   Graph_h_dataFail2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_h_dataFail2);
   
   grae->Draw("zp");
   
   TGraph *graph = new TGraph(87);
   graph->SetName("modelFail_Norm[m]");
   graph->SetTitle("Projection of Model for FAIL sample");
   graph->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,-9.87654321,0);
   graph->SetPoint(1,-9.87654321,0);
   graph->SetPoint(2,0,0);
   graph->SetPoint(3,8.1608e-06,0);
   graph->SetPoint(4,19.99999184,0);
   graph->SetPoint(5,20.00000816,0);
   graph->SetPoint(6,39.99999184,0);
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0);
   graph->SetPoint(10,79.99999184,0);
   graph->SetPoint(11,80.00000816,0);
   graph->SetPoint(12,99.99999184,0);
   graph->SetPoint(13,100.0000082,0);
   graph->SetPoint(14,119.9999918,0);
   graph->SetPoint(15,120.0000082,0);
   graph->SetPoint(16,139.9999918,0);
   graph->SetPoint(17,140.0000082,0);
   graph->SetPoint(18,159.9999918,0);
   graph->SetPoint(19,160.0000082,0);
   graph->SetPoint(20,179.9999918,0);
   graph->SetPoint(21,180.0000082,0);
   graph->SetPoint(22,199.9999918,0);
   graph->SetPoint(23,200.0000082,0);
   graph->SetPoint(24,219.9999918,0);
   graph->SetPoint(25,220.0000082,0);
   graph->SetPoint(26,239.9999918,0);
   graph->SetPoint(27,240.0000082,0);
   graph->SetPoint(28,259.9999918,0);
   graph->SetPoint(29,260.0000082,0);
   graph->SetPoint(30,279.9999918,0);
   graph->SetPoint(31,280.0000082,0);
   graph->SetPoint(32,299.9999918,0);
   graph->SetPoint(33,300.0000082,0);
   graph->SetPoint(34,319.9999918,0);
   graph->SetPoint(35,320.0000082,0);
   graph->SetPoint(36,339.9999918,0);
   graph->SetPoint(37,340.0000082,0);
   graph->SetPoint(38,359.9999918,0);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,0);
   graph->SetPoint(42,399.9999918,0);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0);
   graph->SetPoint(52,499.9999918,0);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,0);
   graph->SetPoint(56,539.9999918,0);
   graph->SetPoint(57,540.0000082,0);
   graph->SetPoint(58,559.9999918,0);
   graph->SetPoint(59,560.0000082,0);
   graph->SetPoint(60,579.9999918,0);
   graph->SetPoint(61,580.0000082,0);
   graph->SetPoint(62,599.9999918,0);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0);
   graph->SetPoint(66,639.9999918,0);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0);
   graph->SetPoint(74,719.9999918,0);
   graph->SetPoint(75,720.0000082,0);
   graph->SetPoint(76,739.9999918,0);
   graph->SetPoint(77,740.0000082,0);
   graph->SetPoint(78,759.9999918,0);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]5 = new TH1F("Graph_modelFail_Norm[m]5","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]5->SetMinimum(0);
   Graph_modelFail_Norm[m]5->SetMaximum(1.1);
   Graph_modelFail_Norm[m]5->SetDirectory(0);
   Graph_modelFail_Norm[m]5->SetStats(0);
   Graph_modelFail_Norm[m]5->SetFillColor(2);
   Graph_modelFail_Norm[m]5->SetFillStyle(0);
   Graph_modelFail_Norm[m]5->SetLineStyle(0);
   Graph_modelFail_Norm[m]5->SetLineWidth(2);
   Graph_modelFail_Norm[m]5->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]5->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]5->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]5->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]5->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]5->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]5->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]5->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]5);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelFail_Norm[m]_Comp[sigModFail]");
   graph->SetTitle("Projection of Model for FAIL sample");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,-9.87654321,0);
   graph->SetPoint(1,-9.87654321,0);
   graph->SetPoint(2,0,0);
   graph->SetPoint(3,8.1608e-06,0);
   graph->SetPoint(4,19.99999184,0);
   graph->SetPoint(5,20.00000816,0);
   graph->SetPoint(6,39.99999184,0);
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0);
   graph->SetPoint(10,79.99999184,0);
   graph->SetPoint(11,80.00000816,0);
   graph->SetPoint(12,99.99999184,0);
   graph->SetPoint(13,100.0000082,0);
   graph->SetPoint(14,119.9999918,0);
   graph->SetPoint(15,120.0000082,0);
   graph->SetPoint(16,139.9999918,0);
   graph->SetPoint(17,140.0000082,0);
   graph->SetPoint(18,159.9999918,0);
   graph->SetPoint(19,160.0000082,0);
   graph->SetPoint(20,179.9999918,0);
   graph->SetPoint(21,180.0000082,0);
   graph->SetPoint(22,199.9999918,0);
   graph->SetPoint(23,200.0000082,0);
   graph->SetPoint(24,219.9999918,0);
   graph->SetPoint(25,220.0000082,0);
   graph->SetPoint(26,239.9999918,0);
   graph->SetPoint(27,240.0000082,0);
   graph->SetPoint(28,259.9999918,0);
   graph->SetPoint(29,260.0000082,0);
   graph->SetPoint(30,279.9999918,0);
   graph->SetPoint(31,280.0000082,0);
   graph->SetPoint(32,299.9999918,0);
   graph->SetPoint(33,300.0000082,0);
   graph->SetPoint(34,319.9999918,0);
   graph->SetPoint(35,320.0000082,0);
   graph->SetPoint(36,339.9999918,0);
   graph->SetPoint(37,340.0000082,0);
   graph->SetPoint(38,359.9999918,0);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,0);
   graph->SetPoint(42,399.9999918,0);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0);
   graph->SetPoint(52,499.9999918,0);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,0);
   graph->SetPoint(56,539.9999918,0);
   graph->SetPoint(57,540.0000082,0);
   graph->SetPoint(58,559.9999918,0);
   graph->SetPoint(59,560.0000082,0);
   graph->SetPoint(60,579.9999918,0);
   graph->SetPoint(61,580.0000082,0);
   graph->SetPoint(62,599.9999918,0);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0);
   graph->SetPoint(66,639.9999918,0);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0);
   graph->SetPoint(74,719.9999918,0);
   graph->SetPoint(75,720.0000082,0);
   graph->SetPoint(76,739.9999918,0);
   graph->SetPoint(77,740.0000082,0);
   graph->SetPoint(78,759.9999918,0);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[sigModFail]6 = new TH1F("Graph_modelFail_Norm[m]_Comp[sigModFail]6","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetMaximum(1.1);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[sigModFail]6->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[sigModFail]6);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelFail_Norm[m]_Comp[bkg1ModFail]");
   graph->SetTitle("Projection of Model for FAIL sample");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,-9.87654321,0);
   graph->SetPoint(1,-9.87654321,0);
   graph->SetPoint(2,0,0);
   graph->SetPoint(3,8.1608e-06,0);
   graph->SetPoint(4,19.99999184,0);
   graph->SetPoint(5,20.00000816,0);
   graph->SetPoint(6,39.99999184,0);
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0);
   graph->SetPoint(10,79.99999184,0);
   graph->SetPoint(11,80.00000816,0);
   graph->SetPoint(12,99.99999184,0);
   graph->SetPoint(13,100.0000082,0);
   graph->SetPoint(14,119.9999918,0);
   graph->SetPoint(15,120.0000082,0);
   graph->SetPoint(16,139.9999918,0);
   graph->SetPoint(17,140.0000082,0);
   graph->SetPoint(18,159.9999918,0);
   graph->SetPoint(19,160.0000082,0);
   graph->SetPoint(20,179.9999918,0);
   graph->SetPoint(21,180.0000082,0);
   graph->SetPoint(22,199.9999918,0);
   graph->SetPoint(23,200.0000082,0);
   graph->SetPoint(24,219.9999918,0);
   graph->SetPoint(25,220.0000082,0);
   graph->SetPoint(26,239.9999918,0);
   graph->SetPoint(27,240.0000082,0);
   graph->SetPoint(28,259.9999918,0);
   graph->SetPoint(29,260.0000082,0);
   graph->SetPoint(30,279.9999918,0);
   graph->SetPoint(31,280.0000082,0);
   graph->SetPoint(32,299.9999918,0);
   graph->SetPoint(33,300.0000082,0);
   graph->SetPoint(34,319.9999918,0);
   graph->SetPoint(35,320.0000082,0);
   graph->SetPoint(36,339.9999918,0);
   graph->SetPoint(37,340.0000082,0);
   graph->SetPoint(38,359.9999918,0);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,0);
   graph->SetPoint(42,399.9999918,0);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0);
   graph->SetPoint(52,499.9999918,0);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,0);
   graph->SetPoint(56,539.9999918,0);
   graph->SetPoint(57,540.0000082,0);
   graph->SetPoint(58,559.9999918,0);
   graph->SetPoint(59,560.0000082,0);
   graph->SetPoint(60,579.9999918,0);
   graph->SetPoint(61,580.0000082,0);
   graph->SetPoint(62,599.9999918,0);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0);
   graph->SetPoint(66,639.9999918,0);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0);
   graph->SetPoint(74,719.9999918,0);
   graph->SetPoint(75,720.0000082,0);
   graph->SetPoint(76,739.9999918,0);
   graph->SetPoint(77,740.0000082,0);
   graph->SetPoint(78,759.9999918,0);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetMaximum(1.1);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[bkg1ModFail]7);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelFail_Norm[m]_Comp[bkg2ModFail]");
   graph->SetTitle("Projection of Model for FAIL sample");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#cccccc");
   graph->SetLineColor(ci);
   graph->SetLineStyle(2);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->SetPoint(0,-9.87654321,0);
   graph->SetPoint(1,-9.87654321,0);
   graph->SetPoint(2,0,0);
   graph->SetPoint(3,8.1608e-06,0);
   graph->SetPoint(4,19.99999184,0);
   graph->SetPoint(5,20.00000816,0);
   graph->SetPoint(6,39.99999184,0);
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0);
   graph->SetPoint(10,79.99999184,0);
   graph->SetPoint(11,80.00000816,0);
   graph->SetPoint(12,99.99999184,0);
   graph->SetPoint(13,100.0000082,0);
   graph->SetPoint(14,119.9999918,0);
   graph->SetPoint(15,120.0000082,0);
   graph->SetPoint(16,139.9999918,0);
   graph->SetPoint(17,140.0000082,0);
   graph->SetPoint(18,159.9999918,0);
   graph->SetPoint(19,160.0000082,0);
   graph->SetPoint(20,179.9999918,0);
   graph->SetPoint(21,180.0000082,0);
   graph->SetPoint(22,199.9999918,0);
   graph->SetPoint(23,200.0000082,0);
   graph->SetPoint(24,219.9999918,0);
   graph->SetPoint(25,220.0000082,0);
   graph->SetPoint(26,239.9999918,0);
   graph->SetPoint(27,240.0000082,0);
   graph->SetPoint(28,259.9999918,0);
   graph->SetPoint(29,260.0000082,0);
   graph->SetPoint(30,279.9999918,0);
   graph->SetPoint(31,280.0000082,0);
   graph->SetPoint(32,299.9999918,0);
   graph->SetPoint(33,300.0000082,0);
   graph->SetPoint(34,319.9999918,0);
   graph->SetPoint(35,320.0000082,0);
   graph->SetPoint(36,339.9999918,0);
   graph->SetPoint(37,340.0000082,0);
   graph->SetPoint(38,359.9999918,0);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,0);
   graph->SetPoint(42,399.9999918,0);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0);
   graph->SetPoint(52,499.9999918,0);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,0);
   graph->SetPoint(56,539.9999918,0);
   graph->SetPoint(57,540.0000082,0);
   graph->SetPoint(58,559.9999918,0);
   graph->SetPoint(59,560.0000082,0);
   graph->SetPoint(60,579.9999918,0);
   graph->SetPoint(61,580.0000082,0);
   graph->SetPoint(62,599.9999918,0);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0);
   graph->SetPoint(66,639.9999918,0);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0);
   graph->SetPoint(74,719.9999918,0);
   graph->SetPoint(75,720.0000082,0);
   graph->SetPoint(76,739.9999918,0);
   graph->SetPoint(77,740.0000082,0);
   graph->SetPoint(78,759.9999918,0);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetMaximum(1.1);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[bkg2ModFail]8);
   
   graph->Draw("l");
   
   TH1D *frame_179c990__4 = new TH1D("frame_179c990__4","",40,0,800);
   frame_179c990__4->SetMaximum(1389.547);
   frame_179c990__4->SetDirectory(0);
   frame_179c990__4->SetStats(0);
   frame_179c990__4->SetFillColor(2);
   frame_179c990__4->SetFillStyle(0);
   frame_179c990__4->SetLineStyle(0);
   frame_179c990__4->SetLineWidth(2);
   frame_179c990__4->SetMarkerStyle(20);
   frame_179c990__4->SetMarkerSize(1.2);
   frame_179c990__4->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_179c990__4->GetXaxis()->SetNdivisions(505);
   frame_179c990__4->GetXaxis()->SetLabelFont(42);
   frame_179c990__4->GetXaxis()->SetLabelSize(0.05);
   frame_179c990__4->GetXaxis()->SetTitleSize(0.055);
   frame_179c990__4->GetXaxis()->SetTitleOffset(1.2);
   frame_179c990__4->GetXaxis()->SetTitleFont(42);
   frame_179c990__4->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_179c990__4->GetYaxis()->SetLabelFont(42);
   frame_179c990__4->GetYaxis()->SetLabelOffset(0.01);
   frame_179c990__4->GetYaxis()->SetLabelSize(0.05);
   frame_179c990__4->GetYaxis()->SetTitleSize(0.055);
   frame_179c990__4->GetYaxis()->SetTitleOffset(1.4);
   frame_179c990__4->GetYaxis()->SetTitleFont(42);
   frame_179c990__4->GetZaxis()->SetLabelFont(42);
   frame_179c990__4->GetZaxis()->SetLabelSize(0.035);
   frame_179c990__4->GetZaxis()->SetTitleSize(0.035);
   frame_179c990__4->GetZaxis()->SetTitleFont(42);
   frame_179c990__4->Draw("AXISSAME");
   
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
   cfail_bf->Modified();
   cfail_bf->cd();
   cfail_bf->SetSelected(cfail_bf);
}
