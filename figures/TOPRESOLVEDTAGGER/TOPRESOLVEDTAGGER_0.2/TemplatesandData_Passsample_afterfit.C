{
//=========Macro generated from canvas: cpass/cpass
//=========  (Mon Oct 26 23:57:35 2015) by ROOT version5.32/00
   TCanvas *cpass = new TCanvas("cpass", "cpass",1,1,800,1176);
   gStyle->SetOptStat(0);
   cpass->SetHighLightColor(2);
   cpass->Range(-187.013,-343.8453,851.9481,1948.457);
   cpass->SetFillColor(0);
   cpass->SetBorderMode(0);
   cpass->SetBorderSize(10);
   cpass->SetTickx(1);
   cpass->SetTicky(1);
   cpass->SetLeftMargin(0.18);
   cpass->SetRightMargin(0.05);
   cpass->SetTopMargin(0.08);
   cpass->SetBottomMargin(0.15);
   cpass->SetFrameFillStyle(0);
   cpass->SetFrameLineStyle(0);
   cpass->SetFrameBorderMode(0);
   cpass->SetFrameBorderSize(10);
   cpass->SetFrameFillStyle(0);
   cpass->SetFrameLineStyle(0);
   cpass->SetFrameBorderMode(0);
   cpass->SetFrameBorderSize(10);
   
   TH1D *frame_29ec0f0__5 = new TH1D("frame_29ec0f0__5","",40,0,800);
   frame_29ec0f0__5->SetMaximum(1765.073);
   frame_29ec0f0__5->SetDirectory(0);
   frame_29ec0f0__5->SetStats(0);
   frame_29ec0f0__5->SetFillColor(2);
   frame_29ec0f0__5->SetFillStyle(0);
   frame_29ec0f0__5->SetLineStyle(0);
   frame_29ec0f0__5->SetLineWidth(2);
   frame_29ec0f0__5->SetMarkerStyle(20);
   frame_29ec0f0__5->SetMarkerSize(1.2);
   frame_29ec0f0__5->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_29ec0f0__5->GetXaxis()->SetNdivisions(505);
   frame_29ec0f0__5->GetXaxis()->SetLabelFont(42);
   frame_29ec0f0__5->GetXaxis()->SetLabelSize(0.05);
   frame_29ec0f0__5->GetXaxis()->SetTitleSize(0.055);
   frame_29ec0f0__5->GetXaxis()->SetTitleOffset(1.2);
   frame_29ec0f0__5->GetXaxis()->SetTitleFont(42);
   frame_29ec0f0__5->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_29ec0f0__5->GetYaxis()->SetLabelFont(42);
   frame_29ec0f0__5->GetYaxis()->SetLabelOffset(0.01);
   frame_29ec0f0__5->GetYaxis()->SetLabelSize(0.05);
   frame_29ec0f0__5->GetYaxis()->SetTitleSize(0.055);
   frame_29ec0f0__5->GetYaxis()->SetTitleOffset(1.4);
   frame_29ec0f0__5->GetYaxis()->SetTitleFont(42);
   frame_29ec0f0__5->GetZaxis()->SetLabelFont(42);
   frame_29ec0f0__5->GetZaxis()->SetLabelSize(0.035);
   frame_29ec0f0__5->GetZaxis()->SetTitleSize(0.035);
   frame_29ec0f0__5->GetZaxis()->SetTitleFont(42);
   frame_29ec0f0__5->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(40);
   grae->SetName("h_dataPass");
   grae->SetTitle("Histogram of dataPass_plot__m");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   grae->SetPoint(0,10,0);
   grae->SetPointError(0,10,10,0,1.147874);
   grae->SetPoint(1,30,0);
   grae->SetPointError(1,10,10,0,1.147874);
   grae->SetPoint(2,50,3);
   grae->SetPointError(2,10,10,1.632705,2.918186);
   grae->SetPoint(3,70,22);
   grae->SetPointError(3,10,10,4.654502,5.761366);
   grae->SetPoint(4,90,115);
   grae->SetPointError(4,10,10,10.23546,11.23546);
   grae->SetPoint(5,110,256);
   grae->SetPointError(5,10,10,15.50781,16.50781);
   grae->SetPoint(6,130,459);
   grae->SetPointError(6,10,10,20.93012,21.93012);
   grae->SetPoint(7,150,760);
   grae->SetPointError(7,10,10,27.07263,28.07263);
   grae->SetPoint(8,170,1618);
   grae->SetPointError(8,10,10,39.72748,40.72748);
   grae->SetPoint(9,190,1414);
   grae->SetPointError(9,10,10,37.10652,38.10652);
   grae->SetPoint(10,210,615);
   grae->SetPointError(10,10,10,24.30423,25.30423);
   grae->SetPoint(11,230,269);
   grae->SetPointError(11,10,10,15.90884,16.90884);
   grae->SetPoint(12,250,152);
   grae->SetPointError(12,10,10,11.83896,12.83896);
   grae->SetPoint(13,270,86);
   grae->SetPointError(13,10,10,9.255555,10.30959);
   grae->SetPoint(14,290,58);
   grae->SetPointError(14,10,10,7.593751,8.659564);
   grae->SetPoint(15,310,37);
   grae->SetPointError(15,10,10,6.055143,7.137555);
   grae->SetPoint(16,330,33);
   grae->SetPointError(16,10,10,5.715302,6.802567);
   grae->SetPoint(17,350,25);
   grae->SetPointError(17,10,10,4.966335,6.066589);
   grae->SetPoint(18,370,25);
   grae->SetPointError(18,10,10,4.966335,6.066589);
   grae->SetPoint(19,390,17);
   grae->SetPointError(19,10,10,4.082184,5.203719);
   grae->SetPoint(20,410,9);
   grae->SetPointError(20,10,10,2.943461,4.110204);
   grae->SetPoint(21,430,11);
   grae->SetPointError(21,10,10,3.265579,4.416521);
   grae->SetPoint(22,450,14);
   grae->SetPointError(22,10,10,3.6965,4.830381);
   grae->SetPoint(23,470,7);
   grae->SetPointError(23,10,10,2.58147,3.770281);
   grae->SetPoint(24,490,5);
   grae->SetPointError(24,10,10,2.159691,3.382473);
   grae->SetPoint(25,510,4);
   grae->SetPointError(25,10,10,1.914339,3.162753);
   grae->SetPoint(26,530,2);
   grae->SetPointError(26,10,10,1.291815,2.63786);
   grae->SetPoint(27,550,0);
   grae->SetPointError(27,10,10,0,1.147874);
   grae->SetPoint(28,570,1);
   grae->SetPointError(28,10,10,0.8272462,2.299527);
   grae->SetPoint(29,590,1);
   grae->SetPointError(29,10,10,0.8272462,2.299527);
   grae->SetPoint(30,610,2);
   grae->SetPointError(30,10,10,1.291815,2.63786);
   grae->SetPoint(31,630,0);
   grae->SetPointError(31,10,10,0,1.147874);
   grae->SetPoint(32,650,1);
   grae->SetPointError(32,10,10,0.8272462,2.299527);
   grae->SetPoint(33,670,2);
   grae->SetPointError(33,10,10,1.291815,2.63786);
   grae->SetPoint(34,690,0);
   grae->SetPointError(34,10,10,0,1.147874);
   grae->SetPoint(35,710,1);
   grae->SetPointError(35,10,10,0.8272462,2.299527);
   grae->SetPoint(36,730,0);
   grae->SetPointError(36,10,10,0,1.147874);
   grae->SetPoint(37,750,0);
   grae->SetPointError(37,10,10,0,1.147874);
   grae->SetPoint(38,770,0);
   grae->SetPointError(38,10,10,0,1.147874);
   grae->SetPoint(39,790,0);
   grae->SetPointError(39,10,10,0,1.147874);
   
   TH1F *Graph_h_dataPass3 = new TH1F("Graph_h_dataPass3","Histogram of dataPass_plot__m",100,0,880);
   Graph_h_dataPass3->SetMinimum(0);
   Graph_h_dataPass3->SetMaximum(1824.6);
   Graph_h_dataPass3->SetDirectory(0);
   Graph_h_dataPass3->SetStats(0);
   Graph_h_dataPass3->SetFillColor(2);
   Graph_h_dataPass3->SetFillStyle(0);
   Graph_h_dataPass3->SetLineStyle(0);
   Graph_h_dataPass3->SetLineWidth(2);
   Graph_h_dataPass3->SetMarkerStyle(20);
   Graph_h_dataPass3->SetMarkerSize(1.2);
   Graph_h_dataPass3->GetXaxis()->SetNdivisions(505);
   Graph_h_dataPass3->GetXaxis()->SetLabelFont(42);
   Graph_h_dataPass3->GetXaxis()->SetLabelSize(0.05);
   Graph_h_dataPass3->GetXaxis()->SetTitleSize(0.055);
   Graph_h_dataPass3->GetXaxis()->SetTitleOffset(1.2);
   Graph_h_dataPass3->GetXaxis()->SetTitleFont(42);
   Graph_h_dataPass3->GetYaxis()->SetLabelFont(42);
   Graph_h_dataPass3->GetYaxis()->SetLabelOffset(0.01);
   Graph_h_dataPass3->GetYaxis()->SetLabelSize(0.05);
   Graph_h_dataPass3->GetYaxis()->SetTitleSize(0.055);
   Graph_h_dataPass3->GetYaxis()->SetTitleOffset(1.4);
   Graph_h_dataPass3->GetYaxis()->SetTitleFont(42);
   Graph_h_dataPass3->GetZaxis()->SetLabelFont(42);
   Graph_h_dataPass3->GetZaxis()->SetLabelSize(0.035);
   Graph_h_dataPass3->GetZaxis()->SetTitleSize(0.035);
   Graph_h_dataPass3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_h_dataPass3);
   
   grae->Draw("zp");
   
   TGraph *graph = new TGraph(87);
   graph->SetName("modelPass_Norm[m]");
   graph->SetTitle("Projection of Model for PASS sample");
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
   graph->SetPoint(7,40.00000816,1.214460925);
   graph->SetPoint(8,59.99999184,1.214460925);
   graph->SetPoint(9,60.00000816,30.59664481);
   graph->SetPoint(10,79.99999184,30.59664481);
   graph->SetPoint(11,80.00000816,121.3015821);
   graph->SetPoint(12,99.99999184,121.3015821);
   graph->SetPoint(13,100.0000082,242.6458432);
   graph->SetPoint(14,119.9999918,242.6458432);
   graph->SetPoint(15,120.0000082,419.3905219);
   graph->SetPoint(16,139.9999918,419.3905219);
   graph->SetPoint(17,140.0000082,782.8555981);
   graph->SetPoint(18,159.9999918,782.8555981);
   graph->SetPoint(19,160.0000082,1681.021716);
   graph->SetPoint(20,179.9999918,1681.021716);
   graph->SetPoint(21,180.0000082,1352.909565);
   graph->SetPoint(22,199.9999918,1352.909565);
   graph->SetPoint(23,200.0000082,596.2569502);
   graph->SetPoint(24,219.9999918,596.2569502);
   graph->SetPoint(25,220.0000082,271.8284163);
   graph->SetPoint(26,239.9999918,271.8284163);
   graph->SetPoint(27,240.0000082,154.8385409);
   graph->SetPoint(28,259.9999918,154.8385409);
   graph->SetPoint(29,260.0000082,92.23634146);
   graph->SetPoint(30,279.9999918,92.23634146);
   graph->SetPoint(31,280.0000082,63.0920442);
   graph->SetPoint(32,299.9999918,63.0920442);
   graph->SetPoint(33,300.0000082,44.10225796);
   graph->SetPoint(34,319.9999918,44.10225796);
   graph->SetPoint(35,320.0000082,31.98818104);
   graph->SetPoint(36,339.9999918,31.98818104);
   graph->SetPoint(37,340.0000082,26.46514545);
   graph->SetPoint(38,359.9999918,26.46514545);
   graph->SetPoint(39,360.0000082,18.18665186);
   graph->SetPoint(40,379.9999918,18.18665186);
   graph->SetPoint(41,380.0000082,18.32478621);
   graph->SetPoint(42,399.9999918,18.32478621);
   graph->SetPoint(43,400.0000082,18.20198828);
   graph->SetPoint(44,419.9999918,18.20198828);
   graph->SetPoint(45,420.0000082,12.81694583);
   graph->SetPoint(46,439.9999918,12.81694583);
   graph->SetPoint(47,440.0000082,8.938214349);
   graph->SetPoint(48,459.9999918,8.938214349);
   graph->SetPoint(49,460.0000082,7.906918979);
   graph->SetPoint(50,479.9999918,7.906918979);
   graph->SetPoint(51,480.0000082,6.231488842);
   graph->SetPoint(52,499.9999918,6.231488842);
   graph->SetPoint(53,500.0000082,4.373613574);
   graph->SetPoint(54,519.9999918,4.373613574);
   graph->SetPoint(55,520.0000082,3.923213867);
   graph->SetPoint(56,539.9999918,3.923213867);
   graph->SetPoint(57,540.0000082,2.929918106);
   graph->SetPoint(58,559.9999918,2.929918106);
   graph->SetPoint(59,560.0000082,2.059304148);
   graph->SetPoint(60,579.9999918,2.059304148);
   graph->SetPoint(61,580.0000082,2.374904011);
   graph->SetPoint(62,599.9999918,2.374904011);
   graph->SetPoint(63,600.0000082,0.9337355315);
   graph->SetPoint(64,619.9999918,0.9337355315);
   graph->SetPoint(65,620.0000082,0.8147026293);
   graph->SetPoint(66,639.9999918,0.8147026293);
   graph->SetPoint(67,640.0000082,1.11442847);
   graph->SetPoint(68,659.9999918,1.11442847);
   graph->SetPoint(69,660.0000082,0.5362324208);
   graph->SetPoint(70,679.9999918,0.5362324208);
   graph->SetPoint(71,680.0000082,0.5792791573);
   graph->SetPoint(72,699.9999918,0.5792791573);
   graph->SetPoint(73,700.0000082,0.1767552225);
   graph->SetPoint(74,719.9999918,0.1767552225);
   graph->SetPoint(75,720.0000082,0.3983722608);
   graph->SetPoint(76,739.9999918,0.3983722608);
   graph->SetPoint(77,740.0000082,0.3013511845);
   graph->SetPoint(78,759.9999918,0.3013511845);
   graph->SetPoint(79,760.0000082,0.04694718121);
   graph->SetPoint(80,779.9999918,0.04694718121);
   graph->SetPoint(81,780.0000082,0.08643842926);
   graph->SetPoint(82,799.9999918,0.08643842926);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]9 = new TH1F("Graph_modelPass_Norm[m]9","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]9->SetMinimum(0);
   Graph_modelPass_Norm[m]9->SetMaximum(1849.124);
   Graph_modelPass_Norm[m]9->SetDirectory(0);
   Graph_modelPass_Norm[m]9->SetStats(0);
   Graph_modelPass_Norm[m]9->SetFillColor(2);
   Graph_modelPass_Norm[m]9->SetFillStyle(0);
   Graph_modelPass_Norm[m]9->SetLineStyle(0);
   Graph_modelPass_Norm[m]9->SetLineWidth(2);
   Graph_modelPass_Norm[m]9->SetMarkerStyle(20);
   Graph_modelPass_Norm[m]9->SetMarkerSize(1.2);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetNdivisions(505);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelPass_Norm[m]9->GetXaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelPass_Norm[m]9->GetYaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]9->GetZaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]9->GetZaxis()->SetLabelSize(0.035);
   Graph_modelPass_Norm[m]9->GetZaxis()->SetTitleSize(0.035);
   Graph_modelPass_Norm[m]9->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelPass_Norm[m]9);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelPass_Norm[m]_Comp[sigModPass]");
   graph->SetTitle("Projection of Model for PASS sample");
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
   graph->SetPoint(11,80.00000816,0.6215446494);
   graph->SetPoint(12,99.99999184,0.6215446494);
   graph->SetPoint(13,100.0000082,2.289633667);
   graph->SetPoint(14,119.9999918,2.289633667);
   graph->SetPoint(15,120.0000082,21.02641383);
   graph->SetPoint(16,139.9999918,21.02641383);
   graph->SetPoint(17,140.0000082,174.7887379);
   graph->SetPoint(18,159.9999918,174.7887379);
   graph->SetPoint(19,160.0000082,719.3609889);
   graph->SetPoint(20,179.9999918,719.3609889);
   graph->SetPoint(21,180.0000082,662.5155188);
   graph->SetPoint(22,199.9999918,662.5155188);
   graph->SetPoint(23,200.0000082,184.7508354);
   graph->SetPoint(24,219.9999918,184.7508354);
   graph->SetPoint(25,220.0000082,36.09530801);
   graph->SetPoint(26,239.9999918,36.09530801);
   graph->SetPoint(27,240.0000082,8.488736501);
   graph->SetPoint(28,259.9999918,8.488736501);
   graph->SetPoint(29,260.0000082,3.554787166);
   graph->SetPoint(30,279.9999918,3.554787166);
   graph->SetPoint(31,280.0000082,1.522189182);
   graph->SetPoint(32,299.9999918,1.522189182);
   graph->SetPoint(33,300.0000082,0.658735643);
   graph->SetPoint(34,319.9999918,0.658735643);
   graph->SetPoint(35,320.0000082,0.3363172507);
   graph->SetPoint(36,339.9999918,0.3363172507);
   graph->SetPoint(37,340.0000082,0.2690919448);
   graph->SetPoint(38,359.9999918,0.2690919448);
   graph->SetPoint(39,360.0000082,0.1639761803);
   graph->SetPoint(40,379.9999918,0.1639761803);
   graph->SetPoint(41,380.0000082,0.14463494);
   graph->SetPoint(42,399.9999918,0.14463494);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0.08226627517);
   graph->SetPoint(48,459.9999918,0.08226627517);
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
   
   TH1F *Graph_modelPass_Norm[m]_Comp[sigModPass]10 = new TH1F("Graph_modelPass_Norm[m]_Comp[sigModPass]10","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMaximum(791.2971);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetDirectory(0);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetStats(0);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetFillColor(2);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetFillStyle(0);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetLineStyle(0);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetLineWidth(2);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMarkerStyle(20);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMarkerSize(1.2);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetNdivisions(505);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetXaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetYaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetZaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetZaxis()->SetLabelSize(0.035);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetZaxis()->SetTitleSize(0.035);
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelPass_Norm[m]_Comp[sigModPass]10);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelPass_Norm[m]_Comp[bkg1ModPass]");
   graph->SetTitle("Projection of Model for PASS sample");
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
   graph->SetPoint(7,40.00000816,0.4474452111);
   graph->SetPoint(8,59.99999184,0.4474452111);
   graph->SetPoint(9,60.00000816,10.77599125);
   graph->SetPoint(10,79.99999184,10.77599125);
   graph->SetPoint(11,80.00000816,55.43572157);
   graph->SetPoint(12,99.99999184,55.43572157);
   graph->SetPoint(13,100.0000082,119.6639812);
   graph->SetPoint(14,119.9999918,119.6639812);
   graph->SetPoint(15,120.0000082,213.1718319);
   graph->SetPoint(16,139.9999918,213.1718319);
   graph->SetPoint(17,140.0000082,350.943591);
   graph->SetPoint(18,159.9999918,350.943591);
   graph->SetPoint(19,160.0000082,555.7176404);
   graph->SetPoint(20,179.9999918,555.7176404);
   graph->SetPoint(21,180.0000082,389.9611555);
   graph->SetPoint(22,199.9999918,389.9611555);
   graph->SetPoint(23,200.0000082,235.8338513);
   graph->SetPoint(24,219.9999918,235.8338513);
   graph->SetPoint(25,220.0000082,134.0754959);
   graph->SetPoint(26,239.9999918,134.0754959);
   graph->SetPoint(27,240.0000082,81.82659975);
   graph->SetPoint(28,259.9999918,81.82659975);
   graph->SetPoint(29,260.0000082,51.58390449);
   graph->SetPoint(30,279.9999918,51.58390449);
   graph->SetPoint(31,280.0000082,34.04107991);
   graph->SetPoint(32,299.9999918,34.04107991);
   graph->SetPoint(33,300.0000082,22.03036888);
   graph->SetPoint(34,319.9999918,22.03036888);
   graph->SetPoint(35,320.0000082,15.76313698);
   graph->SetPoint(36,339.9999918,15.76313698);
   graph->SetPoint(37,340.0000082,13.04165563);
   graph->SetPoint(38,359.9999918,13.04165563);
   graph->SetPoint(39,360.0000082,9.682374615);
   graph->SetPoint(40,379.9999918,9.682374615);
   graph->SetPoint(41,380.0000082,8.514815847);
   graph->SetPoint(42,399.9999918,8.514815847);
   graph->SetPoint(43,400.0000082,7.968155156);
   graph->SetPoint(44,419.9999918,7.968155156);
   graph->SetPoint(45,420.0000082,6.913495195);
   graph->SetPoint(46,439.9999918,6.913495195);
   graph->SetPoint(47,440.0000082,4.055510097);
   graph->SetPoint(48,459.9999918,4.055510097);
   graph->SetPoint(49,460.0000082,4.141993005);
   graph->SetPoint(50,479.9999918,4.141993005);
   graph->SetPoint(51,480.0000082,3.178215585);
   graph->SetPoint(52,499.9999918,3.178215585);
   graph->SetPoint(53,500.0000082,2.467708245);
   graph->SetPoint(54,519.9999918,2.467708245);
   graph->SetPoint(55,520.0000082,2.003299968);
   graph->SetPoint(56,539.9999918,2.003299968);
   graph->SetPoint(57,540.0000082,1.433786649);
   graph->SetPoint(58,559.9999918,1.433786649);
   graph->SetPoint(59,560.0000082,1.051573531);
   graph->SetPoint(60,579.9999918,1.051573531);
   graph->SetPoint(61,580.0000082,1.178559625);
   graph->SetPoint(62,599.9999918,1.178559625);
   graph->SetPoint(63,600.0000082,0.6667895806);
   graph->SetPoint(64,619.9999918,0.6667895806);
   graph->SetPoint(65,620.0000082,0.5067801876);
   graph->SetPoint(66,639.9999918,0.5067801876);
   graph->SetPoint(67,640.0000082,0.4791833735);
   graph->SetPoint(68,659.9999918,0.4791833735);
   graph->SetPoint(69,660.0000082,0.4529201627);
   graph->SetPoint(70,679.9999918,0.4529201627);
   graph->SetPoint(71,680.0000082,0.5414131028);
   graph->SetPoint(72,699.9999918,0.5414131028);
   graph->SetPoint(73,700.0000082,0.105924354);
   graph->SetPoint(74,719.9999918,0.105924354);
   graph->SetPoint(75,720.0000082,0.1218125031);
   graph->SetPoint(76,739.9999918,0.1218125031);
   graph->SetPoint(77,740.0000082,0.1538294045);
   graph->SetPoint(78,759.9999918,0.1538294045);
   graph->SetPoint(79,760.0000082,0.04694718121);
   graph->SetPoint(80,779.9999918,0.04694718121);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMaximum(611.2894);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetDirectory(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetStats(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetFillColor(2);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetFillStyle(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetLineStyle(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetLineWidth(2);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMarkerStyle(20);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMarkerSize(1.2);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetNdivisions(505);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetXaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetYaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetZaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetZaxis()->SetLabelSize(0.035);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetZaxis()->SetTitleSize(0.035);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11);
   
   graph->Draw("l");
   
   graph = new TGraph(87);
   graph->SetName("modelPass_Norm[m]_Comp[bkg2ModPass]");
   graph->SetTitle("Projection of Model for PASS sample");
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
   graph->SetPoint(7,40.00000816,0.7670157138);
   graph->SetPoint(8,59.99999184,0.7670157138);
   graph->SetPoint(9,60.00000816,19.82065355);
   graph->SetPoint(10,79.99999184,19.82065355);
   graph->SetPoint(11,80.00000816,65.24431587);
   graph->SetPoint(12,99.99999184,65.24431587);
   graph->SetPoint(13,100.0000082,120.6922283);
   graph->SetPoint(14,119.9999918,120.6922283);
   graph->SetPoint(15,120.0000082,185.1922762);
   graph->SetPoint(16,139.9999918,185.1922762);
   graph->SetPoint(17,140.0000082,257.1232692);
   graph->SetPoint(18,159.9999918,257.1232692);
   graph->SetPoint(19,160.0000082,405.943087);
   graph->SetPoint(20,179.9999918,405.943087);
   graph->SetPoint(21,180.0000082,300.4328905);
   graph->SetPoint(22,199.9999918,300.4328905);
   graph->SetPoint(23,200.0000082,175.6722635);
   graph->SetPoint(24,219.9999918,175.6722635);
   graph->SetPoint(25,220.0000082,101.6576124);
   graph->SetPoint(26,239.9999918,101.6576124);
   graph->SetPoint(27,240.0000082,64.52320465);
   graph->SetPoint(28,259.9999918,64.52320465);
   graph->SetPoint(29,260.0000082,37.09764981);
   graph->SetPoint(30,279.9999918,37.09764981);
   graph->SetPoint(31,280.0000082,27.52877511);
   graph->SetPoint(32,299.9999918,27.52877511);
   graph->SetPoint(33,300.0000082,21.41315344);
   graph->SetPoint(34,319.9999918,21.41315344);
   graph->SetPoint(35,320.0000082,15.88872681);
   graph->SetPoint(36,339.9999918,15.88872681);
   graph->SetPoint(37,340.0000082,13.15439787);
   graph->SetPoint(38,359.9999918,13.15439787);
   graph->SetPoint(39,360.0000082,8.340301064);
   graph->SetPoint(40,379.9999918,8.340301064);
   graph->SetPoint(41,380.0000082,9.665335422);
   graph->SetPoint(42,399.9999918,9.665335422);
   graph->SetPoint(43,400.0000082,10.23383312);
   graph->SetPoint(44,419.9999918,10.23383312);
   graph->SetPoint(45,420.0000082,5.903450634);
   graph->SetPoint(46,439.9999918,5.903450634);
   graph->SetPoint(47,440.0000082,4.800437977);
   graph->SetPoint(48,459.9999918,4.800437977);
   graph->SetPoint(49,460.0000082,3.764925974);
   graph->SetPoint(50,479.9999918,3.764925974);
   graph->SetPoint(51,480.0000082,3.053273257);
   graph->SetPoint(52,499.9999918,3.053273257);
   graph->SetPoint(53,500.0000082,1.905905329);
   graph->SetPoint(54,519.9999918,1.905905329);
   graph->SetPoint(55,520.0000082,1.9199139);
   graph->SetPoint(56,539.9999918,1.9199139);
   graph->SetPoint(57,540.0000082,1.496131457);
   graph->SetPoint(58,559.9999918,1.496131457);
   graph->SetPoint(59,560.0000082,1.007730617);
   graph->SetPoint(60,579.9999918,1.007730617);
   graph->SetPoint(61,580.0000082,1.196344386);
   graph->SetPoint(62,599.9999918,1.196344386);
   graph->SetPoint(63,600.0000082,0.2669459508);
   graph->SetPoint(64,619.9999918,0.2669459508);
   graph->SetPoint(65,620.0000082,0.3079224417);
   graph->SetPoint(66,639.9999918,0.3079224417);
   graph->SetPoint(67,640.0000082,0.6352450962);
   graph->SetPoint(68,659.9999918,0.6352450962);
   graph->SetPoint(69,660.0000082,0.08331225807);
   graph->SetPoint(70,679.9999918,0.08331225807);
   graph->SetPoint(71,680.0000082,0.03786605453);
   graph->SetPoint(72,699.9999918,0.03786605453);
   graph->SetPoint(73,700.0000082,0.07083086859);
   graph->SetPoint(74,719.9999918,0.07083086859);
   graph->SetPoint(75,720.0000082,0.2765597577);
   graph->SetPoint(76,739.9999918,0.2765597577);
   graph->SetPoint(77,740.0000082,0.14752178);
   graph->SetPoint(78,759.9999918,0.14752178);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0.08643842926);
   graph->SetPoint(82,799.9999918,0.08643842926);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMaximum(446.5374);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetDirectory(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetStats(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetFillColor(2);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetFillStyle(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetLineStyle(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetLineWidth(2);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMarkerStyle(20);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMarkerSize(1.2);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetNdivisions(505);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetXaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetLabelSize(0.05);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetTitleSize(0.055);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetYaxis()->SetTitleFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetZaxis()->SetLabelFont(42);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetZaxis()->SetLabelSize(0.035);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetZaxis()->SetTitleSize(0.035);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12);
   
   graph->Draw("l");
   
   TH1D *frame_29ec0f0__6 = new TH1D("frame_29ec0f0__6","",40,0,800);
   frame_29ec0f0__6->SetMaximum(1765.073);
   frame_29ec0f0__6->SetDirectory(0);
   frame_29ec0f0__6->SetStats(0);
   frame_29ec0f0__6->SetFillColor(2);
   frame_29ec0f0__6->SetFillStyle(0);
   frame_29ec0f0__6->SetLineStyle(0);
   frame_29ec0f0__6->SetLineWidth(2);
   frame_29ec0f0__6->SetMarkerStyle(20);
   frame_29ec0f0__6->SetMarkerSize(1.2);
   frame_29ec0f0__6->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_29ec0f0__6->GetXaxis()->SetNdivisions(505);
   frame_29ec0f0__6->GetXaxis()->SetLabelFont(42);
   frame_29ec0f0__6->GetXaxis()->SetLabelSize(0.05);
   frame_29ec0f0__6->GetXaxis()->SetTitleSize(0.055);
   frame_29ec0f0__6->GetXaxis()->SetTitleOffset(1.2);
   frame_29ec0f0__6->GetXaxis()->SetTitleFont(42);
   frame_29ec0f0__6->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_29ec0f0__6->GetYaxis()->SetLabelFont(42);
   frame_29ec0f0__6->GetYaxis()->SetLabelOffset(0.01);
   frame_29ec0f0__6->GetYaxis()->SetLabelSize(0.05);
   frame_29ec0f0__6->GetYaxis()->SetTitleSize(0.055);
   frame_29ec0f0__6->GetYaxis()->SetTitleOffset(1.4);
   frame_29ec0f0__6->GetYaxis()->SetTitleFont(42);
   frame_29ec0f0__6->GetZaxis()->SetLabelFont(42);
   frame_29ec0f0__6->GetZaxis()->SetLabelSize(0.035);
   frame_29ec0f0__6->GetZaxis()->SetTitleSize(0.035);
   frame_29ec0f0__6->GetZaxis()->SetTitleFont(42);
   frame_29ec0f0__6->Draw("AXISSAME");
   
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
   cpass->Modified();
   cpass->cd();
   cpass->SetSelected(cpass);
}
