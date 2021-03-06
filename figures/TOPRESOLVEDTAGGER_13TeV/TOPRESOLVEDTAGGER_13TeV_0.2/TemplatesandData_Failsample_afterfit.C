{
//=========Macro generated from canvas: cfail/cfail
//=========  (Sat Oct 31 02:00:54 2015) by ROOT version5.32/00
   TCanvas *cfail = new TCanvas("cfail", "cfail",1,1,800,1176);
   gStyle->SetOptStat(0);
   cfail->SetHighLightColor(2);
   cfail->Range(-187.013,-43.4823,851.9481,246.3997);
   cfail->SetFillColor(0);
   cfail->SetBorderMode(0);
   cfail->SetBorderSize(10);
   cfail->SetTickx(1);
   cfail->SetTicky(1);
   cfail->SetLeftMargin(0.18);
   cfail->SetRightMargin(0.05);
   cfail->SetTopMargin(0.08);
   cfail->SetBottomMargin(0.15);
   cfail->SetFrameFillStyle(0);
   cfail->SetFrameLineStyle(0);
   cfail->SetFrameBorderMode(0);
   cfail->SetFrameBorderSize(10);
   cfail->SetFrameFillStyle(0);
   cfail->SetFrameLineStyle(0);
   cfail->SetFrameBorderMode(0);
   cfail->SetFrameBorderSize(10);
   
   TH1D *frame_1f1d040__7 = new TH1D("frame_1f1d040__7","",40,0,800);
   frame_1f1d040__7->SetMaximum(223.2091);
   frame_1f1d040__7->SetDirectory(0);
   frame_1f1d040__7->SetStats(0);
   frame_1f1d040__7->SetFillColor(2);
   frame_1f1d040__7->SetFillStyle(0);
   frame_1f1d040__7->SetLineStyle(0);
   frame_1f1d040__7->SetLineWidth(2);
   frame_1f1d040__7->SetMarkerStyle(20);
   frame_1f1d040__7->SetMarkerSize(1.2);
   frame_1f1d040__7->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_1f1d040__7->GetXaxis()->SetNdivisions(505);
   frame_1f1d040__7->GetXaxis()->SetLabelFont(42);
   frame_1f1d040__7->GetXaxis()->SetLabelSize(0.05);
   frame_1f1d040__7->GetXaxis()->SetTitleSize(0.055);
   frame_1f1d040__7->GetXaxis()->SetTitleOffset(1.2);
   frame_1f1d040__7->GetXaxis()->SetTitleFont(42);
   frame_1f1d040__7->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_1f1d040__7->GetYaxis()->SetLabelFont(42);
   frame_1f1d040__7->GetYaxis()->SetLabelOffset(0.01);
   frame_1f1d040__7->GetYaxis()->SetLabelSize(0.05);
   frame_1f1d040__7->GetYaxis()->SetTitleSize(0.055);
   frame_1f1d040__7->GetYaxis()->SetTitleOffset(1.4);
   frame_1f1d040__7->GetYaxis()->SetTitleFont(42);
   frame_1f1d040__7->GetZaxis()->SetLabelFont(42);
   frame_1f1d040__7->GetZaxis()->SetLabelSize(0.035);
   frame_1f1d040__7->GetZaxis()->SetTitleSize(0.035);
   frame_1f1d040__7->GetZaxis()->SetTitleFont(42);
   frame_1f1d040__7->Draw("");
   
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
   grae->SetPoint(2,50,0);
   grae->SetPointError(2,10,10,0,1.147874);
   grae->SetPoint(3,70,1);
   grae->SetPointError(3,10,10,0.8272462,2.299527);
   grae->SetPoint(4,90,11);
   grae->SetPointError(4,10,10,3.265579,4.416521);
   grae->SetPoint(5,110,21);
   grae->SetPointError(5,10,10,4.545807,5.655182);
   grae->SetPoint(6,130,60);
   grae->SetPointError(6,10,10,7.724317,8.789023);
   grae->SetPoint(7,150,110);
   grae->SetPointError(7,10,10,10,11);
   grae->SetPoint(8,170,129);
   grae->SetPointError(8,10,10,10.86882,11.86882);
   grae->SetPoint(9,190,159);
   grae->SetPointError(9,10,10,12.11943,13.11943);
   grae->SetPoint(10,210,162);
   grae->SetPointError(10,10,10,12.23774,13.23774);
   grae->SetPoint(11,230,198);
   grae->SetPointError(11,10,10,13.58013,14.58013);
   grae->SetPoint(12,250,186);
   grae->SetPointError(12,10,10,13.14734,14.14734);
   grae->SetPoint(13,270,155);
   grae->SetPointError(13,10,10,11.95994,12.95994);
   grae->SetPoint(14,290,151);
   grae->SetPointError(14,10,10,11.79837,12.79837);
   grae->SetPoint(15,310,126);
   grae->SetPointError(15,10,10,10.7361,11.7361);
   grae->SetPoint(16,330,117);
   grae->SetPointError(16,10,10,10.3282,11.3282);
   grae->SetPoint(17,350,123);
   grae->SetPointError(17,10,10,10.6018,11.6018);
   grae->SetPoint(18,370,93);
   grae->SetPointError(18,10,10,9.626284,10.67824);
   grae->SetPoint(19,390,89);
   grae->SetPointError(19,10,10,9.416226,10.46934);
   grae->SetPoint(20,410,76);
   grae->SetPointError(20,10,10,8.698576,9.756061);
   grae->SetPoint(21,430,66);
   grae->SetPointError(21,10,10,8.103403,9.165094);
   grae->SetPoint(22,450,62);
   grae->SetPointError(22,10,10,7.852713,8.916365);
   grae->SetPoint(23,470,57);
   grae->SetPointError(23,10,10,7.527619,8.594007);
   grae->SetPoint(24,490,44);
   grae->SetPointError(24,10,10,6.60794,7.68351);
   grae->SetPoint(25,510,41);
   grae->SetPointError(25,10,10,6.376898,7.455185);
   grae->SetPoint(26,530,33);
   grae->SetPointError(26,10,10,5.715302,6.802567);
   grae->SetPoint(27,550,33);
   grae->SetPointError(27,10,10,5.715302,6.802567);
   grae->SetPoint(28,570,26);
   grae->SetPointError(28,10,10,5.066015,6.164324);
   grae->SetPoint(29,590,23);
   grae->SetPointError(29,10,10,4.760717,5.865235);
   grae->SetPoint(30,610,19);
   grae->SetPointError(30,10,10,4.320219,5.435196);
   grae->SetPoint(31,630,26);
   grae->SetPointError(31,10,10,5.066015,6.164324);
   grae->SetPoint(32,650,12);
   grae->SetPointError(32,10,10,3.415266,4.559819);
   grae->SetPoint(33,670,15);
   grae->SetPointError(33,10,10,3.82938,4.958738);
   grae->SetPoint(34,690,16);
   grae->SetPointError(34,10,10,3.957801,5.083066);
   grae->SetPoint(35,710,21);
   grae->SetPointError(35,10,10,4.545807,5.655182);
   grae->SetPoint(36,730,16);
   grae->SetPointError(36,10,10,3.957801,5.083066);
   grae->SetPoint(37,750,10);
   grae->SetPointError(37,10,10,3.108694,4.26695);
   grae->SetPoint(38,770,15);
   grae->SetPointError(38,10,10,3.82938,4.958738);
   grae->SetPoint(39,790,7);
   grae->SetPointError(39,10,10,2.58147,3.770281);
   
   TH1F *Graph_h_dataFail4 = new TH1F("Graph_h_dataFail4","Histogram of dataFail_plot__m",100,0,880);
   Graph_h_dataFail4->SetMinimum(0);
   Graph_h_dataFail4->SetMaximum(233.8381);
   Graph_h_dataFail4->SetDirectory(0);
   Graph_h_dataFail4->SetStats(0);
   Graph_h_dataFail4->SetFillColor(2);
   Graph_h_dataFail4->SetFillStyle(0);
   Graph_h_dataFail4->SetLineStyle(0);
   Graph_h_dataFail4->SetLineWidth(2);
   Graph_h_dataFail4->SetMarkerStyle(20);
   Graph_h_dataFail4->SetMarkerSize(1.2);
   Graph_h_dataFail4->GetXaxis()->SetNdivisions(505);
   Graph_h_dataFail4->GetXaxis()->SetLabelFont(42);
   Graph_h_dataFail4->GetXaxis()->SetLabelSize(0.05);
   Graph_h_dataFail4->GetXaxis()->SetTitleSize(0.055);
   Graph_h_dataFail4->GetXaxis()->SetTitleOffset(1.2);
   Graph_h_dataFail4->GetXaxis()->SetTitleFont(42);
   Graph_h_dataFail4->GetYaxis()->SetLabelFont(42);
   Graph_h_dataFail4->GetYaxis()->SetLabelOffset(0.01);
   Graph_h_dataFail4->GetYaxis()->SetLabelSize(0.05);
   Graph_h_dataFail4->GetYaxis()->SetTitleSize(0.055);
   Graph_h_dataFail4->GetYaxis()->SetTitleOffset(1.4);
   Graph_h_dataFail4->GetYaxis()->SetTitleFont(42);
   Graph_h_dataFail4->GetZaxis()->SetLabelFont(42);
   Graph_h_dataFail4->GetZaxis()->SetLabelSize(0.035);
   Graph_h_dataFail4->GetZaxis()->SetTitleSize(0.035);
   Graph_h_dataFail4->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_h_dataFail4);
   
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
   graph->SetPoint(7,40.00000816,0.1505410461);
   graph->SetPoint(8,59.99999184,0.1505410461);
   graph->SetPoint(9,60.00000816,3.577985092);
   graph->SetPoint(10,79.99999184,3.577985092);
   graph->SetPoint(11,80.00000816,9.366074396);
   graph->SetPoint(12,99.99999184,9.366074396);
   graph->SetPoint(13,100.0000082,32.31009689);
   graph->SetPoint(14,119.9999918,32.31009689);
   graph->SetPoint(15,120.0000082,75.86025849);
   graph->SetPoint(16,139.9999918,75.86025849);
   graph->SetPoint(17,140.0000082,101.0385432);
   graph->SetPoint(18,159.9999918,101.0385432);
   graph->SetPoint(19,160.0000082,123.4656972);
   graph->SetPoint(20,179.9999918,123.4656972);
   graph->SetPoint(21,180.0000082,168.4847322);
   graph->SetPoint(22,199.9999918,168.4847322);
   graph->SetPoint(23,200.0000082,183.6496413);
   graph->SetPoint(24,219.9999918,183.6496413);
   graph->SetPoint(25,220.0000082,202.1650762);
   graph->SetPoint(26,239.9999918,202.1650762);
   graph->SetPoint(27,240.0000082,172.8766889);
   graph->SetPoint(28,259.9999918,172.8766889);
   graph->SetPoint(29,260.0000082,167.0056017);
   graph->SetPoint(30,279.9999918,167.0056017);
   graph->SetPoint(31,280.0000082,144.0465016);
   graph->SetPoint(32,299.9999918,144.0465016);
   graph->SetPoint(33,300.0000082,124.7095902);
   graph->SetPoint(34,319.9999918,124.7095902);
   graph->SetPoint(35,320.0000082,126.5580874);
   graph->SetPoint(36,339.9999918,126.5580874);
   graph->SetPoint(37,340.0000082,96.28863815);
   graph->SetPoint(38,359.9999918,96.28863815);
   graph->SetPoint(39,360.0000082,90.63701273);
   graph->SetPoint(40,379.9999918,90.63701273);
   graph->SetPoint(41,380.0000082,79.88376068);
   graph->SetPoint(42,399.9999918,79.88376068);
   graph->SetPoint(43,400.0000082,71.18585141);
   graph->SetPoint(44,419.9999918,71.18585141);
   graph->SetPoint(45,420.0000082,66.43001488);
   graph->SetPoint(46,439.9999918,66.43001488);
   graph->SetPoint(47,440.0000082,54.10562027);
   graph->SetPoint(48,459.9999918,54.10562027);
   graph->SetPoint(49,460.0000082,50.70258688);
   graph->SetPoint(50,479.9999918,50.70258688);
   graph->SetPoint(51,480.0000082,49.81436619);
   graph->SetPoint(52,499.9999918,49.81436619);
   graph->SetPoint(53,500.0000082,38.86125999);
   graph->SetPoint(54,519.9999918,38.86125999);
   graph->SetPoint(55,520.0000082,30.47450555);
   graph->SetPoint(56,539.9999918,30.47450555);
   graph->SetPoint(57,540.0000082,35.38685769);
   graph->SetPoint(58,559.9999918,35.38685769);
   graph->SetPoint(59,560.0000082,33.62002243);
   graph->SetPoint(60,579.9999918,33.62002243);
   graph->SetPoint(61,580.0000082,33.40335981);
   graph->SetPoint(62,599.9999918,33.40335981);
   graph->SetPoint(63,600.0000082,16.9974042);
   graph->SetPoint(64,619.9999918,16.9974042);
   graph->SetPoint(65,620.0000082,21.45412276);
   graph->SetPoint(66,639.9999918,21.45412276);
   graph->SetPoint(67,640.0000082,17.26624571);
   graph->SetPoint(68,659.9999918,17.26624571);
   graph->SetPoint(69,660.0000082,19.63938188);
   graph->SetPoint(70,679.9999918,19.63938188);
   graph->SetPoint(71,680.0000082,15.28108316);
   graph->SetPoint(72,699.9999918,15.28108316);
   graph->SetPoint(73,700.0000082,8.594170696);
   graph->SetPoint(74,719.9999918,8.594170696);
   graph->SetPoint(75,720.0000082,13.99378368);
   graph->SetPoint(76,739.9999918,13.99378368);
   graph->SetPoint(77,740.0000082,10.77294773);
   graph->SetPoint(78,759.9999918,10.77294773);
   graph->SetPoint(79,760.0000082,11.46834929);
   graph->SetPoint(80,779.9999918,11.46834929);
   graph->SetPoint(81,780.0000082,7.556911573);
   graph->SetPoint(82,799.9999918,7.556911573);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]13 = new TH1F("Graph_modelFail_Norm[m]13","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]13->SetMinimum(0);
   Graph_modelFail_Norm[m]13->SetMaximum(222.3816);
   Graph_modelFail_Norm[m]13->SetDirectory(0);
   Graph_modelFail_Norm[m]13->SetStats(0);
   Graph_modelFail_Norm[m]13->SetFillColor(2);
   Graph_modelFail_Norm[m]13->SetFillStyle(0);
   Graph_modelFail_Norm[m]13->SetLineStyle(0);
   Graph_modelFail_Norm[m]13->SetLineWidth(2);
   Graph_modelFail_Norm[m]13->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]13->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]13->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]13->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]13->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]13->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]13->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]13->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]13);
   
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
   graph->SetPoint(11,80.00000816,0.06669849766);
   graph->SetPoint(12,99.99999184,0.06669849766);
   graph->SetPoint(13,100.0000082,0.2667939906);
   graph->SetPoint(14,119.9999918,0.2667939906);
   graph->SetPoint(15,120.0000082,0.9171043429);
   graph->SetPoint(16,139.9999918,0.9171043429);
   graph->SetPoint(17,140.0000082,3.28490101);
   graph->SetPoint(18,159.9999918,3.28490101);
   graph->SetPoint(19,160.0000082,2.551217536);
   graph->SetPoint(20,179.9999918,2.551217536);
   graph->SetPoint(21,180.0000082,2.551217536);
   graph->SetPoint(22,199.9999918,2.551217536);
   graph->SetPoint(23,200.0000082,1.333969953);
   graph->SetPoint(24,219.9999918,1.333969953);
   graph->SetPoint(25,220.0000082,0.9337789673);
   graph->SetPoint(26,239.9999918,0.9337789673);
   graph->SetPoint(27,240.0000082,0.7670327231);
   graph->SetPoint(28,259.9999918,0.7670327231);
   graph->SetPoint(29,260.0000082,0.2667939906);
   graph->SetPoint(30,279.9999918,0.2667939906);
   graph->SetPoint(31,280.0000082,0.2334447418);
   graph->SetPoint(32,299.9999918,0.2334447418);
   graph->SetPoint(33,300.0000082,0.1333969953);
   graph->SetPoint(34,319.9999918,0.1333969953);
   graph->SetPoint(35,320.0000082,0.1167223709);
   graph->SetPoint(36,339.9999918,0.1167223709);
   graph->SetPoint(37,340.0000082,0.1667462442);
   graph->SetPoint(38,359.9999918,0.1667462442);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,0.08337312208);
   graph->SetPoint(42,399.9999918,0.08337312208);
   graph->SetPoint(43,400.0000082,0.01667462442);
   graph->SetPoint(44,419.9999918,0.01667462442);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0.01667462442);
   graph->SetPoint(48,459.9999918,0.01667462442);
   graph->SetPoint(49,460.0000082,0.05002387325);
   graph->SetPoint(50,479.9999918,0.05002387325);
   graph->SetPoint(51,480.0000082,0.03334924883);
   graph->SetPoint(52,499.9999918,0.03334924883);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,0.03334924883);
   graph->SetPoint(56,539.9999918,0.03334924883);
   graph->SetPoint(57,540.0000082,0.01667462442);
   graph->SetPoint(58,559.9999918,0.01667462442);
   graph->SetPoint(59,560.0000082,0.01667462442);
   graph->SetPoint(60,579.9999918,0.01667462442);
   graph->SetPoint(61,580.0000082,0);
   graph->SetPoint(62,599.9999918,0);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0.01667462442);
   graph->SetPoint(66,639.9999918,0.01667462442);
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
   
   TH1F *Graph_modelFail_Norm[m]_Comp[sigModFail]14 = new TH1F("Graph_modelFail_Norm[m]_Comp[sigModFail]14","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMaximum(3.613391);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[sigModFail]14);
   
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
   graph->SetPoint(7,40.00000816,0.1505410461);
   graph->SetPoint(8,59.99999184,0.1505410461);
   graph->SetPoint(9,60.00000816,0.6774347074);
   graph->SetPoint(10,79.99999184,0.6774347074);
   graph->SetPoint(11,80.00000816,3.086091445);
   graph->SetPoint(12,99.99999184,3.086091445);
   graph->SetPoint(13,100.0000082,14.97883409);
   graph->SetPoint(14,119.9999918,14.97883409);
   graph->SetPoint(15,120.0000082,36.35566263);
   graph->SetPoint(16,139.9999918,36.35566263);
   graph->SetPoint(17,140.0000082,47.94732318);
   graph->SetPoint(18,159.9999918,47.94732318);
   graph->SetPoint(19,160.0000082,62.32399308);
   graph->SetPoint(20,179.9999918,62.32399308);
   graph->SetPoint(21,180.0000082,84.75460895);
   graph->SetPoint(22,199.9999918,84.75460895);
   graph->SetPoint(23,200.0000082,102.744264);
   graph->SetPoint(24,219.9999918,102.744264);
   graph->SetPoint(25,220.0000082,106.3572491);
   graph->SetPoint(26,239.9999918,106.3572491);
   graph->SetPoint(27,240.0000082,92.95909597);
   graph->SetPoint(28,259.9999918,92.95909597);
   graph->SetPoint(29,260.0000082,77.45336822);
   graph->SetPoint(30,279.9999918,77.45336822);
   graph->SetPoint(31,280.0000082,69.2488812);
   graph->SetPoint(32,299.9999918,69.2488812);
   graph->SetPoint(33,300.0000082,56.00126915);
   graph->SetPoint(34,319.9999918,56.00126915);
   graph->SetPoint(35,320.0000082,59.99060687);
   graph->SetPoint(36,339.9999918,59.99060687);
   graph->SetPoint(37,340.0000082,49.37746312);
   graph->SetPoint(38,359.9999918,49.37746312);
   graph->SetPoint(39,360.0000082,43.05473918);
   graph->SetPoint(40,379.9999918,43.05473918);
   graph->SetPoint(41,380.0000082,38.38796675);
   graph->SetPoint(42,399.9999918,38.38796675);
   graph->SetPoint(43,400.0000082,34.32335851);
   graph->SetPoint(44,419.9999918,34.32335851);
   graph->SetPoint(45,420.0000082,30.78564393);
   graph->SetPoint(46,439.9999918,30.78564393);
   graph->SetPoint(47,440.0000082,22.43061587);
   graph->SetPoint(48,459.9999918,22.43061587);
   graph->SetPoint(49,460.0000082,19.11871285);
   graph->SetPoint(50,479.9999918,19.11871285);
   graph->SetPoint(51,480.0000082,22.12953378);
   graph->SetPoint(52,499.9999918,22.12953378);
   graph->SetPoint(53,500.0000082,18.59181919);
   graph->SetPoint(54,519.9999918,18.59181919);
   graph->SetPoint(55,520.0000082,12.49490683);
   graph->SetPoint(56,539.9999918,12.49490683);
   graph->SetPoint(57,540.0000082,18.51654867);
   graph->SetPoint(58,559.9999918,18.51654867);
   graph->SetPoint(59,560.0000082,13.84977624);
   graph->SetPoint(60,579.9999918,13.84977624);
   graph->SetPoint(61,580.0000082,14.45194043);
   graph->SetPoint(62,599.9999918,14.45194043);
   graph->SetPoint(63,600.0000082,6.397994459);
   graph->SetPoint(64,619.9999918,6.397994459);
   graph->SetPoint(65,620.0000082,6.397994459);
   graph->SetPoint(66,639.9999918,6.397994459);
   graph->SetPoint(67,640.0000082,7.602322828);
   graph->SetPoint(68,659.9999918,7.602322828);
   graph->SetPoint(69,660.0000082,8.88192172);
   graph->SetPoint(70,679.9999918,8.88192172);
   graph->SetPoint(71,680.0000082,6.17218289);
   graph->SetPoint(72,699.9999918,6.17218289);
   graph->SetPoint(73,700.0000082,4.967854521);
   graph->SetPoint(74,719.9999918,4.967854521);
   graph->SetPoint(75,720.0000082,6.623806028);
   graph->SetPoint(76,739.9999918,6.623806028);
   graph->SetPoint(77,740.0000082,3.537714583);
   graph->SetPoint(78,759.9999918,3.537714583);
   graph->SetPoint(79,760.0000082,3.161361968);
   graph->SetPoint(80,779.9999918,3.161361968);
   graph->SetPoint(81,780.0000082,3.914067199);
   graph->SetPoint(82,799.9999918,3.914067199);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMaximum(116.993);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15);
   
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
   graph->SetPoint(9,60.00000816,2.900550384);
   graph->SetPoint(10,79.99999184,2.900550384);
   graph->SetPoint(11,80.00000816,6.213284453);
   graph->SetPoint(12,99.99999184,6.213284453);
   graph->SetPoint(13,100.0000082,17.06446881);
   graph->SetPoint(14,119.9999918,17.06446881);
   graph->SetPoint(15,120.0000082,38.58749152);
   graph->SetPoint(16,139.9999918,38.58749152);
   graph->SetPoint(17,140.0000082,49.80631898);
   graph->SetPoint(18,159.9999918,49.80631898);
   graph->SetPoint(19,160.0000082,58.59048658);
   graph->SetPoint(20,179.9999918,58.59048658);
   graph->SetPoint(21,180.0000082,81.17890571);
   graph->SetPoint(22,199.9999918,81.17890571);
   graph->SetPoint(23,200.0000082,79.57140737);
   graph->SetPoint(24,219.9999918,79.57140737);
   graph->SetPoint(25,220.0000082,94.87404816);
   graph->SetPoint(26,239.9999918,94.87404816);
   graph->SetPoint(27,240.0000082,79.15056022);
   graph->SetPoint(28,259.9999918,79.15056022);
   graph->SetPoint(29,260.0000082,89.28543947);
   graph->SetPoint(30,279.9999918,89.28543947);
   graph->SetPoint(31,280.0000082,74.56417566);
   graph->SetPoint(32,299.9999918,74.56417566);
   graph->SetPoint(33,300.0000082,68.57492405);
   graph->SetPoint(34,319.9999918,68.57492405);
   graph->SetPoint(35,320.0000082,66.45075819);
   graph->SetPoint(36,339.9999918,66.45075819);
   graph->SetPoint(37,340.0000082,46.74442878);
   graph->SetPoint(38,359.9999918,46.74442878);
   graph->SetPoint(39,360.0000082,47.58227355);
   graph->SetPoint(40,379.9999918,47.58227355);
   graph->SetPoint(41,380.0000082,41.4124208);
   graph->SetPoint(42,399.9999918,41.4124208);
   graph->SetPoint(43,400.0000082,36.84581828);
   graph->SetPoint(44,419.9999918,36.84581828);
   graph->SetPoint(45,420.0000082,35.64437095);
   graph->SetPoint(46,439.9999918,35.64437095);
   graph->SetPoint(47,440.0000082,31.65832978);
   graph->SetPoint(48,459.9999918,31.65832978);
   graph->SetPoint(49,460.0000082,31.53385015);
   graph->SetPoint(50,479.9999918,31.53385015);
   graph->SetPoint(51,480.0000082,27.65148316);
   graph->SetPoint(52,499.9999918,27.65148316);
   graph->SetPoint(53,500.0000082,20.2694408);
   graph->SetPoint(54,519.9999918,20.2694408);
   graph->SetPoint(55,520.0000082,17.94624948);
   graph->SetPoint(56,539.9999918,17.94624948);
   graph->SetPoint(57,540.0000082,16.85363439);
   graph->SetPoint(58,559.9999918,16.85363439);
   graph->SetPoint(59,560.0000082,19.75357156);
   graph->SetPoint(60,579.9999918,19.75357156);
   graph->SetPoint(61,580.0000082,18.95141938);
   graph->SetPoint(62,599.9999918,18.95141938);
   graph->SetPoint(63,600.0000082,10.59940974);
   graph->SetPoint(64,619.9999918,10.59940974);
   graph->SetPoint(65,620.0000082,15.03945367);
   graph->SetPoint(66,639.9999918,15.03945367);
   graph->SetPoint(67,640.0000082,9.663922877);
   graph->SetPoint(68,659.9999918,9.663922877);
   graph->SetPoint(69,660.0000082,10.75746016);
   graph->SetPoint(70,679.9999918,10.75746016);
   graph->SetPoint(71,680.0000082,9.108900267);
   graph->SetPoint(72,699.9999918,9.108900267);
   graph->SetPoint(73,700.0000082,3.626316175);
   graph->SetPoint(74,719.9999918,3.626316175);
   graph->SetPoint(75,720.0000082,7.369977653);
   graph->SetPoint(76,739.9999918,7.369977653);
   graph->SetPoint(77,740.0000082,7.235233148);
   graph->SetPoint(78,759.9999918,7.235233148);
   graph->SetPoint(79,760.0000082,8.306987327);
   graph->SetPoint(80,779.9999918,8.306987327);
   graph->SetPoint(81,780.0000082,3.642844374);
   graph->SetPoint(82,799.9999918,3.642844374);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMaximum(104.3615);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetDirectory(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetStats(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetFillColor(2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetFillStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetLineStyle(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetLineWidth(2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMarkerStyle(20);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMarkerSize(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetNdivisions(505);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetTitleOffset(1.2);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetXaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetLabelOffset(0.01);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetLabelSize(0.05);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetTitleSize(0.055);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetTitleOffset(1.4);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetYaxis()->SetTitleFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetZaxis()->SetLabelFont(42);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetZaxis()->SetLabelSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetZaxis()->SetTitleSize(0.035);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16);
   
   graph->Draw("l");
   
   TH1D *frame_1f1d040__8 = new TH1D("frame_1f1d040__8","",40,0,800);
   frame_1f1d040__8->SetMaximum(223.2091);
   frame_1f1d040__8->SetDirectory(0);
   frame_1f1d040__8->SetStats(0);
   frame_1f1d040__8->SetFillColor(2);
   frame_1f1d040__8->SetFillStyle(0);
   frame_1f1d040__8->SetLineStyle(0);
   frame_1f1d040__8->SetLineWidth(2);
   frame_1f1d040__8->SetMarkerStyle(20);
   frame_1f1d040__8->SetMarkerSize(1.2);
   frame_1f1d040__8->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_1f1d040__8->GetXaxis()->SetNdivisions(505);
   frame_1f1d040__8->GetXaxis()->SetLabelFont(42);
   frame_1f1d040__8->GetXaxis()->SetLabelSize(0.05);
   frame_1f1d040__8->GetXaxis()->SetTitleSize(0.055);
   frame_1f1d040__8->GetXaxis()->SetTitleOffset(1.2);
   frame_1f1d040__8->GetXaxis()->SetTitleFont(42);
   frame_1f1d040__8->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_1f1d040__8->GetYaxis()->SetLabelFont(42);
   frame_1f1d040__8->GetYaxis()->SetLabelOffset(0.01);
   frame_1f1d040__8->GetYaxis()->SetLabelSize(0.05);
   frame_1f1d040__8->GetYaxis()->SetTitleSize(0.055);
   frame_1f1d040__8->GetYaxis()->SetTitleOffset(1.4);
   frame_1f1d040__8->GetYaxis()->SetTitleFont(42);
   frame_1f1d040__8->GetZaxis()->SetLabelFont(42);
   frame_1f1d040__8->GetZaxis()->SetLabelSize(0.035);
   frame_1f1d040__8->GetZaxis()->SetTitleSize(0.035);
   frame_1f1d040__8->GetZaxis()->SetTitleFont(42);
   frame_1f1d040__8->Draw("AXISSAME");
   
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
   cfail->Modified();
   cfail->cd();
   cfail->SetSelected(cfail);
}
