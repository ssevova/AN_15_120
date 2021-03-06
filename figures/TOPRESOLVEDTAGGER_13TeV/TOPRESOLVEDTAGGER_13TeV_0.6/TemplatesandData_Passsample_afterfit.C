{
//=========Macro generated from canvas: cpass/cpass
//=========  (Sat Oct 31 02:03:01 2015) by ROOT version5.32/00
   TCanvas *cpass = new TCanvas("cpass", "cpass",1,1,800,1176);
   gStyle->SetOptStat(0);
   cpass->SetHighLightColor(2);
   cpass->Range(-187.013,-121.16,851.9481,686.5736);
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
   
   TH1D *frame_2986f50__5 = new TH1D("frame_2986f50__5","",40,0,800);
   frame_2986f50__5->SetMaximum(621.9549);
   frame_2986f50__5->SetDirectory(0);
   frame_2986f50__5->SetStats(0);
   frame_2986f50__5->SetFillColor(2);
   frame_2986f50__5->SetFillStyle(0);
   frame_2986f50__5->SetLineStyle(0);
   frame_2986f50__5->SetLineWidth(2);
   frame_2986f50__5->SetMarkerStyle(20);
   frame_2986f50__5->SetMarkerSize(1.2);
   frame_2986f50__5->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_2986f50__5->GetXaxis()->SetNdivisions(505);
   frame_2986f50__5->GetXaxis()->SetLabelFont(42);
   frame_2986f50__5->GetXaxis()->SetLabelSize(0.05);
   frame_2986f50__5->GetXaxis()->SetTitleSize(0.055);
   frame_2986f50__5->GetXaxis()->SetTitleOffset(1.2);
   frame_2986f50__5->GetXaxis()->SetTitleFont(42);
   frame_2986f50__5->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_2986f50__5->GetYaxis()->SetLabelFont(42);
   frame_2986f50__5->GetYaxis()->SetLabelOffset(0.01);
   frame_2986f50__5->GetYaxis()->SetLabelSize(0.05);
   frame_2986f50__5->GetYaxis()->SetTitleSize(0.055);
   frame_2986f50__5->GetYaxis()->SetTitleOffset(1.4);
   frame_2986f50__5->GetYaxis()->SetTitleFont(42);
   frame_2986f50__5->GetZaxis()->SetLabelFont(42);
   frame_2986f50__5->GetZaxis()->SetLabelSize(0.035);
   frame_2986f50__5->GetZaxis()->SetTitleSize(0.035);
   frame_2986f50__5->GetZaxis()->SetTitleFont(42);
   frame_2986f50__5->Draw("");
   
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
   grae->SetPoint(2,50,7);
   grae->SetPointError(2,10,10,2.58147,3.770281);
   grae->SetPoint(3,70,29);
   grae->SetPointError(3,10,10,5.353932,6.44702);
   grae->SetPoint(4,90,80);
   grae->SetPointError(4,10,10,8.925539,9.981567);
   grae->SetPoint(5,110,104);
   grae->SetPointError(5,10,10,9.710289,10.71029);
   grae->SetPoint(6,130,165);
   grae->SetPointError(6,10,10,12.35496,13.35496);
   grae->SetPoint(7,150,292);
   grae->SetPointError(7,10,10,16.59532,17.59532);
   grae->SetPoint(8,170,568);
   grae->SetPointError(8,10,10,23.33799,24.33799);
   grae->SetPoint(9,190,271);
   grae->SetPointError(9,10,10,15.96967,16.96967);
   grae->SetPoint(10,210,77);
   grae->SetPointError(10,10,10,8.755868,9.812979);
   grae->SetPoint(11,230,44);
   grae->SetPointError(11,10,10,6.60794,7.68351);
   grae->SetPoint(12,250,19);
   grae->SetPointError(12,10,10,4.320219,5.435196);
   grae->SetPoint(13,270,19);
   grae->SetPointError(13,10,10,4.320219,5.435196);
   grae->SetPoint(14,290,8);
   grae->SetPointError(14,10,10,2.768386,3.945142);
   grae->SetPoint(15,310,6);
   grae->SetPointError(15,10,10,2.379931,3.583642);
   grae->SetPoint(16,330,2);
   grae->SetPointError(16,10,10,1.291815,2.63786);
   grae->SetPoint(17,350,2);
   grae->SetPointError(17,10,10,1.291815,2.63786);
   grae->SetPoint(18,370,1);
   grae->SetPointError(18,10,10,0.8272462,2.299527);
   grae->SetPoint(19,390,2);
   grae->SetPointError(19,10,10,1.291815,2.63786);
   grae->SetPoint(20,410,3);
   grae->SetPointError(20,10,10,1.632705,2.918186);
   grae->SetPoint(21,430,1);
   grae->SetPointError(21,10,10,0.8272462,2.299527);
   grae->SetPoint(22,450,0);
   grae->SetPointError(22,10,10,0,1.147874);
   grae->SetPoint(23,470,3);
   grae->SetPointError(23,10,10,1.632705,2.918186);
   grae->SetPoint(24,490,2);
   grae->SetPointError(24,10,10,1.291815,2.63786);
   grae->SetPoint(25,510,0);
   grae->SetPointError(25,10,10,0,1.147874);
   grae->SetPoint(26,530,2);
   grae->SetPointError(26,10,10,1.291815,2.63786);
   grae->SetPoint(27,550,1);
   grae->SetPointError(27,10,10,0.8272462,2.299527);
   grae->SetPoint(28,570,0);
   grae->SetPointError(28,10,10,0,1.147874);
   grae->SetPoint(29,590,1);
   grae->SetPointError(29,10,10,0.8272462,2.299527);
   grae->SetPoint(30,610,0);
   grae->SetPointError(30,10,10,0,1.147874);
   grae->SetPoint(31,630,0);
   grae->SetPointError(31,10,10,0,1.147874);
   grae->SetPoint(32,650,0);
   grae->SetPointError(32,10,10,0,1.147874);
   grae->SetPoint(33,670,0);
   grae->SetPointError(33,10,10,0,1.147874);
   grae->SetPoint(34,690,0);
   grae->SetPointError(34,10,10,0,1.147874);
   grae->SetPoint(35,710,0);
   grae->SetPointError(35,10,10,0,1.147874);
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
   Graph_h_dataPass3->SetMaximum(651.5718);
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
   graph->SetPoint(5,20.00000816,0.02143075915);
   graph->SetPoint(6,39.99999184,0.02143075915);
   graph->SetPoint(7,40.00000816,8.814229502);
   graph->SetPoint(8,59.99999184,8.814229502);
   graph->SetPoint(9,60.00000816,37.89680618);
   graph->SetPoint(10,79.99999184,37.89680618);
   graph->SetPoint(11,80.00000816,70.09239352);
   graph->SetPoint(12,99.99999184,70.09239352);
   graph->SetPoint(13,100.0000082,105.7753442);
   graph->SetPoint(14,119.9999918,105.7753442);
   graph->SetPoint(15,120.0000082,175.0002464);
   graph->SetPoint(16,139.9999918,175.0002464);
   graph->SetPoint(17,140.0000082,329.1112077);
   graph->SetPoint(18,159.9999918,329.1112077);
   graph->SetPoint(19,160.0000082,558.5372694);
   graph->SetPoint(20,179.9999918,558.5372694);
   graph->SetPoint(21,180.0000082,242.9742112);
   graph->SetPoint(22,199.9999918,242.9742112);
   graph->SetPoint(23,200.0000082,64.67540519);
   graph->SetPoint(24,219.9999918,64.67540519);
   graph->SetPoint(25,220.0000082,47.47927598);
   graph->SetPoint(26,239.9999918,47.47927598);
   graph->SetPoint(27,240.0000082,26.26767514);
   graph->SetPoint(28,259.9999918,26.26767514);
   graph->SetPoint(29,260.0000082,11.65648736);
   graph->SetPoint(30,279.9999918,11.65648736);
   graph->SetPoint(31,280.0000082,8.790844477);
   graph->SetPoint(32,299.9999918,8.790844477);
   graph->SetPoint(33,300.0000082,4.452259728);
   graph->SetPoint(34,319.9999918,4.452259728);
   graph->SetPoint(35,320.0000082,2.743868642);
   graph->SetPoint(36,339.9999918,2.743868642);
   graph->SetPoint(37,340.0000082,3.274846948);
   graph->SetPoint(38,359.9999918,3.274846948);
   graph->SetPoint(39,360.0000082,1.424217501);
   graph->SetPoint(40,379.9999918,1.424217501);
   graph->SetPoint(41,380.0000082,2.882835681);
   graph->SetPoint(42,399.9999918,2.882835681);
   graph->SetPoint(43,400.0000082,2.245701885);
   graph->SetPoint(44,419.9999918,2.245701885);
   graph->SetPoint(45,420.0000082,0.6404012382);
   graph->SetPoint(46,439.9999918,0.6404012382);
   graph->SetPoint(47,440.0000082,0.8264211768);
   graph->SetPoint(48,459.9999918,0.8264211768);
   graph->SetPoint(49,460.0000082,0.5061565533);
   graph->SetPoint(50,479.9999918,0.5061565533);
   graph->SetPoint(51,480.0000082,0.9102728091);
   graph->SetPoint(52,499.9999918,0.9102728091);
   graph->SetPoint(53,500.0000082,0.4254963599);
   graph->SetPoint(54,519.9999918,0.4254963599);
   graph->SetPoint(55,520.0000082,0.1714460732);
   graph->SetPoint(56,539.9999918,0.1714460732);
   graph->SetPoint(57,540.0000082,0.4866698141);
   graph->SetPoint(58,559.9999918,0.4866698141);
   graph->SetPoint(59,560.0000082,0.4596919544);
   graph->SetPoint(60,579.9999918,0.4596919544);
   graph->SetPoint(61,580.0000082,0.5736960486);
   graph->SetPoint(62,599.9999918,0.5736960486);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0.02143075915);
   graph->SetPoint(66,639.9999918,0.02143075915);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0.08572303661);
   graph->SetPoint(74,719.9999918,0.08572303661);
   graph->SetPoint(75,720.0000082,0.02143075915);
   graph->SetPoint(76,739.9999918,0.02143075915);
   graph->SetPoint(77,740.0000082,0.1812536402);
   graph->SetPoint(78,759.9999918,0.1812536402);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0.1812536402);
   graph->SetPoint(82,799.9999918,0.1812536402);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]9 = new TH1F("Graph_modelPass_Norm[m]9","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]9->SetMinimum(0);
   Graph_modelPass_Norm[m]9->SetMaximum(614.391);
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
   graph->SetPoint(9,60.00000816,0.1726591741);
   graph->SetPoint(10,79.99999184,0.1726591741);
   graph->SetPoint(11,80.00000816,2.589887612);
   graph->SetPoint(12,99.99999184,2.589887612);
   graph->SetPoint(13,100.0000082,6.215730269);
   graph->SetPoint(14,119.9999918,6.215730269);
   graph->SetPoint(15,120.0000082,36.94906326);
   graph->SetPoint(16,139.9999918,36.94906326);
   graph->SetPoint(17,140.0000082,175.0764026);
   graph->SetPoint(18,159.9999918,175.0764026);
   graph->SetPoint(19,160.0000082,340.3975618);
   graph->SetPoint(20,179.9999918,340.3975618);
   graph->SetPoint(21,180.0000082,120.5161035);
   graph->SetPoint(22,199.9999918,120.5161035);
   graph->SetPoint(23,200.0000082,8.80561788);
   graph->SetPoint(24,219.9999918,8.80561788);
   graph->SetPoint(25,220.0000082,2.244569264);
   graph->SetPoint(26,239.9999918,2.244569264);
   graph->SetPoint(27,240.0000082,1.122284632);
   graph->SetPoint(28,259.9999918,1.122284632);
   graph->SetPoint(29,260.0000082,1.122284632);
   graph->SetPoint(30,279.9999918,1.122284632);
   graph->SetPoint(31,280.0000082,0.3453183483);
   graph->SetPoint(32,299.9999918,0.3453183483);
   graph->SetPoint(33,300.0000082,0.2589887612);
   graph->SetPoint(34,319.9999918,0.2589887612);
   graph->SetPoint(35,320.0000082,0);
   graph->SetPoint(36,339.9999918,0);
   graph->SetPoint(37,340.0000082,0.08632958706);
   graph->SetPoint(38,359.9999918,0.08632958706);
   graph->SetPoint(39,360.0000082,0.08632958706);
   graph->SetPoint(40,379.9999918,0.08632958706);
   graph->SetPoint(41,380.0000082,0);
   graph->SetPoint(42,399.9999918,0);
   graph->SetPoint(43,400.0000082,0.08632958706);
   graph->SetPoint(44,419.9999918,0.08632958706);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0.08632958706);
   graph->SetPoint(52,499.9999918,0.08632958706);
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
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMaximum(374.4373);
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
   graph->SetPoint(5,20.00000816,0.02143075915);
   graph->SetPoint(6,39.99999184,0.02143075915);
   graph->SetPoint(7,40.00000816,0.9429534028);
   graph->SetPoint(8,59.99999184,0.9429534028);
   graph->SetPoint(9,60.00000816,6.514950783);
   graph->SetPoint(10,79.99999184,6.514950783);
   graph->SetPoint(11,80.00000816,13.26563992);
   graph->SetPoint(12,99.99999184,13.26563992);
   graph->SetPoint(13,100.0000082,21.21645156);
   graph->SetPoint(14,119.9999918,21.21645156);
   graph->SetPoint(15,120.0000082,36.53944436);
   graph->SetPoint(16,139.9999918,36.53944436);
   graph->SetPoint(17,140.0000082,45.60465548);
   graph->SetPoint(18,159.9999918,45.60465548);
   graph->SetPoint(19,160.0000082,59.27747982);
   graph->SetPoint(20,179.9999918,59.27747982);
   graph->SetPoint(21,180.0000082,29.93877054);
   graph->SetPoint(22,199.9999918,29.93877054);
   graph->SetPoint(23,200.0000082,14.46576243);
   graph->SetPoint(24,219.9999918,14.46576243);
   graph->SetPoint(25,220.0000082,8.786611253);
   graph->SetPoint(26,239.9999918,8.786611253);
   graph->SetPoint(27,240.0000082,5.207674474);
   graph->SetPoint(28,259.9999918,5.207674474);
   graph->SetPoint(29,260.0000082,1.864476046);
   graph->SetPoint(30,279.9999918,1.864476046);
   graph->SetPoint(31,280.0000082,2.207368193);
   graph->SetPoint(32,299.9999918,2.207368193);
   graph->SetPoint(33,300.0000082,1.092968717);
   graph->SetPoint(34,319.9999918,1.092968717);
   graph->SetPoint(35,320.0000082,0.5786304971);
   graph->SetPoint(36,339.9999918,0.5786304971);
   graph->SetPoint(37,340.0000082,1.328707068);
   graph->SetPoint(38,359.9999918,1.328707068);
   graph->SetPoint(39,360.0000082,0.3214613873);
   graph->SetPoint(40,379.9999918,0.3214613873);
   graph->SetPoint(41,380.0000082,0.7500765704);
   graph->SetPoint(42,399.9999918,0.7500765704);
   graph->SetPoint(43,400.0000082,0.2143075915);
   graph->SetPoint(44,419.9999918,0.2143075915);
   graph->SetPoint(45,420.0000082,0.08572303661);
   graph->SetPoint(46,439.9999918,0.08572303661);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0.08572303661);
   graph->SetPoint(52,499.9999918,0.08572303661);
   graph->SetPoint(53,500.0000082,0.2571691098);
   graph->SetPoint(54,519.9999918,0.2571691098);
   graph->SetPoint(55,520.0000082,0.1714460732);
   graph->SetPoint(56,539.9999918,0.1714460732);
   graph->SetPoint(57,540.0000082,0.1500153141);
   graph->SetPoint(58,559.9999918,0.1500153141);
   graph->SetPoint(59,560.0000082,0);
   graph->SetPoint(60,579.9999918,0);
   graph->SetPoint(61,580.0000082,0.08572303661);
   graph->SetPoint(62,599.9999918,0.08572303661);
   graph->SetPoint(63,600.0000082,0);
   graph->SetPoint(64,619.9999918,0);
   graph->SetPoint(65,620.0000082,0.02143075915);
   graph->SetPoint(66,639.9999918,0.02143075915);
   graph->SetPoint(67,640.0000082,0);
   graph->SetPoint(68,659.9999918,0);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0);
   graph->SetPoint(72,699.9999918,0);
   graph->SetPoint(73,700.0000082,0.08572303661);
   graph->SetPoint(74,719.9999918,0.08572303661);
   graph->SetPoint(75,720.0000082,0.02143075915);
   graph->SetPoint(76,739.9999918,0.02143075915);
   graph->SetPoint(77,740.0000082,0.02143075915);
   graph->SetPoint(78,759.9999918,0.02143075915);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0.02143075915);
   graph->SetPoint(82,799.9999918,0.02143075915);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMaximum(65.20523);
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
   graph->SetPoint(7,40.00000816,7.8712761);
   graph->SetPoint(8,59.99999184,7.8712761);
   graph->SetPoint(9,60.00000816,31.20919622);
   graph->SetPoint(10,79.99999184,31.20919622);
   graph->SetPoint(11,80.00000816,54.23686599);
   graph->SetPoint(12,99.99999184,54.23686599);
   graph->SetPoint(13,100.0000082,78.34316239);
   graph->SetPoint(14,119.9999918,78.34316239);
   graph->SetPoint(15,120.0000082,101.5117388);
   graph->SetPoint(16,139.9999918,101.5117388);
   graph->SetPoint(17,140.0000082,108.4301497);
   graph->SetPoint(18,159.9999918,108.4301497);
   graph->SetPoint(19,160.0000082,158.8622278);
   graph->SetPoint(20,179.9999918,158.8622278);
   graph->SetPoint(21,180.0000082,92.51933708);
   graph->SetPoint(22,199.9999918,92.51933708);
   graph->SetPoint(23,200.0000082,41.40402488);
   graph->SetPoint(24,219.9999918,41.40402488);
   graph->SetPoint(25,220.0000082,36.44809546);
   graph->SetPoint(26,239.9999918,36.44809546);
   graph->SetPoint(27,240.0000082,19.93771604);
   graph->SetPoint(28,259.9999918,19.93771604);
   graph->SetPoint(29,260.0000082,8.669726685);
   graph->SetPoint(30,279.9999918,8.669726685);
   graph->SetPoint(31,280.0000082,6.238157935);
   graph->SetPoint(32,299.9999918,6.238157935);
   graph->SetPoint(33,300.0000082,3.10030225);
   graph->SetPoint(34,319.9999918,3.10030225);
   graph->SetPoint(35,320.0000082,2.165238144);
   graph->SetPoint(36,339.9999918,2.165238144);
   graph->SetPoint(37,340.0000082,1.859810294);
   graph->SetPoint(38,359.9999918,1.859810294);
   graph->SetPoint(39,360.0000082,1.016426526);
   graph->SetPoint(40,379.9999918,1.016426526);
   graph->SetPoint(41,380.0000082,2.132759111);
   graph->SetPoint(42,399.9999918,2.132759111);
   graph->SetPoint(43,400.0000082,1.945064706);
   graph->SetPoint(44,419.9999918,1.945064706);
   graph->SetPoint(45,420.0000082,0.5546782016);
   graph->SetPoint(46,439.9999918,0.5546782016);
   graph->SetPoint(47,440.0000082,0.8264211768);
   graph->SetPoint(48,459.9999918,0.8264211768);
   graph->SetPoint(49,460.0000082,0.5061565533);
   graph->SetPoint(50,479.9999918,0.5061565533);
   graph->SetPoint(51,480.0000082,0.7382201854);
   graph->SetPoint(52,499.9999918,0.7382201854);
   graph->SetPoint(53,500.0000082,0.16832725);
   graph->SetPoint(54,519.9999918,0.16832725);
   graph->SetPoint(55,520.0000082,0);
   graph->SetPoint(56,539.9999918,0);
   graph->SetPoint(57,540.0000082,0.3366545001);
   graph->SetPoint(58,559.9999918,0.3366545001);
   graph->SetPoint(59,560.0000082,0.4596919544);
   graph->SetPoint(60,579.9999918,0.4596919544);
   graph->SetPoint(61,580.0000082,0.487973012);
   graph->SetPoint(62,599.9999918,0.487973012);
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
   graph->SetPoint(77,740.0000082,0.159822881);
   graph->SetPoint(78,759.9999918,0.159822881);
   graph->SetPoint(79,760.0000082,0);
   graph->SetPoint(80,779.9999918,0);
   graph->SetPoint(81,780.0000082,0.159822881);
   graph->SetPoint(82,799.9999918,0.159822881);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMaximum(174.7485);
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
   
   TH1D *frame_2986f50__6 = new TH1D("frame_2986f50__6","",40,0,800);
   frame_2986f50__6->SetMaximum(621.9549);
   frame_2986f50__6->SetDirectory(0);
   frame_2986f50__6->SetStats(0);
   frame_2986f50__6->SetFillColor(2);
   frame_2986f50__6->SetFillStyle(0);
   frame_2986f50__6->SetLineStyle(0);
   frame_2986f50__6->SetLineWidth(2);
   frame_2986f50__6->SetMarkerStyle(20);
   frame_2986f50__6->SetMarkerSize(1.2);
   frame_2986f50__6->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_2986f50__6->GetXaxis()->SetNdivisions(505);
   frame_2986f50__6->GetXaxis()->SetLabelFont(42);
   frame_2986f50__6->GetXaxis()->SetLabelSize(0.05);
   frame_2986f50__6->GetXaxis()->SetTitleSize(0.055);
   frame_2986f50__6->GetXaxis()->SetTitleOffset(1.2);
   frame_2986f50__6->GetXaxis()->SetTitleFont(42);
   frame_2986f50__6->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_2986f50__6->GetYaxis()->SetLabelFont(42);
   frame_2986f50__6->GetYaxis()->SetLabelOffset(0.01);
   frame_2986f50__6->GetYaxis()->SetLabelSize(0.05);
   frame_2986f50__6->GetYaxis()->SetTitleSize(0.055);
   frame_2986f50__6->GetYaxis()->SetTitleOffset(1.4);
   frame_2986f50__6->GetYaxis()->SetTitleFont(42);
   frame_2986f50__6->GetZaxis()->SetLabelFont(42);
   frame_2986f50__6->GetZaxis()->SetLabelSize(0.035);
   frame_2986f50__6->GetZaxis()->SetTitleSize(0.035);
   frame_2986f50__6->GetZaxis()->SetTitleFont(42);
   frame_2986f50__6->Draw("AXISSAME");
   
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
