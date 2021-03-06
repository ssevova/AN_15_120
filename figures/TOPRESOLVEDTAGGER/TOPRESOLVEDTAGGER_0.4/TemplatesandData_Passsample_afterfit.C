{
//=========Macro generated from canvas: cpass/cpass
//=========  (Mon Oct 26 23:59:18 2015) by ROOT version5.32/00
   TCanvas *cpass = new TCanvas("cpass", "cpass",1,1,800,1176);
   gStyle->SetOptStat(0);
   cpass->SetHighLightColor(2);
   cpass->Range(-187.013,-308.3254,851.9481,1747.177);
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
   
   TH1D *frame_1ea0be0__5 = new TH1D("frame_1ea0be0__5","",40,0,800);
   frame_1ea0be0__5->SetMaximum(1582.737);
   frame_1ea0be0__5->SetDirectory(0);
   frame_1ea0be0__5->SetStats(0);
   frame_1ea0be0__5->SetFillColor(2);
   frame_1ea0be0__5->SetFillStyle(0);
   frame_1ea0be0__5->SetLineStyle(0);
   frame_1ea0be0__5->SetLineWidth(2);
   frame_1ea0be0__5->SetMarkerStyle(20);
   frame_1ea0be0__5->SetMarkerSize(1.2);
   frame_1ea0be0__5->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_1ea0be0__5->GetXaxis()->SetNdivisions(505);
   frame_1ea0be0__5->GetXaxis()->SetLabelFont(42);
   frame_1ea0be0__5->GetXaxis()->SetLabelSize(0.05);
   frame_1ea0be0__5->GetXaxis()->SetTitleSize(0.055);
   frame_1ea0be0__5->GetXaxis()->SetTitleOffset(1.2);
   frame_1ea0be0__5->GetXaxis()->SetTitleFont(42);
   frame_1ea0be0__5->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_1ea0be0__5->GetYaxis()->SetLabelFont(42);
   frame_1ea0be0__5->GetYaxis()->SetLabelOffset(0.01);
   frame_1ea0be0__5->GetYaxis()->SetLabelSize(0.05);
   frame_1ea0be0__5->GetYaxis()->SetTitleSize(0.055);
   frame_1ea0be0__5->GetYaxis()->SetTitleOffset(1.4);
   frame_1ea0be0__5->GetYaxis()->SetTitleFont(42);
   frame_1ea0be0__5->GetZaxis()->SetLabelFont(42);
   frame_1ea0be0__5->GetZaxis()->SetLabelSize(0.035);
   frame_1ea0be0__5->GetZaxis()->SetTitleSize(0.035);
   frame_1ea0be0__5->GetZaxis()->SetTitleFont(42);
   frame_1ea0be0__5->Draw("");
   
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
   grae->SetPoint(3,70,18);
   grae->SetPointError(3,10,10,4.202887,5.321007);
   grae->SetPoint(4,90,101);
   grae->SetPointError(4,10,10,9.562306,10.56231);
   grae->SetPoint(5,110,229);
   grae->SetPointError(5,10,10,14.641,15.641);
   grae->SetPoint(6,130,397);
   grae->SetPointError(6,10,10,19.43113,20.43113);
   grae->SetPoint(7,150,637);
   grae->SetPointError(7,10,10,24.74381,25.74381);
   grae->SetPoint(8,170,1442);
   grae->SetPointError(8,10,10,37.47697,38.47697);
   grae->SetPoint(9,190,1211);
   grae->SetPointError(9,10,10,34.30302,35.30302);
   grae->SetPoint(10,210,455);
   grae->SetPointError(10,10,10,20.83659,21.83659);
   grae->SetPoint(11,230,185);
   grae->SetPointError(11,10,10,13.11066,14.11066);
   grae->SetPoint(12,250,96);
   grae->SetPointError(12,10,10,9.780867,10.83201);
   grae->SetPoint(13,270,61);
   grae->SetPointError(13,10,10,7.788779,8.852952);
   grae->SetPoint(14,290,32);
   grae->SetPointError(14,10,10,5.627135,6.715753);
   grae->SetPoint(15,310,21);
   grae->SetPointError(15,10,10,4.545807,5.655182);
   grae->SetPoint(16,330,21);
   grae->SetPointError(16,10,10,4.545807,5.655182);
   grae->SetPoint(17,350,16);
   grae->SetPointError(17,10,10,3.957801,5.083066);
   grae->SetPoint(18,370,14);
   grae->SetPointError(18,10,10,3.6965,4.830381);
   grae->SetPoint(19,390,8);
   grae->SetPointError(19,10,10,2.768386,3.945142);
   grae->SetPoint(20,410,2);
   grae->SetPointError(20,10,10,1.291815,2.63786);
   grae->SetPoint(21,430,6);
   grae->SetPointError(21,10,10,2.379931,3.583642);
   grae->SetPoint(22,450,7);
   grae->SetPointError(22,10,10,2.58147,3.770281);
   grae->SetPoint(23,470,4);
   grae->SetPointError(23,10,10,1.914339,3.162753);
   grae->SetPoint(24,490,0);
   grae->SetPointError(24,10,10,0,1.147874);
   grae->SetPoint(25,510,2);
   grae->SetPointError(25,10,10,1.291815,2.63786);
   grae->SetPoint(26,530,2);
   grae->SetPointError(26,10,10,1.291815,2.63786);
   grae->SetPoint(27,550,0);
   grae->SetPointError(27,10,10,0,1.147874);
   grae->SetPoint(28,570,0);
   grae->SetPointError(28,10,10,0,1.147874);
   grae->SetPoint(29,590,1);
   grae->SetPointError(29,10,10,0.8272462,2.299527);
   grae->SetPoint(30,610,1);
   grae->SetPointError(30,10,10,0.8272462,2.299527);
   grae->SetPoint(31,630,0);
   grae->SetPointError(31,10,10,0,1.147874);
   grae->SetPoint(32,650,0);
   grae->SetPointError(32,10,10,0,1.147874);
   grae->SetPoint(33,670,1);
   grae->SetPointError(33,10,10,0.8272462,2.299527);
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
   Graph_h_dataPass3->SetMaximum(1628.525);
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
   graph->SetPoint(7,40.00000816,0.8775731429);
   graph->SetPoint(8,59.99999184,0.8775731429);
   graph->SetPoint(9,60.00000816,24.20712018);
   graph->SetPoint(10,79.99999184,24.20712018);
   graph->SetPoint(11,80.00000816,102.7412376);
   graph->SetPoint(12,99.99999184,102.7412376);
   graph->SetPoint(13,100.0000082,208.4095048);
   graph->SetPoint(14,119.9999918,208.4095048);
   graph->SetPoint(15,120.0000082,360.6134348);
   graph->SetPoint(16,139.9999918,360.6134348);
   graph->SetPoint(17,140.0000082,683.3385256);
   graph->SetPoint(18,159.9999918,683.3385256);
   graph->SetPoint(19,160.0000082,1507.368663);
   graph->SetPoint(20,179.9999918,1507.368663);
   graph->SetPoint(21,180.0000082,1141.41951);
   graph->SetPoint(22,199.9999918,1141.41951);
   graph->SetPoint(23,200.0000082,450.4767187);
   graph->SetPoint(24,219.9999918,450.4767187);
   graph->SetPoint(25,220.0000082,187.8653868);
   graph->SetPoint(26,239.9999918,187.8653868);
   graph->SetPoint(27,240.0000082,96.91586714);
   graph->SetPoint(28,259.9999918,96.91586714);
   graph->SetPoint(29,260.0000082,59.13588398);
   graph->SetPoint(30,279.9999918,59.13588398);
   graph->SetPoint(31,280.0000082,40.14289474);
   graph->SetPoint(32,299.9999918,40.14289474);
   graph->SetPoint(33,300.0000082,23.78172257);
   graph->SetPoint(34,319.9999918,23.78172257);
   graph->SetPoint(35,320.0000082,17.0596045);
   graph->SetPoint(36,339.9999918,17.0596045);
   graph->SetPoint(37,340.0000082,15.0423229);
   graph->SetPoint(38,359.9999918,15.0423229);
   graph->SetPoint(39,360.0000082,9.572708764);
   graph->SetPoint(40,379.9999918,9.572708764);
   graph->SetPoint(41,380.0000082,9.080722306);
   graph->SetPoint(42,399.9999918,9.080722306);
   graph->SetPoint(43,400.0000082,9.422568167);
   graph->SetPoint(44,419.9999918,9.422568167);
   graph->SetPoint(45,420.0000082,6.259035335);
   graph->SetPoint(46,439.9999918,6.259035335);
   graph->SetPoint(47,440.0000082,4.395706061);
   graph->SetPoint(48,459.9999918,4.395706061);
   graph->SetPoint(49,460.0000082,3.441269518);
   graph->SetPoint(50,479.9999918,3.441269518);
   graph->SetPoint(51,480.0000082,3.169862739);
   graph->SetPoint(52,499.9999918,3.169862739);
   graph->SetPoint(53,500.0000082,1.842622165);
   graph->SetPoint(54,519.9999918,1.842622165);
   graph->SetPoint(55,520.0000082,1.817891319);
   graph->SetPoint(56,539.9999918,1.817891319);
   graph->SetPoint(57,540.0000082,1.321175915);
   graph->SetPoint(58,559.9999918,1.321175915);
   graph->SetPoint(59,560.0000082,0.931173528);
   graph->SetPoint(60,579.9999918,0.931173528);
   graph->SetPoint(61,580.0000082,1.135192975);
   graph->SetPoint(62,599.9999918,1.135192975);
   graph->SetPoint(63,600.0000082,0.1573062421);
   graph->SetPoint(64,619.9999918,0.1573062421);
   graph->SetPoint(65,620.0000082,0.3729884145);
   graph->SetPoint(66,639.9999918,0.3729884145);
   graph->SetPoint(67,640.0000082,0.6977919225);
   graph->SetPoint(68,659.9999918,0.6977919225);
   graph->SetPoint(69,660.0000082,0.3069619953);
   graph->SetPoint(70,679.9999918,0.3069619953);
   graph->SetPoint(71,680.0000082,0.3141608497);
   graph->SetPoint(72,699.9999918,0.3141608497);
   graph->SetPoint(73,700.0000082,0.0626548458);
   graph->SetPoint(74,719.9999918,0.0626548458);
   graph->SetPoint(75,720.0000082,0.1628234099);
   graph->SetPoint(76,739.9999918,0.1628234099);
   graph->SetPoint(77,740.0000082,0.0767583907);
   graph->SetPoint(78,759.9999918,0.0767583907);
   graph->SetPoint(79,760.0000082,0.0626548458);
   graph->SetPoint(80,779.9999918,0.0626548458);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]9 = new TH1F("Graph_modelPass_Norm[m]9","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]9->SetMinimum(0);
   Graph_modelPass_Norm[m]9->SetMaximum(1658.106);
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
   graph->SetPoint(11,80.00000816,0.5825002505);
   graph->SetPoint(12,99.99999184,0.5825002505);
   graph->SetPoint(13,100.0000082,2.065111721);
   graph->SetPoint(14,119.9999918,2.065111721);
   graph->SetPoint(15,120.0000082,17.82080798);
   graph->SetPoint(16,139.9999918,17.82080798);
   graph->SetPoint(17,140.0000082,151.4832615);
   graph->SetPoint(18,159.9999918,151.4832615);
   graph->SetPoint(19,160.0000082,648.4131358);
   graph->SetPoint(20,179.9999918,648.4131358);
   graph->SetPoint(21,180.0000082,582.8082638);
   graph->SetPoint(22,199.9999918,582.8082638);
   graph->SetPoint(23,200.0000082,150.2610829);
   graph->SetPoint(24,219.9999918,150.2610829);
   graph->SetPoint(25,220.0000082,27.05898725);
   graph->SetPoint(26,239.9999918,27.05898725);
   graph->SetPoint(27,240.0000082,5.265726107);
   graph->SetPoint(28,259.9999918,5.265726107);
   graph->SetPoint(29,260.0000082,2.467116203);
   graph->SetPoint(30,279.9999918,2.467116203);
   graph->SetPoint(31,280.0000082,0.905105198);
   graph->SetPoint(32,299.9999918,0.905105198);
   graph->SetPoint(33,300.0000082,0.341373917);
   graph->SetPoint(34,319.9999918,0.341373917);
   graph->SetPoint(35,320.0000082,0.06038689355);
   graph->SetPoint(36,339.9999918,0.06038689355);
   graph->SetPoint(37,340.0000082,0.1850745298);
   graph->SetPoint(38,359.9999918,0.1850745298);
   graph->SetPoint(39,360.0000082,0.15367547);
   graph->SetPoint(40,379.9999918,0.15367547);
   graph->SetPoint(41,380.0000082,0.1355492142);
   graph->SetPoint(42,399.9999918,0.1355492142);
   graph->SetPoint(43,400.0000082,0);
   graph->SetPoint(44,419.9999918,0);
   graph->SetPoint(45,420.0000082,0);
   graph->SetPoint(46,439.9999918,0);
   graph->SetPoint(47,440.0000082,0.07709844488);
   graph->SetPoint(48,459.9999918,0.07709844488);
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
   Graph_modelPass_Norm[m]_Comp[sigModPass]10->SetMaximum(713.2544);
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
   graph->SetPoint(7,40.00000816,0.4994814026);
   graph->SetPoint(8,59.99999184,0.4994814026);
   graph->SetPoint(9,60.00000816,13.57862234);
   graph->SetPoint(10,79.99999184,13.57862234);
   graph->SetPoint(11,80.00000816,68.55764996);
   graph->SetPoint(12,99.99999184,68.55764996);
   graph->SetPoint(13,100.0000082,147.1469095);
   graph->SetPoint(14,119.9999918,147.1469095);
   graph->SetPoint(15,120.0000082,253.8987491);
   graph->SetPoint(16,139.9999918,253.8987491);
   graph->SetPoint(17,140.0000082,413.1194663);
   graph->SetPoint(18,159.9999918,413.1194663);
   graph->SetPoint(19,160.0000082,659.9017839);
   graph->SetPoint(20,179.9999918,659.9017839);
   graph->SetPoint(21,180.0000082,425.6351006);
   graph->SetPoint(22,199.9999918,425.6351006);
   graph->SetPoint(23,200.0000082,231.556189);
   graph->SetPoint(24,219.9999918,231.556189);
   graph->SetPoint(25,220.0000082,122.0853459);
   graph->SetPoint(26,239.9999918,122.0853459);
   graph->SetPoint(27,240.0000082,70.47975252);
   graph->SetPoint(28,259.9999918,70.47975252);
   graph->SetPoint(29,260.0000082,45.46473091);
   graph->SetPoint(30,279.9999918,45.46473091);
   graph->SetPoint(31,280.0000082,28.73735761);
   graph->SetPoint(32,299.9999918,28.73735761);
   graph->SetPoint(33,300.0000082,16.37211045);
   graph->SetPoint(34,319.9999918,16.37211045);
   graph->SetPoint(35,320.0000082,12.35599848);
   graph->SetPoint(36,339.9999918,12.35599848);
   graph->SetPoint(37,340.0000082,11.08044867);
   graph->SetPoint(38,359.9999918,11.08044867);
   graph->SetPoint(39,360.0000082,7.041707638);
   graph->SetPoint(40,379.9999918,7.041707638);
   graph->SetPoint(41,380.0000082,6.185887472);
   graph->SetPoint(42,399.9999918,6.185887472);
   graph->SetPoint(43,400.0000082,5.946082099);
   graph->SetPoint(44,419.9999918,5.946082099);
   graph->SetPoint(45,420.0000082,4.758143016);
   graph->SetPoint(46,439.9999918,4.758143016);
   graph->SetPoint(47,440.0000082,3.399389801);
   graph->SetPoint(48,459.9999918,3.399389801);
   graph->SetPoint(49,460.0000082,2.768248304);
   graph->SetPoint(50,479.9999918,2.768248304);
   graph->SetPoint(51,480.0000082,2.742196252);
   graph->SetPoint(52,499.9999918,2.742196252);
   graph->SetPoint(53,500.0000082,1.615561044);
   graph->SetPoint(54,519.9999918,1.615561044);
   graph->SetPoint(55,520.0000082,1.071901067);
   graph->SetPoint(56,539.9999918,1.071901067);
   graph->SetPoint(57,540.0000082,0.8602777445);
   graph->SetPoint(58,559.9999918,0.8602777445);
   graph->SetPoint(59,560.0000082,0.7370027534);
   graph->SetPoint(60,579.9999918,0.7370027534);
   graph->SetPoint(61,580.0000082,0.5980922345);
   graph->SetPoint(62,599.9999918,0.5980922345);
   graph->SetPoint(63,600.0000082,0.07623432762);
   graph->SetPoint(64,619.9999918,0.07623432762);
   graph->SetPoint(65,620.0000082,0.3143862481);
   graph->SetPoint(66,639.9999918,0.3143862481);
   graph->SetPoint(67,640.0000082,0.3582676806);
   graph->SetPoint(68,659.9999918,0.3582676806);
   graph->SetPoint(69,660.0000082,0.3069619953);
   graph->SetPoint(70,679.9999918,0.3069619953);
   graph->SetPoint(71,680.0000082,0.2920814494);
   graph->SetPoint(72,699.9999918,0.2920814494);
   graph->SetPoint(73,700.0000082,0.0626548458);
   graph->SetPoint(74,719.9999918,0.0626548458);
   graph->SetPoint(75,720.0000082,0.08588908402);
   graph->SetPoint(76,739.9999918,0.08588908402);
   graph->SetPoint(77,740.0000082,0.0767583907);
   graph->SetPoint(78,759.9999918,0.0767583907);
   graph->SetPoint(79,760.0000082,0.0626548458);
   graph->SetPoint(80,779.9999918,0.0626548458);
   graph->SetPoint(81,780.0000082,0);
   graph->SetPoint(82,799.9999918,0);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg1ModPass]11->SetMaximum(725.892);
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
   graph->SetPoint(7,40.00000816,0.3780917403);
   graph->SetPoint(8,59.99999184,0.3780917403);
   graph->SetPoint(9,60.00000816,10.62849784);
   graph->SetPoint(10,79.99999184,10.62849784);
   graph->SetPoint(11,80.00000816,33.60108734);
   graph->SetPoint(12,99.99999184,33.60108734);
   graph->SetPoint(13,100.0000082,59.19748357);
   graph->SetPoint(14,119.9999918,59.19748357);
   graph->SetPoint(15,120.0000082,88.89387772);
   graph->SetPoint(16,139.9999918,88.89387772);
   graph->SetPoint(17,140.0000082,118.7357978);
   graph->SetPoint(18,159.9999918,118.7357978);
   graph->SetPoint(19,160.0000082,199.0537433);
   graph->SetPoint(20,179.9999918,199.0537433);
   graph->SetPoint(21,180.0000082,132.9761454);
   graph->SetPoint(22,199.9999918,132.9761454);
   graph->SetPoint(23,200.0000082,68.65944684);
   graph->SetPoint(24,219.9999918,68.65944684);
   graph->SetPoint(25,220.0000082,38.72105369);
   graph->SetPoint(26,239.9999918,38.72105369);
   graph->SetPoint(27,240.0000082,21.17038851);
   graph->SetPoint(28,259.9999918,21.17038851);
   graph->SetPoint(29,260.0000082,11.20403687);
   graph->SetPoint(30,279.9999918,11.20403687);
   graph->SetPoint(31,280.0000082,10.50043194);
   graph->SetPoint(32,299.9999918,10.50043194);
   graph->SetPoint(33,300.0000082,7.0682382);
   graph->SetPoint(34,319.9999918,7.0682382);
   graph->SetPoint(35,320.0000082,4.643219123);
   graph->SetPoint(36,339.9999918,4.643219123);
   graph->SetPoint(37,340.0000082,3.7767997);
   graph->SetPoint(38,359.9999918,3.7767997);
   graph->SetPoint(39,360.0000082,2.377325656);
   graph->SetPoint(40,379.9999918,2.377325656);
   graph->SetPoint(41,380.0000082,2.75928562);
   graph->SetPoint(42,399.9999918,2.75928562);
   graph->SetPoint(43,400.0000082,3.476486068);
   graph->SetPoint(44,419.9999918,3.476486068);
   graph->SetPoint(45,420.0000082,1.500892319);
   graph->SetPoint(46,439.9999918,1.500892319);
   graph->SetPoint(47,440.0000082,0.9192178155);
   graph->SetPoint(48,459.9999918,0.9192178155);
   graph->SetPoint(49,460.0000082,0.6730212136);
   graph->SetPoint(50,479.9999918,0.6730212136);
   graph->SetPoint(51,480.0000082,0.4276664871);
   graph->SetPoint(52,499.9999918,0.4276664871);
   graph->SetPoint(53,500.0000082,0.2270611211);
   graph->SetPoint(54,519.9999918,0.2270611211);
   graph->SetPoint(55,520.0000082,0.7459902524);
   graph->SetPoint(56,539.9999918,0.7459902524);
   graph->SetPoint(57,540.0000082,0.4608981703);
   graph->SetPoint(58,559.9999918,0.4608981703);
   graph->SetPoint(59,560.0000082,0.1941707746);
   graph->SetPoint(60,579.9999918,0.1941707746);
   graph->SetPoint(61,580.0000082,0.5371007402);
   graph->SetPoint(62,599.9999918,0.5371007402);
   graph->SetPoint(63,600.0000082,0.08107191447);
   graph->SetPoint(64,619.9999918,0.08107191447);
   graph->SetPoint(65,620.0000082,0.05860216645);
   graph->SetPoint(66,639.9999918,0.05860216645);
   graph->SetPoint(67,640.0000082,0.3395242418);
   graph->SetPoint(68,659.9999918,0.3395242418);
   graph->SetPoint(69,660.0000082,0);
   graph->SetPoint(70,679.9999918,0);
   graph->SetPoint(71,680.0000082,0.02207940033);
   graph->SetPoint(72,699.9999918,0.02207940033);
   graph->SetPoint(73,700.0000082,0);
   graph->SetPoint(74,719.9999918,0);
   graph->SetPoint(75,720.0000082,0.07693432587);
   graph->SetPoint(76,739.9999918,0.07693432587);
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
   
   TH1F *Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12 = new TH1F("Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12","Projection of Model for PASS sample",100,-91.85185,891.8519);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMinimum(0);
   Graph_modelPass_Norm[m]_Comp[bkg2ModPass]12->SetMaximum(218.9591);
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
   
   TH1D *frame_1ea0be0__6 = new TH1D("frame_1ea0be0__6","",40,0,800);
   frame_1ea0be0__6->SetMaximum(1582.737);
   frame_1ea0be0__6->SetDirectory(0);
   frame_1ea0be0__6->SetStats(0);
   frame_1ea0be0__6->SetFillColor(2);
   frame_1ea0be0__6->SetFillStyle(0);
   frame_1ea0be0__6->SetLineStyle(0);
   frame_1ea0be0__6->SetLineWidth(2);
   frame_1ea0be0__6->SetMarkerStyle(20);
   frame_1ea0be0__6->SetMarkerSize(1.2);
   frame_1ea0be0__6->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_1ea0be0__6->GetXaxis()->SetNdivisions(505);
   frame_1ea0be0__6->GetXaxis()->SetLabelFont(42);
   frame_1ea0be0__6->GetXaxis()->SetLabelSize(0.05);
   frame_1ea0be0__6->GetXaxis()->SetTitleSize(0.055);
   frame_1ea0be0__6->GetXaxis()->SetTitleOffset(1.2);
   frame_1ea0be0__6->GetXaxis()->SetTitleFont(42);
   frame_1ea0be0__6->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_1ea0be0__6->GetYaxis()->SetLabelFont(42);
   frame_1ea0be0__6->GetYaxis()->SetLabelOffset(0.01);
   frame_1ea0be0__6->GetYaxis()->SetLabelSize(0.05);
   frame_1ea0be0__6->GetYaxis()->SetTitleSize(0.055);
   frame_1ea0be0__6->GetYaxis()->SetTitleOffset(1.4);
   frame_1ea0be0__6->GetYaxis()->SetTitleFont(42);
   frame_1ea0be0__6->GetZaxis()->SetLabelFont(42);
   frame_1ea0be0__6->GetZaxis()->SetLabelSize(0.035);
   frame_1ea0be0__6->GetZaxis()->SetTitleSize(0.035);
   frame_1ea0be0__6->GetZaxis()->SetTitleFont(42);
   frame_1ea0be0__6->Draw("AXISSAME");
   
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
