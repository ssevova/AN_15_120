{
//=========Macro generated from canvas: cfail/cfail
//=========  (Sat Oct 31 01:56:56 2015) by ROOT version5.32/00
   TCanvas *cfail = new TCanvas("cfail", "cfail",1,1,800,1176);
   gStyle->SetOptStat(0);
   cfail->SetHighLightColor(2);
   cfail->Range(-187.013,-18.62122,851.9481,105.5202);
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
   
   TH1D *frame_31ed210__7 = new TH1D("frame_31ed210__7","",40,0,800);
   frame_31ed210__7->SetMaximum(95.58892);
   frame_31ed210__7->SetDirectory(0);
   frame_31ed210__7->SetStats(0);
   frame_31ed210__7->SetFillColor(2);
   frame_31ed210__7->SetFillStyle(0);
   frame_31ed210__7->SetLineStyle(0);
   frame_31ed210__7->SetLineWidth(2);
   frame_31ed210__7->SetMarkerStyle(20);
   frame_31ed210__7->SetMarkerSize(1.2);
   frame_31ed210__7->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_31ed210__7->GetXaxis()->SetNdivisions(505);
   frame_31ed210__7->GetXaxis()->SetLabelFont(42);
   frame_31ed210__7->GetXaxis()->SetLabelSize(0.05);
   frame_31ed210__7->GetXaxis()->SetTitleSize(0.055);
   frame_31ed210__7->GetXaxis()->SetTitleOffset(1.2);
   frame_31ed210__7->GetXaxis()->SetTitleFont(42);
   frame_31ed210__7->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_31ed210__7->GetYaxis()->SetLabelFont(42);
   frame_31ed210__7->GetYaxis()->SetLabelOffset(0.01);
   frame_31ed210__7->GetYaxis()->SetLabelSize(0.05);
   frame_31ed210__7->GetYaxis()->SetTitleSize(0.055);
   frame_31ed210__7->GetYaxis()->SetTitleOffset(1.4);
   frame_31ed210__7->GetYaxis()->SetTitleFont(42);
   frame_31ed210__7->GetZaxis()->SetLabelFont(42);
   frame_31ed210__7->GetZaxis()->SetLabelSize(0.035);
   frame_31ed210__7->GetZaxis()->SetTitleSize(0.035);
   frame_31ed210__7->GetZaxis()->SetTitleFont(42);
   frame_31ed210__7->Draw("");
   
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
   grae->SetPoint(3,70,0);
   grae->SetPointError(3,10,10,0,1.147874);
   grae->SetPoint(4,90,1);
   grae->SetPointError(4,10,10,0.8272462,2.299527);
   grae->SetPoint(5,110,1);
   grae->SetPointError(5,10,10,0.8272462,2.299527);
   grae->SetPoint(6,130,6);
   grae->SetPointError(6,10,10,2.379931,3.583642);
   grae->SetPoint(7,150,10);
   grae->SetPointError(7,10,10,3.108694,4.26695);
   grae->SetPoint(8,170,10);
   grae->SetPointError(8,10,10,3.108694,4.26695);
   grae->SetPoint(9,190,25);
   grae->SetPointError(9,10,10,4.966335,6.066589);
   grae->SetPoint(10,210,35);
   grae->SetPointError(10,10,10,5.887675,6.97241);
   grae->SetPoint(11,230,56);
   grae->SetPointError(11,10,10,7.4609,8.527879);
   grae->SetPoint(12,250,81);
   grae->SetPointError(12,10,10,8.981384,10.03706);
   grae->SetPoint(13,270,59);
   grae->SetPointError(13,10,10,7.659312,8.724565);
   grae->SetPoint(14,290,54);
   grae->SetPointError(14,10,10,7.32564,8.39385);
   grae->SetPoint(15,310,56);
   grae->SetPointError(15,10,10,7.4609,8.527879);
   grae->SetPoint(16,330,53);
   grae->SetPointError(16,10,10,7.257065,8.325916);
   grae->SetPoint(17,350,50);
   grae->SetPointError(17,10,10,7.047337,8.118225);
   grae->SetPoint(18,370,41);
   grae->SetPointError(18,10,10,6.376898,7.455185);
   grae->SetPoint(19,390,36);
   grae->SetPointError(19,10,10,5.971996,7.055545);
   grae->SetPoint(20,410,35);
   grae->SetPointError(20,10,10,5.887675,6.97241);
   grae->SetPoint(21,430,30);
   grae->SetPointError(21,10,10,5.446522,6.538046);
   grae->SetPoint(22,450,30);
   grae->SetPointError(22,10,10,5.446522,6.538046);
   grae->SetPoint(23,470,34);
   grae->SetPointError(23,10,10,5.802128,6.888101);
   grae->SetPoint(24,490,24);
   grae->SetPointError(24,10,10,4.864612,5.966932);
   grae->SetPoint(25,510,20);
   grae->SetPointError(25,10,10,4.434448,5.546519);
   grae->SetPoint(26,530,16);
   grae->SetPointError(26,10,10,3.957801,5.083066);
   grae->SetPoint(27,550,18);
   grae->SetPointError(27,10,10,4.202887,5.321007);
   grae->SetPoint(28,570,13);
   grae->SetPointError(28,10,10,3.558662,4.697573);
   grae->SetPoint(29,590,15);
   grae->SetPointError(29,10,10,3.82938,4.958738);
   grae->SetPoint(30,610,14);
   grae->SetPointError(30,10,10,3.6965,4.830381);
   grae->SetPoint(31,630,12);
   grae->SetPointError(31,10,10,3.415266,4.559819);
   grae->SetPoint(32,650,8);
   grae->SetPointError(32,10,10,2.768386,3.945142);
   grae->SetPoint(33,670,9);
   grae->SetPointError(33,10,10,2.943461,4.110204);
   grae->SetPoint(34,690,8);
   grae->SetPointError(34,10,10,2.768386,3.945142);
   grae->SetPoint(35,710,9);
   grae->SetPointError(35,10,10,2.943461,4.110204);
   grae->SetPoint(36,730,6);
   grae->SetPointError(36,10,10,2.379931,3.583642);
   grae->SetPoint(37,750,6);
   grae->SetPointError(37,10,10,2.379931,3.583642);
   grae->SetPoint(38,770,7);
   grae->SetPointError(38,10,10,2.58147,3.770281);
   grae->SetPoint(39,790,4);
   grae->SetPointError(39,10,10,1.914339,3.162753);
   
   TH1F *Graph_h_dataFail4 = new TH1F("Graph_h_dataFail4","Histogram of dataFail_plot__m",100,0,880);
   Graph_h_dataFail4->SetMinimum(0);
   Graph_h_dataFail4->SetMaximum(100.1408);
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
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0.2866595086);
   graph->SetPoint(10,79.99999184,0.2866595086);
   graph->SetPoint(11,80.00000816,1.116908799);
   graph->SetPoint(12,99.99999184,1.116908799);
   graph->SetPoint(13,100.0000082,4.264716418);
   graph->SetPoint(14,119.9999918,4.264716418);
   graph->SetPoint(15,120.0000082,11.98710458);
   graph->SetPoint(16,139.9999918,11.98710458);
   graph->SetPoint(17,140.0000082,15.53654191);
   graph->SetPoint(18,159.9999918,15.53654191);
   graph->SetPoint(19,160.0000082,14.29749879);
   graph->SetPoint(20,179.9999918,14.29749879);
   graph->SetPoint(21,180.0000082,27.62399962);
   graph->SetPoint(22,199.9999918,27.62399962);
   graph->SetPoint(23,200.0000082,39.39614471);
   graph->SetPoint(24,219.9999918,39.39614471);
   graph->SetPoint(25,220.0000082,61.78883203);
   graph->SetPoint(26,239.9999918,61.78883203);
   graph->SetPoint(27,240.0000082,61.6166895);
   graph->SetPoint(28,259.9999918,61.6166895);
   graph->SetPoint(29,260.0000082,69.69361505);
   graph->SetPoint(30,279.9999918,69.69361505);
   graph->SetPoint(31,280.0000082,61.01180288);
   graph->SetPoint(32,299.9999918,61.01180288);
   graph->SetPoint(33,300.0000082,59.4678729);
   graph->SetPoint(34,319.9999918,59.4678729);
   graph->SetPoint(35,320.0000082,56.41033334);
   graph->SetPoint(36,339.9999918,56.41033334);
   graph->SetPoint(37,340.0000082,39.48451728);
   graph->SetPoint(38,359.9999918,39.48451728);
   graph->SetPoint(39,360.0000082,42.74387643);
   graph->SetPoint(40,379.9999918,42.74387643);
   graph->SetPoint(41,380.0000082,37.74194191);
   graph->SetPoint(42,399.9999918,37.74194191);
   graph->SetPoint(43,400.0000082,31.83657762);
   graph->SetPoint(44,419.9999918,31.83657762);
   graph->SetPoint(45,420.0000082,26.81269213);
   graph->SetPoint(46,439.9999918,26.81269213);
   graph->SetPoint(47,440.0000082,24.66092316);
   graph->SetPoint(48,459.9999918,24.66092316);
   graph->SetPoint(49,460.0000082,26.92184141);
   graph->SetPoint(50,479.9999918,26.92184141);
   graph->SetPoint(51,480.0000082,23.90231128);
   graph->SetPoint(52,499.9999918,23.90231128);
   graph->SetPoint(53,500.0000082,17.93495433);
   graph->SetPoint(54,519.9999918,17.93495433);
   graph->SetPoint(55,520.0000082,14.20238559);
   graph->SetPoint(56,539.9999918,14.20238559);
   graph->SetPoint(57,540.0000082,16.52401908);
   graph->SetPoint(58,559.9999918,16.52401908);
   graph->SetPoint(59,560.0000082,18.65836939);
   graph->SetPoint(60,579.9999918,18.65836939);
   graph->SetPoint(61,580.0000082,15.74411931);
   graph->SetPoint(62,599.9999918,15.74411931);
   graph->SetPoint(63,600.0000082,6.114505129);
   graph->SetPoint(64,619.9999918,6.114505129);
   graph->SetPoint(65,620.0000082,11.22537796);
   graph->SetPoint(66,639.9999918,11.22537796);
   graph->SetPoint(67,640.0000082,9.290455037);
   graph->SetPoint(68,659.9999918,9.290455037);
   graph->SetPoint(69,660.0000082,9.136603555);
   graph->SetPoint(70,679.9999918,9.136603555);
   graph->SetPoint(71,680.0000082,6.576257849);
   graph->SetPoint(72,699.9999918,6.576257849);
   graph->SetPoint(73,700.0000082,4.76868515);
   graph->SetPoint(74,719.9999918,4.76868515);
   graph->SetPoint(75,720.0000082,7.65375881);
   graph->SetPoint(76,739.9999918,7.65375881);
   graph->SetPoint(77,740.0000082,5.036859954);
   graph->SetPoint(78,759.9999918,5.036859954);
   graph->SetPoint(79,760.0000082,5.814454597);
   graph->SetPoint(80,779.9999918,5.814454597);
   graph->SetPoint(81,780.0000082,4.71582731);
   graph->SetPoint(82,799.9999918,4.71582731);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]13 = new TH1F("Graph_modelFail_Norm[m]13","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]13->SetMinimum(0);
   Graph_modelFail_Norm[m]13->SetMaximum(76.66298);
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
   graph->SetPoint(11,80.00000816,0);
   graph->SetPoint(12,99.99999184,0);
   graph->SetPoint(13,100.0000082,1.143253104e-05);
   graph->SetPoint(14,119.9999918,1.143253104e-05);
   graph->SetPoint(15,120.0000082,3.429759311e-05);
   graph->SetPoint(16,139.9999918,3.429759311e-05);
   graph->SetPoint(17,140.0000082,0.0001028927793);
   graph->SetPoint(18,159.9999918,0.0001028927793);
   graph->SetPoint(19,160.0000082,4.573012414e-05);
   graph->SetPoint(20,179.9999918,4.573012414e-05);
   graph->SetPoint(21,180.0000082,0.0001028927793);
   graph->SetPoint(22,199.9999918,0.0001028927793);
   graph->SetPoint(23,200.0000082,2.286506207e-05);
   graph->SetPoint(24,219.9999918,2.286506207e-05);
   graph->SetPoint(25,220.0000082,0.0001257578414);
   graph->SetPoint(26,239.9999918,0.0001257578414);
   graph->SetPoint(27,240.0000082,0.0001257578414);
   graph->SetPoint(28,259.9999918,0.0001257578414);
   graph->SetPoint(29,260.0000082,5.716265518e-05);
   graph->SetPoint(30,279.9999918,5.716265518e-05);
   graph->SetPoint(31,280.0000082,1.143253104e-05);
   graph->SetPoint(32,299.9999918,1.143253104e-05);
   graph->SetPoint(33,300.0000082,1.143253104e-05);
   graph->SetPoint(34,319.9999918,1.143253104e-05);
   graph->SetPoint(35,320.0000082,3.429759311e-05);
   graph->SetPoint(36,339.9999918,3.429759311e-05);
   graph->SetPoint(37,340.0000082,4.573012414e-05);
   graph->SetPoint(38,359.9999918,4.573012414e-05);
   graph->SetPoint(39,360.0000082,0);
   graph->SetPoint(40,379.9999918,0);
   graph->SetPoint(41,380.0000082,1.143253104e-05);
   graph->SetPoint(42,399.9999918,1.143253104e-05);
   graph->SetPoint(43,400.0000082,2.286506207e-05);
   graph->SetPoint(44,419.9999918,2.286506207e-05);
   graph->SetPoint(45,420.0000082,1.143253104e-05);
   graph->SetPoint(46,439.9999918,1.143253104e-05);
   graph->SetPoint(47,440.0000082,0);
   graph->SetPoint(48,459.9999918,0);
   graph->SetPoint(49,460.0000082,0);
   graph->SetPoint(50,479.9999918,0);
   graph->SetPoint(51,480.0000082,0);
   graph->SetPoint(52,499.9999918,0);
   graph->SetPoint(53,500.0000082,0);
   graph->SetPoint(54,519.9999918,0);
   graph->SetPoint(55,520.0000082,1.143253104e-05);
   graph->SetPoint(56,539.9999918,1.143253104e-05);
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
   
   TH1F *Graph_modelFail_Norm[m]_Comp[sigModFail]14 = new TH1F("Graph_modelFail_Norm[m]_Comp[sigModFail]14","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[sigModFail]14->SetMaximum(0.0001383336);
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
   graph->SetPoint(7,40.00000816,0);
   graph->SetPoint(8,59.99999184,0);
   graph->SetPoint(9,60.00000816,0);
   graph->SetPoint(10,79.99999184,0);
   graph->SetPoint(11,80.00000816,0.1613771035);
   graph->SetPoint(12,99.99999184,0.1613771035);
   graph->SetPoint(13,100.0000082,0.6455084142);
   graph->SetPoint(14,119.9999918,0.6455084142);
   graph->SetPoint(15,120.0000082,1.936525243);
   graph->SetPoint(16,139.9999918,1.936525243);
   graph->SetPoint(17,140.0000082,3.442711542);
   graph->SetPoint(18,159.9999918,3.442711542);
   graph->SetPoint(19,160.0000082,4.948897842);
   graph->SetPoint(20,179.9999918,4.948897842);
   graph->SetPoint(21,180.0000082,8.552986488);
   graph->SetPoint(22,199.9999918,8.552986488);
   graph->SetPoint(23,200.0000082,14.20118511);
   graph->SetPoint(24,219.9999918,14.20118511);
   graph->SetPoint(25,220.0000082,19.41904479);
   graph->SetPoint(26,239.9999918,19.41904479);
   graph->SetPoint(27,240.0000082,21.30177767);
   graph->SetPoint(28,259.9999918,21.30177767);
   graph->SetPoint(29,260.0000082,21.89349371);
   graph->SetPoint(30,279.9999918,21.89349371);
   graph->SetPoint(31,280.0000082,17.05218061);
   graph->SetPoint(32,299.9999918,17.05218061);
   graph->SetPoint(33,300.0000082,15.76116378);
   graph->SetPoint(34,319.9999918,15.76116378);
   graph->SetPoint(35,320.0000082,15.54599431);
   graph->SetPoint(36,339.9999918,15.54599431);
   graph->SetPoint(37,340.0000082,12.69499881);
   graph->SetPoint(38,359.9999918,12.69499881);
   graph->SetPoint(39,360.0000082,13.39429959);
   graph->SetPoint(40,379.9999918,13.39429959);
   graph->SetPoint(41,380.0000082,10.43571936);
   graph->SetPoint(42,399.9999918,10.43571936);
   graph->SetPoint(43,400.0000082,7.584723867);
   graph->SetPoint(44,419.9999918,7.584723867);
   graph->SetPoint(45,420.0000082,7.315762027);
   graph->SetPoint(46,439.9999918,7.315762027);
   graph->SetPoint(47,440.0000082,5.379236785);
   graph->SetPoint(48,459.9999918,5.379236785);
   graph->SetPoint(49,460.0000082,6.293707038);
   graph->SetPoint(50,479.9999918,6.293707038);
   graph->SetPoint(51,480.0000082,5.594406256);
   graph->SetPoint(52,499.9999918,5.594406256);
   graph->SetPoint(53,500.0000082,3.980635221);
   graph->SetPoint(54,519.9999918,3.980635221);
   graph->SetPoint(55,520.0000082,3.0123726);
   graph->SetPoint(56,539.9999918,3.0123726);
   graph->SetPoint(57,540.0000082,5.110274946);
   graph->SetPoint(58,559.9999918,5.110274946);
   graph->SetPoint(59,560.0000082,4.518558899);
   graph->SetPoint(60,579.9999918,4.518558899);
   graph->SetPoint(61,580.0000082,3.49650391);
   graph->SetPoint(62,599.9999918,3.49650391);
   graph->SetPoint(63,600.0000082,1.291016828);
   graph->SetPoint(64,619.9999918,1.291016828);
   graph->SetPoint(65,620.0000082,1.882732875);
   graph->SetPoint(66,639.9999918,1.882732875);
   graph->SetPoint(67,640.0000082,2.097902346);
   graph->SetPoint(68,659.9999918,2.097902346);
   graph->SetPoint(69,660.0000082,3.066164967);
   graph->SetPoint(70,679.9999918,3.066164967);
   graph->SetPoint(71,680.0000082,1.183432093);
   graph->SetPoint(72,699.9999918,1.183432093);
   graph->SetPoint(73,700.0000082,1.237224461);
   graph->SetPoint(74,719.9999918,1.237224461);
   graph->SetPoint(75,720.0000082,2.635826025);
   graph->SetPoint(76,739.9999918,2.635826025);
   graph->SetPoint(77,740.0000082,0.9682626213);
   graph->SetPoint(78,759.9999918,0.9682626213);
   graph->SetPoint(79,760.0000082,0.699300782);
   graph->SetPoint(80,779.9999918,0.699300782);
   graph->SetPoint(81,780.0000082,1.721355771);
   graph->SetPoint(82,799.9999918,1.721355771);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg1ModFail]15->SetMaximum(24.08284);
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
   graph->SetPoint(9,60.00000816,0.2866595086);
   graph->SetPoint(10,79.99999184,0.2866595086);
   graph->SetPoint(11,80.00000816,0.9555316955);
   graph->SetPoint(12,99.99999184,0.9555316955);
   graph->SetPoint(13,100.0000082,3.619196571);
   graph->SetPoint(14,119.9999918,3.619196571);
   graph->SetPoint(15,120.0000082,10.05054504);
   graph->SetPoint(16,139.9999918,10.05054504);
   graph->SetPoint(17,140.0000082,12.09372747);
   graph->SetPoint(18,159.9999918,12.09372747);
   graph->SetPoint(19,160.0000082,9.348555216);
   graph->SetPoint(20,179.9999918,9.348555216);
   graph->SetPoint(21,180.0000082,19.07091024);
   graph->SetPoint(22,199.9999918,19.07091024);
   graph->SetPoint(23,200.0000082,25.19493673);
   graph->SetPoint(24,219.9999918,25.19493673);
   graph->SetPoint(25,220.0000082,42.36966148);
   graph->SetPoint(26,239.9999918,42.36966148);
   graph->SetPoint(27,240.0000082,40.31478607);
   graph->SetPoint(28,259.9999918,40.31478607);
   graph->SetPoint(29,260.0000082,47.80006417);
   graph->SetPoint(30,279.9999918,47.80006417);
   graph->SetPoint(31,280.0000082,43.95961084);
   graph->SetPoint(32,299.9999918,43.95961084);
   graph->SetPoint(33,300.0000082,43.70669768);
   graph->SetPoint(34,319.9999918,43.70669768);
   graph->SetPoint(35,320.0000082,40.86430473);
   graph->SetPoint(36,339.9999918,40.86430473);
   graph->SetPoint(37,340.0000082,26.78947273);
   graph->SetPoint(38,359.9999918,26.78947273);
   graph->SetPoint(39,360.0000082,29.34957684);
   graph->SetPoint(40,379.9999918,29.34957684);
   graph->SetPoint(41,380.0000082,27.30621112);
   graph->SetPoint(42,399.9999918,27.30621112);
   graph->SetPoint(43,400.0000082,24.25183089);
   graph->SetPoint(44,419.9999918,24.25183089);
   graph->SetPoint(45,420.0000082,19.49691867);
   graph->SetPoint(46,439.9999918,19.49691867);
   graph->SetPoint(47,440.0000082,19.28168637);
   graph->SetPoint(48,459.9999918,19.28168637);
   graph->SetPoint(49,460.0000082,20.62813437);
   graph->SetPoint(50,479.9999918,20.62813437);
   graph->SetPoint(51,480.0000082,18.30790503);
   graph->SetPoint(52,499.9999918,18.30790503);
   graph->SetPoint(53,500.0000082,13.95431911);
   graph->SetPoint(54,519.9999918,13.95431911);
   graph->SetPoint(55,520.0000082,11.19000156);
   graph->SetPoint(56,539.9999918,11.19000156);
   graph->SetPoint(57,540.0000082,11.41374414);
   graph->SetPoint(58,559.9999918,11.41374414);
   graph->SetPoint(59,560.0000082,14.13981049);
   graph->SetPoint(60,579.9999918,14.13981049);
   graph->SetPoint(61,580.0000082,12.2476154);
   graph->SetPoint(62,599.9999918,12.2476154);
   graph->SetPoint(63,600.0000082,4.823488301);
   graph->SetPoint(64,619.9999918,4.823488301);
   graph->SetPoint(65,620.0000082,9.342645085);
   graph->SetPoint(66,639.9999918,9.342645085);
   graph->SetPoint(67,640.0000082,7.192552691);
   graph->SetPoint(68,659.9999918,7.192552691);
   graph->SetPoint(69,660.0000082,6.070438588);
   graph->SetPoint(70,679.9999918,6.070438588);
   graph->SetPoint(71,680.0000082,5.392825756);
   graph->SetPoint(72,699.9999918,5.392825756);
   graph->SetPoint(73,700.0000082,3.53146069);
   graph->SetPoint(74,719.9999918,3.53146069);
   graph->SetPoint(75,720.0000082,5.017932785);
   graph->SetPoint(76,739.9999918,5.017932785);
   graph->SetPoint(77,740.0000082,4.068597332);
   graph->SetPoint(78,759.9999918,4.068597332);
   graph->SetPoint(79,760.0000082,5.115153815);
   graph->SetPoint(80,779.9999918,5.115153815);
   graph->SetPoint(81,780.0000082,2.994471539);
   graph->SetPoint(82,799.9999918,2.994471539);
   graph->SetPoint(83,800.0000082,0);
   graph->SetPoint(84,800,0);
   graph->SetPoint(85,809.8765432,0);
   graph->SetPoint(86,809.8765432,0);
   
   TH1F *Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16 = new TH1F("Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16","Projection of Model for FAIL sample",100,-91.85185,891.8519);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMinimum(0);
   Graph_modelFail_Norm[m]_Comp[bkg2ModFail]16->SetMaximum(52.58007);
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
   
   TH1D *frame_31ed210__8 = new TH1D("frame_31ed210__8","",40,0,800);
   frame_31ed210__8->SetMaximum(95.58892);
   frame_31ed210__8->SetDirectory(0);
   frame_31ed210__8->SetStats(0);
   frame_31ed210__8->SetFillColor(2);
   frame_31ed210__8->SetFillStyle(0);
   frame_31ed210__8->SetLineStyle(0);
   frame_31ed210__8->SetLineWidth(2);
   frame_31ed210__8->SetMarkerStyle(20);
   frame_31ed210__8->SetMarkerSize(1.2);
   frame_31ed210__8->GetXaxis()->SetTitle("Top Mass [GeV/c^2]");
   frame_31ed210__8->GetXaxis()->SetNdivisions(505);
   frame_31ed210__8->GetXaxis()->SetLabelFont(42);
   frame_31ed210__8->GetXaxis()->SetLabelSize(0.05);
   frame_31ed210__8->GetXaxis()->SetTitleSize(0.055);
   frame_31ed210__8->GetXaxis()->SetTitleOffset(1.2);
   frame_31ed210__8->GetXaxis()->SetTitleFont(42);
   frame_31ed210__8->GetYaxis()->SetTitle("Events / ( 20 )");
   frame_31ed210__8->GetYaxis()->SetLabelFont(42);
   frame_31ed210__8->GetYaxis()->SetLabelOffset(0.01);
   frame_31ed210__8->GetYaxis()->SetLabelSize(0.05);
   frame_31ed210__8->GetYaxis()->SetTitleSize(0.055);
   frame_31ed210__8->GetYaxis()->SetTitleOffset(1.4);
   frame_31ed210__8->GetYaxis()->SetTitleFont(42);
   frame_31ed210__8->GetZaxis()->SetLabelFont(42);
   frame_31ed210__8->GetZaxis()->SetLabelSize(0.035);
   frame_31ed210__8->GetZaxis()->SetTitleSize(0.035);
   frame_31ed210__8->GetZaxis()->SetTitleFont(42);
   frame_31ed210__8->Draw("AXISSAME");
   
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
