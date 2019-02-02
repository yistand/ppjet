/*
 * Draw particle level comparison plots for:
 *
 * Cases
 * 1: Pi Plus
 * 2: Pi Minus
 *
 */

TH1D *hP0_0118, *hP0_0130;
TH1D *hP0_0118_lowpt, *hP0_0130_lowpt;
TGraphErrors *gP0_0118 = new TGraphErrors();
TGraphErrors *gP0_0130 = new TGraphErrors();

TGraphErrors *gP0_0118_lowpt = new TGraphErrors();
TGraphErrors *gP0_0130_lowpt = new TGraphErrors();

Double_t pilowptbins[16]={0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
Double_t klowptbins[11]={0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6};
Double_t plowptbins[15]={0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
Double_t ptbins[15] = {3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,10.0,12.0,15.0};
Double_t pi = TMath::Pi();

void makeInvYieldRatios_Pt0Study(const Int_t iCase, const Int_t iPrint)
{
  TLatex plotText;
  plotText.SetTextFont(43);
  plotText.SetTextSize(16);
  plotText.SetNDC();

  TFile *file0_0118 = new TFile("./pt0files/Pythia.Perugia0.NominalPt0.NNPDF0118.root");
  TFile *file0_0130 = new TFile("./pt0files/Pythia.Perugia0.NominalPt0.NNPDF0130.root");

  gStyle->SetOptDate(0);
  gStyle->SetOptStat("");

  switch(iCase){
  case 1:{ // Pi Plus
    hP0_0118 = (TH1D*)file0_0118->Get("piPlusPt");
    hP0_0130 = (TH1D*)file0_0130->Get("piPlusPt");

    hP0_0118_lowpt = (TH1D*)file0_0118->Get("piPlusPtLowPt");
    hP0_0130_lowpt = (TH1D*)file0_0130->Get("piPlusPtLowPt");

    ifstream infile2012("./PionYields2012.txt");
    ifstream infile2005("./PionYields2005.txt");
    if (infile2012.is_open() == NULL || infile2005.is_open()==NULL){
      cout << "input file is no good, exiting..." << endl;
      exit(1);
    }
    
    Double_t binWidth, xVal, yVal;
    Double_t xValTrue, yPlus, yStatPlus, ySystPlus;
    Double_t yMinus, yStatMinus, ySystMinus;
    for (Int_t iBin = 0; iBin < hP0_0118->GetNbinsX(); ++iBin){
      xValTrue = yPlus = yStatPlus = ySystPlus = 0.;
      yMinus = yStatMinus = ySystMinus = 0.;
      infile2012 >> xValTrue >> yPlus >> yStatPlus >> ySystPlus >> yMinus >> yStatMinus >> ySystMinus;

      binWidth = xVal = yVal = 0.;
      binWidth = ptbins[iBin+1]-ptbins[iBin];
      xVal = (ptbins[iBin]+ptbins[iBin+1])/2.;

      yVal = hP0_0118->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*1.0*1e7);
      gP0_0118->SetPoint(iBin,xVal,yVal/yPlus);
      Double_t yErr = sqrt(hP0_0118->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*1.0*1e7);
      Double_t rErr = (yVal/yPlus)*sqrt(pow((yStatPlus/yPlus),2.)+pow((yErr/yVal),2.));
      gP0_0118->SetPointError(iBin,0.,rErr);

      yVal = hP0_0130->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*1.0*1e7);
      gP0_0130->SetPoint(iBin,xVal,yVal/yPlus);
      Double_t yErr = sqrt(hP0_0130->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*1.0*1e7);
      Double_t rErr = (yVal/yPlus)*sqrt(pow((yStatPlus/yPlus),2.)+pow((yErr/yVal),2.));
      gP0_0130->SetPointError(iBin,0.,rErr);
    }
    infile2012.close();
    gP0_0118->SetPoint(14,0.,-1.);

    for (Int_t iBin = 0; iBin < hP0_0118_lowpt->GetNbinsX(); ++iBin){
      xValTrue = yPlus = yStatPlus = ySystPlus = 0.;
      yMinus = yStatMinus = ySystMinus = 0.;
      infile2005 >> xValTrue >> yPlus >> yStatPlus >> ySystPlus >> yMinus >> yStatMinus >> ySystMinus;

      binWidth = xVal = yVal = 0.;
      binWidth = pilowptbins[iBin+1]-pilowptbins[iBin];
      xVal = (pilowptbins[iBin]+pilowptbins[iBin+1])/2.;

      yVal = hP0_0118_lowpt->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*0.5*1e7);
      gP0_0118_lowpt->SetPoint(iBin,xVal,yVal/yPlus);
      Double_t yErr = sqrt(hP0_0118_lowpt->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*0.5*1e7);
      Double_t rErr = (yVal/yPlus)*sqrt(pow((yStatPlus/yPlus),2.)+pow((yErr/yVal),2.));
      gP0_0118_lowpt->SetPointError(iBin,0.,rErr);

      yVal = hP0_0130_lowpt->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*0.5*1e7);
      gP0_0130_lowpt->SetPoint(iBin,xVal,yVal/yPlus);
      Double_t yErr = sqrt(hP0_0130_lowpt->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*0.5*1e7);
      Double_t rErr = (yVal/yPlus)*sqrt(pow((yStatPlus/yPlus),2.)+pow((yErr/yVal),2.));
      gP0_0130_lowpt->SetPointError(iBin,0.,rErr);
    }
    infile2005.close();

    gP0_0118->SetMarkerStyle(20);
    gP0_0118->SetMarkerColor(2);
    gP0_0130->SetMarkerStyle(22);
    gP0_0130->SetMarkerColor(4);

    gP0_0118_lowpt->SetMarkerStyle(20);
    gP0_0118_lowpt->SetMarkerColor(2);
    gP0_0130_lowpt->SetMarkerStyle(22);
    gP0_0130_lowpt->SetMarkerColor(4);

    gP0_0118->GetXaxis()->SetRangeUser(0.,15.);
    gP0_0118->GetXaxis()->SetTitleFont(43);
    gP0_0118->GetXaxis()->SetTitle("p_{T}");
    gP0_0118->GetXaxis()->SetTitleSize(22);
    gP0_0118->GetXaxis()->SetTitleOffset(1.1);
    gP0_0118->GetXaxis()->CenterTitle();
    gP0_0118->GetXaxis()->SetLabelFont(43);
    gP0_0118->GetXaxis()->SetLabelSize(15);
    gP0_0118->GetXaxis()->SetLabelOffset(0.01);

    gP0_0118->GetYaxis()->SetTitleFont(43);
    gP0_0118->GetYaxis()->SetTitle("Simulation/Data Ratio");
    gP0_0118->GetYaxis()->SetTitleSize(22);
    gP0_0118->GetYaxis()->SetTitleOffset(1.1);
    gP0_0118->GetYaxis()->CenterTitle();

    gP0_0118->GetYaxis()->SetRangeUser(0.,5.);
    gP0_0118->GetYaxis()->SetLabelFont(43);
    gP0_0118->GetYaxis()->SetLabelSize(18);
    gP0_0118->GetYaxis()->SetLabelOffset(0.0025);
    gP0_0118->GetYaxis()->SetTickLength(0.02);

    gP0_0118->SetName("nnpdf0118");
    gP0_0130->SetName("nnpdf0130");

    TCanvas *canvas = new TCanvas("c","c1",800,600);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gP0_0118->Draw("AP");
    gP0_0130->Draw("P");

    gP0_0118_lowpt->Draw("P");
    gP0_0130_lowpt->Draw("P");

    plotText.DrawLatex(0.14,0.95,"1e7 thrown events");
    plotText.DrawLatex(0.14,0.9,"Charge State: +1");
    plotText.DrawLatex(0.14,0.85,"PYTHIA 6.4.28");
    plotText.DrawLatex(0.14,0.8,"Primordial k_{T} = 0.5 GeV/c");
    plotText.DrawLatex(0.14,0.75,"UE p_{T,0} = 2.0 GeV/c");
    plotText.DrawLatex(0.14,0.7,"Tune: Perugia 0");

    TLine *hline = new TLine(0.,1.,15.,1.);
    hline->SetLineStyle(2);
    hline->SetLineWidth(2);
    hline->SetLineColor(13);
    hline->Draw();

    TLegend *leg = new TLegend(0.65,0.85,0.99,0.99);
    leg->SetFillColor(0);
    leg->SetTextFont(43);
    leg->SetTextSize(16);
    leg->AddEntry("nnpdf0118","NNPDF 3.0 LO #alpha_{S} = 0.118","P");
    leg->AddEntry("nnpdf0130","NNPDF 3.0 LO #alpha_{S} = 0.130","P");
    leg->Draw();

    if (iPrint == 1)
      canvas->Print("InvariantYieldRatios.PiPlus.Perugia0.ReducedKt.NominalPt0.NNPDF.pdf");
    
    break;
  }

  case 2:{ // Pi Minus
    hP0_0118 = (TH1D*)file0_0118->Get("piMinusPt");
    hP0_0130 = (TH1D*)file0_0130->Get("piMinusPt");

    hP0_0118_lowpt = (TH1D*)file0_0118->Get("piMinusPtLowPt");
    hP0_0130_lowpt = (TH1D*)file0_0130->Get("piMinusPtLowPt");

    ifstream infile2012("./PionYields2012.txt");
    ifstream infile2005("./PionYields2005.txt");
    if (infile2012.is_open() == NULL || infile2005.is_open()==NULL){
      cout << "input file is no good, exiting..." << endl;
      exit(1);
    }
    
    Double_t binWidth, xVal, yVal;
    Double_t xValTrue, yPlus, yStatPlus, ySystPlus;
    Double_t yMinus, yStatMinus, ySystMinus;
    for (Int_t iBin = 0; iBin < hP0_0118->GetNbinsX(); ++iBin){
      xValTrue = yPlus = yStatPlus = ySystPlus = 0.;
      yMinus = yStatMinus = ySystMinus = 0.;
      infile2012 >> xValTrue >> yPlus >> yStatPlus >> ySystPlus >> yMinus >> yStatMinus >> ySystMinus;

      binWidth = xVal = yVal = 0.;
      binWidth = ptbins[iBin+1]-ptbins[iBin];
      xVal = (ptbins[iBin]+ptbins[iBin+1])/2.;

      yVal = hP0_0118->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*1.0*1e7);
      gP0_0118->SetPoint(iBin,xVal,yVal/yMinus);
      Double_t yErr = sqrt(hP0_0118->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*1.0*1e7);
      Double_t rErr = (yVal/yMinus)*sqrt(pow((yStatMinus/yMinus),2.)+pow((yErr/yVal),2.));
      gP0_0118->SetPointError(iBin,0.,rErr);

      yVal = hP0_0130->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*1.0*1e7);
      gP0_0130->SetPoint(iBin,xVal,yVal/yMinus);
      Double_t yErr = sqrt(hP0_0130->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*1.0*1e7);
      Double_t rErr = (yVal/yMinus)*sqrt(pow((yStatMinus/yMinus),2.)+pow((yErr/yVal),2.));
      gP0_0130->SetPointError(iBin,0.,rErr);
    }
    infile2012.close();
    gP0_0118->SetPoint(14,0.,-1.);

    for (Int_t iBin = 0; iBin < hP0_0118_lowpt->GetNbinsX(); ++iBin){
      xValTrue = yPlus = yStatPlus = ySystPlus = 0.;
      yMinus = yStatMinus = ySystMinus = 0.;
      infile2005 >> xValTrue >> yPlus >> yStatPlus >> ySystPlus >> yMinus >> yStatMinus >> ySystMinus;

      binWidth = xVal = yVal = 0.;
      binWidth = pilowptbins[iBin+1]-pilowptbins[iBin];
      xVal = (pilowptbins[iBin]+pilowptbins[iBin+1])/2.;

      yVal = hP0_0118_lowpt->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*0.5*1e7);
      gP0_0118_lowpt->SetPoint(iBin,xVal,yVal/yMinus);
      Double_t yErr = sqrt(hP0_0118_lowpt->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*0.5*1e7);
      Double_t rErr = (yVal/yMinus)*sqrt(pow((yStatMinus/yMinus),2.)+pow((yErr/yVal),2.));
      gP0_0118_lowpt->SetPointError(iBin,0.,rErr);

      yVal = hP0_0130_lowpt->GetBinContent(iBin+1)/(2*pi*xVal*binWidth*0.5*1e7);
      gP0_0130_lowpt->SetPoint(iBin,xVal,yVal/yMinus);
      Double_t yErr = sqrt(hP0_0130_lowpt->GetBinContent(iBin+1))/(2*pi*xVal*binWidth*0.5*1e7);
      Double_t rErr = (yVal/yMinus)*sqrt(pow((yStatMinus/yMinus),2.)+pow((yErr/yVal),2.));
      gP0_0130_lowpt->SetPointError(iBin,0.,rErr);
    }
    infile2005.close();

    gP0_0118->SetMarkerStyle(20);
    gP0_0118->SetMarkerColor(2);
    gP0_0130->SetMarkerStyle(22);
    gP0_0130->SetMarkerColor(4);

    gP0_0118_lowpt->SetMarkerStyle(20);
    gP0_0118_lowpt->SetMarkerColor(2);
    gP0_0130_lowpt->SetMarkerStyle(22);
    gP0_0130_lowpt->SetMarkerColor(4);

    gP0_0118->GetXaxis()->SetRangeUser(0.,15.);
    gP0_0118->GetXaxis()->SetTitleFont(43);
    gP0_0118->GetXaxis()->SetTitle("p_{T}");
    gP0_0118->GetXaxis()->SetTitleSize(22);
    gP0_0118->GetXaxis()->SetTitleOffset(1.1);
    gP0_0118->GetXaxis()->CenterTitle();
    gP0_0118->GetXaxis()->SetLabelFont(43);
    gP0_0118->GetXaxis()->SetLabelSize(15);
    gP0_0118->GetXaxis()->SetLabelOffset(0.01);

    gP0_0118->GetYaxis()->SetTitleFont(43);
    gP0_0118->GetYaxis()->SetTitle("Simulation/Data Ratio");
    gP0_0118->GetYaxis()->SetTitleSize(22);
    gP0_0118->GetYaxis()->SetTitleOffset(1.1);
    gP0_0118->GetYaxis()->CenterTitle();

    gP0_0118->GetYaxis()->SetRangeUser(0.,5.);
    gP0_0118->GetYaxis()->SetLabelFont(43);
    gP0_0118->GetYaxis()->SetLabelSize(18);
    gP0_0118->GetYaxis()->SetLabelOffset(0.0025);
    gP0_0118->GetYaxis()->SetTickLength(0.02);

    gP0_0118->SetName("nnpdf0118");
    gP0_0130->SetName("nnpdf0130");

    TCanvas *canvas = new TCanvas("c","c1",800,600);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.1);
    gP0_0118->Draw("AP");
    gP0_0130->Draw("P");

    gP0_0118_lowpt->Draw("P");
    gP0_0130_lowpt->Draw("P");

    plotText.DrawLatex(0.14,0.95,"1e7 thrown events");
    plotText.DrawLatex(0.14,0.9,"Charge State: -1");
    plotText.DrawLatex(0.14,0.85,"PYTHIA 6.4.28");
    plotText.DrawLatex(0.14,0.8,"Primordial k_{T} = 0.5 GeV/c");
    plotText.DrawLatex(0.14,0.75,"UE p_{T,0} = 2.0 GeV/c");
    plotText.DrawLatex(0.14,0.7,"Tune: Perugia 0");

    TLine *hline = new TLine(0.,1.,15.,1.);
    hline->SetLineStyle(2);
    hline->SetLineWidth(2);
    hline->SetLineColor(13);
    hline->Draw();

    TLegend *leg = new TLegend(0.65,0.85,0.99,0.99);
    leg->SetFillColor(0);
    leg->SetTextFont(43);
    leg->SetTextSize(16);
    leg->AddEntry("nnpdf0118","NNPDF 3.0 LO #alpha_{S} = 0.118","P");
    leg->AddEntry("nnpdf0130","NNPDF 3.0 LO #alpha_{S} = 0.130","P");
    leg->Draw();

    if (iPrint == 1)
      canvas->Print("InvariantYieldRatios.PiMinus.Perugia0.ReducedKt.NominalPt0.NNPDF.pdf");
    break;
  }

  default:
    cout << "Bad choice, try gain..." << endl;
  }
}
