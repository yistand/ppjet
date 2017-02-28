// From Steven Horvat: Fit pT spectra
// 2017.02.24
//



TF1* GetFunc(TH1D* hist1) {// ,Int_t msqrts,Int_t cent,Int_t species){
  TF1* func = new TF1("func","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))",0.20,20.0);
  //TF1* func = new TF1("func","[0]*pow(1.-[1]*x*x,[2])",0.10,100.0);
  func->SetParameters(10.*hist1->GetBinContent(10),1.5,1.1);
  //func->SetParLimits(1,.001,10.);
  func->SetParLimits(2,1.0001,10.);
  //TCanvas *c = new TCanvas("c", "canvas", 800, 600);
  Int_t fitFlag =-1;
  hist1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hist1->GetYaxis()->SetTitle("#frac{d^{2}N}{2#pip_{T}d#etadp_{T}} (GeV/c)^{-2}");
  //hist1->DrawCopy("AP");
  fitFlag = hist1->Fit(func,"QRsame","",0.10,100.0);
  for(int jk = 1; jk < 100; jk++){
    if((TString)gMinuit->fCstatu == TString("SUCCESSFUL"))break;
    for(int kj=1; kj < 100; kj++){
      if((TString)gMinuit->fCstatu == TString("SUCCESSFUL")){
    cout << "{jk,kj} = {" <<  .00099*(double)jk << "," << 10*(double)kj << "}" << endl;
    break;
      }
func->SetParameters(10.*hist1->GetBinContent(10),1.+.00099*(double)jk,1.+.01*(double)kj);
      //func->SetParLimits(1,.05,10.);
      //func->SetParLimits(2,1.0001,10.);
      fitFlag=hist1->Fit(func,"QRsame","",0.10,100.0);
    }
  }
  TFitResultPtr fitptr=hist1->Fit(func,"QRSEIsame","",0.10,100.0);
  cout << "fit flag = " << fitFlag << endl;
  gPad->SetLogy();
  cout << "gMinuit->fCstatu = " << gMinuit->fCstatu << "\tchisquare/NDF = " << func->GetChisquare() << "/" << func->GetNDF() << endl;
  cout << "[0] = " << func->GetParameter(0) << " and [1] = " << func->GetParameter(1) << " and [2] = " << func->GetParameter(2) << endl;
//  Char_t buffer[100];
//  sprintf(buffer,"fit_spc%d_cent%d_%d",species,cent,msqrts);
//  TString pName = TString("effQA/spectra/")+TString(buffer)+TString(".png");
//  c1->Print(pName);
  return (TF1*)func;
}
