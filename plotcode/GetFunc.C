// From Steven Horvat: Fit pT spectra
// 2017.02.24
//
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

// LY get mean and its error
double GetMean(TF1 *func, double xmin = 0.2, double xmax = 20) {
	return func->Mean(xmin, xmax);
}

double GetMeanError(TF1 *func, double xmin = 0.2, double xmax = 20) {
	return sqrt(func->Variance(xmin,xmax)/func->Integral(xmin,xmax));// mean error = StdDev / sqrt( Neff )
	
}

//Double_t FuncMean(Double_t *x, Double_t *par) {	// 
//	//double xmin = 0.2, xmax = 20;
//	//double xmin = par[0];
//	//double xmax = par[1];
//	double xmin = x[0];
//	double xmax = x[1];
//	TF1* func = new TF1("func","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))");	
//	TF1 *xfunc = new TF1("xfunc",Form("x*%s",func->GetName()));
//	func->SetParameters(par);
//	xfunc->SetParameters(par);
//	double y = xfunc->Integral(xmin, xmax)/func->Integral(xmin,xmax);
//	return y;
//}
struct FuncMean {
	FuncMean(TF1 *f1):func(f1) {
		cout<<f1->GetName()<<endl; 
		xfunc = new TF1("xfunc",Form("x*%s",func->GetName()));
	}
	double operator() (double *x, double *par)  const {
		double xmin = x[0];
		double xmax = x[1];
		//func = new TF1("func","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))");	
		func->SetParameters(par);
		xfunc->SetParameters(par);
		double y = xfunc->Integral(xmin, xmax)/func->Integral(xmin,xmax);
		return y;
	}
	TF1 *func;
	TF1 *xfunc;
};

double GetMeanFitError(TF1 *func, double xmin = 0.2, double xmax = 20, const Double_t *params = 0, const Double_t *covmat = 0) {
	int Npar = func->GetNpar();
	//vector<double> si;		// (Delta par_i)^2
	//for(int i = 0; i<Npar; i++) {
	//	si.push_back(covmat[Npar*i+i]);
	//}
	//TF1 *xfunc = new TF1("xfunc",Form("x*%s",func->GetName()),xmin, xmax);
	double epsilon = 0.001;	
	double xminxmax[2] = {xmin, xmax};
	double *pxminxmax = xminxmax;
	FuncMean *funmean  = new FuncMean(func);
	TF2 *fmean = new TF2("fmean",funmean,0.2,20,0.2,20,3);
	//TF3 *fmean = new TF3("fmean",FuncMean,0.2,20,0.2,20,3);
	//fmean->SetParameters(xmin,xmax);
	vector<double> der;	// derivative
	for(int i = 0; i<Npar ;i++) {
		der.push_back(fmean->GradientPar(i, pxminxmax, epsilon*covmat[i*Npar+i]));
	}

	double varsum = 0;
	for(int i = 0; i<Npar ;i++) {
		for(int j = 0; j<Npar ;j++) {
			varsum+=der.at(i)*der.at(j)*covmat[i*Npar+j];
		}
	}
	
	return sqrt(varsum);
}




TF1* GetFunc(TH1D* hist1, double &mean, double &meanerr, double &meanfiterr, double xmin = 0.2, double xmax = 20) {// ,Int_t msqrts,Int_t cent,Int_t species){
  TF1* func = new TF1("func","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))",xmin,xmax);
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
  mean = GetMean(func, xmin, xmax);
  meanerr = GetMeanError(func, xmin, xmax);
  double *cov = fitptr->GetCovarianceMatrix().GetMatrixArray();
  meanfiterr = GetMeanFitError(func,xmin, xmax, func->GetParameters(),cov);
  return (TF1*)func;
}


