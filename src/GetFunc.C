// From Steven Horvat: Fit pT spectra
// 2017.02.24
//
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TProfile.h"

#include <iomanip> 

TH1D *ProjAndScale(TFile *f, int ip, TString histname = "hreco") {
	if(!f->IsOpen()) {cout<<"TH1D *ProjAndScale(): TFile not open"<<endl; return NULL;}
	TH2D *h = (TH2D*)f->Get(histname);
	if(!h) { cout<<"TH1D *ProjAndScale(): Cannot find "<<histname<<endl; return NULL;}
	TH1D *h1 = (TH1D*)h->ProjectionY(Form("%s%d",h->GetName(),ip),ip,ip);
	h1->Scale(1,"width");
	return h1;
}


// LY get mean and its error
double GetMean(TF1 *func, double xmin = 0.2, double xmax = 20) {
	return func->Mean(xmin, xmax);
}

double GetMeanError(TF1 *func, double xmin = 0.2, double xmax = 20) {
	return sqrt(func->Variance(xmin,xmax)/func->Integral(xmin,xmax));// mean error = StdDev / sqrt( Neff )
	
}

struct FuncMean {
	FuncMean(TF1 *f1):func(f1) {
		xfunc = new TF1("xfunc",Form("x*%s",func->GetName()));
	}
	double operator() (double *x, double *par)  const {
		//double xmin = x[0];
		//double xmax = x[1];
		//func = new TF1("func","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))");	
		double xmax = x[0];		// Somehow 2D does not work... use 1D only.. 
		double xmin = 0.2;
		func->SetParameters(par);
		xfunc->SetParameters(par);
		double y = xfunc->Integral(xmin, xmax)/func->Integral(xmin,xmax);
		//double y = par[0]*x[0]+par[1];
		cout<<"In FuncMean:"<<endl;
		cout<<"xmin = "<<xmin<<" xmax = "<<xmax<<endl;
		cout<<"p0 = "<<par[0]<<" p1 = "<<par[1]<<" p2 = "<<par[2]<<" 	x[0] = "<<x[0]<<" y = "<<std::setprecision(9)<<y<<endl;
		return y;
	}
	TF1 *func;
	TF1 *xfunc;
};



//double GetMeanFitError(TF1 *func, double xmin = 0.2, double xmax = 20, const Double_t *params = 0, const Double_t *covmat = 0) {
double GetMeanFitError(TF1 *func, double xmin = 0.2, double xmax = 20, const Double_t *covmat = 0) {
	int Npar = func->GetNpar();
	double *params = func->GetParameters();
	//cout<<"pars: "<<endl;
	//for(int i = 0; i<Npar; i++) {
	//	cout<<params[i]<<endl;
	//}
	//vector<double> si;		// (Delta par_i)^2
	//for(int i = 0; i<Npar; i++) {
	//	si.push_back(covmat[Npar*i+i]);
	//}
	//TF1 *xfunc = new TF1("xfunc",Form("x*%s",func->GetName()),xmin, xmax);
	double epsilon = 0.001;	
	double xminxmax[2] = {xmin, xmax};
	double *pxminxmax = xminxmax;
	FuncMean *funmean  = new FuncMean(func);
	TF1 *fmean = new TF1("fmean",funmean,0.2,20, 3);
	fmean->SetParameters(params);
	vector<double> der;	// derivative
	cout<<"Set der[0] = 0. As it should be cancelled in mean calculations"<<endl;
	der.push_back(0);
	for(int i = 1; i<Npar ;i++) {
		//der.push_back(fmean->GradientPar(i, pxminxmax, epsilon));
		der.push_back(fmean->GradientPar(i, &xmax, epsilon));
		cout<<epsilon*covmat[i*Npar+i]<<endl;
		cout<<"der["<<i<<"] = "<<der.at(i)<<endl;
	}


	double varsum = 0;
	for(int i = 0; i<Npar ;i++) {
		for(int j = 0; j<Npar ;j++) {
			varsum+=der.at(i)*der.at(j)*covmat[i*Npar+j];
		}
	}
	
	return sqrt(varsum);
}




TF1* GetFunc4Hist(TH1D* hist1, double &mean, double &meanerr, double &meanfiterr, double xmin = 0.2, double xmax = 20) {// ,Int_t msqrts,Int_t cent,Int_t species){


  TF1* func = new TF1("func","[0]*x*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))",xmin,xmax);
  //TF1* func = new TF1("func","[0]*pow(1.-[1]*x*x,[2])",0.10,100.0);
  func->SetParameters(10.*hist1->GetBinContent(10),1.5,1.1);
  //func->SetParLimits(1,.001,10.);
  func->SetParLimits(2,1.0001,10.);
  //TCanvas *c = new TCanvas("c", "canvas", 800, 600);
  Int_t fitFlag =-1;
  hist1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hist1->GetYaxis()->SetTitle("#frac{d^{2}N}{2#pip_{T}d#etadp_{T}} (GeV/c)^{-2}");
  //hist1->DrawCopy("AP");
  fitFlag = hist1->Fit(func,"QRsame","",xmin,xmax);
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
      fitFlag=hist1->Fit(func,"QRsame","",xmin,xmax);
    }
  }
  TFitResultPtr fitptr=hist1->Fit(func,"QRSEIsame","",xmin, xmax);
  cout << "fit flag = " << fitFlag << endl;
  gPad->SetLogy();
  cout << "gMinuit->fCstatu = " << gMinuit->fCstatu << "\tchisquare/NDF = " << func->GetChisquare() << "/" << func->GetNDF() << endl;
  cout << "[0] = " << func->GetParameter(0) << " +/- "<< func->GetParError(0) << " and [1] = " << func->GetParameter(1) << " +/- "<< func->GetParError(1) << " and [2] = " << func->GetParameter(2) << " +/- "<< func->GetParError(2) << endl;
//  Char_t buffer[100];
//  sprintf(buffer,"fit_spc%d_cent%d_%d",species,cent,msqrts);
//  TString pName = TString("effQA/spectra/")+TString(buffer)+TString(".png");
//  c1->Print(pName);
  mean = GetMean(func, xmin, xmax);
  meanerr = GetMeanError(func, xmin, xmax);
  double *cov = fitptr->GetCovarianceMatrix().GetMatrixArray();
  //meanfiterr = GetMeanFitError(func,xmin, xmax, func->GetParameters(),cov);
  meanfiterr = GetMeanFitError(func,xmin, xmax, cov);
  cout<<"mean = "<<mean<<" meanerr = "<<meanerr<<" meanfiterr = "<<meanfiterr<<endl;
  return (TF1*)func;
}




TH1D* GetFunc(vector<double> &Amean, vector<double> &Ameanerr, vector<double> &Ameanfiterr, TString filename = "Unfolding_TranPtAveJPCharged_NFweight_Baye5_McPt02.root", TString histname = "hreco", double xmin = 0.2, double xmax = 10) {		// Don't use TProfile, SetBinContent, SetBinError does not work

  TFile *file = new TFile(filename);
  TH2D *hreco = (TH2D*)file->Get(histname);
  TH1D *pfx = (TH1D*)hreco->ProjectionX(histname+"pfx");
  pfx->Reset();

  for(int ibin = 1; ibin<hreco->GetNbinsX()+1; ibin++) { 
  	TH1D *hist1 = (TH1D*)hreco->ProjectionY(Form("%s%d",hreco->GetName(),ibin),ibin,ibin);
	if(!hist1->GetEntries()) continue;
  	hist1->Scale(1,"width");
  	double mean = 0, meanerr = 0, meanfiterr = 0;
	TF1 *f1 = GetFunc4Hist(hist1, mean, meanerr, meanfiterr, xmin, xmax); 
	Amean.push_back(mean);
	Ameanerr.push_back(meanerr);
	Ameanfiterr.push_back(meanfiterr);

	pfx->SetBinContent(ibin, mean);
	pfx->SetBinError(ibin, meanerr);
  }

  return pfx;

}


