//====================================================================================================
//
//	2017.12.05	Li YI
//	Fit 1Ds with Gaussians
//
//====================================================================================================

#include "FitPid2D.h"

#include "TF2.h"

#include <iostream>

using std::cout;
using std::endl;

FitPid2D::FitPid2D() {}

FitPid2D::~FitPid2D() { 
	//if(func2d) {cout<<func2d<<endl; delete func2d; }
	//if(func1dx) delete func1dx;
	//if(func1dy) delete func1dy;
	//cout<<"delete func"<<endl;
	//for(int i = 0; i<NSpecies; i++) {
	//	if(func1x_enhanced[i]) delete func1x_enhanced[i];
	//	if(func1y_enhanced[i]) delete func1y_enhanced[i];
	//}
	//cout<<"delete array func"<<endl;
	if(h2d) delete h2d;
	if(h1dx) delete h1dx;
	if(h1dy) delete h1dy;
	cout<<"delete th1d"<<endl;
	for(int i = 0; i<NSpecies; i++) {
		if(h1x_enhanced[i]) delete h1x_enhanced[i];
		if(h1y_enhanced[i]) delete h1y_enhanced[i];
	}
	cout<<"delete array th1d"<<endl;
}

//================================================ input ====================================================
bool FitPid2D::Input2D(TH2D* h) {
	if(h) {
		h2d = (TH2D*)h->Clone(Form("%s_copy",h->GetName()));
		h2d->Sumw2();
		return true;
	}
	else return false;
}

//================================================ Rebin ====================================================
bool FitPid2D::Rebin2D(int rebin){
	cout<<"Info: FitPid2D::Rebin2D() needs to be called at the begining, before Project1Ds()"<<endl;
	return h2d->Rebin2D(rebin,rebin);
}

//================================================ Projection ====================================================
bool FitPid2D::Project1Ds(double *ProjEnhanceRangeX, double *ProjEnhanceRangeY) {	 // project to 1D histograms. arguments are projection ranges for enhanced species
// ProjEnhanceRangeX and ProjEnhanceRangeY length is 2*NSpecies each
// for enhanced 1Ds: each dimension has NSpecies range pairs: pion_min, pion_max, kaon_min, kaon_max, proton_min, proton_max
	if(!h2d) return false;

	h1dx = (TH1D*)h2d->ProjectionX();
	h1dy = (TH1D*)h2d->ProjectionY();

	h1dx->Sumw2();
	h1dy->Sumw2();

	// Normalized by bin width
	h1dx->Scale(1./h1dx->GetBinWidth(0));
	h1dy->Scale(1./h1dy->GetBinWidth(0));

	for(int i = 0; i<NSpecies; i++) {
		int ymin = h1dy->FindBin(ProjEnhanceRangeY[i*2]);
		int ymax = h1dy->FindBin(ProjEnhanceRangeY[i*2+1]);
		if(ymin>ymax) {cout<<"WARNING!! Project1Ds-x: ymin="<<ymin<<" > ymax="<<ymax<<"!!! I'm going to swap them"<<endl; std::swap(ymin,ymax);}
		cout<<"ProjectX"<<i<<": "<<ProjEnhanceRangeY[i*2]<<", "<<ProjEnhanceRangeY[i*2+1]<<" -> "<<ymin<<", "<<ymax<<endl;
		h1x_enhanced[i] = (TH1D*)h2d->ProjectionX(Form("%sh1x_enhanced%d",h2d->GetName(),i),ymin, ymax);
		h1x_enhanced[i]->Sumw2();
		h1x_enhanced[i]->Scale(1./h1x_enhanced[i]->GetBinWidth(0));

		int xmin = h1dx->FindBin(ProjEnhanceRangeX[i*2]);
		int xmax = h1dx->FindBin(ProjEnhanceRangeX[i*2+1]);
		if(xmin>xmax) {cout<<"WARNING!! Project1Ds-y: xmin="<<xmin<<" > xmax="<<xmax<<"!!! I'm going to swap them"<<endl; std::swap(xmin,xmax);}
		cout<<"ProjectY"<<i<<": "<<ProjEnhanceRangeX[i*2]<<", "<<ProjEnhanceRangeX[i*2+1]<<" -> "<<xmin<<", "<<xmax<<endl;
		h1y_enhanced[i] = (TH1D*)h2d->ProjectionY(Form("%sh1y_enhanced%d",h2d->GetName(),i),xmin, xmax);
		h1y_enhanced[i]->Sumw2();
		h1y_enhanced[i]->Scale(1./h1y_enhanced[i]->GetBinWidth(0));
	}

	return true;
}

//================================================ Initial Fit Parameters ====================================================
bool FitPid2D::InitialPars(double *par){  // Initial fitting parameters
	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = par[i];
		meanx[i] = par[i+NSpecies];
		sigmax[i] = par[i+NSpecies*2];
		meany[i] = par[i+NSpecies*3];
		sigmay[i] = par[i+NSpecies*4];
		for(int j = 0; j<NSpecies; j++) {
			Yieldx_enhanced[i][j] = par[j+NSpecies*i+NSpecies*5];
			Yieldy_enhanced[i][j] = par[j+NSpecies*i+NSpecies*5+NSpecies*NSpecies];
		}
	}
	return true;
}

//================================================ GuessYieldReInit ====================================================
bool FitPid2D::GuessYieldReInit(double *rangex, double *rangey) {
	int rangexbin[NSpecies*2] = {0};
	int rangeybin[NSpecies*2] = {0};
	for(int i = 0; i<NSpecies*2; i++) {
		rangexbin[i] = h1dx->FindBin(rangex[i]);
		rangeybin[i] = h1dy->FindBin(rangey[i]);
	}

	double xbinwidth = h1dx->GetBinWidth(0);
	double ybinwidth = h1dy->GetBinWidth(0);

	for(int i = 0; i<NSpecies; i++) {
		double y2 = h2d->Integral(rangexbin[i*2],rangexbin[i*2+1],rangeybin[i*2],rangeybin[i*2+1]);	// h2d is NOT normalized by bin width, Don't integral with 'width'
		//if(y2>0) Yield[i] = y2*xbinwidth*ybinwidth;
		if(y2>0) Yield[i] = y2;
		for(int j = 0; j<NSpecies; j++) {
			double yx = h1x_enhanced[i]->Integral(rangexbin[j*2],rangexbin[j*2+1],"width");		// h1x_enhanced is normalized by bin width
			double yy = h1y_enhanced[i]->Integral(rangeybin[j*2],rangeybin[j*2+1],"width");
			//if(yx>0) Yieldx_enhanced[i][j] = yx*xbinwidth;
			//if(yy>0) Yieldy_enhanced[i][j] = yy*ybinwidth;
			if(yx>0) Yieldx_enhanced[i][j] = yx;
			if(yy>0) Yieldy_enhanced[i][j] = yy;
		}
	}
	return true;
}

//================================================ Print ====================================================
void FitPid2D::PrintPars(){  // Print fitting parameters
	for(int i = 0; i<NSpecies; i++) {
		cout<<"Yield["<<i<<"] = "<<Yield[i]<<endl;
		cout<<"meanx["<<i<<"] = "<<meanx[i]<<endl;
		cout<<"sigmax["<<i<<"] = "<<sigmax[i]<<endl;
		cout<<"meany["<<i<<"] = "<<meany[i]<<endl;
		cout<<"sigmay["<<i<<"] = "<<sigmay[i]<<endl;
		for(int j = 0; j<NSpecies; j++) {
			cout<<"Yieldx_enhanced["<<i<<"]["<<j<<"] = "<<Yieldx_enhanced[i][j]<<endl;
		}
		for(int j = 0; j<NSpecies; j++) {
			cout<<"Yieldy_enhanced["<<i<<"]["<<j<<"] = "<<Yieldy_enhanced[i][j]<<endl;
		}
	}
}

//================================================ Fit (the Main) ====================================================
bool FitPid2D::CombinedFit1Ds(bool useInitPar) {
	const int Nfitfunc = NSpecies*2+2;
	
	TH1D *h1[Nfitfunc];
	h1[0] = h1dx;
	h1[1] = h1dy;	
	for(int i = 0; i<NSpecies; i++) {
		h1[2+i] = h1x_enhanced[i];
	}
	for(int i = 0; i<NSpecies; i++) {
		h1[2+NSpecies+i] = h1y_enhanced[i];
	}
	
	TF1 *f1[Nfitfunc];	// Nfitfunc 3-Gaussian
	for(int i = 0; i<Nfitfunc; i++) {
		f1[i] = new TF1(Form("fun%d",i),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
	}

	ROOT::Math::WrappedMultiTF1 *wf1[Nfitfunc];
	for(int i = 0; i<Nfitfunc; i++) {
		wf1[i] = new ROOT::Math::WrappedMultiTF1(*f1[i],1);
	}
	
	ROOT::Fit::DataOptions opt;
	
	ROOT::Fit::DataRange range;

	ROOT::Fit::BinData *data[Nfitfunc];

	ROOT::Fit::Chi2Function *chi2[Nfitfunc];

	for(int i = 0; i<Nfitfunc; i++) {
		data[i] = new ROOT::Fit::BinData(opt,range);
		ROOT::Fit::FillData(*(data[i]),h1[i]);
		chi2[i] = new ROOT::Fit::Chi2Function(*(data[i]),*(wf1[i]));
	}
	
	GlobalChi2 globalChi2(*chi2[0],*chi2[1],*chi2[2],*chi2[3],*chi2[4],*chi2[5],*chi2[6],*chi2[7]);

	ROOT::Fit::Fitter fitter;
	
	const int Npar = 33;
	double par0[Npar] = {1e4,1e4,1e4,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1e3};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
			for(int j = 0; j<NSpecies; j++) {
				par0[j+NSpecies*i+NSpecies*5] = Yieldx_enhanced[i][j];
				par0[j+NSpecies*i+NSpecies*5+NSpecies*NSpecies] = Yieldy_enhanced[i][j];
			}
		}
	}
	
	fitter.Config().SetParamsSettings(Npar,par0);
	fitter.Config().ParSettings(0).SetLimits(1,1e7);
	fitter.Config().ParSettings(1).SetLimits(1,1e7);
	fitter.Config().ParSettings(2).SetLimits(1,1e7);
	fitter.Config().ParSettings(3).SetLimits(-2,2);
	fitter.Config().ParSettings(4).SetLimits(-2,2);
	fitter.Config().ParSettings(5).SetLimits(-2,2);
	fitter.Config().ParSettings(6).SetLimits(0,1);
	fitter.Config().ParSettings(7).SetLimits(0,1);
	fitter.Config().ParSettings(8).SetLimits(0,1);
	fitter.Config().ParSettings(9).SetLimits(-2,2);
	fitter.Config().ParSettings(10).SetLimits(-2,2);
	fitter.Config().ParSettings(11).SetLimits(-2,2);
	fitter.Config().ParSettings(12).SetLimits(0,1);
	fitter.Config().ParSettings(13).SetLimits(0,1);
	fitter.Config().ParSettings(14).SetLimits(0,1);

	fitter.Config().MinimizerOptions().SetPrintLevel(10);
	fitter.Config().SetMinimizer("Minuit2","Migrad");

	// fit FCN function directly
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	int totalsize = 0;
	for(int i = 0; i<Nfitfunc; i++) {
		totalsize+=data[i]->Size();
	}
	fitter.FitFCN(Npar,globalChi2,0,totalsize,true);
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);

	if(!result.IsValid()) {
		cout<<"CombinedFit Failed!!"<<endl;
		return false;
	}

	// each one of the 8 3-gaus has 9 parameters
	// location of pars for each fit function
	//int pars[8][9] = {{0,3,6,1,4,7,2,5,8},		// projection X
	//		{0,9,12,1,10,13,2,11,14},	// projection Y
	//		{15,3,6,16,4,7,17,5,8},		// projection enhanced pion X
	//		{18,3,6,19,4,7,20,5,8},		// projection enhanced kaon X
	//		{21,3,6,22,4,7,23,5,8},		// projection enhanced proton X
	//		{24,9,12,25,10,13,26,11,14},	// projection enhanced pion Y
	//		{27,9,12,28,10,13,29,11,14},	// projection enhanced kaon Y
	//		{30,9,12,31,10,13,32,11,14}};	// projection enhanced proton Y
	int pars[8][9] = {{0,1,2,3,4,5,6,7,8},		// projection X
			{0,1,2,9,10,11,12,13,14},	// projection Y
			{15,16,17,3,4,5,6,7,8},		// projection enhanced pion X
			{18,19,20,3,4,5,6,7,8},		// projection enhanced kaon X
			{21,22,23,3,4,5,6,7,8},		// projection enhanced proton X
			{24,25,26,9,10,11,12,13,14},	// projection enhanced pion Y
			{27,28,29,9,10,11,12,13,14},	// projection enhanced kaon Y
			{30,31,32,9,10,11,12,13,14}};	// projection enhanced proton Y


	for(int i = 0; i<Nfitfunc; i++) {
		cout<<"set "<<i<<endl;
		f1[i]->SetFitResult(result,pars[i]);
		//f1[i]->SetRange(range().first,range().second);
	}
	cout<<"set done"<<endl;


	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f1[0]->GetParameter(i);
		eYield[i] = f1[0]->GetParError(i);
		meanx[i] = f1[0]->GetParameter(i+NSpecies*1);
		sigmax[i] = f1[0]->GetParameter(i+NSpecies*2);
		meany[i] = f1[1]->GetParameter(i+NSpecies*1);
		sigmay[i] = f1[1]->GetParameter(i+NSpecies*2);
		for(int j = 0; j<NSpecies; j++) {
			Yieldx_enhanced[j][i] = f1[2+j]->GetParameter(i+NSpecies);
			Yieldy_enhanced[j][i] = f1[2+NSpecies+j]->GetParameter(i+NSpecies);
		}
	}

	for(int i = 0; i<Nfitfunc; i++) {
		for(int j = 0; j<f1[i]->GetNpar(); j++) {
			cout<<"f1["<<i<<"] par"<<j<<" = "<<f1[i]->GetParameter(j)<<endl;
		}
	}



	//h1dx->GetListOfFunctions()->Add(f1[0]);
	//h1dy->GetListOfFunctions()->Add(f1[1]);
	func1dx = (TF1*)f1[0]->Clone("func1dx");
	for(int i = 0; i<f1[0]->GetNpar(); i++) {
		func1dx->SetParameter(i,f1[0]->GetParameter(i));
	}
	func1dy = (TF1*)f1[1]->Clone("func1dy");
	for(int i = 0; i<f1[0]->GetNpar(); i++) {
		func1dy->SetParameter(i,f1[0]->GetParameter(i));
	}
	for(int i = 0; i<NSpecies; i++) {
		//h1x_enhanced[i]->GetListOfFunctions()->Add(f1[i+2]);
		//h1y_enhanced[i]->GetListOfFunctions()->Add(f1[i+2+NSpecies]);
		func1x_enhanced[i] = (TF1*)f1[i+2]->Clone(Form("funcx1_enhanced%d",i));
		func1y_enhanced[i] = (TF1*)f1[i+2+NSpecies]->Clone(Form("funcx1_enhanced%d",i));
		for(int j = 0; j<f1[i+2]->GetNpar(); j++) {
			func1x_enhanced[i]->SetParameter(j,f1[i+2]->GetParameter(j));
		}
		for(int j = 0; j<f1[i+2]->GetNpar(); j++) {
			func1y_enhanced[i]->SetParameter(j,f1[i+2+NSpecies]->GetParameter(j));
		}
	}	

	return true;
}


//================================================ Fit: 1gaus for enhanced instead 3gaus for 1Ds, 2-D 3gaus for 2D  ====================================================
bool FitPid2D::CombinedFit1Ds_2(bool useInitPar) {
	const int Nfitfunc = NSpecies*2+1;	// 1gaus fit for each 1Ds, 2-D 3gaus fit for 2D
	
	TH2D *h2d_copy = (TH2D*)h2d->Clone("h2d_copy");
	h2d_copy->Sumw2();
	h2d_copy->Scale(1./(h2d_copy->GetXaxis()->GetBinWidth(0)*h2d_copy->GetYaxis()->GetBinWidth(0)));

	TH1D *h1[NSpecies*2];
	for(int i = 0; i<NSpecies; i++) {
		h1[i] = h1x_enhanced[i];
	}
	for(int i = 0; i<NSpecies; i++) {
		h1[NSpecies+i] = h1y_enhanced[i];
	}
	
	TF2 *f2 = new TF2("fun0","[0]*TMath::Gaus(x,[3],[6],1)*TMath::Gaus(y,[9],[12],1)+[1]*TMath::Gaus(x,[4],[7],1)*TMath::Gaus(y,[2],[13],1)+[2]*TMath::Gaus(x,[5],[8],1)*TMath::Gaus(y,[11],[14],1)",-2,2,-2,2);
	f2->SetNpx(100);
	f2->SetNpy(100);
	//TF2 *f2 = new TF2("fun0","[0]/((2*pi)*[6]*[12])*exp(-0.5*(pow((x-[3])/[6],2)+pow((y-[9])/[12],2)))+[1]/((2*pi)*[7]*[13])*exp(-0.5*(pow((x-[4])/[7],2)+pow((y-[10])/[13],2)))+[2]/((2*pi)*[8]*[14])*exp(-0.5*(pow((x-[5])/[8],2)+pow((y-[11])/[14],2)))",-10,10,-10,10);
	TF1 *f1[NSpecies*2];	// 1-Gaussian
	for(int i = 0; i<NSpecies*2; i++) {
		//f1[i] = new TF1(Form("fun%d",i+1),"[0]/(pow(2*pi,0.5)*[2])*exp(-0.5*(pow((x-[1])/[2],2)))",-10,10);
		f1[i] = new TF1(Form("fun%d",i+1),"[0]*TMath::Gaus(x,[1],[2],1)",-2,2);
		f1[i]->SetNpx(1000);
	}

	ROOT::Math::WrappedMultiTF1 *wf2;
	ROOT::Math::WrappedMultiTF1 *wf1[NSpecies*2];
	wf2 = new ROOT::Math::WrappedMultiTF1(*f2,1);
	for(int i = 0; i<NSpecies*2; i++) {
		wf1[i] = new ROOT::Math::WrappedMultiTF1(*f1[i],1);
	}
	
	ROOT::Fit::DataOptions opt;
	
	ROOT::Fit::DataRange range;

	ROOT::Fit::BinData *data[Nfitfunc];

	ROOT::Fit::Chi2Function *chi2[Nfitfunc];

	for(int i = 0; i<Nfitfunc; i++) {
		data[i] = new ROOT::Fit::BinData(opt,range);
	}
	ROOT::Fit::FillData(*(data[0]),h2d_copy);
	chi2[0] = new ROOT::Fit::Chi2Function(*(data[0]),*(wf2));
	for(int i = 0; i<NSpecies*2; i++) {
		ROOT::Fit::FillData(*(data[i+1]),h1[i]);
		chi2[i+1] = new ROOT::Fit::Chi2Function(*(data[i+1]),*(wf1[i]));
	}
	
	GlobalChi2_2 globalChi2_2(*chi2[0],*chi2[1],*chi2[2],*chi2[3],*chi2[4],*chi2[5],*chi2[6]);

	ROOT::Fit::Fitter fitter;
	
	const int Npar = 21;
	double par0[Npar] = {1e4,1e4,1e4,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01,1e3,1e3,1e3,1e3,1e3,1e3};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
			par0[i+NSpecies*5] = Yieldx_enhanced[i][i];
			par0[i+NSpecies*6] = Yieldy_enhanced[i][i];
		}
	}
	
	fitter.Config().SetParamsSettings(Npar,par0);

	fitter.Config().MinimizerOptions().SetPrintLevel(10);
	fitter.Config().SetMinimizer("Minuit2","Migrad");

	// fit FCN function directly
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	int totalsize = 0;
	for(int i = 0; i<Nfitfunc; i++) {
		totalsize+=data[i]->Size();
	}
	fitter.FitFCN(Npar,globalChi2_2,0,totalsize,true);
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);

	if(!result.IsValid()) {
		cout<<"CombinedFit2 Failed!!"<<endl;
		return false;
	}

	// 1 2-D 3-gaus (15 parameters) + 6 1-D 1-gaus (each 3 parameters)
	// location of pars for each fit function
	int par2d[15] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
	int pars1d[6][3] = {
			{15,3,6},		// projection enhanced pion X
			{16,4,7},		// projection enhanced kaon X
			{17,5,8},		// projection enhanced proton X
			{18,9,12},	// projection enhanced pion Y
			{19,10,13},	// projection enhanced kaon Y
			{20,11,14}};	// projection enhanced proton Y


	f2->SetFitResult(result,par2d);
	for(int i = 0; i<NSpecies*2; i++) {
		cout<<"set "<<i<<endl;
		f1[i]->SetFitResult(result,pars1d[i]);
		//f1[i]->SetRange(range().first,range().second);
	}
	cout<<"set done"<<endl;

	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f2->GetParameter(i);
		eYield[i] = f2->GetParError(i);
		meanx[i] = f2->GetParameter(i+NSpecies);
		sigmax[i] = f2->GetParameter(i+NSpecies*2);
		meany[i] = f2->GetParameter(i+NSpecies*3);
		sigmay[i] = f2->GetParameter(i+NSpecies*4);
	}
	for(int i = 0; i<NSpecies; i++) {
		Yieldx_enhanced[i][i] = f1[i]->GetParameter(0);
	}

	h2d_copy->GetListOfFunctions()->Add(f2);

	func2d = (TF2*)f2->Clone("f2_h2d");
	func2d->SetParameters(f2->GetParameters());
	//func2d->SetParameter(0,f2->GetParameter(0)*(h2d_copy->GetXaxis()->GetBinWidth(0)*h2d_copy->GetYaxis()->GetBinWidth(0)));
	h2d->GetListOfFunctions()->Add(func2d);



	for(int i = 0; i<NSpecies; i++) {
		h1x_enhanced[i]->GetListOfFunctions()->Add(f1[i]);
		h1y_enhanced[i]->GetListOfFunctions()->Add(f1[i+NSpecies]);

		func1x_enhanced[i] = (TF1*)f1[i]->Clone(Form("func1x_enhanced%d",i));
		func1y_enhanced[i] = (TF1*)f1[i]->Clone(Form("func1y_enhanced%d",i));
	}	
	TCanvas *c = new TCanvas();
	gStyle->SetOptFit(1);
	c->SetLogz();
	h2d_copy->Draw("surf2");
	c->SaveAs("h2d.png");
	c->SetLogy();
	c->Clear();
	h1dx->Draw();
	//TF1 *f3gausx = new TF1("f3gausx","[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
	TF1 *f3gausx = new TF1("f3gausx","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	f3gausx->SetNpx(1000);
	cout<<Yield[0]<<", "<<meanx[0]<<", "<<sigmax[0]<<", "<<Yield[1]<<", "<<meanx[1]<<", "<<sigmax[1]<<", "<<Yield[2]<<", "<<meanx[2]<<", "<<sigmax[2]<<endl;
	f3gausx->SetParameters(Yield[0],meanx[0],sigmax[0],Yield[1],meanx[1],sigmax[1],Yield[2],meanx[2],sigmax[2]);
	f3gausx->Draw("same");
	c->SaveAs("h2d_px.png");
	h1dx->GetListOfFunctions()->Add(f3gausx);
	c->Clear();
	h1dy->Draw();
	//TF1 *f3gausy = new TF1("f3gausy","[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
	TF1 *f3gausy = new TF1("f3gausy","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	f3gausy->SetNpx(10000);
	f3gausy->SetParameters(Yield[0],meany[0],sigmay[0],Yield[1],meany[1],sigmay[1],Yield[2],meany[2],sigmay[2]);
	f3gausy->Draw("same");
	c->SaveAs("h2d_py.png");
	h1dy->GetListOfFunctions()->Add(f3gausy);
	for(int i = 0; i<NSpecies; i++) {
		c->Clear();
		h1x_enhanced[i]->Draw();
		c->SaveAs(Form("h1dx_%d.png",i));
		c->Clear();
		h1y_enhanced[i]->Draw();
		c->SaveAs(Form("h1dy_%d.png",i));
	}

	h2d_copy->Delete();

	return true;
}


//================================================ Fit: 3gaus for 2 non-enhanced 1Ds  ====================================================
bool FitPid2D::CombinedFit1Ds_3(bool useInitPar) {
	const int Nfitfunc = 2;	// 3gaus for 2 non-enhanced 1Ds
	
	TH1D *h1[Nfitfunc];
	if(!h1dx || !h1dy) {cout<<"ERR!! no projection. Call bool Project1Ds() first!"<<endl; return false;}
	h1[0] = h1dx;
	h1[1] = h1dy;

	TF1 *f1[Nfitfunc];	// 3-Gaussians
	for(int i = 0; i<Nfitfunc; i++) {
		f1[i] = new TF1(Form("fun%d",i+1),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
	}

	ROOT::Math::WrappedMultiTF1 *wf1[Nfitfunc];
	for(int i = 0; i<Nfitfunc; i++) {
		wf1[i] = new ROOT::Math::WrappedMultiTF1(*f1[i],1);
	}
	
	ROOT::Fit::DataOptions opt;
	
	ROOT::Fit::DataRange range;

	ROOT::Fit::BinData *data[Nfitfunc];

	ROOT::Fit::Chi2Function *chi2[Nfitfunc];

	for(int i = 0; i<Nfitfunc; i++) {
		data[i] = new ROOT::Fit::BinData(opt,range);
	}
	for(int i = 0; i<Nfitfunc; i++) {
		ROOT::Fit::FillData(*(data[i]),h1[i]);
		chi2[i] = new ROOT::Fit::Chi2Function(*(data[i]),*(wf1[i]));
	}
	
	GlobalChi2_3 globalChi2_3(*chi2[0],*chi2[1]);

	ROOT::Fit::Fitter fitter;
	
	const int Npar = 15;
	double par0[Npar] = {1e4,1e4,1e4,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
		}
	}
	
	fitter.Config().SetParamsSettings(Npar,par0);

	for(int i = 3; i<Npar; i++) {
		fitter.Config().ParSettings(i).Fix();
	}
	fitter.Config().ParSettings(0).SetLimits(1e5,1e7);
	fitter.Config().ParSettings(0).SetStepSize(5);

	fitter.Config().MinimizerOptions().SetPrintLevel(10);
	fitter.Config().SetMinimizer("Minuit2","Migrad");

	// fit FCN function directly
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	int totalsize = 0;
	for(int i = 0; i<Nfitfunc; i++) {
		totalsize+=data[i]->Size();
	}
	fitter.FitFCN(Npar,globalChi2_3,0,totalsize,true);
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);

	if(!result.IsValid()) {
		cout<<"CombinedFit3 Failed!!"<<endl;
		return false;
	}

	// 2 1-D 3-gaus (each 9 parameters) 
	// location of pars for each fit function
	int pars1d[Nfitfunc][9] = {
			{0,1,2,3,4,5,6,7,8},	
			{0,1,2,9,10,11,12,13,14}
			   };

	for(int i = 0; i<Nfitfunc; i++) {
		cout<<"set "<<i<<endl;
		f1[i]->SetFitResult(result,pars1d[i]);
		//f1[i]->SetRange(range().first,range().second);
	}
	cout<<"set done"<<endl;

	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f1[0]->GetParameter(i);
		eYield[i] = f1[0]->GetParError(i);
		meanx[i] = f1[0]->GetParameter(i+NSpecies);
		sigmax[i] = f1[0]->GetParameter(i+NSpecies*2);
		meany[i] = f1[1]->GetParameter(i+NSpecies);
		sigmay[i] = f1[1]->GetParameter(i+NSpecies*2);
	}

	h1dx->GetListOfFunctions()->Add(f1[0]);
	h1dy->GetListOfFunctions()->Add(f1[1]);

	func1dx = (TF1*)f1[0]->Clone("func1dx");
	func1dy = (TF1*)f1[1]->Clone("func1dx");

	TCanvas *c = new TCanvas();
	gStyle->SetOptFit(1);
	h1dx->Draw();
	cout<<Yield[0]<<", "<<meanx[0]<<", "<<sigmax[0]<<", "<<Yield[1]<<", "<<meanx[1]<<", "<<sigmax[1]<<", "<<Yield[2]<<", "<<meanx[2]<<", "<<sigmax[2]<<endl;
	c->SaveAs("h2d_px.png");
	c->Clear();
	h1dy->Draw();
	c->SaveAs("h2d_py.png");

	return true;
}


//================================================ 2D fit 3 gaus ====================================================
bool FitPid2D::Fit2D_3gaus(bool useInitPar) {
	if(!h2d) {
		cout<<"Err: bool FitPid2D::Fit2D_3gaus(bool useInitPar). Need call Project1Ds() first"<<endl;
		return false;
	}

	TH2D *h2d_copy = (TH2D*)h2d->Clone("h2d_copy");
	h2d_copy->Scale(1./(h2d_copy->GetXaxis()->GetBinWidth(0)*h2d_copy->GetYaxis()->GetBinWidth(0)));

	//TF2 *f2 = new TF2("fun0","[0]/((2*pi)*[6]*[12])*exp(-0.5*(pow((x-[3])/[6],2)+pow((y-[9])/[12],2)))+[1]/((2*pi)*[7]*[13])*exp(-0.5*(pow((x-[4])/[7],2)+pow((y-[10])/[13],2)))+[2]/((2*pi)*[8]*[14])*exp(-0.5*(pow((x-[5])/[8],2)+pow((y-[11])/[14],2)))",-10,10,-10,10);
	TF2 *f2 = new TF2("fun0","[0]*TMath::Gaus(x,[3],[6],1)*TMath::Gaus(y,[9],[12],1)+[1]*TMath::Gaus(x,[4],[7],1)*TMath::Gaus(y,[10],[13],1)+[2]*TMath::Gaus(x,[5],[8],1)*TMath::Gaus(y,[11],[14],1)",-2,2,-2,2);
			//[0]/((2*pi)*[6]*[12])*exp(-0.5*(pow((x-[3])/[6],2)+pow((y-[9])/[12],2)))+[1]/((2*pi)*[7]*[13])*exp(-0.5*(pow((x-[4])/[7],2)+pow((y-[10])/[13],2)))+[2]/((2*pi)*[8]*[14])*exp(-0.5*(pow((x-[5])/[8],2)+pow((y-[11])/[14],2)))",-10,10,-10,10);
	f2->SetNpx(100);
	f2->SetNpy(100);

	const int Npar = 15;
	double par0[Npar] = {1e4,1e4,1e4,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
		}
	}
	
	f2->SetParameters(par0);

	int fitStatus = h2d_copy->Fit(f2,"R");
	if(fitStatus!=0) {cout<<" Fit2D_3gaus failed!!!"<<endl; return false;}
	for(int i = 0; i<15; i++) {
		cout<<f2->GetParameter(i)<<", ";
	}
	cout<<endl;
	h2d_copy->GetListOfFunctions()->Add(f2);

	func2d = (TF2*)f2->Clone("func2d");
	func2d->SetParameters(f2->GetParameters());
	func2d->SetParameter(0,f2->GetParameter(0)*(h2d_copy->GetXaxis()->GetBinWidth(0)*h2d_copy->GetYaxis()->GetBinWidth(0)));
	h2d->GetListOfFunctions()->Add(func2d);


	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f2->GetParameter(i);
		eYield[i] = f2->GetParError(i);
		meanx[i] = f2->GetParameter(i+NSpecies);
		sigmax[i] = f2->GetParameter(i+NSpecies*2);
		meany[i] = f2->GetParameter(i+NSpecies*3);
		sigmay[i] = f2->GetParameter(i+NSpecies*4);
	}
	chi = f2->GetChisquare();
	ndf = f2->GetNDF();

	TCanvas *c = new TCanvas();
	gStyle->SetOptFit(1);

	c->SetLogz();
	h2d_copy->Draw("surf2");
	c->SaveAs("h2d.png");

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	c->SetLogy();
	h1dx->Draw();
	//TF1 *f3gausx = new TF1("f3gausx","[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-2,2);
	//TF1 *f3gausx = new TF1("f3gausx","gaus+gaus(3)+gaus(6)",-2,2);
	TF1 *f3gausx = new TF1("f3gausx","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	cout<<Yield[0]<<", "<<meanx[0]<<", "<<sigmax[0]<<", "<<Yield[1]<<", "<<meanx[1]<<", "<<sigmax[1]<<", "<<Yield[2]<<", "<<meanx[2]<<", "<<sigmax[2]<<endl;
	f3gausx->SetParameters(Yield[0],meanx[0],sigmax[0],Yield[1],meanx[1],sigmax[1],Yield[2],meanx[2],sigmax[2]);
	f3gausx->Draw("same");
	h1dx->GetListOfFunctions()->Add(f3gausx);
	func1dx = (TF1*)f3gausx->Clone("func1dx");
	c->SaveAs("h2d_px.png");
	c->Clear();
	h1dy->Draw();
	//TF1 *f3gausy = new TF1("f3gausy","[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-2,2);
	//TF1 *f3gausy = new TF1("f3gausy","gaus+gaus(3)+gaus(6)",-2,2);
	TF1 *f3gausy = new TF1("f3gausy","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	f3gausy->SetParameters(Yield[0],meany[0],sigmay[0],Yield[1],meany[1],sigmay[1],Yield[2],meany[2],sigmay[2]);
	f3gausy->Draw("same");
	h1dy->GetListOfFunctions()->Add(f3gausy);
	func1dy = (TF1*)f3gausy->Clone("func1dy");
	f3gausy->SetNpx(10000);
	c->SaveAs("h2d_py.png");


	h2d_copy->Delete();

	return true;
}

//================================================ individual Fit: 1gaus ====================================================
bool FitPid2D::SeperateFit1Ds(bool useInitPar) {
	const int Nfitfunc = NSpecies*2+2;

	TF1 *f1[Nfitfunc];	// Nfitfunc 3-Gaussian
	for(int i = 0; i<Nfitfunc; i++) {
		//f1[i] = new TF1(Form("fun%d",i),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
		//f1[i] = new TF1(Form("fun%d",i),"gaus+gaus(3)+gaus(6)",-10,10);
		f1[i] = new TF1("f3gausy","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
		f1[i]->SetNpx(1000);
	}

	const int Npar = 33;
	double par0[Npar] = {1e5,1e5,1e5,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
			for(int j = 0; j<NSpecies; j++) {
				par0[j+NSpecies*i+NSpecies*5] = Yieldx_enhanced[i][j];
				par0[j+NSpecies*i+NSpecies*5+NSpecies*NSpecies] = Yieldy_enhanced[i][j];
			}
		}
	}
	
	int pars[8][9] = {{0,3,6,1,4,7,2,5,8},		// projection X
			{0,9,12,1,10,13,2,11,14},	// projection Y
			{15,3,6,16,4,7,17,5,8},		// projection enhanced pion X
			{18,3,6,19,4,7,20,5,8},		// projection enhanced kaon X
			{21,3,6,22,4,7,23,5,8},		// projection enhanced proton X
			{24,9,12,25,10,13,26,11,14},	// projection enhanced pion Y
			{27,9,12,28,10,13,29,11,14},	// projection enhanced kaon Y
			{30,9,12,31,10,13,32,11,14}};	// projection enhanced proton Y


	for(int i = 0; i<Nfitfunc; i++) {
		//cout<<"set "<<i<<endl;
		for(int j = 0; j<9; j++) {
			//if( j==0 || j==3 || j==6 ) {			// ROOT gaus: [0] is the constant, instead of yield
			//	f1[i]->SetParameter(j,par0[pars[i][j]]*pow(2*pi,0.5)*par0[pars[i][j+2]]);
			//}
			//else {
			//	f1[i]->SetParameter(j,par0[pars[i][j]]);
			//}

			f1[i]->SetParameter(j,par0[pars[i][j]]);
		}
		//f1[i]->SetRange(range().first,range().second);
	}

	//for(int i = 0; i<NSpecies; i++) {
	//	int fitstatus1 = h1x_enhanced[i]->Fit(f1[i+NSpecies]);
	//	int fitstatus2 = h1y_enhanced[i]->Fit(f1[i+NSpecies*2]);
	//	if(fitstatus1!=0 || fitstatus2!=0) { cout<<"Species "<<i<<" fit failed"<<endl; return false;}
	//}
	//for(int i = 0; i<NSpecies; i++) {
	//	meanx[i] = f1[i+NSpecies]->GetParameter(1+i*3);
	//	meany[i] = f1[i+NSpecies*2]->GetParameter(1+i*3);
	//	sigmax[i] = f1[i+NSpecies]->GetParameter(2+i*3);
	//	sigmay[i] = f1[i+NSpecies*2]->GetParameter(2+i*3);
	//	for(int j = 0; j<NSpecies; j++) {
	//		Yieldx_enhanced[i][j] = f1[i+NSpecies]->GetParameter(0+j*3);
	//		Yieldy_enhanced[i][j] = f1[i+NSpecies*2]->GetParameter(0+j*3);;
	//	}
	//}

	int pars1gausx[NSpecies][3] = {
					{15,3,6},
					{19,4,7},
					{23,5,8}
					};

	int pars1gausy[NSpecies][3] = {
					{24,9,12},
					{28,10,13},
					{32,11,14}
					};

	//TF1 *fgaus = new TF1("fgaus","[0]/(pow(2*pi,0.5)*[2])*exp(-0.5*(pow((x-[1])/[2],2)))",-10,10);
	//TF1 *fgaus = new TF1("fgaus","gaus",-10,10);
	TF1 *fgaus = new TF1("fgaus","[0]*TMath::Gaus(x,[1],[2],1)",-2,2);
	fgaus->SetNpx(1000);
	TCanvas *c1gaus = new TCanvas();
	gStyle->SetOptFit(1);
	c1gaus->SetLogy();
	for(int i = 0; i<NSpecies; i++) {
		for(int j=0; j<3; j++) {
			cout<<par0[pars1gausx[i][j]]<<endl;
			fgaus->SetParameter(j,par0[pars1gausx[i][j]]);
		}
		// Set Parameter Limits
		fgaus->SetParLimits(0,0,1e6);
		fgaus->SetParLimits(1,-2,2);
		fgaus->SetParLimits(2,0,1);
		cout<<"Fit "<<i<<"-x:"<<endl;
		int fitstatus1 = h1x_enhanced[i]->Fit(fgaus);
		if(fitstatus1!=0) { cout<<"Species "<<i<<"-x fit failed"<<endl; return false;}
		Yieldx_enhanced[i][i] = fgaus->GetParameter(0);
		//Yieldx_enhanced[i][i] = fgaus->GetParameter(0)*pow(2*pi,0.5)*fgaus->GetParameter(2);
		meanx[i] = fgaus->GetParameter(1);
		sigmax[i] = fgaus->GetParameter(2);

		c1gaus->Clear();
		h1x_enhanced[i]->Draw();
		c1gaus->SaveAs(Form("h1x_enhanced_1g_%d.png",i));

		for(int j=0; j<3; j++) {
			cout<<par0[pars1gausy[i][j]]<<endl;
			fgaus->SetParameter(j,par0[pars1gausy[i][j]]);
		}
		// Set Parameter Limits
		fgaus->SetParLimits(0,0,1e6);
		fgaus->SetParLimits(1,-2,2);
		fgaus->SetParLimits(2,0,1);
		cout<<"Fit "<<i<<"-y:"<<endl;
		int fitstatus2 = h1y_enhanced[i]->Fit(fgaus);
		if(fitstatus2!=0) { cout<<"Species "<<i<<"-y fit failed"<<endl; return false;}
		Yieldy_enhanced[i][i] = fgaus->GetParameter(0);
		//Yieldy_enhanced[i][i] = fgaus->GetParameter(0)*pow(2*pi,0.5)*fgaus->GetParameter(2);
		meany[i] = fgaus->GetParameter(1);
		sigmay[i] = fgaus->GetParameter(2);

		c1gaus->Clear();
		h1y_enhanced[i]->Draw();
		c1gaus->SaveAs(Form("h1y_enhanced_1g_%d.png",i));
	}

	for(int i = 0; i<NSpecies; i++) {
		f1[0]->FixParameter(1+i*3,meanx[i]);
		f1[1]->FixParameter(1+i*3,meany[i]);
		f1[0]->FixParameter(2+i*3,sigmax[i]);
		f1[1]->FixParameter(2+i*3,sigmay[i]);
		f1[0]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmax[i]));
		//f1[0]->SetParameter(1+i*3,meanx[i]);
		//f1[0]->SetParameter(2+i*3,sigmax[i]);
		//f1[0]->SetParLimits(2+i*3,0,1);

		f1[1]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmay[i]));
		//f1[1]->SetParameter(1+i*3,meany[i]);
		//f1[1]->SetParameter(2+i*3,sigmay[i]);
		//f1[1]->SetParLimits(2+i*3,0,1);
	}
	
	TCanvas *c1dx = new TCanvas();
	c1dx->SetLogy();
	gStyle->SetOptFit(1);
	h1dx->Fit(f1[0]);	
	c1dx->SaveAs("h1dx.png");
	for(int i = 0; i<NSpecies; i++) {
		//Yield[i] = f1[0]->GetParameter(0+i*3);
		meanx[i] = f1[0]->GetParameter(1+i*NSpecies);
		sigmax[i] = f1[0]->GetParameter(2+i*NSpecies);
		Yield[i] = f1[0]->GetParameter(0+i*NSpecies)*pow(2*pi,0.5)*sigmax[i];
		cout<<"x Species "<<i<<" Yield = "<<Yield[i]<<endl;
	}
	TCanvas *c1dy = new TCanvas();
	c1dy->SetLogy();
	h1dy->Fit(f1[1]);	
	c1dy->SaveAs("h1dy.png");
	for(int i = 0; i<NSpecies; i++) {
		meany[i] = f1[1]->GetParameter(1+i*NSpecies);
		sigmay[i] = f1[1]->GetParameter(2+i*NSpecies);
		//Yield[i] += f1[1]->GetParameter(0+i*3)*pow(2*pi,0.5)*sigmay[i];
		Yield[i] = f1[1]->GetParameter(0+i*3)*pow(2*pi,0.5)*sigmay[i];
		cout<<"y Species "<<i<<" Yield = "<<f1[1]->GetParameter(0+i*3)<<endl;
	}
	//cout<<"Average:"<<endl;
	//for(int i = 0; i<NSpecies; i++) {
	//	Yield[i] /= 2.;
	//	cout<<"Species "<<i<<" Yield = "<<Yield[i]<<endl;
	//}

	cout<<endl;

	return true;
	

}


//================================================ individual Fit: 3gaus ====================================================
bool FitPid2D::SeperateFit1Ds_3gaus(bool useInitPar) {
	const int Nfitfunc = NSpecies*2+2;

	TF1 *f1[Nfitfunc];	// Nfitfunc 3-Gaussian
	for(int i = 0; i<Nfitfunc; i++) {
		//f1[i] = new TF1(Form("fun%d",i),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
		//f1[i] = new TF1(Form("fun%d",i),"gaus+gaus(3)+gaus(6)",-10,10);
		f1[i] = new TF1("f3gaus","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	}

	const int Npar = 33;
	double par0[Npar] = {1e5,1e5,1e5,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
			for(int j = 0; j<NSpecies; j++) {
				par0[j+NSpecies*i+NSpecies*5] = Yieldx_enhanced[i][j];
				par0[j+NSpecies*i+NSpecies*5+NSpecies*NSpecies] = Yieldy_enhanced[i][j];
			}
		}
	}
	
	int pars[8][9] = {{0,3,6,1,4,7,2,5,8},		// projection X
			{0,9,12,1,10,13,2,11,14},	// projection Y
			{15,3,6,16,4,7,17,5,8},		// projection enhanced pion X
			{18,3,6,19,4,7,20,5,8},		// projection enhanced kaon X
			{21,3,6,22,4,7,23,5,8},		// projection enhanced proton X
			{24,9,12,25,10,13,26,11,14},	// projection enhanced pion Y
			{27,9,12,28,10,13,29,11,14},	// projection enhanced kaon Y
			{30,9,12,31,10,13,32,11,14}};	// projection enhanced proton Y


	for(int i = 0; i<Nfitfunc; i++) {
		//cout<<"set "<<i<<endl;
		for(int j = 0; j<9; j++) {
			//if( j==0 || j==3 || j==6 ) {			// ROOT gaus: [0] is the constant, instead of yield
			//	f1[i]->SetParameter(j,par0[pars[i][j]]*pow(2*pi,0.5)*par0[pars[i][j+2]]);
			//}
			//else {
			//	f1[i]->SetParameter(j,par0[pars[i][j]]);
			//}

			f1[i]->SetParameter(j,par0[pars[i][j]]);
		}
		//f1[i]->SetRange(range().first,range().second);
		f1[i]->SetParLimits(0,0,1e8);
		f1[i]->SetParLimits(1,-2,0);
		f1[i]->SetParLimits(2,0,1);
		f1[i]->SetParLimits(3,0,1e8);
		f1[i]->SetParLimits(4,-1,1);
		f1[i]->SetParLimits(5,0,1);
		f1[i]->SetParLimits(6,0,1e8);
		f1[i]->SetParLimits(7,0,2);
		f1[i]->SetParLimits(8,0,1);
	}

	for(int i = 0; i<NSpecies; i++) {
		int fitstatus1 = h1x_enhanced[i]->Fit(f1[i+2],"R");
		int fitstatus2 = h1y_enhanced[i]->Fit(f1[i+NSpecies+2],"R");
		if(fitstatus1!=0 || fitstatus2!=0) { cout<<"Species "<<i<<" fit failed"<<endl; return false;}
	}
	for(int i = 0; i<NSpecies; i++) {
		meanx[i] = f1[i+2]->GetParameter(1+i*3);
		meany[i] = f1[i+NSpecies+2]->GetParameter(1+i*3);
		sigmax[i] = f1[i+2]->GetParameter(2+i*3);
		sigmay[i] = f1[i+NSpecies+2]->GetParameter(2+i*3);
		for(int j = 0; j<NSpecies; j++) {
			Yieldx_enhanced[i][j] = f1[i+2]->GetParameter(0+j*3);
			Yieldy_enhanced[i][j] = f1[i+NSpecies+2]->GetParameter(0+j*3);;
		}
	}
	TCanvas *c1 = new TCanvas();
	c1->SetLogy();
	for(int i = 0; i<NSpecies; i++) {
		h1x_enhanced[i]->Draw();
		f1[i+2]->Draw("same");
		c1->SaveAs(Form("h1x_enhanced_3g_%d.png",i));
		h1x_enhanced[i]->Draw();

		c1->Clear();
		h1y_enhanced[i]->Draw();
		f1[i+NSpecies+2]->Draw("same");
		c1->SaveAs(Form("h1y_enhanced_3g_%d.png",i));
	}


	// Now fit non-enhanced one
	for(int i = 0; i<NSpecies; i++) {
		f1[0]->FixParameter(1+i*3,meanx[i]);
		f1[1]->FixParameter(1+i*3,meany[i]);
		f1[0]->FixParameter(2+i*3,sigmax[i]);
		f1[1]->FixParameter(2+i*3,sigmay[i]);
		f1[0]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmax[i]));
		//f1[0]->SetParameter(1+i*3,meanx[i]);
		//f1[0]->SetParameter(2+i*3,sigmax[i]);
		//f1[0]->SetParLimits(2+i*3,0,1);

		f1[1]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmay[i]));
		//f1[1]->SetParameter(1+i*3,meany[i]);
		//f1[1]->SetParameter(2+i*3,sigmay[i]);
		//f1[1]->SetParLimits(2+i*3,0,1);
	}
	
	TCanvas *c1dx = new TCanvas();
	c1dx->SetLogy();
	gStyle->SetOptFit(1);
	h1dx->Fit(f1[0]);	
	c1dx->SaveAs("h1dx.png");
	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f1[0]->GetParameter(0+i*3);
		meanx[i] = f1[0]->GetParameter(1+i*NSpecies);
		sigmax[i] = f1[0]->GetParameter(2+i*NSpecies);
		//Yield[i] = f1[0]->GetParameter(0+i*NSpecies)*pow(2*pi,0.5)*sigmax[i];
		cout<<"x Species "<<i<<" Yield = "<<Yield[i]<<endl;
	}
	chi = f1[0]->GetChisquare();
	ndf = f1[0]->GetNDF();
	//cout<<"chi2: "<<chi<<" ndf:"<<ndf<<endl;
	TCanvas *c1dy = new TCanvas();
	c1dy->SetLogy();
	h1dy->Fit(f1[1]);	
	c1dy->SaveAs("h1dy.png");
	for(int i = 0; i<NSpecies; i++) {
		meany[i] = f1[1]->GetParameter(1+i*NSpecies);
		sigmay[i] = f1[1]->GetParameter(2+i*NSpecies);
		//Yield[i] += f1[1]->GetParameter(0+i*3)*pow(2*pi,0.5)*sigmay[i];
		Yield[i] = f1[1]->GetParameter(0+i*3)*pow(2*pi,0.5)*sigmay[i];
		cout<<"y Species "<<i<<" Yield = "<<f1[1]->GetParameter(0+i*3)<<endl;
	}
	//cout<<"Average:"<<endl;
	//for(int i = 0; i<NSpecies; i++) {
	//	Yield[i] /= 2.;
	//	cout<<"Species "<<i<<" Yield = "<<Yield[i]<<endl;
	//}

	cout<<endl;

	return true;
	

}



//================================================ individual Fit: 3gaus or 1gaus for x (zd), 1gaus for y(zb) ====================================================
bool FitPid2D::SeperateFit1Ds_x3g_y1g(bool useInitPar) {
	const int Nfitfunc = NSpecies*3;

	TF1 *f3gausx[NSpecies];	// Nfitfunc 3-Gaussian for x (zd)
	TF1 *f1gausx[NSpecies];	// Nfitfunc 1-Gaussian for x (zd)
	TF1 *f1gausy[NSpecies];	// Nfitfunc 1-Gaussian for y (zb)
	for(int i = 0; i<NSpecies; i++) {
		//f1[i] = new TF1(Form("fun%d",i),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
		//f1[i] = new TF1(Form("fun%d",i),"gaus+gaus(3)+gaus(6)",-10,10);
		f3gausx[i] = new TF1("f3gausx","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
		f1gausx[i] = new TF1("f1gausx","[0]*TMath::Gaus(x,[1],[2],1)",-2,2);
		f1gausy[i] = new TF1("f1gausy","[0]*TMath::Gaus(x,[1],[2],1)",-2,2);

		f3gausx[i]->SetNpx(1000);
		f1gausx[i]->SetNpx(1000);
		f1gausy[i]->SetNpx(1000);
	}

	const int Npar = 33;
	double par0[Npar] = {1e5,1e5,1e5,-0.5,0,0.5,0.1,0.1,0.1,-0.5,0,0.5,0.01,0.01,0.01,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5};
	if(useInitPar) {
		for(int i = 0; i<NSpecies; i++) {
			par0[i] = Yield[i];
			par0[i+NSpecies] = meanx[i];
			par0[i+NSpecies*2] = sigmax[i];
			par0[i+NSpecies*3] = meany[i];
			par0[i+NSpecies*4] = sigmay[i];
			for(int j = 0; j<NSpecies; j++) {
				par0[j+NSpecies*i+NSpecies*5] = Yieldx_enhanced[i][j];
				par0[j+NSpecies*i+NSpecies*5+NSpecies*NSpecies] = Yieldy_enhanced[i][j];
			}
		}
	}
	
	int pars[Nfitfunc][9] = { 
				{15,3,6,16,4,7,17,5,8},		// projection enhanced pion X (3gaus)
				{18,3,6,19,4,7,20,5,8},		// projection enhanced kaon X (3gaus)
				{21,3,6,22,4,7,23,5,8},		// projection enhanced proton X (3gaus)
				{15,3,6},		// projection enhanced pion X (1gaus)
				{19,4,7},		// projection enhanced kaon X (1gaus)
				{23,5,8},		// projection enhanced proton X (1gaus)
				{24,9,12},	// projection enhanced pion Y
				{28,10,13},	// projection enhanced kaon Y
				{32,11,14}};	// projection enhanced proton Y



	TCanvas *c1gaus = new TCanvas();
	gStyle->SetOptFit(1);
	c1gaus->SetLogy();
	for(int i=0; i<NSpecies; i++) {
		// Set Parameter Limits
		cout<<"f1gauxs"<<i<<" Yield range<"<<h1x_enhanced[i]->Integral("width")*2<<endl;
		//f1gausx[i]->SetParLimits(0,0,h1x_enhanced[i]->Integral("width")*2);
		//f1gausx[i]->SetParLimits(1,-2,2);
		//f1gausx[i]->SetParLimits(2,0,0.5);
		cout<<"Fit "<<i<<"-x (1gaus):"<<endl;
		for(int j = 0; j<3; j++) {
			f1gausx[i]->SetParameter(j,par0[pars[i+3][j]]);
			cout<<"f1gaux"<<i<<" par"<<j<<" = "<<par0[pars[i+3][j]]<<endl;
		}
		int fitstatusx = h1x_enhanced[i]->Fit(f1gausx[i],"R");
		if(fitstatusx==0) {
			Yieldx_enhanced[i][i]=f1gausx[i]->GetParameter(0);
			meanx[i]=f1gausx[i]->GetParameter(1);
			sigmax[i]=f1gausx[i]->GetParameter(2);
		}
		else {
			cout<<"Species "<<i<<"-x 1-gaus fit failed. Switch to 3-gaus"<<endl;
			for(int j = 0; j<9; j++) {
				f3gausx[i]->SetParameter(j,par0[pars[i][j]]);
				cout<<"f3gaux"<<i<<" par"<<j<<" = "<<par0[pars[i][j]]<<endl;
			}
			double maxy = h1x_enhanced[i]->Integral("width")*2;
			f3gausx[i]->SetParLimits(0,0,maxy);
			f3gausx[i]->SetParLimits(1,-2,0);
			f3gausx[i]->SetParLimits(2,0,0.5);
			f3gausx[i]->SetParLimits(3,0,maxy);
			f3gausx[i]->SetParLimits(4,-1,0.5);
			f3gausx[i]->SetParLimits(5,0,0.5);
			f3gausx[i]->SetParLimits(6,0,maxy);
			f3gausx[i]->SetParLimits(7,0,2);
			f3gausx[i]->SetParLimits(8,0,0.5);

			cout<<"Fit "<<i<<"-x:"<<endl;
			fitstatusx = h1x_enhanced[i]->Fit(f3gausx[i],"R");
			if(fitstatusx!=0) { cout<<"Species "<<i<<"-x fit (3gaus) failed!!"<<endl; }//return false;}
			else {
				Yieldx_enhanced[i][i] = f3gausx[i]->GetParameter(i*3);
				//Yieldx_enhanced[i][i] = fgaus->GetParameter(0)*pow(2*pi,0.5)*fgaus->GetParameter(2);
				meanx[i] = f3gausx[i]->GetParameter(i*3+1);
				sigmax[i] = f3gausx[i]->GetParameter(i*3+2);
			}
		}

		if(fitstatusx==0) {
			c1gaus->Clear();
			h1x_enhanced[i]->Draw();
			c1gaus->SaveAs(Form("h1x_enhanced_1g_%d.png",i));
		}

		//for(int j=0; j<3; j++) {
		//	cout<<par0[pars1gausy[i][j]]<<endl;
		//	fgaus->SetParameter(j,par0[pars1gausy[i][j]]);
		//}
		// Set Parameter Limits
		//f1gausy[i]->SetParLimits(0,0,h1y_enhanced[i]->Integral("width")*2);
		//f1gausy[i]->SetParLimits(1,-2,2);
		//f1gausy[i]->SetParLimits(2,0,1);


		cout<<"Fit "<<i<<"-y:"<<endl;
		for(int j = 0; j<3; j++) {
			f1gausy[i]->SetParameter(j,par0[pars[i+6][j]]);
			cout<<"f1gauy"<<i<<" par"<<j<<" = "<<par0[pars[i+6][j]]<<endl;
		}
		int fitstatusy = h1y_enhanced[i]->Fit(f1gausy[i],"R");
		if(fitstatusy!=0) { cout<<"Species "<<i<<"-y fit failed"<<endl; }//return false;}
		else {
			Yieldy_enhanced[i][i] = f1gausy[i]->GetParameter(0);
			//Yieldy_enhanced[i][i] = fgaus->GetParameter(0)*pow(2*pi,0.5)*fgaus->GetParameter(2);
			meany[i] = f1gausy[i]->GetParameter(1);
			sigmay[i] = f1gausy[i]->GetParameter(2);

			c1gaus->Clear();
			h1y_enhanced[i]->Draw();
			c1gaus->SaveAs(Form("h1y_enhanced_1g_%d.png",i));
		}
	}


	cout<<endl;

	return true;
	

}


//================================================ 1D Fit for non-enhancec X 3gaus ====================================================
bool FitPid2D::Fit1Dx(bool useInitPar) {

	const int Nfitfunc = NSpecies;

	TF1 *f1[Nfitfunc];	// Nfitfunc 3-Gaussian
	for(int i = 0; i<Nfitfunc; i++) {
		//f1[i] = new TF1(Form("fun%d",i),"[0]/(pow(2*pi,0.5)*[6])*exp(-0.5*(pow((x-[3])/[6],2)))+[1]/(pow(2*pi,0.5)*[7])*exp(-0.5*(pow((x-[4])/[7],2)))+[2]/(pow(2*pi,0.5)*[8])*exp(-0.5*(pow((x-[5])/[8],2)))",-10,10);
		//f1[i] = new TF1(Form("fun%d",i),"gaus+gaus(3)+gaus(6)",-10,10);
		f1[i] = new TF1("f3gausx","[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1)+[6]*TMath::Gaus(x,[7],[8],1)",-2,2);
	}

	// Now fit non-enhanced one
	for(int i = 0; i<NSpecies; i++) {
		f1[0]->FixParameter(1+i*3,meanx[i]);
		f1[1]->FixParameter(1+i*3,meany[i]);
		f1[0]->FixParameter(2+i*3,sigmax[i]);
		f1[1]->FixParameter(2+i*3,sigmay[i]);
		f1[0]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmax[i]));
		//f1[0]->SetParameter(1+i*3,meanx[i]);
		//f1[0]->SetParameter(2+i*3,sigmax[i]);
		//f1[0]->SetParLimits(2+i*3,0,1);

		f1[1]->SetParameter(0+i*3,Yield[i]/(pow(2*pi,0.5)*sigmay[i]));
		//f1[1]->SetParameter(1+i*3,meany[i]);
		//f1[1]->SetParameter(2+i*3,sigmay[i]);
		//f1[1]->SetParLimits(2+i*3,0,1);
	}
	
	TCanvas *c1dx = new TCanvas();
	c1dx->SetLogy();
	gStyle->SetOptFit(1);
	h1dx->Fit(f1[0]);	
	c1dx->SaveAs("h1dx.png");
	for(int i = 0; i<NSpecies; i++) {
		Yield[i] = f1[0]->GetParameter(0+i*3);
		meanx[i] = f1[0]->GetParameter(1+i*NSpecies);
		sigmax[i] = f1[0]->GetParameter(2+i*NSpecies);
		eYield[i] = f1[0]->GetParError(0+i*3);
		//Yield[i] = f1[0]->GetParameter(0+i*NSpecies)*pow(2*pi,0.5)*sigmax[i];
		cout<<"x Species "<<i<<" Yield = "<<Yield[i]<<endl;
	}
	chi = f1[0]->GetChisquare();
	ndf = f1[0]->GetNDF();
	cout<<"chi2: "<<chi<<" ndf:"<<ndf<<endl;
	
	cout<<endl;

	return true;
}

//================================================ Get Fit Output ====================================================
double* FitPid2D::GetYields() {
	return Yield;
}
double* FitPid2D::GetYieldsError() {
	return eYield;
}
double* FitPid2D::GetParameters(){
	for(int i = 0; i<NSpecies; i++) {	
		Parameters[i] = Yield[i];
		Parameters[NSpecies+i] = meanx[i];
		Parameters[2*NSpecies+i] = sigmax[i];
		Parameters[3*NSpecies+i] = meany[i];
		Parameters[4*NSpecies+i] = sigmay[i];
			for(int j = 0; j<NSpecies; j++) {	
				Parameters[5*NSpecies+i*NSpecies+j] = Yieldx_enhanced[i][j];
				Parameters[5*NSpecies+NSpecies*NSpecies+i*NSpecies+j] = Yieldy_enhanced[i][j];
			}
	}
	return Parameters;
}


//================================================ Write Fit Output ====================================================
bool FitPid2D::Write2Root(const char* outputfilename) {
	TFile *fout = new TFile(outputfilename,"RECREATE");
	if(!fout->IsOpen()) {cout<<"Cannot write to "<<outputfilename<<endl; return false;}
	h2d->Write();
	h1dx->Write();		
	h1dy->Write();		
	for(int i = 0; i<NSpecies; i++) {
		h1x_enhanced[i]->Write();
		h1y_enhanced[i]->Write();
	}
	if(func2d) func2d->Write();
	if(func1dx) func1dx->Write();
	if(func1dy) func1dy->Write();
	for(int i = 0; i<NSpecies; i++) {
		if(func1x_enhanced[i]) func1x_enhanced[i]->Write();
		if(func1y_enhanced[i]) func1y_enhanced[i]->Write();
	}
	fout->Close();
	return true;
}
