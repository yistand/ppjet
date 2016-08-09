//======================================================================================================================================================
//
//		2016.07.27		Li Yi
//		estimate Jet resolution per jet pt. so that we will know what jet pt bin shall be used
//
//
//======================================================================================================================================================

#if !(defined(__CINT__) || defined(__CLING__))
#include <iostream>
using std::cout;
using std::endl;

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"


#endif

void JetResolution() {

	TFile *fin = new TFile("UnfoldMatrxLead.root");
	TH2D *cov = (TH2D*)fin->Get("CovMatrix");
	int Nbins = cov->GetNbinsY();
	float ymax = cov->ProjectionY()->GetBinLowEdge(Nbins+1);
	TH1D *ptmp;

	TH1D *hjrs = new TH1D("hjrs","jet resolution vs pT",Nbins,0,ymax);		// output histogram
	TH1D *hA = new TH1D("hA","A vs pT",Nbins,0,ymax);		// test histogram

	TF1 *f1 = new TF1("f1","[0]*exp(-pow(x-[1],2)/(2*[2]*[2]))",0,ymax);		// gaussian

	for(int i = 0; i<Nbins; i++) {
		ptmp = (TH1D*)cov->ProjectionX(Form("McJet%d",i),i+1,i+1);
		if(ptmp->Integral()>0) {
			f1->SetParameters(1,i,1);
			ptmp->Fit(f1);
		}
		float sigma = f1->GetParameter(2);
		float esigma = f1->GetParError(2);
		hjrs->SetBinContent(i+1,fabs(sigma));	
		hjrs->SetBinError(i+1,esigma);	
		hA->SetBinContent(i+1,f1->GetParameter(0));	
		hA->SetBinError(i+1,f1->GetParError(0));	
	}

	hjrs->Draw();
	hA->Draw();
}



#ifndef __CINT__
int main () { JetResolution(); return 0; }  // Main program when run stand-alone
#endif

