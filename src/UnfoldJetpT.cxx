//======================================================================================================================================================
//
//		2016.07.25		Li YI
//		use RooUnfold to unfold leading jet pT using Mc - Rc from embedding
//
//======================================================================================================================================================
#if !(defined(__CINT__) || defined(__CLING__))
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#endif


void UnfoldJetpT() {

	cout << "==================================== TRAIN ====================================" << endl;
	TFile *fin = new TFile("UnfoldMatrxLead.root");
	TH1D *measured = (TH1D*)fin->Get("Rc");
	TH1D *truth = (TH1D*)fin->Get("Mc");
	TH2D *cov = (TH2D*)fin->Get("CovMatrix");

	RooUnfoldResponse response (measured, truth, cov);

	cout << "==================================== self TEST =====================================" << endl;
	RooUnfoldBayes   unfold (&response, measured, 4);  
	TH1D* reco= (TH1D*) unfold.Hreco();


	TCanvas *c1 = new TCanvas("ctest");
	reco->SetLineColor(1);
	reco->Draw("h");
	measured->SetLineColor(kBlue);
	measured->Draw("psame");
	truth->SetLineColor(kRed);
	truth->Draw("psame");


	cout << "==================================== UNFOLDING  =====================================" << endl;
	TFile *freal = new TFile("/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP2_151030_P12id.root");
	TH1D *LeadJetPt_data = (TH1D*)freal->Get("LeadingJetPt");
	RooUnfoldBayes   unfold_data (&response, LeadJetPt_data , 4);
	TH1D *UnfoldedLeadJetPt_data = (TH1D*)unfold_data.Hreco();
	

	TCanvas *c2 = new TCanvas("cunfold");
	UnfoldedLeadJetPt_data->SetLineColor(1);
	UnfoldedLeadJetPt_data->Draw("h");
	LeadJetPt_data->SetLineColor(kBlue);
	LeadJetPt_data->Draw("hpsame");
	truth->SetLineColor(kRed);
	truth->Draw("psame");
	measured->SetLineColor(kBlue-2);
	measured->Draw("psame");
}



#ifndef __CINT__
int main () { UnfoldJetpT(); return 0; }  // Main program when run stand-alone
#endif
