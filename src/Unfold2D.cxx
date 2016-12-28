//=====================================================================-*-C++-*-
//
//	2016.08.23		Li YI
//	try to unfold underlying event activity vs leading jet pt
//
//==============================================================================
//
//
// File and Version Information:
//      $Id: RooUnfoldTestHarness2D.icc 336 2012-06-12 19:39:44Z T.J.Adye $
//
// Description:
//      Test Harness class for the RooUnfold package using 2D toy MC.
//      Inherits from RooUnfoldTestHarness.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef __UNFOLD2D_CXX
#define __UNFOLD2D_CXX

#include "Unfold2D.hh"
#include "CrossSectionPerpT.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <cmath>
#include <iostream>
#include <algorithm>

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TLine.h"
#include "TStyle.h"

#include "RooUnfoldErrors.h"
#include "RooUnfoldParms.h"
#include "RooUnfoldResponse.h"
#include "RooUnfold.h"
#endif


#if !defined(__CINT__) || defined(__MAKECINT__)
using std::cerr;
using std::cout;
using std::endl;
using std::sin;
using std::cos;
#endif



//==============================================================================
// parameters
//==============================================================================
//double Unfold2D::Wptbins[] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};
double Unfold2D::Wptbins[] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,100};

void Unfold2D::Help() {

	cout<<"For parameters set up:"<<endl;
	cout<<"method = args[0];			// use 1 if not sure"<<endl;
	cout<<"// 0	  no unfolding (output copied from measured input)"<<endl;
	cout<<"// 1	  Bayes"<<endl;
	cout<<"// 2	  SVD"<<endl;
	cout<<"// 3	  bin-by-bin"<<endl;
	cout<<"// 4	  TUnfold"<<endl;
	cout<<"// 5	  matrix inversion"<<endl<<endl;

	cout<<"ntx = args[1];				// x-axis bininng for Mc"<<endl;
	cout<<"nty = args[2];				// y-axis bininng for Mc"<<endl;
	cout<<"nmx = args[3];				// x-axis bininng for Rc"<<endl;
	cout<<"nmy = args[4];				// y-axis bininng for Rc"<<endl;
	cout<<"xlo = args[5];				// x-axis minimum"<<endl;
	cout<<"xhi = args[6];				// x-axis maximum"<<endl;
	cout<<"ylo = args[7];				// y-axis minimum"<<endl;
	cout<<"yhi = args[8];				// y-axis maximum"<<endl;
	cout<<"overflow = args[9];			// use histogram under/overflows if 1 (set from RooUnfoldResponse) "<<endl;
	cout<<"verbose = args[10];			// debug print out level"<<endl<<endl;
	cout<<"doerror = args[11];			// use 1 if not sure"<<endl;
	cout<<"// RooUnfold: "<<endl;
	cout<<"// enum ErrorTreatment {  // Error treatment:"<<endl;
	cout<<"//    kNoError,            //   no error treatment: returns sqrt(N)"<<endl;
	cout<<"//    kErrors,             //   bin-by-bin errors (diagonal covariance matrix)"<<endl;
	cout<<"//    kCovariance,         //   covariance matrix from unfolding"<<endl;
	cout<<"//    kCovToy,             //   covariance matrix from toy MC"<<endl;
	cout<<"//    kDefault=-1          //   not specified"<<endl<<endl;
	cout<<"regparm = args[12];			// regularisation parameter (If not sure use: Bayes niter=3, SVD kterm=ntx/2)"<<endl;
	cout<<"WIDEBIN = args[13];			// 1: for use wide pT bins (various bin size). 0: for use fine fixed pT bins"<<endl;

	cout<<endl<<endl;

	cout<<"XvariableName = "<<XvariableName<<endl;
	cout<<"YvariableName = "<<YvariableName<<endl;

	cout<<endl<<endl;
}

void Unfold2D::SetParms (const char* const* argv) 
{
	const int Npar = 14;
	float args[Npar] = {0};
	for(int i = 0; i<Npar; i++) 
	{
		args[i] = atof(argv[i]);
	}
	SetParms(args);

}


int Unfold2D::Float2Int(float aFloat) 
{
	int aInt = 0;
	if(aFloat>=0)	aInt = (int) (aFloat + 0.5);
	else		aInt = (int) (aFloat - 0.5);	
	return aInt;
}


void Unfold2D::SetParms (float* args)
{
	method = args[0];			// use 1 if not sure
	// 0	  no unfolding (output copied from measured input)
	// 1	  Bayes
	// 2	  SVD
	// 3	  bin-by-bin
	// 4	  TUnfold
	// 5	  matrix inversion

	ntx = Float2Int(args[1]);				// x-axis bininng for Mc
	nty = Float2Int(args[2]);				// y-axis bininng for Mc
	nmx = Float2Int(args[3]);				// x-axis bininng for Rc
	nmy = Float2Int(args[4]);				// y-axis bininng for Rc
	xlo = args[5];				// x-axis minimum
	xhi = args[6];				// x-axis maximum
	ylo = args[7];				// y-axis minimum
	yhi = args[8];				// y-axis maximum
	overflow = Float2Int(args[9]);			// use histogram under/overflows if 1 (set from RooUnfoldResponse) 
	verbose = Float2Int(args[10]);			// debug print out level

	doerror = Float2Int(args[11]);			// use 1 if not sure
	// RooUnfold: 
	// enum ErrorTreatment {  // Error treatment:
	//    kNoError,            //   no error treatment: returns sqrt(N)
	//    kErrors,             //   bin-by-bin errors (diagonal covariance matrix)
	//    kCovariance,         //   covariance matrix from unfolding
	//    kCovToy,             //   covariance matrix from toy MC
	//    kDefault=-1          //   not specified
	//  };

	regparm = args[12];			// regularisation parameter (If not sure use: Bayes niter=3, SVD kterm=ntx/2)

	if(args[13]>0) WIDEBIN = true;		// Use wide and various pT bins or not
	else WIDEBIN = false;
}

void Unfold2D::SetDefaultParms() 		// use default parameterization
{
	float args[14] = {1, 60, 100, 60, 100, 0, 60, 0, 1, 0, 1, 1, 4, 0};	// neutralfrac
	//float args[14] = {3, 50, 40, 50, 40, 0, 100, -0.5, 39.5, 0, 1, 1, 4, 0};	// TranNtrk
	//float args[14] = {1, 50, 500, 50, 500, 0, 100, 0, 50, 0, 1, 1, 4, 1};		// PtSum
	//		method, ntx, nty, nmx, nmy, xlo, xhi, ylo, yhi, overflow, verbose, doerror, regparm, WIDEBIN
	SetParms(args);	
}

void Unfold2D::PrintParms() 
{
	cout<<"method = "<<method<<endl;
	cout<<"ntx = "<<ntx<<endl;
	cout<<"nty = "<<nty<<endl;
	cout<<"nmx = "<<nmx<<endl;
	cout<<"nmy = "<<nmy<<endl;
	cout<<"xlo = "<<xlo<<endl;
	cout<<"xhi = "<<xhi<<endl;
	cout<<"ylo = "<<ylo<<endl;
	cout<<"yhi = "<<yhi<<endl;
	cout<<"overflow = "<<overflow<<endl;
	cout<<"verbose = "<<verbose<<endl;
	cout<<"doerror = "<<doerror<<endl;
	cout<<"regparm = "<<regparm<<endl;
	cout<<"WIDEBIN = "<<WIDEBIN<<endl;
}

//==============================================================================
Int_t Unfold2D::FillbyXsec4Train ()
{
	int Nevents[NUMBEROFPT];
	std::fill_n(Nevents, NUMBEROFPT, 1e8);
	return FillbyXsec4Train (Nevents);
}

Int_t Unfold2D::FillbyXsec4Train (int *Nevents)
{
	Bool_t flagMatch2Lead;
	Bool_t flagMatch2Sub;
	Bool_t flagIsTrigger;
	Bool_t flagtrigmatch;

	Float_t Mcjneutralfrac;
	Int_t  McLeadAreaNtrk;
	Int_t  McSubAreaNtrk;
	Int_t  McTranMaxNtrk;
	Int_t  McTranMinNtrk;
	Float_t McLeadAreaPtSum;
	Float_t McSubAreaPtSum;
	Float_t McTranMaxPtSum;
	Float_t McTranMinPtSum;

	Float_t Rcjneutralfrac;
	Int_t  RcTranMaxNtrk;
	Int_t  RcTranMinNtrk;
	Int_t  RcLeadAreaNtrk;
	Int_t  RcSubAreaNtrk;
	Float_t RcLeadAreaPtSum;
	Float_t RcSubAreaPtSum;
	Float_t RcTranMaxPtSum;
	Float_t RcTranMinPtSum;


	if(WIDEBIN) pfxTrain= new TProfile("pfxtrain", "Training", WNbins, Wptbins);
	else pfxTrain= new TProfile ("pfxtrain", "Training", nmx, xlo, xhi);
	pfxTrain->SetLineColor(kRed);
	pfxTrain->GetXaxis()->SetTitle(XvariableName);
	pfxTrain->GetYaxis()->SetTitle(YvariableName);

	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename;
		if(ExcludeOpt&&TrigName.Contains("JP")) {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_excluded.root",TrigName.Data(),PTBINS[i],TranCharge.Data());
		}
		else {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data());
		}
		cout<<"Read in "<<ifilename<<endl;
		ftrain = new TFile(ifilename);
		tree = (TTree*)ftrain->Get("ResultTree");
		if(!tree) {cout<<"Cannot find tree from input in "<<__PRETTY_FUNCTION__<<endl;exit;}
		tree->SetBranchAddress("runid",&runid);
		tree->SetBranchAddress("eventid",&eventid);
		tree->SetBranchAddress("flagMatch2Lead",&flagMatch2Lead);
		tree->SetBranchAddress("flagMatch2Sub",&flagMatch2Sub);
		tree->SetBranchAddress("flagIsTrigger",&flagIsTrigger);
		tree->SetBranchAddress("trigmatch",&flagtrigmatch);

		tree->SetBranchAddress("Mcj1pt",&McJet);
		tree->SetBranchAddress("Mcj1neutralfrac",&Mcjneutralfrac);
		tree->SetBranchAddress("McLeadAreaNtrk",&McLeadAreaNtrk);
		tree->SetBranchAddress("McSubAreaNtrk",&McSubAreaNtrk);
		tree->SetBranchAddress("McTranMaxNtrk",&McTranMaxNtrk);
		tree->SetBranchAddress("McTranMinNtrk",&McTranMinNtrk);
		tree->SetBranchAddress("McLeadAreaPtSum",&McLeadAreaPtSum);
		tree->SetBranchAddress("McSubLeadAreaPtSum",&McSubAreaPtSum);
		tree->SetBranchAddress("McTranMaxPtSum",&McTranMaxPtSum);
		tree->SetBranchAddress("McTranMinPtSum",&McTranMinPtSum);
		//tree->SetBranchAddress("McPart",&McPart);
		
		tree->SetBranchAddress("Rcj1pt",&RcJet);
		tree->SetBranchAddress("Rcj1neutralfrac",&Rcjneutralfrac);
		tree->SetBranchAddress("RcLeadAreaNtrk",&RcLeadAreaNtrk);
		tree->SetBranchAddress("RcSubAreaNtrk",&RcSubAreaNtrk);
		tree->SetBranchAddress("RcTranMaxNtrk",&RcTranMaxNtrk);
		tree->SetBranchAddress("RcTranMinNtrk",&RcTranMinNtrk);
		tree->SetBranchAddress("RcLeadAreaPtSum",&RcLeadAreaPtSum);
		tree->SetBranchAddress("RcSubLeadAreaPtSum",&RcSubAreaPtSum);
		tree->SetBranchAddress("RcTranMaxPtSum",&RcTranMaxPtSum);
		tree->SetBranchAddress("RcTranMinPtSum",&RcTranMinPtSum);
		//tree->SetBranchAddress("RcPart",&RcPart);

		if(verbose && Nevents[i]>tree->GetEntries()) {cout<<__PRETTY_FUNCTION__<<": pT bin "<<PTBINS[i]<<" Read # of events["<<i<<"] = "<<tree->GetEntries()<<endl; }

		// Weight per pT bin
		weight = XSEC[i]/Nevents[i];
		if(Nevents[i]>tree->GetEntries()) {		// use all events, then NUMBEROFEVENT from include/CrossSectionPerpT.h
			weight = XSEC[i]/NUMBEROFEVENT[i];
			if(verbose) { cout<<__PRETTY_FUNCTION__<<": pT bin "<<PTBINS[i]<<" Use "<<NUMBEROFEVENT[i]<<" for NUMBEROFEVENT in weight = XSEC["<<i<<"]/NUMBEROFEVENT["<<i<<"]"<<endl;}
		}

		for (Int_t j= 0; j<Nevents[i] && j<tree->GetEntries(); j++) {		// loop over entries
			tree->GetEntry(j);

			flag = -999;
			//////JPs
			if(TrigName.Contains("JP",TString::kIgnoreCase)) {
				if(flagIsTrigger && flagtrigmatch && (flagMatch2Lead||flagMatch2Sub) && RcJet>0 && Rcjneutralfrac<0.9) flag = 1; 		
				else if(McJet>0 && (RcJet<=0 || !flagtrigmatch) ) flag = 0;		// non-zero Mc, zero Rc
				else if((RcJet>0&&flagtrigmatch&&flagIsTrigger&&Rcjneutralfrac<0.9) && McJet<=0) flag = -1;		// non-zero Rc, zero Mc		additional neutral fraction here. It is also included in flagMatch2LeadGood and flagMatch2SubGood
			}
			else {
			//if(!flagIsTrigger) continue;			// no correction for trig now..		
			////MB only need the following 3 lines
			if(flagIsTrigger && (flagMatch2Lead||flagMatch2Sub) && RcJet>0 && Rcjneutralfrac<0.9)     flag = 1; 				//2016.11.15 add naive simulated trigger for VPDMB
			else if(McJet>0 && RcJet<=0) flag = 0;		// non-zero Mc, zero Rc
			else if((flagIsTrigger && RcJet>0 && Rcjneutralfrac<0.9) && McJet<=0) flag = -1;		// non-zero Rc, zero Mc			//2016.11.15 add naive simulated trigger for VPDMB
			//
			//if((flagMatch2Lead||flagMatch2Sub) && RcJet>0 && Rcjneutralfrac<0.9)      flag = 1; 				
			//else if(McJet>0 && RcJet<=0) flag = 0;		// non-zero Mc, zero Rc
			//else if(RcJet>0 && Rcjneutralfrac<0.9 && McJet<=0) flag = -1;		// non-zero Rc, zero Mc			
			}

			if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
			{
				if(McTranMaxNtrk>McTranMinNtrk) McPart = (float)McTranMaxNtrk;
				else McPart = (float)McTranMinNtrk;
				if(RcTranMaxNtrk>RcTranMinNtrk) RcPart = (float)RcTranMaxNtrk;
				else RcPart = (float)RcTranMinNtrk;
			}
			else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
			{
				if(McTranMaxNtrk<McTranMinNtrk) McPart = (float)McTranMaxNtrk;
				else McPart = (float)McTranMinNtrk;
				if(RcTranMaxNtrk>RcTranMinNtrk) RcPart = (float)RcTranMaxNtrk;
				else RcPart = (float)RcTranMinNtrk;
			}
			else if(YvariableName.Contains("TranNtrk",TString::kIgnoreCase))
			{
				McPart = 0.5*(McTranMaxNtrk+McTranMinNtrk);
				RcPart = 0.5*(RcTranMaxNtrk+RcTranMinNtrk);
			}	
			else if(YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase))
			{
				McPart = McTranMaxNtrk+McTranMinNtrk;
				RcPart = RcTranMaxNtrk+RcTranMinNtrk;
			}	
			else if(YvariableName.Contains("TranMaxPtSum",TString::kIgnoreCase))
			{
				if(McTranMaxPtSum>McTranMinPtSum) McPart = McTranMaxPtSum;
				else McPart = McTranMinPtSum;
				if(RcTranMaxPtSum>RcTranMinPtSum) RcPart = RcTranMaxPtSum;
				else RcPart = RcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranMinPtSum",TString::kIgnoreCase)) 
			{
				if(McTranMaxPtSum<McTranMinPtSum) McPart = McTranMaxPtSum;
				else McPart = McTranMinPtSum;
				if(RcTranMaxPtSum>RcTranMinPtSum) RcPart = RcTranMaxPtSum;
				else RcPart = RcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranPtSum",TString::kIgnoreCase))
			{
				McPart = 0.5*(McTranMaxPtSum+McTranMinPtSum);
				RcPart = 0.5*(RcTranMaxPtSum+RcTranMinPtSum);
			}
			else if(YvariableName.Contains("TranTotPtSum",TString::kIgnoreCase))
			{
				McPart = McTranMaxPtSum+McTranMinPtSum;
				RcPart = RcTranMaxPtSum+RcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranMaxPtAve",TString::kIgnoreCase))
			{
				float McTranMaxPtAve = 1.*McTranMaxPtSum/McTranMaxNtrk;
				float McTranMinPtAve = 1.*McTranMinPtSum/McTranMinNtrk;
				if(McTranMaxPtAve>McTranMinPtAve) McPart = McTranMaxPtAve;
				else McPart = McTranMinPtAve;
				float RcTranMaxPtAve = 1.*RcTranMaxPtSum/RcTranMaxNtrk;
				float RcTranMinPtAve = 1.*RcTranMinPtSum/RcTranMinNtrk;
				if(RcTranMaxPtAve>RcTranMinPtAve) RcPart = RcTranMaxPtAve;
				else RcPart = RcTranMinPtAve;
			}
			else if(YvariableName.Contains("TranMinPtAve",TString::kIgnoreCase)) 
			{
				float McTranMaxPtAve = 1.*McTranMaxPtSum/McTranMaxNtrk;
				float McTranMinPtAve = 1.*McTranMinPtSum/McTranMinNtrk;
				if(McTranMaxPtAve<McTranMinPtAve) McPart = McTranMaxPtAve;
				else McPart = McTranMinPtAve;
				float RcTranMaxPtAve = 1.*RcTranMaxPtSum/RcTranMaxNtrk;
				float RcTranMinPtAve = 1.*RcTranMinPtSum/RcTranMinNtrk;
				if(RcTranMaxPtAve>RcTranMinPtAve) RcPart = RcTranMaxPtAve;
				else RcPart = RcTranMinPtAve;
			}
			else if(YvariableName.Contains("TranPtAve",TString::kIgnoreCase))
			{
				McPart = (McTranMaxPtSum+McTranMinPtSum)/(McTranMaxNtrk+McTranMinNtrk);
				RcPart = (RcTranMaxPtSum+RcTranMinPtSum)/(RcTranMaxNtrk+RcTranMinNtrk);
			}
			else if(YvariableName.Contains("LeadAreaNtrk",TString::kIgnoreCase)) 
			{
				McPart = McLeadAreaNtrk;
				RcPart = RcLeadAreaNtrk;
			}
			else if(YvariableName.Contains("SubAreaNtrk",TString::kIgnoreCase)) 
			{
				McPart = McSubAreaNtrk;
				RcPart = RcSubAreaNtrk;
			}
			else if(YvariableName.Contains("LeadAreaPtSum",TString::kIgnoreCase)) 
			{
				McPart = McLeadAreaPtSum;
				RcPart = RcLeadAreaPtSum;
			}
			else if(YvariableName.Contains("SubAreaPtSum",TString::kIgnoreCase)) 
			{
				McPart = McSubAreaPtSum;
				RcPart = RcSubAreaPtSum;
			}
			else if(YvariableName.Contains("LeadAreaPtAve",TString::kIgnoreCase)) 
			{
				McPart = McLeadAreaPtSum/McLeadAreaNtrk;
				RcPart = RcLeadAreaPtSum/RcLeadAreaNtrk;
			}
			else if(YvariableName.Contains("SubAreaPtAve",TString::kIgnoreCase)) 
			{
				McPart = McSubAreaPtSum/McSubAreaNtrk;
				RcPart = RcSubAreaPtSum/RcSubAreaNtrk;
			}
			else if(YvariableName.Contains("NeutralFrac",TString::kIgnoreCase)) 
			{
				McPart = Mcjneutralfrac;
				RcPart = Rcjneutralfrac;
			}
			else // default: TranNtrk
			{
				McPart = 0.5*(McTranMaxNtrk+McTranMinNtrk);
				RcPart = 0.5*(RcTranMaxNtrk+RcTranMinNtrk);
			}

			if(flag==1) {			// one-to-one Mc to Rc matched
				hTrainTrue->Fill(McJet, McPart, weight);	
				hTrain->Fill(RcJet, RcPart, weight);	
				response->Fill(RcJet, RcPart,McJet, McPart, weight);
				pfxTrain->Fill(RcJet, RcPart, weight);	
			}
			else if(flag==0) {		// true Mc, no Rc
				hTrainTrue->Fill(McJet, McPart, weight);
				response->Miss(McJet,McPart, weight);
			}
			else if(flag==-1){				// Fake Rc, no Mc
				hTrain->Fill(RcJet,RcPart, weight);
				hTrainFake->Fill(RcJet,RcPart, weight);
				response->Fake(RcJet,RcPart, weight);
				pfxTrain->Fill(RcJet,RcPart, weight);
			}
		}
		ftrain->Close();
	}

	return 1;
}


//==============================================================================
Int_t Unfold2D::Fill4Train ()
{
	//ftrain = new TFile("EmbedTree.root");		// MC embedding files with response matrix and so on
	ftrain = new TFile(inputname);		// MC embedding files with response matrix and so on
	if(!ftrain->IsOpen()) {
		std::cout<<"cannot open "<<inputname<<std::endl;
		return 0;
	}

	tree = (TTree*)ftrain->Get("Tree");
	if(!tree) {cout<<"Cannot find tree from input in "<<__PRETTY_FUNCTION__<<endl;exit;}
	tree->SetBranchAddress("runid",&runid);
	tree->SetBranchAddress("eventid",&eventid);
	tree->SetBranchAddress("flag",&flag);
	tree->SetBranchAddress("McJet",&McJet);
	tree->SetBranchAddress("McPart",&McPart);
	tree->SetBranchAddress("RcJet",&RcJet);
	tree->SetBranchAddress("RcPart",&RcPart);
	//tree->SetBranchAddress("weight",&weight);

	for (Int_t i= 0; i<tree->GetEntries(); i++) {		// loop over entries
		tree->GetEntry(i);
		if(flag==1) {			// one-to-one Mc to Rc matched
			hTrainTrue->Fill(McJet, McPart, weight);	
			hTrain->Fill(RcJet, RcPart, weight);	
			response->Fill(RcJet, RcPart,McJet, McPart, weight);
		}
		else if(flag==0) {		// true Mc, no Rc
			hTrainTrue->Fill(McJet, McPart, weight);
			response->Miss(McJet,McPart, weight);
		}
		else {				// Fake Rc, no Mc
			hTrain->Fill(RcJet,RcPart, weight);
			hTrainFake->Fill(RcJet,RcPart, weight);
			response->Fake(RcJet,RcPart, weight);
		}
	}

	return 1;
}

//==============================================================================
// Train: create response matrix
//==============================================================================

Int_t Unfold2D::Train()
{
	if(WIDEBIN) hTrainTrue = InitWideX2DHisto("htraintrue", "Training Truth", nty, ylo, yhi);
	else hTrainTrue= new TH2D ("htraintrue", "Training Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrainTrue->SetLineColor(kBlue);
	hTrainTrue->GetXaxis()->SetTitle(XvariableName);
	hTrainTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrain= InitWideX2DHisto("htrain", "Training Measured", nmy, ylo, yhi);
	else hTrain= new TH2D ("htrain", "Training Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrain->SetLineColor(kRed);
	hTrain->GetXaxis()->SetTitle(XvariableName);
	hTrain->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrainFake= InitWideX2DHisto("htrainfake", "Training Fakes", nmy, ylo, yhi);
	else hTrainFake= new TH2D ("htrainfake", "Training Fakes", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrainFake->SetLineColor(93);
	hTrainFake->GetXaxis()->SetTitle(XvariableName);
	hTrainFake->GetYaxis()->SetTitle(YvariableName);


	response = new RooUnfoldResponse ("response", "");

	response->Setup (hTrain, hTrainTrue);

	FillbyXsec4Train();
	//Fill4Train();

	hTrainFakeX= ProjectionX (hTrainFake, "hTrainFakeX", "Training Fakes X");
	hTrainFakeY= ProjectionY (hTrainFake, "hTrainFakeY", "Training Fakes Y");

	hTrainTrueX= ProjectionX (hTrainTrue, "hTrainTrueX", "Training X");
	hTrainTrueY= ProjectionY (hTrainTrue, "hTrainTrueY", "Training Y");
	hTrainX=     ProjectionX (hTrain,     "hTrainX",     "Training Measured X");
	hTrainY=     ProjectionY (hTrain,     "hTrainY",     "Training Measured Y");

	return 1;
}


Int_t Unfold2D::ReadResponseMatrix() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	TString ifilename = Form("ResponseMatrix%s_%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data());	
	cout<<__PRETTY_FUNCTION__<<" Read in "<<ifilename<<endl;
	ftrain = new TFile(ifilename);
	hTrainTrue = (TH2D*)ftrain->Get("htraintrue");
	hTrain = (TH2D*)ftrain->Get("htrain");
	hTrainFake = (TH2D*)ftrain->Get("htrainfake");
	response = (RooUnfoldResponse*)ftrain->Get("response");
	
	return 1;
}

void Unfold2D::TrainResults()
{
	setmax (hTrainTrueX, hTrainX, hTrainFakeX);
	setmax (hTrainTrueY, hTrainY, hTrainFakeY);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TLegend *lTrain;
	TCanvas *ctrainx = new TCanvas();
	hTrainTrueX->Draw();
	hTrainX    ->Draw("SAME");
	if (hTrainFakeX) hTrainFakeX->Draw("SAME");
	Legend (lTrain, hTrainTrueX, hTrainFakeX, hTrainX);
	ctrainx->Update();

	TCanvas *ctrainy = new TCanvas();
	hTrainTrueY->Draw();
	hTrainY    ->Draw("SAME");
	if (hTrainFakeY) hTrainFakeY->Draw("SAME");
	lTrain->Draw();
	ctrainy->Update();

}


//==============================================================================
// Read and fill histogram for unfolding
//==============================================================================
Int_t Unfold2D::Fill4Unfold() {

	if(WIDEBIN) hMeas= InitWideX2DHisto("hmeas", "Measured", nmy, ylo, yhi);
	else hMeas= new TH2D ("hmeas", "Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hMeas->SetLineColor(kRed);
	hMeas->GetXaxis()->SetTitle(XvariableName);
	hMeas->GetYaxis()->SetTitle(YvariableName);

	Float_t InputJet;
	Float_t InputPart;	

	Int_t  InputLeadAreaNtrk;
	Int_t  InputSubAreaNtrk;
	Int_t  InputTranMaxNtrk;
	Int_t  InputTranMinNtrk;
	Float_t InputLeadAreaPtSum;
	Float_t InputSubAreaPtSum;
	Float_t InputTranMaxPtSum;
	Float_t InputTranMinPtSum;

	TString ifilename;
	if(TrigName.Contains("MB")) {
		ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_160829.root");
	}
	else {
		ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_160829.root");
	}
	cout<<__PRETTY_FUNCTION__<<__func__<<ifilename<<endl;
	TFile *fin = new TFile(ifilename);
	if(!fin) {cout<<"Cannot open file in "<<__PRETTY_FUNCTION__<<endl; exit;}
	TTree *ttree = (TTree*)fin->Get("ResultTree");
	if(!ttree) {cout<<"Cannot find tree from input in "<<__PRETTY_FUNCTION__<<endl;exit;}

	ttree->SetBranchAddress("j1pt",&InputJet);
	ttree->SetBranchAddress("LeadAreaNtrk",&InputLeadAreaNtrk);
	ttree->SetBranchAddress("SubAreaNtrk",&InputSubAreaNtrk);
	ttree->SetBranchAddress("TranMaxNtrk",&InputTranMaxNtrk);
	ttree->SetBranchAddress("TranMinNtrk",&InputTranMinNtrk);
	ttree->SetBranchAddress("LeadAreaPtSum",&InputLeadAreaPtSum);
	ttree->SetBranchAddress("SubLeadAreaPtSum",&InputSubAreaPtSum);
	ttree->SetBranchAddress("TranMaxPtSum",&InputTranMaxPtSum);
	ttree->SetBranchAddress("TranMinPtSum",&InputTranMinPtSum);

	for (Int_t j=0;j<ttree->GetEntries() ; j++) {			// Caution!!! Need to update loop if needed
		ttree->GetEntry(j);

		if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
		{
			if(InputTranMaxNtrk>InputTranMinNtrk) InputPart = (float)InputTranMaxNtrk;
			else InputPart = (float)InputTranMinNtrk;
		}
		else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
		{
			if(InputTranMaxNtrk<InputTranMinNtrk) InputPart = (float)InputTranMaxNtrk;
			else InputPart = (float)InputTranMinNtrk;
		}
		else if(YvariableName.Contains("TranNtrk",TString::kIgnoreCase))
		{
			InputPart = 0.5*(InputTranMaxNtrk+InputTranMinNtrk);
		}	
		else if(YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase))
		{
			InputPart = InputTranMaxNtrk+InputTranMinNtrk;
		}	
		else if(YvariableName.Contains("TranMaxPtSum",TString::kIgnoreCase))
		{
			if(InputTranMaxPtSum>InputTranMinPtSum) InputPart = InputTranMaxPtSum;
			else InputPart = InputTranMinPtSum;
		}
		else if(YvariableName.Contains("TranMinPtSum",TString::kIgnoreCase)) 
		{
			if(InputTranMaxPtSum<InputTranMinPtSum) InputPart = InputTranMaxPtSum;
			else InputPart = InputTranMinPtSum;
		}
		else if(YvariableName.Contains("TranPtSum",TString::kIgnoreCase))
		{
			InputPart = 0.5*(InputTranMaxPtSum+InputTranMinPtSum);
		}
		else if(YvariableName.Contains("TranTotPtSum",TString::kIgnoreCase))
		{
			InputPart = InputTranMaxPtSum+InputTranMinPtSum;
		}
		else if(YvariableName.Contains("TranMaxPtAve",TString::kIgnoreCase))
		{
			float InputTranMaxPtAve = 1.*InputTranMaxPtSum/InputTranMaxNtrk;
			float InputTranMinPtAve = 1.*InputTranMinPtSum/InputTranMinNtrk;
			if(InputTranMaxPtAve>InputTranMinPtAve) InputPart = InputTranMaxPtAve;
			else InputPart = InputTranMinPtAve;
		}
		else if(YvariableName.Contains("TranMinPtAve",TString::kIgnoreCase)) 
		{
			float InputTranMaxPtAve = 1.*InputTranMaxPtSum/InputTranMaxNtrk;
			float InputTranMinPtAve = 1.*InputTranMinPtSum/InputTranMinNtrk;
			if(InputTranMaxPtAve<InputTranMinPtAve) InputPart = InputTranMaxPtAve;
			else InputPart = InputTranMinPtAve;
		}
		else if(YvariableName.Contains("TranPtAve",TString::kIgnoreCase))
		{
			InputPart = (InputTranMaxPtSum+InputTranMinPtSum)/(InputTranMaxNtrk+InputTranMinNtrk);
		}
		else if(YvariableName.Contains("LeadAreaNtrk",TString::kIgnoreCase)) 
		{
			InputPart = InputLeadAreaNtrk;
		}
		else if(YvariableName.Contains("SubAreaNtrk",TString::kIgnoreCase)) 
		{
			InputPart = InputSubAreaNtrk;
		}
		else if(YvariableName.Contains("LeadAreaPtSum",TString::kIgnoreCase)) 
		{
			InputPart = InputLeadAreaPtSum;
		}
		else if(YvariableName.Contains("SubAreaPtSum",TString::kIgnoreCase)) 
		{
			InputPart = InputSubAreaPtSum;
		}
		else if(YvariableName.Contains("LeadAreaPtAve",TString::kIgnoreCase)) 
		{
			InputPart = InputLeadAreaPtSum/InputLeadAreaNtrk;
		}
		else if(YvariableName.Contains("SubAreaPtAve",TString::kIgnoreCase)) 
		{
			InputPart = InputSubAreaPtSum/InputSubAreaNtrk;
		}
		else // default: TranNtrk
		{
			InputPart = 0.5*(InputTranMaxNtrk+InputTranMinNtrk);
		}

		// Fill histogram
		hMeas->Fill(InputJet, InputPart);	

	}// End Loop of events
	fin->Close();

	hMeas->Sumw2();

	hMeasX= ProjectionX (hMeas, "hMeasX", "Measured X");
	hMeasY= ProjectionY (hMeas, "hMeasY", "Measured Y");


	return 1;
}


//==============================================================================
// Read hMeas histogram from input root file to unfold
//==============================================================================
Int_t Unfold2D::ReadMeasHist4Unfold() {

	TString ifilename = TString("~/Scratch/pp200Y12_jetunderlying/leadjetpthist4NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_160829_NoEffCorr.root");
	cout<<__PRETTY_FUNCTION__<<" Read in "<<ifilename<<endl;
	TFile *fin = new TFile(ifilename);
	if(!fin) {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<" Cannot open "<< ifilename<<endl; 
		exit; 
	}
	TString ihistname = "h"+ToLower(YvariableName)+"vsleadjetpt";			// Default is Tran Multiplicity 
	TH2D *htmp = (TH2D*)fin->Get(ihistname);
	if(hTrain && htmp) hMeas = RebinAs(htmp,hTrain);
	else {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<": hTrain is zero."<<endl; 
		exit; 
	}
	if(!hMeas) {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<": Cannot get hmeas from "<< ifilename<<endl; 
		exit; 
	}

	hMeasX= ProjectionX (hMeas, "hMeasX", "Test Measured X");
	hMeasY= ProjectionY (hMeas, "hMeasY", "Test Measured Y");

	return 1;

}

//==============================================================================
// Unfold 
//==============================================================================
Int_t Unfold2D::Unfold() {
	if (verbose>=0) cout << "Create RooUnfold object for method " << method << endl;
	unfold= RooUnfold::New ((RooUnfold::Algorithm)method, response, hMeas, regparm, "unfold");
	if (!unfold) return 0;
	unfold->SetVerbose (verbose);
	if (verbose>=0) {cout << "Created "; unfold->Print();}
	hReco= unfold->Hreco((RooUnfold::ErrorTreatment)doerror);
	if (!hReco) return 0;
	hReco->SetName("hreco");
	hReco->SetLineColor(kBlack);  // otherwise inherits style from hTrainTrue
	if (verbose>=0) unfold->PrintTable (cout, hTrue, (RooUnfold::ErrorTreatment)doerror);
	if (verbose>=2 && doerror>=RooUnfold::kCovariance) {
		TMatrixD covmat= unfold->Ereco((RooUnfold::ErrorTreatment)doerror);
		int ntbins = ntx;
		TMatrixD errmat(ntbins,ntbins);
		for (Int_t i=0; i<ntbins; i++) {
			for (Int_t j=0; j<ntbins; j++) {
				errmat(i,j)= covmat(i,j)>=0 ? sqrt(covmat(i,j)) : -sqrt(-covmat(i,j));
			}
			RooUnfoldResponse::PrintMatrix(errmat,"covariance matrix");
		}
	}

	hCorr= CorrelationHist (unfold->Ereco((RooUnfold::ErrorTreatment)doerror),
			"corr", "Unfolded correlation matrix",
			response->Hresponse()->GetYaxis()->GetXmin(),
			response->Hresponse()->GetYaxis()->GetXmax());




	return 1;
}


//==============================================================================
// Test unfolding algorithm
//==============================================================================

Int_t Unfold2D::TrainAndTest ()
{

	if(WIDEBIN) hTrainTrue = InitWideX2DHisto("htraintrue", "Training Truth", nty, ylo, yhi);
	else hTrainTrue= new TH2D ("htraintrue", "Training Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrainTrue->SetLineColor(kBlue);
	hTrainTrue->GetXaxis()->SetTitle(XvariableName);
	hTrainTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrain= InitWideX2DHisto("htrain", "Training Measured", nmy, ylo, yhi);
	else hTrain= new TH2D ("htrain", "Training Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrain->SetLineColor(kRed);
	hTrain->GetXaxis()->SetTitle(XvariableName);
	hTrain->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrainFake= InitWideX2DHisto("htrainfake", "Training Fakes", nmy, ylo, yhi);
	else hTrainFake= new TH2D ("htrainfake", "Training Fakes", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrainFake->SetLineColor(93);
	hTrainFake->GetXaxis()->SetTitle(XvariableName);
	hTrainFake->GetYaxis()->SetTitle(YvariableName);

	cout<<"Set response"<<endl;
	response = new RooUnfoldResponse ("response", "");
	response->Setup (hTrain, hTrainTrue);

	int Nevents[NUMBEROFPT] = {600000,200000,300000,150000,150000,150000,80000,50000,40000,38000,12000};
	//for(int i=0; i<NUMBEROFPT; i++) {
	//	cout<<Nevents[i]<<endl;
	//}
	FillbyXsec4Train(Nevents);

	hTrainFakeX= ProjectionX (hTrainFake, "hTrainFakeX", "Training Fakes X");
	hTrainFakeY= ProjectionY (hTrainFake, "hTrainFakeY", "Training Fakes Y");

	hTrainTrueX= ProjectionX (hTrainTrue, "hTrainTrueX", "Training X");
	hTrainTrueY= ProjectionY (hTrainTrue, "hTrainTrueY", "Training Y");
	hTrainX=     ProjectionX (hTrain,     "hTrainX",     "Training Measured X");
	hTrainY=     ProjectionY (hTrain,     "hTrainY",     "Training Measured Y");

	TrainResults();

	Fill4Test(Nevents);

	Unfold();

	return 1;
}

Int_t Unfold2D::Fill4Test (int *Nevents)
{
	//TFile *fin = new TFile("UnfoldMatrix_part1.root");			// measured result
	cout<<"I am using parts of simulations for train, rest for test"<<endl;
	cout<<"If this is not what you want, change loop in "<< __PRETTY_FUNCTION__<<endl;

	if(WIDEBIN) hTrue= InitWideX2DHisto("htrue", "Test Truth", nty, ylo, yhi);
	else hTrue= new TH2D ("htrue", "Test Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrue->SetLineColor(kBlue);
	hTrue->GetXaxis()->SetTitle(XvariableName);
	hTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hMeas= InitWideX2DHisto("hmeas", "Test Measured", nmy, ylo, yhi);
	else hMeas= new TH2D ("hmeas", "Test Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hMeas->SetLineColor(kRed);
	hMeas->GetXaxis()->SetTitle(XvariableName);
	hMeas->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hFake= InitWideX2DHisto("hfake", "Test Fake", nmy, ylo, yhi);
	else hFake= new TH2D ("hfake", "Test Fake", nmx, xlo, xhi, nmy, ylo, yhi);
	hFake->SetLineColor(93);
	hFake->GetXaxis()->SetTitle(XvariableName);
	hFake->GetYaxis()->SetTitle(YvariableName);

	TProfile *pfxMeas;
	if(WIDEBIN) pfxMeas= new TProfile("pfxmeas", "Test Measured", WNbins, Wptbins);
	else pfxMeas= new TProfile ("pfxmeas", "Test Measured", nmx, xlo, xhi);
	pfxMeas->SetLineColor(kRed);
	pfxMeas->GetXaxis()->SetTitle(XvariableName);
	pfxMeas->GetYaxis()->SetTitle(YvariableName);


	int tflag;
	float tMcJet, tMcPart, tRcJet, tRcPart;
	double tweight;

	Bool_t tflagMatch2Lead;
	Bool_t tflagMatch2Sub;
	Bool_t tflagIsTrigger;
	Bool_t tflagtrigmatch;

	Int_t  tMcLeadAreaNtrk;
	Int_t  tMcSubAreaNtrk;
	Int_t  tMcTranMaxNtrk;
	Int_t  tMcTranMinNtrk;
	Float_t tMcLeadAreaPtSum;
	Float_t tMcSubAreaPtSum;
	Float_t tMcTranMaxPtSum;
	Float_t tMcTranMinPtSum;

	Int_t  tRcTranMaxNtrk;
	Int_t  tRcTranMinNtrk;
	Int_t  tRcLeadAreaNtrk;
	Int_t  tRcSubAreaNtrk;
	Float_t tRcLeadAreaPtSum;
	Float_t tRcSubAreaPtSum;
	Float_t tRcTranMaxPtSum;
	Float_t tRcTranMinPtSum;




	double totalxsec = 0;
	double totalnevents = 0;
	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename;
		if(ExcludeOpt&&TrigName.Contains("JP")) {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_excluded.root",TrigName.Data(),PTBINS[i],TranCharge.Data());
		}
		else {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data());
		}

		cout<<"Read in "<<ifilename<<endl;
		TFile *fin = new TFile(ifilename);
		TTree *ttree = (TTree*)fin->Get("ResultTree");
		if(!ttree) {cout<<"Cannot find tree from input in "<<__PRETTY_FUNCTION__<<endl;exit;}

		ttree->SetBranchAddress("flagMatch2Lead",&tflagMatch2Lead);
		ttree->SetBranchAddress("flagMatch2Sub",&tflagMatch2Sub);
		ttree->SetBranchAddress("flagIsTrigger",&tflagIsTrigger);
		ttree->SetBranchAddress("trigmatch",&tflagtrigmatch);

		ttree->SetBranchAddress("Mcj1pt",&tMcJet);
		ttree->SetBranchAddress("McLeadAreaNtrk",&tMcLeadAreaNtrk);
		ttree->SetBranchAddress("McSubAreaNtrk",&tMcSubAreaNtrk);
		ttree->SetBranchAddress("McTranMaxNtrk",&tMcTranMaxNtrk);
		ttree->SetBranchAddress("McTranMinNtrk",&tMcTranMinNtrk);
		ttree->SetBranchAddress("McLeadAreaPtSum",&tMcLeadAreaPtSum);
		ttree->SetBranchAddress("McSubLeadAreaPtSum",&tMcSubAreaPtSum);
		ttree->SetBranchAddress("McTranMaxPtSum",&tMcTranMaxPtSum);
		ttree->SetBranchAddress("McTranMinPtSum",&tMcTranMinPtSum);

		ttree->SetBranchAddress("Rcj1pt",&tRcJet);
		ttree->SetBranchAddress("RcLeadAreaNtrk",&tRcLeadAreaNtrk);
		ttree->SetBranchAddress("RcSubAreaNtrk",&tRcSubAreaNtrk);
		ttree->SetBranchAddress("RcTranMaxNtrk",&tRcTranMaxNtrk);
		ttree->SetBranchAddress("RcTranMinNtrk",&tRcTranMinNtrk);
		ttree->SetBranchAddress("RcLeadAreaPtSum",&tRcLeadAreaPtSum);
		ttree->SetBranchAddress("RcSubLeadAreaPtSum",&tRcSubAreaPtSum);
		ttree->SetBranchAddress("RcTranMaxPtSum",&tRcTranMaxPtSum);
		ttree->SetBranchAddress("RcTranMinPtSum",&tRcTranMinPtSum);

		if(Nevents[i]>ttree->GetEntries()) {cout<<"Error!! "<<__PRETTY_FUNCTION__<<" \n pT bin "<<PTBINS[i]<<" (Nevents["<<i<<"] = "<<Nevents[i]<<" >ttree->GetEntries() = "<<ttree->GetEntries()<<"\nPlease adjust it and rerun the code."<<endl; exit; }

		// Weight per pT bin
		tweight = XSEC[i]/(ttree->GetEntries()-Nevents[i]);
		totalxsec+=XSEC[i];
		totalnevents+=(ttree->GetEntries()-Nevents[i]);

		cout<<"From "<<Nevents[i]<<" to "<<ttree->GetEntries()<<endl;
		for (Int_t j=Nevents[i];j<ttree->GetEntries() ; j++) {			// Caution!!! Need to update loop if needed
			ttree->GetEntry(j);

			if(!tflagIsTrigger) continue;			// no correction for trig now..
			tflag = -999;
			if(tflagtrigmatch && (tflagMatch2Lead||tflagMatch2Sub) ) tflag = 1; 
			else if(tMcJet>0 && (tRcJet<=0 || !tflagtrigmatch) ) tflag = 0;		// non-zero tMc, zero tRc
			else if((tRcJet>0&&tflagtrigmatch) && tMcJet<=0) tflag = -1;		// non-zero tRc, zero tMc
			//if(tflagMatch2Lead||tflagMatch2Sub) tflag = 1; 
			//else if(tMcJet>=0 && tRcJet<=0) tflag = 0;		// non-zero Mc, zero Rc
			//else if(tRcJet>=0 && tMcJet<=0) tflag = -1;		// non-zero Rc, zero Mc


			// Read tMcPart and tRcPart
			if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
			{
				if(tMcTranMaxNtrk>tMcTranMinNtrk) tMcPart = (float)tMcTranMaxNtrk;
				else tMcPart = (float)tMcTranMinNtrk;
				if(tRcTranMaxNtrk>tRcTranMinNtrk) tRcPart = (float)tRcTranMaxNtrk;
				else tRcPart = (float)tRcTranMinNtrk;
			}
			else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
			{
				if(tMcTranMaxNtrk<tMcTranMinNtrk) tMcPart = (float)tMcTranMaxNtrk;
				else tMcPart = (float)tMcTranMinNtrk;
				if(tRcTranMaxNtrk>tRcTranMinNtrk) tRcPart = (float)tRcTranMaxNtrk;
				else tRcPart = (float)tRcTranMinNtrk;
			}
			else if(YvariableName.Contains("TranNtrk",TString::kIgnoreCase))
			{
				tMcPart = 0.5*(tMcTranMaxNtrk+tMcTranMinNtrk);
				tRcPart = 0.5*(tRcTranMaxNtrk+tRcTranMinNtrk);
			}	
			else if(YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase))
			{
				tMcPart = tMcTranMaxNtrk+tMcTranMinNtrk;
				tRcPart = tRcTranMaxNtrk+tRcTranMinNtrk;
			}	
			else if(YvariableName.Contains("TranMaxPtSum",TString::kIgnoreCase))
			{
				if(tMcTranMaxPtSum>tMcTranMinPtSum) tMcPart = tMcTranMaxPtSum;
				else tMcPart = tMcTranMinPtSum;
				if(tRcTranMaxPtSum>tRcTranMinPtSum) tRcPart = tRcTranMaxPtSum;
				else tRcPart = tRcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranMinPtSum",TString::kIgnoreCase)) 
			{
				if(tMcTranMaxPtSum<tMcTranMinPtSum) tMcPart = tMcTranMaxPtSum;
				else tMcPart = tMcTranMinPtSum;
				if(tRcTranMaxPtSum>tRcTranMinPtSum) tRcPart = tRcTranMaxPtSum;
				else tRcPart = tRcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranPtSum",TString::kIgnoreCase))
			{
				tMcPart = 0.5*(tMcTranMaxPtSum+tMcTranMinPtSum);
				tRcPart = 0.5*(tRcTranMaxPtSum+tRcTranMinPtSum);
			}
			else if(YvariableName.Contains("TranTotPtSum",TString::kIgnoreCase))
			{
				tMcPart = tMcTranMaxPtSum+tMcTranMinPtSum;
				tRcPart = tRcTranMaxPtSum+tRcTranMinPtSum;
			}
			else if(YvariableName.Contains("TranMaxPtAve",TString::kIgnoreCase))
			{
				float tMcTranMaxPtAve = 1.*tMcTranMaxPtSum/tMcTranMaxNtrk;
				float tMcTranMinPtAve = 1.*tMcTranMinPtSum/tMcTranMinNtrk;
				if(tMcTranMaxPtAve>tMcTranMinPtAve) tMcPart = tMcTranMaxPtAve;
				else tMcPart = tMcTranMinPtAve;
				float tRcTranMaxPtAve = 1.*tRcTranMaxPtSum/tRcTranMaxNtrk;
				float tRcTranMinPtAve = 1.*tRcTranMinPtSum/tRcTranMinNtrk;
				if(tRcTranMaxPtAve>tRcTranMinPtAve) tRcPart = tRcTranMaxPtAve;
				else tRcPart = tRcTranMinPtAve;
			}
			else if(YvariableName.Contains("TranMinPtAve",TString::kIgnoreCase)) 
			{
				float tMcTranMaxPtAve = 1.*tMcTranMaxPtSum/tMcTranMaxNtrk;
				float tMcTranMinPtAve = 1.*tMcTranMinPtSum/tMcTranMinNtrk;
				if(tMcTranMaxPtAve<tMcTranMinPtAve) tMcPart = tMcTranMaxPtAve;
				else tMcPart = tMcTranMinPtAve;
				float tRcTranMaxPtAve = 1.*tRcTranMaxPtSum/tRcTranMaxNtrk;
				float tRcTranMinPtAve = 1.*tRcTranMinPtSum/tRcTranMinNtrk;
				if(tRcTranMaxPtAve>tRcTranMinPtAve) tRcPart = tRcTranMaxPtAve;
				else tRcPart = tRcTranMinPtAve;
			}
			else if(YvariableName.Contains("TranPtAve",TString::kIgnoreCase))
			{
				tMcPart = (tMcTranMaxPtSum+tMcTranMinPtSum)/(tMcTranMaxNtrk+tMcTranMinNtrk);
				tRcPart = (tRcTranMaxPtSum+tRcTranMinPtSum)/(tRcTranMaxNtrk+tRcTranMinNtrk);
			}
			else if(YvariableName.Contains("LeadAreaNtrk",TString::kIgnoreCase)) 
			{
				tMcPart = tMcLeadAreaNtrk;
				tRcPart = tRcLeadAreaNtrk;
			}
			else if(YvariableName.Contains("SubAreaNtrk",TString::kIgnoreCase)) 
			{
				tMcPart = tMcSubAreaNtrk;
				tRcPart = tRcSubAreaNtrk;
			}
			else if(YvariableName.Contains("LeadAreaPtSum",TString::kIgnoreCase)) 
			{
				tMcPart = tMcLeadAreaPtSum;
				tRcPart = tRcLeadAreaPtSum;
			}
			else if(YvariableName.Contains("SubAreaPtSum",TString::kIgnoreCase)) 
			{
				tMcPart = tMcSubAreaPtSum;
				tRcPart = tRcSubAreaPtSum;
			}
			else if(YvariableName.Contains("LeadAreaPtAve",TString::kIgnoreCase)) 
			{
				tMcPart = tMcLeadAreaPtSum/tMcLeadAreaNtrk;
				tRcPart = tRcLeadAreaPtSum/tRcLeadAreaNtrk;
			}
			else if(YvariableName.Contains("SubAreaPtAve",TString::kIgnoreCase)) 
			{
				tMcPart = tMcSubAreaPtSum/tMcSubAreaNtrk;
				tRcPart = tRcSubAreaPtSum/tRcSubAreaNtrk;
			}
			else // default: TranNtrk
			{
				tMcPart = 0.5*(tMcTranMaxNtrk+tMcTranMinNtrk);
				tRcPart = 0.5*(tRcTranMaxNtrk+tRcTranMinNtrk);
			}

			// Fill histograms according to flag
			if(tflag==1) {			// one-to-one tMc to tRc matched
				hTrue->Fill(tMcJet, tMcPart, tweight);	
				hMeas->Fill(tRcJet, tRcPart, tweight);	
				pfxMeas->Fill(tRcJet, tRcPart, tweight);	
			}
			else if(tflag==0) {		// true tMc, no tRc
				hTrue->Fill(tMcJet, tMcPart, tweight);
			}
			else if(tflag==-1){				// Fake tRc, no tMc
				hMeas->Fill(tRcJet,tRcPart, tweight);
				hFake->Fill(tRcJet,tRcPart, tweight);
				pfxMeas->Fill(tRcJet,tRcPart, tweight);
			}

		}// End Loop of events
		fin->Close();
	}// End Loop of pT

	//Scale the histogram to simulate statistic of how many events. Note: sumw2() is called after wards, otherwise it would serve the purpose of simulate the statistic . 
	//double toscale = totalnevents/totalxsec;
	double toscale = 100e6/totalxsec;		// say if we have 100M MB data. how does the real unfolding probably looks like
	hMeas->Scale(toscale);
	hTrue->Scale(toscale);
	hFake->Scale(toscale);

	hMeas->Sumw2();
	hTrue->Sumw2();
	hFake->Sumw2();

	hFakeX= ProjectionX (hFake, "hFakeX", "Test Fakes X");
	hFakeY= ProjectionY (hFake, "hFakeY", "Test Fakes Y");

	hTrueX= ProjectionX (hTrue, "hTrueX", "Test X");
	hTrueY= ProjectionY (hTrue, "hTrueY", "Test Y");

	hMeasX= ProjectionX (hMeas, "hMeasX", "Test Measured X");
	hMeasY= ProjectionY (hMeas, "hMeasY", "Test Measured Y");


	return 1;
}

//==============================================================================
// Show results
//==============================================================================

void Unfold2D::Results()
{

	double xmin = 11;
	double xmax = 60;

	if (hReco) {
		hRecoX= ProjectionX (hReco, "hRecoX", "Reconstructed X", "E");
		hRecoY= ProjectionY (hReco, "hRecoY", "Reconstructed Y", "E");

		hRecoX->SetMarkerStyle(kFullDotLarge);
		hRecoY->SetMarkerStyle(kFullDotLarge);
		setmax (hMeasX, hRecoX);
		setmax (hMeasY, hRecoY);

		if(hFakeX&&hTrueX) setmax(hMeasX,hFakeX,hTrueX,hRecoX);
		if(hFakeY&&hTrueY) setmax(hMeasY,hFakeY,hTrueY,hRecoY);
	}

	TLegend *lTest;
	TLine *line = new TLine();
	TCanvas *ctestx = new TCanvas();
	setmax(hMeasX,hTrueX,hRecoX);
	if(1) {
		ctestx->SetLogy();
		hMeasX->SetMinimum(1e-3);
		hMeasX->SetMaximum(1e8);
	}
	hMeasX->GetXaxis()->SetRangeUser(xmin,xmax);
	hMeasX   ->Draw();
	if (hTrueX) hTrueX->Draw("SAME");
	//if (hFakeX) hFakeX->Draw("SAME");
	if (hRecoX) hRecoX->Draw("SAME");
	//Legend (lTest, hTrueX, hFakeX, hMeasX, hRecoX);
	Legend (lTest, hTrueX, NULL, hMeasX, hRecoX);
	ctestx->SaveAs(Form("fig/ctestx.png"));
	ctestx->Update();

	if(hRecoX && hTrueX) {
		TCanvas *ctestxratio = new TCanvas();
		TH1D *hXRecoRTrue = (TH1D*)hRecoX->Clone("hXRecoRTrue");
		hXRecoRTrue->Divide(hTrueX);
		hXRecoRTrue->GetYaxis()->SetTitle(Form("%s Unfolded / Mc particle-level",XvariableName.Data()));
		hXRecoRTrue->GetYaxis()->SetRangeUser(0.7,1.4);
		hXRecoRTrue->GetXaxis()->SetRangeUser(xmin,xmax);
		hXRecoRTrue->Draw();
		//line->DrawLine(hXRecoRTrue->GetXaxis()->GetXmin(),1,hXRecoRTrue->GetXaxis()->GetXmax(),1);
		line->DrawLine(xmin,1,xmax+hXRecoRTrue->GetBinWidth(hXRecoRTrue->FindBin(xmax))/2.,1);
		ctestxratio->SaveAs(Form("fig/jetunfold_ratio.png"));
	}

	TCanvas *ctesty = new TCanvas();
	hMeasY   ->Draw();
	if (hTrueY) hTrueY->Draw("SAME");
	//if (hFakeY) hFakeY->Draw("SAME");
	if (hRecoY) hRecoY->Draw("SAME");
	lTest->Draw();
	ctesty->SaveAs(Form("fig/ctesty.png"));
	ctesty->Update();

	if(hTrue) {
		TCanvas *ctesttrue = new TCanvas();
		hTrue->Draw("colz");
		ctesttrue->SaveAs(Form("fig/ctesttrue.png"));
	}

	TCanvas *ctestmeas = new TCanvas();
	hMeas->Draw("colz");
	ctestmeas->SaveAs(Form("fig/ctestmeas.png"));

	if (!hReco) return;

	TCanvas *ctestreco = new TCanvas();
	hReco->Draw("colz");
	ctestreco->SaveAs(Form("fig/ctestreco.png"));

	//if(hCorr) {
	//	gStyle->SetPalette(1,0);
	//	TCanvas *ccorr = new TCanvas();
	//	hCorr->Draw("COLZ");
	//	ccorr->SaveAs(Form("fig/ccorr.png"));
	//}


	TCanvas *cprofile = new TCanvas();
	TProfile *hMeas_pfx = (TProfile*)hMeas->ProfileX("hmeas_pfx");
	hMeas_pfx->GetYaxis()->SetTitle(YvariableName);
	hMeas_pfx->SetMinimum(0);	
	hMeas_pfx->GetXaxis()->SetRangeUser(xmin,xmax);
	hMeas_pfx->Draw();
	TProfile *True_pfx;
	//if(True_pfx) cout<<"Getname1 = "<<True_pfx->GetName()<<endl;

	if(hTrue) {
		True_pfx = ProfileX(hTrue,"true_pfx","Truth Y vs X");
		True_pfx->Draw("same");
	}
	//if(hFake) hFake->ProfileX()->Draw("same");
	TProfile *Reco_pfx;
	if(hReco) {
		Reco_pfx = ProfileX(hReco,"reco_pfx","Reconstruction Y vs X");
		Reco_pfx->SetMarkerStyle(8);
		Reco_pfx->GetYaxis()->SetTitle(YvariableName);
		Reco_pfx->Draw("same");
	}
	if(hTrue && hReco) setmax(hMeas_pfx,True_pfx,Reco_pfx);
	if(hMeas_pfx->GetMaximum()>1000) hMeas_pfx->SetMaximum(6);
	lTest->Draw();
	cprofile->SaveAs(Form("fig/cprofile.png"));


	if(hTrue && hReco) {
		TCanvas *cprofileratio = new TCanvas();
		TH1D *hTrue_pfx= TProfile2TH1D(True_pfx);
		TH1D *hReco_pfx= TProfile2TH1D(Reco_pfx);
		TH1D *hpfxRecoRTrue = (TH1D*)hReco_pfx->Clone("hpfxRecoRTrue");
		hpfxRecoRTrue->Divide(hTrue_pfx);
		hpfxRecoRTrue->GetYaxis()->SetTitle(Form("%s Unfolded / Mc particle-level",YvariableName.Data()));
		hpfxRecoRTrue->GetXaxis()->SetTitle(XvariableName);
		hpfxRecoRTrue->GetYaxis()->SetRangeUser(0.7,1.4);
		hpfxRecoRTrue->GetXaxis()->SetRangeUser(xmin,xmax);
		hpfxRecoRTrue->SetMarkerStyle(8);
		hpfxRecoRTrue->Draw();
		//line->DrawLine(hpfxRecoRTrue->GetXaxis()->GetXmin(),1,hpfxRecoRTrue->GetXaxis()->GetXmax(),1);
		line->DrawLine(xmin,1,xmax+hpfxRecoRTrue->GetBinWidth(hpfxRecoRTrue->FindBin(xmax))/2.,1);
		cprofileratio->SaveAs(Form("fig/%sVsjpt_ratio.png",XvariableName.Data()));
	}

}

//==============================================================================
// Write Test Results
//==============================================================================

Int_t Unfold2D::WriteTest() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	ftout = new TFile(Form("ResponseMatrix%s_%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data()),"RECREATE");
	if(ftout) cout<<"Write to "<<ftout->GetName()<<endl;
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	response->Write();

	if (hTrue) hTrue->Write();
	if (hMeas) hMeas->Write();
	if (hFake) hFake->Write();
	if (hReco) hReco->Write();


	return 1;
}


//==============================================================================
// Write: record response matrix from training for future unfolding
//==============================================================================

Int_t Unfold2D::WriteTrain() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	cout<<"Write to "<<Form("ResponseMatrix%s_%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data())<<endl;
	ftout = new TFile(Form("ResponseMatrix%s_%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data()),"RECREATE");
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	response->Write();
	if(pfxTrain) pfxTrain->Write();

	cout<<"WriteTrain() done"<<endl;

	return 1;
}


//==============================================================================
// Write: record unfolding result
//==============================================================================

Int_t Unfold2D::WriteUnfoldResult() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	ftout = new TFile(Form("Unfolding%s_%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data()),"RECREATE");
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	response->Write();
	hMeas->Write();
	if(hReco) hReco->Write();

	return 1;
}

//==============================================================================
// Constructors and destructor
//==============================================================================

Unfold2D::Unfold2D (TString name): inputname(name), XvariableName("X"), YvariableName("Y"), TrigName("JP2"), TranCharge("Charged"), ExcludeOpt(0)
{
	Reset();
	SetDefaultParms();
}

	Unfold2D::Unfold2D (const char* name, int argc, const char* const* argv)
: inputname(name), XvariableName("X"), YvariableName("Y"), TrigName("JP2"), TranCharge("Charged"), ExcludeOpt(0)
{
	Reset();
	if(argc<14){ Help(); return; }
	SetParms(argv);
}


Unfold2D::~Unfold2D()
{
	delete response; response= 0;
	delete unfold;   unfold=   0;
	if(ftrain->IsOpen()) ftrain->Close();
	if(ftout->IsOpen()) ftout->Close();
}


//==============================================================================
// Utility routines
//==============================================================================

void Unfold2D::Reset()
{
	response= 0;
	unfold= 0;
	hTrain= hTrainTrue= hTrainFake= hTrue= hMeas= hFake=0 ;
	hReco= 0;
	hCorr= 0;

	hTrainX= hTrainTrueX= hTrueX= hTrainFakeX= hFakeX= hMeasX= hRecoX= 
		hTrainY= hTrainTrueY= hTrueY= hTrainFakeY= hFakeY= hMeasY= hRecoY= 0;
}

void Unfold2D::Init()
{
	//*nothing for now
}

Int_t Unfold2D::CheckParms()
{
	int error = 0;
	if (verbose>=0) PrintParms ();
	if (ntx<=0)     {cerr << "Error: ntx ("    << ntx    << ") <= 0"                  << endl; error = 2;}
	if (nmx<=0)     {cerr << "Error: nmx ("    << nmx    << ") <= 0"                  << endl; error = 2;}
	if (xlo >= xhi) {cerr << "Error: xlo ("    << xlo    << ") >= xhi(" << xhi << ")" << endl; error = 2;}
	if (nty<=0)     {cerr << "Error: nty ("    << nty    << ") <= 0"                  << endl; error = 2;}
	if (nmy<=0)     {cerr << "Error: nmy ("    << nmy    << ") <= 0"                  << endl; error = 2;}
	if (ylo >= yhi) {cerr << "Error: ylo ("    << ylo    << ") >= yhi(" << yhi << ")" << endl; error = 2;}
	return error;
}

TH1D* Unfold2D::ProjectionX (const TH1* h, const char* name, const char* title, Option_t* opt)
{
	const TH2* h2= dynamic_cast<const TH2*>(h);
	TH1D* h1= h2->ProjectionX (name, 1, h->GetNbinsY(), opt);
	if (title) h1->SetTitle (title);
	return h1;
}

TH1D* Unfold2D::ProjectionY (const TH1* h, const char* name, const char* title, Option_t* opt)
{
	const TH2* h2= dynamic_cast<const TH2*>(h);
	TH1D* h1= h2->ProjectionY (name, 1, h->GetNbinsX(), opt);
	if (title) h1->SetTitle (title);
	return h1;
}

TProfile* Unfold2D::ProfileX (const TH1* h, const char* name, const char* title, Option_t* opt)
{
	const TH2* h2= dynamic_cast<const TH2*>(h);
	TProfile* h1= h2->ProfileX (name, 1, h->GetNbinsY(), opt);
	if (title) h1->SetTitle (title);
	return h1;
}

TH1D* Unfold2D::TProfile2TH1D(TProfile *pf) 
{
	int nbins = pf->GetNbinsX();
	TH1D *hnew = new TH1D("hnew","", nbins, pf->GetXaxis()->GetXbins()->GetArray());  
	hnew->SetName(Form("h%s",pf->GetName()));
	hnew->SetTitle(Form("%s",pf->GetTitle()));

	for(int i = 0; i<nbins; i++) {
		hnew->SetBinContent(i+1,pf->GetBinContent(i+1));
		hnew->SetBinError(i+1,pf->GetBinError(i+1));
	}
	
	return hnew;
}


//==============================================================================
// Set histogram Y-axis display range
//==============================================================================

void Unfold2D::setmax (TH1* h,
		const TH1* h1, const TH1* h2, const TH1* h3,
		const TH1* h4, const TH1* h5, const TH1* h6)
{
	// Get the maximum y value of up to 7 histograms
	// Add 10% to match behaviour of ROOT's automatic scaling
	Double_t maxval= h1 ? h1->GetMaximum() : -DBL_MAX;
	if (h2 && h2->GetMaximum() > maxval) maxval= h2->GetMaximum();
	if (h3 && h3->GetMaximum() > maxval) maxval= h3->GetMaximum();
	if (h4 && h4->GetMaximum() > maxval) maxval= h4->GetMaximum();
	if (h5 && h5->GetMaximum() > maxval) maxval= h5->GetMaximum();
	if (h6 && h6->GetMaximum() > maxval) maxval= h6->GetMaximum();
	h->SetMinimum (0.0);
	if (maxval > h->GetMaximum()) h->SetMaximum (1.1*maxval);
}


//==============================================================================
// Set X and Y names for all histograms 
// This will determine what variables to be used in the analysis
//==============================================================================
void Unfold2D::SetXYname(TString xname, TString yname) {
	XvariableName = xname;
	YvariableName = yname;
}

void Unfold2D::SetTrigName(TString tname) {
	TrigName = tname;
}

void Unfold2D::SetTranCharge(TString tname) {
	TranCharge = tname;
}

void Unfold2D::SetExcludeJPTrig(Int_t exclude) {
	ExcludeOpt = exclude;
}

void Unfold2D::Legend (TLegend*& legend, TH1* truth, TH1* fake, TH1* meas, TH1* reco, TF1* ff, TF1* tf)
{
	legend= new TLegend (0.7, (tf ? 0.62 : reco ? 0.72 : 0.75), 0.894, 0.89);
	legend->SetTextSize(0.035);
	legend->SetTextFont(42);
	if (truth) legend->AddEntry (truth, "Particle-level",         "L");
	if (fake)
		legend->AddEntry (fake,  "fakes",         "L");
	legend->AddEntry (meas,  "Detector-level",      "L");
	if (reco)
		legend->AddEntry (reco,  "unfolded", "P");
	if (ff)
		legend->AddEntry (ff->GetHistogram(), "measured fit", "L");
	if (tf)
		legend->AddEntry (tf,    "truth fit",     "L");
	legend->Draw();
}


TH2D* Unfold2D::CorrelationHist (const TMatrixD& cov,
		const char* name, const char* title,
		Double_t lo, Double_t hi)
{
	Int_t nb= cov.GetNrows();
	TH2D* h= new TH2D (name, title, nb, lo, hi, nb, lo, hi);
	h->SetAxisRange (-1.0, 1.0, "Z");
	for(int i=0; i < nb; i++)
		for(int j=0; j < nb; j++) {
			Double_t Viijj= cov(i,i)*cov(j,j);
			if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
		}
	return h;
}

//==============================================================================
// Use wide bin
// somehow the rebin to wide bin does not work... Need to use initwidebin instead
//==============================================================================
TH2D* Unfold2D::RebinAs(TH2D *old, TH2D *model)		// rebin old as the same binining as model
{
	int Nx = model->GetNbinsX();
	double *xbins = new double[Nx+1];
	for(int ibin = 0; ibin<Nx; ibin++) {
		*(xbins+ibin) = model->GetXaxis()->GetBinLowEdge(ibin+1);
	}
	*(xbins+Nx) = model->GetXaxis()->GetBinLowEdge(Nx) + model->GetXaxis()->GetBinWidth(Nx);
	int Ny = model->GetNbinsY();
	double *ybins = new double[Ny+1];
	for(int ibin = 0; ibin<Ny; ibin++) {
		*(ybins+ibin) = model->GetYaxis()->GetBinLowEdge(ibin+1);
	}
	*(ybins+Ny) = model->GetYaxis()->GetBinLowEdge(Ny) + model->GetYaxis()->GetBinWidth(Ny);

	return Rebin2DHisto(old,Nx, xbins, Ny, ybins);
}

TH2D* Unfold2D::RebinX2DHisto(TH2D *old)	
{
	//const int WNbins = 15;
        //double Wptbins[WNbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};

	int Ny = old->GetNbinsY();
	double *ybins = new double[Ny+1];
	for(int ibin = 0; ibin<Ny; ibin++) {
		*(ybins+ibin) = old->GetYaxis()->GetBinLowEdge(ibin+1);
	}
	*(ybins+Ny) = old->GetYaxis()->GetBinLowEdge(Ny) + old->GetYaxis()->GetBinWidth(Ny);

	TH2D *hnew = Rebin2DHisto(old, WNbins, Wptbins, Ny, ybins);
	return hnew;
}


TH2D* Unfold2D::Rebin2DHisto(TH2D *old, int Nx, double *xbins, int Ny, double *ybins)
{
	TH2D *h = new TH2D(Form("Rebin%s",old->GetName()),old->GetTitle(),Nx,xbins,Ny,ybins);	
	TAxis *xaxis = old->GetXaxis();
	TAxis *yaxis = old->GetYaxis();
	for (int j=1; j<=yaxis->GetNbins();j++) {
		for (int i=1; i<=xaxis->GetNbins();i++) {
			h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j));
		}
	}
	h->Sumw2();
	return h;
}

TH2D* Unfold2D::InitWideX2DHisto(TString name, TString title, int Ny, double ylo, double yhi) 
{
	//const int WNbins = 15;
        //double Wptbins[WNbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};

	if(Ny<=0) {cout<<"Err: "<<__PRETTY_FUNCTION__<<", Ny = "<<Ny<<"<=0"<<endl;return NULL;}
	double *ybins = new double[Ny+1];
	for(int ibin = 0; ibin<Ny; ibin++) {
		*(ybins+ibin) = ylo + ibin*(yhi-ylo)/Ny;	
	}
	*(ybins+Ny) = yhi;

	TH2D *h = new TH2D(name, title, WNbins, Wptbins, Ny,ybins);
	return h;
}



#endif
