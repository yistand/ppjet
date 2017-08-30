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
#include "ReWeightByNF.h"
#include "ScaleByY.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <cmath>
#include <iostream>
#include <algorithm>

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2F.h"
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
//double Unfold2D::Wptbins[] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,100};
//double Unfold2D::Wptbins[] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55};	//12JetBin
double Unfold2D::Wptbins[] = {0,1,3,5,7,9,11,15,20,25,35,45,55};// 12JetBinv2

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
	cout<<"flagjetweight = args[14];		// 1: if JPs (mix JP0, JP1, JP2), reweight 2D histograms with Neutral Fraction vs jet pT distribution same as MB. 0: no reweighting procedure"<<endl;
	cout<<"						// if flagjetweight==1, WIDEBIN will be reset to 1, no mattter the input. Because, the reweighting procedure used wide jet pt bin, here needs to keep consistence"<<endl;
	cout<<"flagscaley = args[15];		// 1: if JPs (mix JP0, JP1, JP2), scale each JP Y down to have same mean value as MB vs jet pT . 0: no scaling procedure"<<endl;

	cout<<endl<<endl;

	cout<<"XvariableName = "<<XvariableName<<endl;
	cout<<"YvariableName = "<<YvariableName<<endl;

	cout<<endl<<endl;
}

//==============================================================================
// Constructors and destructor
//==============================================================================

Unfold2D::Unfold2D (TString name): inputname(name), XvariableName("X"), YvariableName("Y"), TrigName("JP2"), TranCharge("Charged"), ExcludeOpt(0), NoFakeOpt(false), NoLossOpt(false), CHANGEPRIOR(0), OPT_TPCSYS(false), OPT_TPCSYSPM(false), OPT_TPCSYSABS(false), OPT_RC02MC05(false)
{
	Reset();
	SetDefaultParms();
	FlagDefaultCalled = 1;
}

Unfold2D::Unfold2D (const char* name, int argc, const char* const* argv)
: inputname(name), XvariableName("X"), YvariableName("Y"), TrigName("JP2"), TranCharge("Charged"), ExcludeOpt(0), NoFakeOpt(false), NoLossOpt(false), CHANGEPRIOR(0), OPT_TPCSYS(false), OPT_TPCSYSPM(false), OPT_TPCSYSABS(false), OPT_RC02MC05(false)
{
	Reset();
	if(argc<14){ Help(); return; }
	SetParms(argv);
	FlagDefaultCalled = 0;
}


Unfold2D::~Unfold2D()
{
	delete response; response= 0;
	delete unfold;   unfold=   0;
	if(ftrain->IsOpen()) ftrain->Close();
	if(ftout->IsOpen()) ftout->Close();
}

//============================================================
//
//	Set Parameters
//
//============================================================

void Unfold2D::SetParms (const char* const* argv) 
{
	const int Npar = 16;
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

	if(args[14]>0) {
		flagjetweight = true;		// Use Neutral Fraction ratio weight for JPs, so that it has the same NF distribution per jet pt as MB
		WIDEBIN = true;			// When NF reweighting is used, NEED to use WIDEBIN = TRUE as this is what used for Reweighting procedure
	}
	else flagjetweight = false;

	if(args[14]<=0 && args[15]>0) {			// When NF reweighting is used, NEED to use flagscaley == FALSE so that we won't over do correction
		flagscaley = true;		// Scale <Y> in JP to be same as MB 
	}
	else flagscaley = false;
}

void Unfold2D::SetDefaultParms() 		// use default parameterization
{
	FlagDefaultCalled = 1;

	float argsNF[16] = {1, 60, 100, 60, 100, 0, 60, 0, 1, 0, 1, 1, 4, 0, 0, 0};	// neutralfrac
	float argsNtrk[16] 	  = {1, 50, 40, 50, 40, 0, 100, -0.5, 39.5, 0, 1, 1, 4, 1, 1, 0};	// TranNtrk, TranTotNtrk
	float argsPtSum[16] = {1, 50, 500, 50, 500, 0, 100, 0, 50, 0, 1, 1, 4, 1, 1, 0};		// PtSum
	float argsPtAve[16] = {1, 50, 1000, 50, 1000, 0, 100, 0, 10, 0, 1, 1, 4, 1, 1, 0};		// PtAve
	//		method, ntx, nty, nmx, nmy, xlo, xhi, ylo, yhi, overflow, verbose, doerror, regparm, WIDEBIN, flagjetweight, flagscaley

	if(YvariableName.Contains("PtAve")) SetParms(argsPtAve);
	else if(YvariableName.Contains("PtSum")) SetParms(argsPtSum);
	else if(YvariableName.Contains("NeutralFrac")) SetParms(argsNF);
	else SetParms(argsNtrk);
	//SetParms(args);	
}

void Unfold2D::PrintParms() 
{
	cout<<"method = "<<method<<endl;
	if(!WIDEBIN) {			// WIDEBIN will change those binning values
		cout<<"ntx = "<<ntx<<endl;
		cout<<"nty = "<<nty<<endl;
		cout<<"nmx = "<<nmx<<endl;
		cout<<"nmy = "<<nmy<<endl;
		cout<<"xlo = "<<xlo<<endl;
		cout<<"xhi = "<<xhi<<endl;
		cout<<"ylo = "<<ylo<<endl;
		cout<<"yhi = "<<yhi<<endl;
	}
	cout<<"overflow = "<<overflow<<endl;
	cout<<"verbose = "<<verbose<<endl;
	cout<<"doerror = "<<doerror<<endl;
	cout<<"regparm = "<<regparm<<endl;
	cout<<"WIDEBIN = "<<WIDEBIN<<endl;
	cout<<"flagjetweight = "<<flagjetweight<<endl;
	cout<<"flagscaley = "<<flagscaley<<endl;
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

	Float_t Mcjpt;
	Float_t Mcjaspt;
	Float_t Mcjneutralfrac;
	Int_t Mcjconstntrk;
	Int_t Mcjconstchargentrk;
	Int_t  McLeadAreaNtrk;
	Int_t  McSubAreaNtrk;
	Int_t  McTranMaxNtrk;
	Int_t  McTranMinNtrk;
	Float_t McLeadAreaPtSum;
	Float_t McSubAreaPtSum;
	Float_t McTranMaxPtSum;
	Float_t McTranMinPtSum;
	Float_t McTrkLeadAreaPt[MAXARRAY];
	Float_t McTrkSubAreaPt[MAXARRAY];
	Float_t McTrkTranMaxPt[MAXARRAY];
	Float_t McTrkTranMinPt[MAXARRAY];
	Float_t McTrkLeadAreaPhi[MAXARRAY];
	Float_t McTrkSubAreaPhi[MAXARRAY];
	Float_t McTrkTranMaxPhi[MAXARRAY];
	Float_t McTrkTranMinPhi[MAXARRAY];
	Float_t McTrkLeadAreaEta[MAXARRAY];
	Float_t McTrkSubAreaEta[MAXARRAY];
	Float_t McTrkTranMaxEta[MAXARRAY];
	Float_t McTrkTranMinEta[MAXARRAY];
	Int_t McTrkLeadAreaId[MAXARRAY];
	Int_t McTrkSubAreaId[MAXARRAY];
	Int_t McTrkTranMaxId[MAXARRAY];
	Int_t McTrkTranMinId[MAXARRAY];

	Float_t Rcjneutralfrac;
	Int_t Rcjconstntrk;
	Int_t Rcjconstchargentrk;
	Int_t  RcTranMaxNtrk;
	Int_t  RcTranMinNtrk;
	Int_t  RcLeadAreaNtrk;
	Int_t  RcSubAreaNtrk;
	Float_t RcLeadAreaPtSum;
	Float_t RcSubAreaPtSum;
	Float_t RcTranMaxPtSum;
	Float_t RcTranMinPtSum;
	Float_t RcTrkLeadAreaPt[MAXARRAY];
	Float_t RcTrkSubAreaPt[MAXARRAY];
	Float_t RcTrkTranMaxPt[MAXARRAY];
	Float_t RcTrkTranMinPt[MAXARRAY];
	Float_t RcTrkLeadAreaPhi[MAXARRAY];
	Float_t RcTrkSubAreaPhi[MAXARRAY];
	Float_t RcTrkTranMaxPhi[MAXARRAY];
	Float_t RcTrkTranMinPhi[MAXARRAY];
	Float_t RcTrkLeadAreaEta[MAXARRAY];
	Float_t RcTrkSubAreaEta[MAXARRAY];
	Float_t RcTrkTranMaxEta[MAXARRAY];
	Float_t RcTrkTranMinEta[MAXARRAY];
	Int_t RcTrkLeadAreaId[MAXARRAY];
	Int_t RcTrkSubAreaId[MAXARRAY];
	Int_t RcTrkTranMaxId[MAXARRAY];
	Int_t RcTrkTranMinId[MAXARRAY];


	if(WIDEBIN) pfxTrain= new TProfile("pfxtrain", "Training", WNbins, Wptbins);
	else pfxTrain= new TProfile ("pfxtrain", "Training", nmx, xlo, xhi);
	pfxTrain->SetLineColor(hTrain->GetLineColor());
	pfxTrain->GetXaxis()->SetTitle(XvariableName);
	pfxTrain->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) pfxTrainTrue= new TProfile("pfxtraintrue", "TrainTrue", WNbins, Wptbins);
	else pfxTrainTrue= new TProfile ("pfxtraintrue", "TrainTrue", nmx, xlo, xhi);
	pfxTrainTrue->SetLineColor(hTrainTrue->GetLineColor());
	pfxTrainTrue->GetXaxis()->SetTitle(XvariableName);
	pfxTrainTrue->GetYaxis()->SetTitle(YvariableName);


	//// fortest to remove 8 events cause a large fluctuation 
	//int countremove = 0;

	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}

	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename;
		if(TrigName.Contains("JP")&&(flagjetweight||flagscaley)) {
		// if flagjetweight==1, JP hist is weighted by MB NF distribution. Then we shall use MB embeding data (no trigger effect). This means we don't unfold trigger effect on jet pt distribution, but the neutral jet trigger effect should be already corrected by weighting procedure. However, if we don't do weighting by MB Neutral Fraction dist, we could also use JP embedding to unfold all trigger effects together
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_McPt02%s.root",  "MB",   PTBINS[i],TranCharge.Data(),STpcSys.Data());
		}
		else if(ExcludeOpt&&TrigName.Contains("JP")) {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_excluded_McPt02%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data(),STpcSys.Data());
		}
		else {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_McPt02%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data(),STpcSys.Data());
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

		//tree->SetBranchAddress("Mcj1pt",&McJet);
		tree->SetBranchAddress("Mcj1pt",&Mcjpt);
		tree->SetBranchAddress("Mcjaspt",&Mcjaspt);
		tree->SetBranchAddress("Mcj1neutralfrac",&Mcjneutralfrac);
		tree->SetBranchAddress("Mcj1constntrk",&Mcjconstntrk);
		tree->SetBranchAddress("Mcj1constchargentrk",&Mcjconstchargentrk);
		tree->SetBranchAddress("McLeadAreaNtrk",&McLeadAreaNtrk);
		tree->SetBranchAddress("McSubAreaNtrk",&McSubAreaNtrk);
		tree->SetBranchAddress("McTranMaxNtrk",&McTranMaxNtrk);
		tree->SetBranchAddress("McTranMinNtrk",&McTranMinNtrk);
		tree->SetBranchAddress("McLeadAreaPtSum",&McLeadAreaPtSum);
		tree->SetBranchAddress("McSubLeadAreaPtSum",&McSubAreaPtSum);
		tree->SetBranchAddress("McTranMaxPtSum",&McTranMaxPtSum);
		tree->SetBranchAddress("McTranMinPtSum",&McTranMinPtSum);
		tree->SetBranchAddress("McTrkLeadAreaPt",McTrkLeadAreaPt);
		tree->SetBranchAddress("McTrkSubAreaPt",McTrkSubAreaPt);
		tree->SetBranchAddress("McTrkTranMaxPt",McTrkTranMaxPt);
		tree->SetBranchAddress("McTrkTranMinPt",McTrkTranMinPt);
		tree->SetBranchAddress("McTrkLeadAreaPhi",McTrkLeadAreaPhi);
		tree->SetBranchAddress("McTrkSubAreaPhi",McTrkSubAreaPhi);
		tree->SetBranchAddress("McTrkTranMaxPhi",McTrkTranMaxPhi);
		tree->SetBranchAddress("McTrkTranMinPhi",McTrkTranMinPhi);
		tree->SetBranchAddress("McTrkLeadAreaEta",McTrkLeadAreaEta);
		tree->SetBranchAddress("McTrkSubAreaEta",McTrkSubAreaEta);
		tree->SetBranchAddress("McTrkTranMaxEta",McTrkTranMaxEta);
		tree->SetBranchAddress("McTrkTranMinEta",McTrkTranMinEta);
		tree->SetBranchAddress("McTrkLeadAreaId",McTrkLeadAreaId);
		tree->SetBranchAddress("McTrkSubAreaId",McTrkSubAreaId);
		tree->SetBranchAddress("McTrkTranMaxId",McTrkTranMaxId);
		tree->SetBranchAddress("McTrkTranMinId",McTrkTranMinId);
		//tree->SetBranchAddress("McPart",&McPart);
		
		tree->SetBranchAddress("Rcj1pt",&RcJet);
		tree->SetBranchAddress("Rcj1neutralfrac",&Rcjneutralfrac);
		tree->SetBranchAddress("Rcj1constntrk",&Rcjconstntrk);
		tree->SetBranchAddress("Rcj1constchargentrk",&Rcjconstchargentrk);
		tree->SetBranchAddress("RcLeadAreaNtrk",&RcLeadAreaNtrk);
		tree->SetBranchAddress("RcSubAreaNtrk",&RcSubAreaNtrk);
		tree->SetBranchAddress("RcTranMaxNtrk",&RcTranMaxNtrk);
		tree->SetBranchAddress("RcTranMinNtrk",&RcTranMinNtrk);
		tree->SetBranchAddress("RcLeadAreaPtSum",&RcLeadAreaPtSum);
		tree->SetBranchAddress("RcSubLeadAreaPtSum",&RcSubAreaPtSum);
		tree->SetBranchAddress("RcTranMaxPtSum",&RcTranMaxPtSum);
		tree->SetBranchAddress("RcTranMinPtSum",&RcTranMinPtSum);
		tree->SetBranchAddress("RcTrkLeadAreaPt",RcTrkLeadAreaPt);
		tree->SetBranchAddress("RcTrkSubAreaPt",RcTrkSubAreaPt);
		tree->SetBranchAddress("RcTrkTranMaxPt",RcTrkTranMaxPt);
		tree->SetBranchAddress("RcTrkTranMinPt",RcTrkTranMinPt);
		tree->SetBranchAddress("RcTrkLeadAreaPhi",RcTrkLeadAreaPhi);
		tree->SetBranchAddress("RcTrkSubAreaPhi",RcTrkSubAreaPhi);
		tree->SetBranchAddress("RcTrkTranMaxPhi",RcTrkTranMaxPhi);
		tree->SetBranchAddress("RcTrkTranMinPhi",RcTrkTranMinPhi);
		tree->SetBranchAddress("RcTrkLeadAreaEta",RcTrkLeadAreaEta);
		tree->SetBranchAddress("RcTrkSubAreaEta",RcTrkSubAreaEta);
		tree->SetBranchAddress("RcTrkTranMaxEta",RcTrkTranMaxEta);
		tree->SetBranchAddress("RcTrkTranMinEta",RcTrkTranMinEta);
		tree->SetBranchAddress("RcTrkLeadAreaMcId",RcTrkLeadAreaId);
		tree->SetBranchAddress("RcTrkSubAreaMcId",RcTrkSubAreaId);
		tree->SetBranchAddress("RcTrkTranMaxMcId",RcTrkTranMaxId);
		tree->SetBranchAddress("RcTrkTranMinMcId",RcTrkTranMinId);
		//tree->SetBranchAddress("RcPart",&RcPart);


		//// fortest to remove 8 events cause a large fluctuation 
		//// pythia pt bins: {2,3,4,5,7,9,11,15,20,25,35,inf}
		//if( (i==3 && runid==13049080 && eventid ==42) 
		//  &&(i==5 && runid==13051020 && eventid ==216) 
		//  &&(i==5 && runid==13055075 && eventid ==323) 
		//  &&(i==5 && runid==13057057 && eventid ==422) 
		//  &&(i==5 && runid==13059076 && eventid ==16) 
		//  &&(i==5 && runid==13059076 && eventid ==16) 
		//  &&(i==8 && runid==13047126 && eventid ==23) 
		//  &&(i==9 && runid==13051006 && eventid ==43) 
		//  &&(i==9 && runid==13055008 && eventid ==70) 
		//  ) {
		//	countremove ++;
		//	continue;	
		//}

		if(verbose && Nevents[i]>tree->GetEntries()) {cout<<__PRETTY_FUNCTION__<<": pT bin "<<PTBINS[i]<<" Read # of events["<<i<<"] = "<<tree->GetEntries()<<endl; }

		// Weight per pT bin
		double xsecweight = XSEC[i]/Nevents[i];
		if(Nevents[i]>tree->GetEntries()) {		// use all events, then NUMBEROFEVENT from include/CrossSectionPerpT.h
			xsecweight = XSEC[i]/NUMBEROFEVENT[i];
			if(verbose) { cout<<__PRETTY_FUNCTION__<<": pT bin "<<PTBINS[i]<<" Use "<<NUMBEROFEVENT[i]<<" for NUMBEROFEVENT in weight = XSEC["<<i<<"]/NUMBEROFEVENT["<<i<<"]"<<endl;}
		}
		weight = xsecweight;
		//weight = XSEC[i]/tree->GetEntries();			// 	in case we want to do half for train half for test, weight need to be same
		//cout<<"train weight = "<<weight<<endl;			// 

		//// if Dijest 
		//TRandom3 *rndm = new TRandom3();
		//rndm->SetSeed(0);

		for (Int_t j= 0; j<Nevents[i] && j<tree->GetEntries(); j++) {		// loop over entries
		//for (Int_t j= 0; j<tree->GetEntries(); j+=2) {		// in case we want to do half for train half for test, *random* select (even and odd for train, test)
			tree->GetEntry(j);


			if(CHANGEPRIOR) {		// if need to change prio, use weight to change. In JP unfolded, we used (170412) MB pythia, but applied NF weight to JP data. However, as the jet shape dramatic different from truth, it pulled the unfolding. In order to avoid this situation, we reweighted the pythia MB particle-level to have same jet pt as measured data. So we don't start too far away
				weight = xsecweight*Weight2ChargePrior(Mcjpt);
				//If use a function form to change prior
				//TF1 *func1 = new TF1("func1","exp(-[0]*x)");
				//func1->SetParameter(0,0.1);
				//weight = xsecweight*Weight2ChargePrior(func1, Mcjpt);	
				//func1->Delete();	// need to delete *new*
				//Or use the default function: 
				//weight = xsecweight*Weight2ChargePrior(NULL, Mcjpt);	
			}

			//// if Dijet 
			//if(Mcjaspt<=0) {
			//	continue;	
			//	//if(rndm->Uniform()>0.5) {
			//	//	McJet = Mcjpt;
			//	//}
			//	//else {
			//	//	McJet = Mcjaspt;
			//	//}
			//}
			McJet = Mcjpt;
			

			flag = -999;
			//////JPs
			if(TrigName.Contains("JP",TString::kIgnoreCase)&&(!flagjetweight) && (!flagscaley)) {
				if(flagIsTrigger && flagtrigmatch && (flagMatch2Lead||flagMatch2Sub) && (RcJet>0) && (Rcjneutralfrac<0.9) && (Rcjneutralfrac>0)) flag = 1; 		
				//else if(McJet>0 && (RcJet<=0 || (!flagtrigmatch)) ) flag = 0;		// non-zero Mc, zero Rc
				else if(McJet>0 && (RcJet<=0 || (!flagtrigmatch) || !(Rcjneutralfrac<0.9 && Rcjneutralfrac>0)) ) flag = 0;		// non-zero Mc, zero Rc ; or  non-zero Mc, non-zero Rc But Rc does not pass Rc jneutralfrac cuts or trigmatch cut
				//else if((RcJet>0&&flagtrigmatch&&flagIsTrigger&&Rcjneutralfrac<0.9&&Rcjneutralfrac>0) && McJet<=0) flag = -1;		// non-zero Rc, zero Mc		additional neutral fraction here. It is also included in flagMatch2LeadGood and flagMatch2SubGood
				else if((RcJet>0&&flagtrigmatch&&flagIsTrigger&&Rcjneutralfrac<0.9&&Rcjneutralfrac>0) && (McJet<=0 || !(flagMatch2Lead||flagMatch2Sub)) ) flag = -1;		// non-zero Rc, zero Mc;  or non-zero Rc, non-zero Mc but Rc not matched to Mc leading/sublead 		additional neutral fraction here. It is also included in flagMatch2LeadGood and flagMatch2SubGood
				//cout<<"flag = "<<flag<<endl;	
			}
			else {
				//if(!flagIsTrigger) continue;			// no correction for trig now..		
				////MB only need the following 3 lines
				//if(flagIsTrigger && (flagMatch2Lead||flagMatch2Sub) && RcJet>0 && Rcjneutralfrac<0.9)     flag = 1; 				//2016.11.15 add naive simulated trigger for VPDMB	// 2017.01.04 don't use this naive simulated trigger: it doesn't describe data
				//else if(McJet>0 && RcJet<=0) flag = 0;		// non-zero Mc, zero Rc
				//else if((flagIsTrigger && RcJet>0 && Rcjneutralfrac<0.9) && McJet<=0) flag = -1;		// non-zero Rc, zero Mc			//2016.11.15 add naive simulated trigger for VPDMB	// 2017.01.04 don't use this naive simulated trigger: it doesn't describe data
				//
				// use ReWeightByNF which requires Rcjneutralfrac!=0 or 0.9	2017.01.04
				if((flagMatch2Lead||flagMatch2Sub) && RcJet>0 && Rcjneutralfrac<0.9 && Rcjneutralfrac>0)      flag = 1; 				
				else if(McJet>0 && (RcJet<=0 || !(Rcjneutralfrac<0.9 && Rcjneutralfrac>0)) ) flag = 0;		// inefficiency events: non-zero Mc, zero Rc; or  non-zero Mc, non-zero Rc But Rc does not pass Rc jneutralfrac cuts
				else if(RcJet>0 && Rcjneutralfrac<0.9 && Rcjneutralfrac>0 && (McJet<=0 || !(flagMatch2Lead||flagMatch2Sub)) ) flag = -1;		// fake events: non-zero Rc, zero Mc; or non-zero Rc, non-zero Mc but Rc not matched to Mc leading/sublead

				// simply case 
				//if(flagMatch2Lead||flagMatch2Sub) flag = 1; 								
				//else if(McJet>=0 && RcJet<=0) flag = 0;		// non-zero Mc, zero Rc			
				//else if(RcJet>=0 && McJet<=0) flag = -1;		// non-zero Rc, zero Mc		
				
			}

			if(OPT_RC02MC05) {		// MC pt>0.5 only
				if(!(YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase)||(YvariableName.Contains("TranPtAve",TString::kIgnoreCase)&&!(YvariableName.Contains("EventWise",TString::kIgnoreCase))))) {
					cout<<"WARNING!!!!!!!!!! in "<<__PRETTY_FUNCTION__<<" OPT_RC02MC05=true is only implemented for TranTotNtrk or particle-wise TranPtAve only. Will reset OPT_RC02MC05 = false instead for "<<YvariableName<<endl;
					OPT_RC02MC05 = false;
				}
				else if(YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase)) {
					RcPart = RcTranMaxNtrk+RcTranMinNtrk;		// RC part not change, still pT>0.2
					
					McPart = 0;					// MC part only count pT>0.5 particles
					for(int im = 0; im<McTranMaxNtrk; im++) {
						if(McTrkTranMaxPt[im]>0.5) McPart+=1;
					}
					for(int im = 0; im<McTranMinNtrk; im++) {
						if(McTrkTranMinPt[im]>0.5) McPart+=1;
					}
				}
				// we will deal with PtAve later in the loop
			}
			{	// as we need to deal with PtAve for OPT_RC02MC05==true case, we just continue here 
			//else {

				if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
				{
					if(McTranMaxNtrk>McTranMinNtrk) McPart = static_cast<float>(McTranMaxNtrk);
					else McPart = static_cast<float>(McTranMinNtrk);
					if(RcTranMaxNtrk>RcTranMinNtrk) RcPart = static_cast<float>(RcTranMaxNtrk);
					else RcPart = static_cast<float>(RcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
				{
					if(McTranMaxNtrk<McTranMinNtrk) McPart = static_cast<float>(McTranMaxNtrk);
					else McPart = static_cast<float>(McTranMinNtrk);
					if(RcTranMaxNtrk<RcTranMinNtrk) RcPart = static_cast<float>(RcTranMaxNtrk);
					else RcPart = static_cast<float>(RcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranDiffNtrk",TString::kIgnoreCase)) 
				{
					McPart = fabs(McTranMaxNtrk-McTranMinNtrk);
					RcPart = fabs(RcTranMaxNtrk-RcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranNtrk",TString::kIgnoreCase))
				{
					McPart = 0.5*(McTranMaxNtrk+McTranMinNtrk);
					RcPart = 0.5*(RcTranMaxNtrk+RcTranMinNtrk);
				}	
				else if(!OPT_RC02MC05 && YvariableName.Contains("TranTotNtrk",TString::kIgnoreCase))	// as we change OPT_RC02MC05 routine, now OPT_RC02MC05 can be true here...
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
				else if(YvariableName.Contains("TranMaxPtAveEventWise",TString::kIgnoreCase))
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
				else if(YvariableName.Contains("TranMinPtAveEventWise",TString::kIgnoreCase)) 
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
				else if(YvariableName.Contains("TranPtAveEventWise",TString::kIgnoreCase))
				{
					McPart = (McTranMaxPtSum+McTranMinPtSum)/(McTranMaxNtrk+McTranMinNtrk);
					RcPart = (RcTranMaxPtSum+RcTranMinPtSum)/(RcTranMaxNtrk+RcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranTotPtAveEventWise",TString::kIgnoreCase))
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
				else if(YvariableName.Contains("LeadAreaPtAveEventWise",TString::kIgnoreCase)) 
				{
					McPart = McLeadAreaPtSum/McLeadAreaNtrk;
					RcPart = RcLeadAreaPtSum/RcLeadAreaNtrk;
				}
				else if(YvariableName.Contains("SubAreaPtAveEventWise",TString::kIgnoreCase)) 
				{
					McPart = McSubAreaPtSum/McSubAreaNtrk;
					RcPart = RcSubAreaPtSum/RcSubAreaNtrk;
				}
				else if(YvariableName.Contains("NeutralFrac",TString::kIgnoreCase)) 
				{
					McPart = Mcjneutralfrac;
					RcPart = Rcjneutralfrac;
				}
				else if(YvariableName.Contains("PtAve",TString::kIgnoreCase)&&!(YvariableName.Contains("EventWise",TString::kIgnoreCase))) 		// Particle-wise variable, need to loop over Mc and Rc particles
				{
					std::vector<float> vMcTrkTotPt;
					std::vector<float> vMcTrkTotPhi;
					std::vector<float> vMcTrkTotEta;
					std::vector<int> vMcTrkTotId;
					std::vector<float> vRcTrkTotPt;
					std::vector<float> vRcTrkTotPhi;
					std::vector<float> vRcTrkTotEta;
					std::vector<int> vRcTrkTotId;

					std::vector<int> vMcTrkTotMatchedId;

					//cout<<endl;
					//cout<<"McJet = "<<McJet<<endl;
					//cout<<"RcJet = "<<RcJet<<endl;
					//cout<<" "<<McTranMaxNtrk+McTranMinNtrk<<" McTran"<<endl;
					//cout<<" "<<RcTranMaxNtrk+RcTranMinNtrk<<" RcTran"<<endl;

					// Merge TranMax and TranMin
					// Mc
					if(flag!=-1) 	// one-to-one Mc jet to Rc jet Or no Rc jet, only Mc jet
					{
						if(YvariableName.Contains("LeadPtAve",TString::kIgnoreCase)) {
							for(int im = 0; im<McLeadAreaNtrk; im++) 
							{
								vMcTrkTotPt.push_back(McTrkLeadAreaPt[im]);
								vMcTrkTotPhi.push_back(McTrkLeadAreaPhi[im]);
								vMcTrkTotEta.push_back(McTrkLeadAreaEta[im]);
								vMcTrkTotId.push_back(McTrkLeadAreaId[im]);
								vMcTrkTotMatchedId.push_back(0);			// assuming not matched
							}
						}
						else if(YvariableName.Contains("SubPtAve",TString::kIgnoreCase)) {
							for(int im = 0; im<McSubAreaNtrk; im++) 
							{
								vMcTrkTotPt.push_back(McTrkSubAreaPt[im]);
								vMcTrkTotPhi.push_back(McTrkSubAreaPhi[im]);
								vMcTrkTotEta.push_back(McTrkSubAreaEta[im]);
								vMcTrkTotId.push_back(McTrkSubAreaId[im]);
								vMcTrkTotMatchedId.push_back(0);			// assuming not matched
							}
						}
						else {	// default is TranPtAve
							for(int im = 0; im<McTranMaxNtrk; im++) 
							{
								if(OPT_RC02MC05 && McTrkTranMaxPt[im]<0.5) continue;		// deal with TranPtAve for OPT_RC02MC05 case

								vMcTrkTotPt.push_back(McTrkTranMaxPt[im]);
								vMcTrkTotPhi.push_back(McTrkTranMaxPhi[im]);
								vMcTrkTotEta.push_back(McTrkTranMaxEta[im]);
								vMcTrkTotId.push_back(McTrkTranMaxId[im]);
								vMcTrkTotMatchedId.push_back(0);			// assuming not matched
							}
							for(int im = 0; im<McTranMinNtrk; im++) 
							{
								if(OPT_RC02MC05 && McTrkTranMinPt[im]<0.5) continue;		// deal with TranPtAve for OPT_RC02MC05 case

								vMcTrkTotPt.push_back(McTrkTranMinPt[im]);
								vMcTrkTotPhi.push_back(McTrkTranMinPhi[im]);
								vMcTrkTotEta.push_back(McTrkTranMinEta[im]);
								vMcTrkTotId.push_back(McTrkTranMinId[im]);
								vMcTrkTotMatchedId.push_back(0);			// assuming not matched
							}
						}
					}



					// Rc
					if(flag!=0) 	// ! No RcJet: one-to-one Rc jet to Mc jet Or Fake Rc jet, no Mc jet
					{
						if(YvariableName.Contains("LeadPtAve",TString::kIgnoreCase)) {
							for(int ir = 0; ir<RcLeadAreaNtrk; ir++) 
							{
								vRcTrkTotPt.push_back(RcTrkLeadAreaPt[ir]);
								vRcTrkTotPhi.push_back(RcTrkLeadAreaPhi[ir]);
								vRcTrkTotEta.push_back(RcTrkLeadAreaEta[ir]);
								vRcTrkTotId.push_back(RcTrkLeadAreaId[ir]);
							}
						}
						else if(YvariableName.Contains("SubPtAve",TString::kIgnoreCase)) {
							for(int ir = 0; ir<RcSubAreaNtrk; ir++) 
							{
								vRcTrkTotPt.push_back(RcTrkSubAreaPt[ir]);
								vRcTrkTotPhi.push_back(RcTrkSubAreaPhi[ir]);
								vRcTrkTotEta.push_back(RcTrkSubAreaEta[ir]);
								vRcTrkTotId.push_back(RcTrkSubAreaId[ir]);
							}
						}
						else {  // default: TranPtAve
							for(int ir = 0; ir<RcTranMaxNtrk; ir++) 
							{
								vRcTrkTotPt.push_back(RcTrkTranMaxPt[ir]);
								vRcTrkTotPhi.push_back(RcTrkTranMaxPhi[ir]);
								vRcTrkTotEta.push_back(RcTrkTranMaxEta[ir]);
								vRcTrkTotId.push_back(RcTrkTranMaxId[ir]);
							}
							for(int ir = 0; ir<RcTranMinNtrk; ir++) 
							{
								vRcTrkTotPt.push_back(RcTrkTranMinPt[ir]);
								vRcTrkTotPhi.push_back(RcTrkTranMinPhi[ir]);
								vRcTrkTotEta.push_back(RcTrkTranMinEta[ir]);
								vRcTrkTotId.push_back(RcTrkTranMinId[ir]);
							}
						}
					}


					
					// Loop over Rc
					if(flag!=0) 	// ! No RcJet: one-to-one Mc jet to Mc jet Or Fake Rc jet, no Mc jet
					{
						for(int ir = 0; ir<vRcTrkTotPt.size(); ir++)
						{
							RcPart = vRcTrkTotPt.at(ir);

							hTrain->Fill(RcJet, RcPart, weight);		
							pfxTrain->Fill(RcJet, RcPart, weight);
							
							// Find the Mc id matched to Rc id
							int iter_matched = -1;
							iter_matched = Unfold2D::Find(vMcTrkTotId, vRcTrkTotId.at(ir));
							
							if(iter_matched!=-1 && flag==1) 	// match succesfully Rc Id to Mc Id	 && Rc jet Matched to Mc jet
							{
								McPart = vMcTrkTotPt.at(iter_matched);
								vMcTrkTotMatchedId.at(iter_matched) = 1;		// Mark the matched MC
								response->Fill(RcJet, RcPart, McJet, McPart, weight);
								hXRcVsMc->Fill(McJet, RcJet, weight);
								hYRcVsMc->Fill(McPart, RcPart, weight);
								
								//cout<<"Match RcJet = "<<RcJet<<" RcPart = "<<RcPart<<" RcPhi = "<< vRcTrkTotPhi.at(ir) << " RcEta = "<< vRcTrkTotEta.at(ir) <<" McJet = "<<McJet<<" McPart = "<<McPart<<" McPhi = "<< vMcTrkTotPhi.at(iter_matched)<<" McEta = "<<vMcTrkTotEta.at(iter_matched)<<" weight = "<<weight<<endl;
							}
							else if(!NoFakeOpt) 		// No Mc Id or No matched Mc jet
							{
								hTrainFake->Fill(RcJet, RcPart, weight);
								response->Fake(RcJet, RcPart, weight);
								hXRcVsMc->Fill(-1, RcJet, weight);
								hYRcVsMc->Fill(-1, RcPart, weight);
								//cout<<"Fake RcJet = "<<RcJet<<" RcPart = "<<RcPart<<" RcPhi = "<< vRcTrkTotPhi.at(ir) << " RcEta = "<< vRcTrkTotEta.at(ir)<<" weight = "<<weight<<endl;
							}
						}
					}

					// Loop over Mc 
					if(flag!=-1) 
					{
						for(int im = 0; im<vMcTrkTotPt.size(); im++)
						{
							McPart = vMcTrkTotPt.at(im);

							hTrainTrue->Fill(McJet, McPart, weight);
							pfxTrainTrue->Fill(McJet, McPart, weight);

							if((!vMcTrkTotMatchedId.at(im) || flag==0) && !NoLossOpt) 		// Not matched to Rc Id or No good Rc Jet
							{
								response->Miss(McJet, McPart, weight);
								hXRcVsMc->Fill(McJet, -1, weight);
								hYRcVsMc->Fill(McPart, -1, weight);
								//cout<<"Miss McJet = "<<McJet<<" McPart = "<<McPart<<" McPhi = "<< vMcTrkTotPhi.at(im)<<" McEta = "<<vMcTrkTotEta.at(im)<<" weight = "<<weight<<endl;
							}

						}
					}

				}
				else if(YvariableName.Contains("JetNtrk",TString::kIgnoreCase))
				{
					McPart = Mcjconstntrk;
					RcPart = Rcjconstntrk;
				}
				else if(YvariableName.Contains("JetChargeNtrk",TString::kIgnoreCase))
				{
					McPart = Mcjconstchargentrk;
					RcPart = Rcjconstchargentrk;
				}
				else // default: TranNtrk
				{
					McPart = 0.5*(McTranMaxNtrk+McTranMinNtrk);
					RcPart = 0.5*(RcTranMaxNtrk+RcTranMinNtrk);
					
					cout<<"WARNING!!!!!!!!!! in "<<__PRETTY_FUNCTION__<<" "<<YvariableName<<" is not valid YvariableName name. Will use TranNtrk instead!"<<endl;
					YvariableName = "TranNtrk";
				}
			}	//end if else OPT_RC02MC05


			// Fill histograms
			if( !(YvariableName.Contains("PtAve",TString::kIgnoreCase)&&(!YvariableName.Contains("EventWise",TString::kIgnoreCase))) ) {
				if(flag==1) {			// one-to-one Mc to Rc matched
					hTrainTrue->Fill(McJet, McPart, weight);	
					hTrain->Fill(RcJet, RcPart, weight);	
					response->Fill(RcJet, RcPart,McJet, McPart, weight);
					pfxTrainTrue->Fill(McJet, McPart, weight);	
					pfxTrain->Fill(RcJet, RcPart, weight);	
					hXRcVsMc->Fill(McJet, RcJet, weight);
					hYRcVsMc->Fill(McPart, RcPart, weight);
				}
				else if(flag==0 && !NoLossOpt) {		// true Mc, no Rc
					hTrainTrue->Fill(McJet, McPart, weight);
					response->Miss(McJet, McPart, weight);
					pfxTrainTrue->Fill(McJet, McPart, weight);
					hXRcVsMc->Fill(McJet, -1, weight);
					hYRcVsMc->Fill(McPart, -1, weight);
				}
				else if(flag==-1 && !NoFakeOpt){				// Fake Rc, no Mc
					hTrain->Fill(RcJet,RcPart, weight);
					hTrainFake->Fill(RcJet,RcPart, weight);
					response->Fake(RcJet,RcPart, weight);
					pfxTrain->Fill(RcJet,RcPart, weight);
					hXRcVsMc->Fill(-1, RcJet, weight);
					hYRcVsMc->Fill(-1, RcPart, weight);
				}
			}
		}	// Event loop end
		//// if Dijet 
		//rndm->Delete(); 

		ftrain->Close();
	}
	//// fortest to remove 3 events cause a large fluctuation 
	//cout<<" -- Removed "<<countremove<<" events"<<endl;

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
			hXRcVsMc->Fill(McJet, RcJet, weight);
			hYRcVsMc->Fill(McPart, RcPart, weight);
		}
		else if(flag==0 && !NoLossOpt) {		// true Mc, no Rc
			hTrainTrue->Fill(McJet, McPart, weight);
			response->Miss(McJet,McPart, weight);
			hXRcVsMc->Fill(McJet, -1, weight);
			hYRcVsMc->Fill(McPart, -1, weight);
		}
		else if(!NoFakeOpt){				// Fake Rc, no Mc
			hTrain->Fill(RcJet,RcPart, weight);
			hTrainFake->Fill(RcJet,RcPart, weight);
			response->Fake(RcJet,RcPart, weight);
			hXRcVsMc->Fill(-1, RcJet, weight);
			hYRcVsMc->Fill(-1, RcPart, weight);
		}
	}

	return 1;
}

//==============================================================================
// Train: create response matrix
//==============================================================================

Int_t Unfold2D::Train()
{
	if(WIDEBIN) hTrainTrue = InitWideXY2DHisto("htraintrue", "Training Truth");
	else hTrainTrue= new TH2F ("htraintrue", "Training Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrainTrue->SetLineColor(kBlue);
	hTrainTrue->GetXaxis()->SetTitle(XvariableName);
	hTrainTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN ) hTrain= InitWideXY2DHisto("htrain", "Training Measured");
	else hTrain= new TH2F ("htrain", "Training Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrain->SetLineColor(kRed);
	hTrain->GetXaxis()->SetTitle(XvariableName);
	hTrain->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN ) hTrainFake= InitWideXY2DHisto("htrainfake", "Training Fakes");
	else hTrainFake= new TH2F ("htrainfake", "Training Fakes", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrainFake->SetLineColor(93);
	hTrainFake->GetXaxis()->SetTitle(XvariableName);
	hTrainFake->GetYaxis()->SetTitle(YvariableName);


	if(WIDEBIN) hXRcVsMc= InitWideXX2DHisto("hXRcVsMc", XvariableName+" - Training Measured vs Truth");
	else hXRcVsMc= new TH2F ("hXRcVsMc", XvariableName+" - Training Measured vs Truth", nmx, xlo, xhi, nmx, xlo, xhi);
	hXRcVsMc->SetLineColor(kRed);
	hXRcVsMc->GetXaxis()->SetTitle(Form("%s Mc",XvariableName.Data()));
	hXRcVsMc->GetYaxis()->SetTitle(Form("%s Rc",XvariableName.Data()));

	if(WIDEBIN) hYRcVsMc= InitWideYY2DHisto("hYRcVsMc", YvariableName+" - Training Measured vs Truth");
	else hYRcVsMc= new TH2F ("hYRcVsMc", YvariableName+" - Training Measured vs Truth", nmy, ylo, yhi, nmy, ylo, yhi);
	hYRcVsMc->SetLineColor(kRed);
	hYRcVsMc->GetXaxis()->SetTitle(Form("%s Mc",YvariableName.Data()));
	hYRcVsMc->GetYaxis()->SetTitle(Form("%s Rc",YvariableName.Data()));


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
	TString SFineBin="";
	if(!WIDEBIN) SFineBin = "_FineBin";
	TString SNoFake="";
	if(NoFakeOpt) SNoFake="NoFake";
	TString SNoLoss="";
	if(NoLossOpt) SNoLoss="NoLoss";
	TString SPrior="";
	if(CHANGEPRIOR) SPrior="mPrior";
	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}
	TString Spt="Pt02";
	if(OPT_RC02MC05) {
		Spt="PtRC02MC05";
	}
	TString ifilename = Form("ResponseMatrix%s_%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SFineBin.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data());	
	if(flagjetweight) {
		cout<<"INFO -------  Int_t Unfold2D::ReadResponseMatrix():  flagjetweight (args[14]) == 1. Use MB embedding for JP unfolding."<<endl;
		ifilename = Form("ResponseMatrix%s_%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),"MB",TranCharge.Data(),SExclude.Data(),SFineBin.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data());	
	}
	cout<<__PRETTY_FUNCTION__<<" Read in "<<ifilename<<endl;
	ftrain = new TFile(ifilename);
	hTrainTrue = (TH2F*)ftrain->Get("htraintrue");
	hTrain = (TH2F*)ftrain->Get("htrain");
	hTrainFake = (TH2F*)ftrain->Get("htrainfake");
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

	if(WIDEBIN) hMeas= InitWideXY2DHisto("hmeas", "Measured");
	else hMeas= new TH2F ("hmeas", "Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hMeas->SetLineColor(kRed);
	hMeas->GetXaxis()->SetTitle(XvariableName);
	hMeas->GetYaxis()->SetTitle(YvariableName);

	Float_t InputJet;
	Float_t InputPart;	

	Int_t  InputJetNtrk;
	Int_t  InputJetChargeNtrk;
	Int_t  InputLeadAreaNtrk;
	Int_t  InputSubAreaNtrk;
	Int_t  InputTranMaxNtrk;
	Int_t  InputTranMinNtrk;
	Float_t InputLeadAreaPtSum;
	Float_t InputSubAreaPtSum;
	Float_t InputTranMaxPtSum;
	Float_t InputTranMinPtSum;

	// some optional parameters for cuts and reweighting
	Float_t InputJetNF;		// Neutral fraction
	Int_t InputRunid;	

	// Tracks
	Float_t InputTrkLeadAreaPt[MAXARRAY];
	Float_t InputTrkSubAreaPt[MAXARRAY];
	Float_t InputTrkTranMaxPt[MAXARRAY];
	Float_t InputTrkTranMinPt[MAXARRAY];

	if(WIDEBIN) pfxMeas= new TProfile("pfxmeas", "Measure", WNbins, Wptbins);
	else pfxMeas= new TProfile ("pfxmeas", "Measure", nmx, xlo, xhi);
	pfxMeas->SetLineColor(hMeas->GetLineColor());
	pfxMeas->GetXaxis()->SetTitle(XvariableName);
	pfxMeas->GetYaxis()->SetTitle(YvariableName);


	TString ifilename;
	if(TrigName.Contains("MB")) {
		//ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_161209.root");
		ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_170418.root");
	}
	else {
		//ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_161209.root");
		//ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_170127.root");
		ifilename = TString("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_170418.root");
	}
	cout<<__PRETTY_FUNCTION__<<__func__<<" Reading in "<<ifilename<<endl;
	TFile *fin = new TFile(ifilename);
	if(!fin) {cout<<"Cannot open file in "<<__PRETTY_FUNCTION__<<endl; exit;}
	TTree *ttree = (TTree*)fin->Get("ResultTree");
	if(!ttree) {cout<<"Cannot find tree from input in "<<__PRETTY_FUNCTION__<<endl;exit;}

	ttree->SetBranchAddress("j1pt",&InputJet);
	ttree->SetBranchAddress("j1constntrk",&InputJetNtrk);
	ttree->SetBranchAddress("j1constchargentrk",&InputJetChargeNtrk);
	ttree->SetBranchAddress("LeadAreaNtrk",&InputLeadAreaNtrk);
	ttree->SetBranchAddress("SubAreaNtrk",&InputSubAreaNtrk);
	ttree->SetBranchAddress("TranMaxNtrk",&InputTranMaxNtrk);
	ttree->SetBranchAddress("TranMinNtrk",&InputTranMinNtrk);
	ttree->SetBranchAddress("LeadAreaPtSum",&InputLeadAreaPtSum);
	ttree->SetBranchAddress("SubLeadAreaPtSum",&InputSubAreaPtSum);
	ttree->SetBranchAddress("TranMaxPtSum",&InputTranMaxPtSum);
	ttree->SetBranchAddress("TranMinPtSum",&InputTranMinPtSum);
	ttree->SetBranchAddress("runid",&InputRunid);
	ttree->SetBranchAddress("j1neutralfrac",&InputJetNF);
	
	ttree->SetBranchAddress("TrkLeadAreaPt",InputTrkLeadAreaPt);
	ttree->SetBranchAddress("TrkSubAreaPt",InputTrkSubAreaPt);
	ttree->SetBranchAddress("TrkTranMaxPt",InputTrkTranMaxPt);
	ttree->SetBranchAddress("TrkTranMinPt",InputTrkTranMinPt);


	// Reweight JPs to have same Neutral Fraction distribution per jet pt bin as MB 
	ReWeightByNF *rwc;
	bool WriteNFfile = false;
	if( flagjetweight && ifilename.Contains("ppJP_") && ifilename.Contains("FullJet",TString::kIgnoreCase) ) {
		cout<<"Do per jet reweight";
		rwc = new ReWeightByNF();
		if(!rwc->ReadNFfile(Form("NFWeight_NoTofMatch_FullJet_TransCharged_160811P12id_R06_HadrCorr_VPDcut_170418.root"))) {
			//rwc->Init4Read("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209.root","~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209.root");
			rwc->Init4Read("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_170418.root","~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_170418.root");
			rwc->FillNFRatios();
			WriteNFfile = true;
		}
		cout<<"."<<endl;
	}

	// Scale JPs to have same *Y* distribution per jet pt bin as MB 
	ScaleByY *sby;
	bool DoScaleByY = !flagjetweight && flagscaley && ifilename.Contains("ppJP_") && ifilename.Contains("FullJet",TString::kIgnoreCase) ; 
	bool WriteYfile = false;
	if( DoScaleByY ) {
		cout<<"Scale "<<TrigName<<" "<<YvariableName<<" to be same mean value as MB"<<endl;
		sby = new ScaleByY(YvariableName);
		if(!sby->ReadYfile(Form("Scale%s_NoTofMatch_FullJet_TransCharged_160811P12id_R06_HadrCorr_VPDcut_170418.root",YvariableName.Data()))) {	// if scaleY root file exist, no need to read again
			//sby->Init4Read("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209.root","~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209.root");
			sby->Init4Read("~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_170418.root","~/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_170418.root");
			if(!sby->FillYRatios()) { cout<<"cannot fill ScaleByY:FillYRatios() in Unfold2D.cxx "<<endl; return 0;  }
			WriteYfile = true;
		}
		cout<<"."<<endl;
	}


	// Jet weight or event weight		for now (2016.12.21), this weight is used to get JPs jet neutral fraction to be same as MB trigger one
	TH1D *htmpx = (TH1D*)hMeas->ProjectionX("htmpx");		// htmpx is the InputJet distribution without reweighting; need to use same X-binning as hMeas. htmpx is used by ReWeightByNF class in order to insure the Nevts distribution to be same as before reweighting

	for (Int_t j=0;j<ttree->GetEntries() ; j++) {			// Caution!!! Need to update loop if needed
		ttree->GetEntry(j);

		if(InputJetNF<=0 || InputJetNF>=0.9) continue;			// this cut is also used by ReWeightByNF 	2017.01.04
		if(InputJet<=0) continue;					// if no jet, continue;
		if(InputRunid<=13048000) continue;				// begining of MB run looks different 

		// Jet weight or event weight		for now (2016.12.21), this weight is used to get JPs jet neutral fraction to be same as MB trigger one
		float xweight = 1;
		if( flagjetweight && ifilename.Contains("ppJP_") && ifilename.Contains("FullJet",TString::kIgnoreCase) ) {
			xweight = rwc->GetNFWeight(InputJet, InputJetNF);	
		}

		// brutally scale JP result to have same <Y> as MB 
		float scaley = 1;
		if( DoScaleByY ) {
			scaley = sby->GetYScale(InputJet, InputJetNF);	
		}


		if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
		{
			if(InputTranMaxNtrk>InputTranMinNtrk) InputPart = static_cast<float>(InputTranMaxNtrk);
			else InputPart = static_cast<float>(InputTranMinNtrk);
		}
		else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
		{
			if(InputTranMaxNtrk<InputTranMinNtrk) InputPart = static_cast<float>(InputTranMaxNtrk);
			else InputPart = static_cast<float>(InputTranMinNtrk);
		}
		else if(YvariableName.Contains("TranDiffNtrk",TString::kIgnoreCase)) 
		{
			InputPart = fabs(InputTranMaxNtrk-InputTranMinNtrk);
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
		else if(YvariableName.Contains("TranMaxPtAveEventWise",TString::kIgnoreCase))
		{
			float InputTranMaxPtAve = 1.*InputTranMaxPtSum/InputTranMaxNtrk;
			float InputTranMinPtAve = 1.*InputTranMinPtSum/InputTranMinNtrk;
			if(InputTranMaxPtAve>InputTranMinPtAve) InputPart = InputTranMaxPtAve;
			else InputPart = InputTranMinPtAve;
		}
		else if(YvariableName.Contains("TranMinPtAveEventWise",TString::kIgnoreCase)) 
		{
			float InputTranMaxPtAve = 1.*InputTranMaxPtSum/InputTranMaxNtrk;
			float InputTranMinPtAve = 1.*InputTranMinPtSum/InputTranMinNtrk;
			if(InputTranMaxPtAve<InputTranMinPtAve) InputPart = InputTranMaxPtAve;
			else InputPart = InputTranMinPtAve;
		}
		else if(YvariableName.Contains("TranPtAveEventWise",TString::kIgnoreCase))
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
		else if(YvariableName.Contains("LeadAreaPtAveEventWise",TString::kIgnoreCase)) 
		{
			InputPart = InputLeadAreaPtSum/InputLeadAreaNtrk;
		}
		else if(YvariableName.Contains("SubAreaPtAveEventWise",TString::kIgnoreCase)) 
		{
			InputPart = InputSubAreaPtSum/InputSubAreaNtrk;
		}
		else if(YvariableName.Contains("LeadPtAve",TString::kIgnoreCase)&&!(YvariableName.Contains("EventWise",TString::kIgnoreCase))) 		// Particle-wise variable, need to loop over Mc and Rc particles
		{

			if( DoScaleByY )  InputPart = InputPart*scaley;

			for(int it = 0; it<InputLeadAreaNtrk; it++) 
			{
				InputPart = InputTrkLeadAreaPt[it];
				hMeas->Fill(InputJet, InputPart, xweight);	
				htmpx->Fill(InputJet);		// without weighting; will later be used by ReWeightByNF for adjusting to have same Nevts distribution before weighting
				pfxMeas->Fill(InputJet, InputPart, xweight);	
			}
		}
		else if(YvariableName.Contains("SubPtAve",TString::kIgnoreCase)&&!(YvariableName.Contains("EventWise",TString::kIgnoreCase))) 		// Particle-wise variable, need to loop over Mc and Rc particles
		{

			if( DoScaleByY )  InputPart = InputPart*scaley;

			for(int it = 0; it<InputSubAreaNtrk; it++) 
			{
				InputPart = InputTrkSubAreaPt[it];
				hMeas->Fill(InputJet, InputPart, xweight);	
				htmpx->Fill(InputJet);		// without weighting; will later be used by ReWeightByNF for adjusting to have same Nevts distribution before weighting
				pfxMeas->Fill(InputJet, InputPart, xweight);	
			}
		}
		else if(YvariableName.Contains("TranPtAve",TString::kIgnoreCase)&&!(YvariableName.Contains("EventWise",TString::kIgnoreCase))) 		// Particle-wise variable, need to loop over Mc and Rc particles
		{

			if( DoScaleByY )  InputPart = InputPart*scaley;

			for(int it = 0; it<InputTranMaxNtrk; it++) 
			{
				InputPart = InputTrkTranMaxPt[it];
				hMeas->Fill(InputJet, InputPart, xweight);	
				htmpx->Fill(InputJet);		// without weighting; will later be used by ReWeightByNF for adjusting to have same Nevts distribution before weighting
				pfxMeas->Fill(InputJet, InputPart, xweight);	
			}
			for(int it = 0; it<InputTranMinNtrk; it++) 
			{
				InputPart = InputTrkTranMinPt[it];
				hMeas->Fill(InputJet, InputPart, xweight);	
				htmpx->Fill(InputJet);		// without weighting; will later be used by ReWeightByNF for adjusting to have same Nevts distribution before weighting
				pfxMeas->Fill(InputJet, InputPart, xweight);	
			}
		}
		else if(YvariableName.Contains("JetNtrk",TString::kIgnoreCase))
		{
			InputPart = InputJetNtrk;
		}
		else if(YvariableName.Contains("JetChargeNtrk",TString::kIgnoreCase))
		{
			InputPart = InputJetChargeNtrk;
		}
		else // default: TranNtrk
		{
			InputPart = 0.5*(InputTranMaxNtrk+InputTranMinNtrk);
			cout<<"WARNING!!!!!!!!!! in "<<__PRETTY_FUNCTION__<<" "<<YvariableName<<" is not valid YvariableName name. Will use TranNtrk instead!"<<endl;
			YvariableName = "TranNtrk";
		}



		// Fill histogram
		if( !(YvariableName.Contains("PtAve",TString::kIgnoreCase)&&(!YvariableName.Contains("EventWise",TString::kIgnoreCase))) ) {

			if( DoScaleByY )  InputPart = InputPart*scaley;

			hMeas->Fill(InputJet, InputPart, xweight);	
			htmpx->Fill(InputJet);		// without weighting; will later be used by ReWeightByNF for adjusting to have same Nevts dist. before weighting
			pfxMeas->Fill(InputJet, InputPart, xweight);	
		}

	}// End Loop of events
	fin->Close();

	// Adjust jet pt Nevts distribution if ReWeightByNF was used
	if( flagjetweight && ifilename.Contains("ppJP_") && ifilename.Contains("FullJet",TString::kIgnoreCase) ) {
		rwc->AdjustNevts(hMeas, htmpx);
	}

	hMeas->Sumw2();

	hMeasX= ProjectionX (hMeas, "hMeasX", "Measured X");
	hMeasY= ProjectionY (hMeas, "hMeasY", "Measured Y");

	if(WriteNFfile && rwc!=NULL) {
		rwc->WriteNF("NFWeight_NoTofMatch_FullJet_TransCharged_160811P12id_R06_HadrCorr_VPDcut_170418.root");
	}

	if(WriteYfile && sby!=NULL) {
		sby->WriteY(Form("Scale%s_NoTofMatch_FullJet_TransCharged_160811P12id_R06_HadrCorr_VPDcut_170418.root",YvariableName.Data()));
	}

	return 1;
}


//==============================================================================
// Read hMeas histogram from input root file to unfold
//==============================================================================
Int_t Unfold2D::ReadMeasHist4Unfold() {

	//TString ifilename = TString("~/Scratch/pp200Y12_jetunderlying/leadjetpthist4NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_161209_NoLossCorr.root");
	TString ifilename = TString("~/Scratch/pp200Y12_jetunderlying/leadjetpthist4NoTofMatch_FullJet_Trans")+TranCharge+TString("_MatchTrig_pp")+TrigName+TString("_160811P12id_R06_HadrCorr_170418_NoLossCorr.root");
	cout<<__PRETTY_FUNCTION__<<" Read in "<<ifilename<<endl;
	TFile *fin = new TFile(ifilename);
	if(!fin) {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<" Cannot open "<< ifilename<<endl; 
		exit; 
	}
	TString ihistname = "h"+ToLower(YvariableName)+"vsleadjetpt";			// Default is Tran Multiplicity 
	TH2F *htmp = (TH2F*)fin->Get(ihistname);
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

	if(WIDEBIN) hTrainTrue = InitWideXY2DHisto("htraintrue", "Training Truth");
	else hTrainTrue= new TH2F ("htraintrue", "Training Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrainTrue->SetLineColor(kBlue);
	hTrainTrue->GetXaxis()->SetTitle(XvariableName);
	hTrainTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrain= InitWideXY2DHisto("htrain", "Training Measured");
	else hTrain= new TH2F ("htrain", "Training Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrain->SetLineColor(kRed);
	hTrain->GetXaxis()->SetTitle(XvariableName);
	hTrain->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hTrainFake= InitWideXY2DHisto("htrainfake", "Training Fakes");
	else hTrainFake= new TH2F ("htrainfake", "Training Fakes", nmx, xlo, xhi, nmy, ylo, yhi);
	hTrainFake->SetLineColor(93);
	hTrainFake->GetXaxis()->SetTitle(XvariableName);
	hTrainFake->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hXRcVsMc= InitWideXX2DHisto("hXRcVsMc", XvariableName+" - Training Measured vs Truth");
	else hXRcVsMc= new TH2F ("hXRcVsMc", XvariableName+" - Training Measured vs Truth", nmx, xlo, xhi, nmx, xlo, xhi);
	hXRcVsMc->SetLineColor(kRed);
	hXRcVsMc->GetXaxis()->SetTitle(Form("%s Mc",XvariableName.Data()));
	hXRcVsMc->GetYaxis()->SetTitle(Form("%s Rc",XvariableName.Data()));

	if(WIDEBIN) hYRcVsMc= InitWideYY2DHisto("hYRcVsMc", YvariableName+" - Training Measured vs Truth");
	else hYRcVsMc= new TH2F ("hYRcVsMc", YvariableName+" - Training Measured vs Truth", nmy, ylo, yhi, nmy, ylo, yhi);
	hYRcVsMc->SetLineColor(kRed);
	hYRcVsMc->GetXaxis()->SetTitle(Form("%s Mc",YvariableName.Data()));
	hYRcVsMc->GetYaxis()->SetTitle(Form("%s Rc",YvariableName.Data()));

	cout<<"Set response"<<endl;
	response = new RooUnfoldResponse ("response", "");
	response->Setup (hTrain, hTrainTrue);

	//int Nevents[NUMBEROFPT] = {600000,200000,300000,150000,150000,150000,80000,50000,40000,38000,12000};
	//int Nevents[NUMBEROFPT] = {1050000,300000,300000,150000,150000,150000,80000,50000,40000,38000,12000};
	int Nevents[NUMBEROFPT] = {1550000,450000,450000,220000,220000,220000,120000,75000,60000,58000,17000};
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

	if(WIDEBIN) hTrue= InitWideXY2DHisto("htrue", "Test Truth");	
	else hTrue= new TH2F ("htrue", "Test Truth", ntx, xlo, xhi, nty, ylo, yhi);
	hTrue->SetLineColor(kBlue);
	hTrue->GetXaxis()->SetTitle(XvariableName);
	hTrue->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hMeas= InitWideXY2DHisto("hmeas", "Test Measured");
	else hMeas= new TH2F ("hmeas", "Test Measured", nmx, xlo, xhi, nmy, ylo, yhi);
	hMeas->SetLineColor(kRed);
	hMeas->GetXaxis()->SetTitle(XvariableName);
	hMeas->GetYaxis()->SetTitle(YvariableName);

	if(WIDEBIN) hFake= InitWideXY2DHisto("hfake", "Test Fake");
	else hFake= new TH2F ("hfake", "Test Fake", nmx, xlo, xhi, nmy, ylo, yhi);
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

	Float_t tRcjneutralfrac;
	Int_t  tRcTranMaxNtrk;
	Int_t  tRcTranMinNtrk;
	Int_t  tRcLeadAreaNtrk;
	Int_t  tRcSubAreaNtrk;
	Float_t tRcLeadAreaPtSum;
	Float_t tRcSubAreaPtSum;
	Float_t tRcTranMaxPtSum;
	Float_t tRcTranMinPtSum;


	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}


	double totalxsec = 0;
	double totalnevents = 0;
	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename;
		if(TrigName.Contains("JP")&&(flagjetweight||flagscaley)) {
		// if flagjetweight==true, JP hist is weighted by MB NF distribution. Then we shall use MB embeding data (no trigger effect). This means we don't unfold trigger effect on jet pt distribution, but the neutral jet trigger effect should be already corrected by weighting procedure. However, if we don't do weighting by MB Neutral Fraction dist, we could also use JP embedding to unfold all trigger effects together
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_McPt02%s.root",  "MB",   PTBINS[i],TranCharge.Data(),STpcSys.Data());
		}
		else if(ExcludeOpt&&TrigName.Contains("JP")) {		
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_excluded_McPt02%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data(),STpcSys.Data());
		}
		else {
			ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTrans%s_McPt02%s.root",TrigName.Data(),PTBINS[i],TranCharge.Data(),STpcSys.Data());
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
		ttree->SetBranchAddress("Rcj1neutralfrac",&tRcjneutralfrac);
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
		//tweight = XSEC[i]/(ttree->GetEntries());				//	in case we want to do half for train half for test, weight need to be same
		//cout<<"test weight = "<<tweight<<endl;				
		totalxsec+=XSEC[i];
		totalnevents+=(ttree->GetEntries()-Nevents[i]);

		cout<<"From "<<Nevents[i]<<" to "<<ttree->GetEntries()<<endl;
		for (Int_t j=Nevents[i];j<ttree->GetEntries() ; j++) {			// Caution!!! Need to update loop if needed
		//for (Int_t j=1;j<ttree->GetEntries() ; j+=2) {			// in case we want to do half for train half for test, *random* select (even and odd for train, test)	
			ttree->GetEntry(j);

			//as it is not in use in Train filling, we also comment it out here..  if(!tflagIsTrigger) continue;			// uncomment for no correction for trig now..
			tflag = -999;
			if(TrigName.Contains("JP")&&!(flagjetweight||flagscaley)) {
				if(tflagtrigmatch && (tflagMatch2Lead||tflagMatch2Sub) ) tflag = 1; 
				else if(tMcJet>0 && (tRcJet<=0 || !tflagtrigmatch) ) tflag = 0;		// non-zero tMc, zero tRc
				else if((tRcJet>0&&tflagtrigmatch) && tMcJet<=0) tflag = -1;		// non-zero tRc, zero tMc
			}
			else {
				//// simply case. check with Train Fill to see what is in use
				//if(tflagMatch2Lead||tflagMatch2Sub) tflag = 1; 							
				//else if(tMcJet>=0 && tRcJet<=0) tflag = 0;		// non-zero Mc, zero Rc
				//else if(tRcJet>=0 && tMcJet<=0) tflag = -1;		// non-zero Rc, zero Mc
				
				//check with Train Fill to see what is in use
				if((tflagMatch2Lead||tflagMatch2Sub) && tRcJet>0 && tRcjneutralfrac<0.9 && tRcjneutralfrac>0) tflag = 1; 
				else if(tMcJet>=0 && (tRcJet<=0 || !(tRcjneutralfrac<0.9 && tRcjneutralfrac>0))) tflag = 0;		// non-zero Mc, zero Rc
				else if(tRcJet>0 && tRcjneutralfrac<0.9 && tRcjneutralfrac>0 && (tMcJet<=0 || !(tflagMatch2Lead||tflagMatch2Sub)) ) tflag = -1;		// fake events: non-zero Rc, zero Mc; or non-zero Rc, non-zero Mc but Rc not matched to Mc leading/sublead
			}


			// Read tMcPart and tRcPart
			if(OPT_RC02MC05) {              // MC pt>0.5 only
				cout<<"ERR!! in "<<__PRETTY_FUNCTION__<<" OPT_RC02MC05=true is only implemented for Real data, not Test data. Will reset OPT_RC02MC05 = false instead"<<endl;
				OPT_RC02MC05 = false;
			}
			else {	
				if(YvariableName.Contains("TranMaxNtrk",TString::kIgnoreCase))
				{
					if(tMcTranMaxNtrk>tMcTranMinNtrk) tMcPart = static_cast<float>(tMcTranMaxNtrk);
					else tMcPart = static_cast<float>(tMcTranMinNtrk);
					if(tRcTranMaxNtrk>tRcTranMinNtrk) tRcPart = static_cast<float>(tRcTranMaxNtrk);
					else tRcPart = static_cast<float>(tRcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranMinNtrk",TString::kIgnoreCase)) 
				{
					if(tMcTranMaxNtrk<tMcTranMinNtrk) tMcPart = static_cast<float>(tMcTranMaxNtrk);
					else tMcPart = static_cast<float>(tMcTranMinNtrk);
					if(tRcTranMaxNtrk<tRcTranMinNtrk) tRcPart = static_cast<float>(tRcTranMaxNtrk);
					else tRcPart = static_cast<float>(tRcTranMinNtrk);
				}
				else if(YvariableName.Contains("TranDiffNtrk",TString::kIgnoreCase)) 
				{
					tMcPart = fabs(tMcTranMaxNtrk-tMcTranMinNtrk);
					tRcPart = fabs(tRcTranMaxNtrk-tRcTranMinNtrk);
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
				else if(YvariableName.Contains("LeadAreaPtAve",TString::kIgnoreCase)||YvariableName.Contains("LeadPtAve",TString::kIgnoreCase)) 
				{
					tMcPart = tMcLeadAreaPtSum/tMcLeadAreaNtrk;
					tRcPart = tRcLeadAreaPtSum/tRcLeadAreaNtrk;
				}
				else if(YvariableName.Contains("SubAreaPtAve",TString::kIgnoreCase)||YvariableName.Contains("SubPtAve",TString::kIgnoreCase)) 
				{
					tMcPart = tMcSubAreaPtSum/tMcSubAreaNtrk;
					tRcPart = tRcSubAreaPtSum/tRcSubAreaNtrk;
				}
				else // default: TranNtrk
				{
					tMcPart = 0.5*(tMcTranMaxNtrk+tMcTranMinNtrk);
					tRcPart = 0.5*(tRcTranMaxNtrk+tRcTranMinNtrk);
					cout<<"WARNING!!!!!!!!!! in "<<__PRETTY_FUNCTION__<<" "<<YvariableName<<" is not valid YvariableName name. Will use TranNtrk instead!"<<endl;
					YvariableName = "TranNtrk";
				}
			}// end if else OPT_RC02MC05


			// Fill histograms according to flag
			if(tflag==1) {			// one-to-one tMc to tRc matched
				hTrue->Fill(tMcJet, tMcPart, tweight);	
				hMeas->Fill(tRcJet, tRcPart, tweight);	
				pfxMeas->Fill(tRcJet, tRcPart, tweight);	
			}
			else if(tflag==0) {		// true tMc, no tRc
				hTrue->Fill(tMcJet, tMcPart, tweight);
			}
			else if(tflag==-1 && !NoFakeOpt){				// Fake tRc, no tMc
				hMeas->Fill(tRcJet,tRcPart, tweight);
				hFake->Fill(tRcJet,tRcPart, tweight);
				pfxMeas->Fill(tRcJet,tRcPart, tweight);
			}

		}// End Loop of events
		fin->Close();
	}// End Loop of pT

	//Scale the histogram to simulate statistic of how many events. Note: sumw2() is called after wards, otherwise it would serve the purpose of simulate the statistic . 
	//double toscale = totalnevents/totalxsec;
	double toscale = 1; // 100e6/totalxsec;		// say if we have 100M MB data. how does the real unfolding probably looks like
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
	TString SFineBin="";
	if(!WIDEBIN) SFineBin = "_FineBin";
	TString SNoFake="";
	if(NoFakeOpt) SNoFake="NoFake";
	TString SNoLoss="";
	if(NoLossOpt) SNoLoss="NoLoss";
	TString SPrior="";
	if(CHANGEPRIOR) SPrior="mPrior";
	TString SMethod="";
	if(method==3) SMethod="_BinByBin";
	else if(method==2) SMethod="_SVD";
	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}
	ftout = new TFile(Form("TrainTest%s_%s%s%s%s%s%s%s_Mc%s%sPt02%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SFineBin.Data(),SMethod.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),STpcSys.Data()),"RECREATE");
	if(ftout) cout<<"Write to "<<ftout->GetName()<<endl;
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	response->Write();
	hXRcVsMc->Write();
	hYRcVsMc->Write();

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
	TString SFineBin="";
	if(!WIDEBIN) SFineBin = "_FineBin";
	TString SIter="";
	TString SNoFake="";
	if(NoFakeOpt) SNoFake="NoFake";
	TString SNoLoss="";
	if(NoLossOpt) SNoLoss="NoLoss";
	TString SPrior="";
        if(CHANGEPRIOR) SPrior="mPrior";
	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}
	TString Spt="Pt02";
	if(OPT_RC02MC05) {
		Spt="PtRC02MC05";
	}
	TString ifilename = Form("ResponseMatrix%s_%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SFineBin.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data());
	if(flagjetweight) {
		cout<<"INFO -------  Int_t Unfold2D::WriteTrain():  flagjetweight (args[14]) == 1. Use MB embedding for JP unfolding, which will be written"<<endl;
		ifilename = Form("ResponseMatrix%s_%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),"MB",TranCharge.Data(),SExclude.Data(),SFineBin.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data());	
	}


	cout<<"Write to "<<ifilename<<endl;
	ftout = new TFile(ifilename,"RECREATE");
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	response->Write();
	hXRcVsMc->Write();
	hYRcVsMc->Write();
	if(pfxTrain) pfxTrain->Write();

	cout<<"WriteTrain() done"<<endl;

	return 1;
}


//========================================================================================================================================
// Write: record histograms for unfolding 	This is an intermediate step for debug
// No unfolding result. Can be called after Fill4Unfold. Only write out hmeas, htrain, htraintrue ..
//========================================================================================================================================

Int_t Unfold2D::WriteHist4Unfold() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	TString SNFWeight="";
	if(flagjetweight) SNFWeight = "_NFweight";
	TString SYScale="";
	if(!flagjetweight && flagscaley) SYScale = "_Yscale";	
	TString SFineBin="";
	if(!WIDEBIN) SFineBin = "_FineBin";
	TString SNoFake="";
	if(NoFakeOpt) SNoFake="NoFake";
	TString SNoLoss="";
	if(NoLossOpt) SNoLoss="NoLoss";
        TString SPrior="";
        if(CHANGEPRIOR) SPrior="mPrior";
	TString SIter="";
	if( method==1 && (regparm!=4) ) {
		SIter = Form("_Baye%d",regparm);
	}
	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}
	TString Spt="Pt02";
	if(OPT_RC02MC05) {
		Spt="PtRC02MC05";
	}
	ftout = new TFile(Form("Hist%s_%s%s%s%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SIter.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data()),"RECREATE");
	hMeas->Write();
	if(hTrainTrue) hTrainTrue->Write();
	if(hTrain)hTrain->Write();
	if(hTrainFake) hTrainFake->Write();
	if(hXRcVsMc) hXRcVsMc->Write();
	if(hYRcVsMc) hYRcVsMc->Write();
	if(response) response->Write();
	if(pfxTrain) pfxTrain->Write();
	if(pfxTrainTrue) pfxTrainTrue->Write();
	if(pfxMeas) pfxMeas->Write();

	return 1;
}


//==============================================================================
// Write: record unfolding result
//==============================================================================

Int_t Unfold2D::WriteUnfoldResult() 
{
	TString SExclude="";
	if(ExcludeOpt) SExclude = "_excluded";
	TString SNFWeight="";
	if(flagjetweight) SNFWeight = "_NFweight";
	TString SYScale="";
	if(!flagjetweight && flagscaley) SYScale = "_Yscale";	
	TString SFineBin="";
	if(!WIDEBIN) SFineBin = "_FineBin";
	TString SNoFake="";
	if(NoFakeOpt) SNoFake="NoFake";
	TString SNoLoss="";
	if(NoLossOpt) SNoLoss="NoLoss";
        TString SPrior="";
        if(CHANGEPRIOR) SPrior="mPrior";
	TString SIter="";
	if( method==1 && (regparm!=4) ) {
		SIter = Form("_Baye%d",regparm);
	}
	else if(method==3) {
		SIter = "_BinByBin";
	}
	else if(method==2) {
		SIter = "_SVD";
	}
	TString STpcSys="";
	if(OPT_TPCSYS) {
		if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
		else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
		else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
		else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
	}
	TString Spt="Pt02";
	if(OPT_RC02MC05) {
		Spt="PtRC02MC05";
	}
	//ftout = new TFile(Form("Unfolding%s_%s%s%s%s%s%s%s%s%s_Mc%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SIter.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data()),"RECREATE");
	ftout = new TFile(Form("Unfolding%s_%s%s%s%s%s%s%s%s%s_12JetBinv2_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SIter.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data()),"RECREATE");	// 12JetBinv2
	//ftout = new TFile(Form("Unfolding%s_%s%s%s%s%s%s%s%s%s_12JetBinv2_pTBinv2_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SIter.Data(),SPrior.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data()),"RECREATE");	// 12JetBinv2 and can profile to 0.5 instead of 0.48 GeV/c for tranptave pt binnning
	hTrainTrue->Write();
	hTrain->Write();
	hTrainFake->Write();
	hXRcVsMc->Write();
	hYRcVsMc->Write();
	response->Write();
	hMeas->Write();
	if(hReco) hReco->Write();
	if(pfxTrain) pfxTrain->Write();
	if(pfxTrainTrue) pfxTrainTrue->Write();
	if(pfxMeas) pfxMeas->Write();
	if(hpriorratio) hpriorratio->Write();

	return 1;
}

//==============================================================================
// Apply event weight to change prior for unfolding
// use the ratio of measured jet pt / pythia MC jet pt
//==============================================================================

float Unfold2D::Weight2ChargePrior(float jetpt) 
{
	int pbin = -1;
	if(!hpriorratio)  {
		TString SExclude="";
		if(ExcludeOpt) SExclude = "_excluded";
		TString SNFWeight="";
		if(flagjetweight) SNFWeight = "_NFweight";
		TString SYScale="";
		if(!flagjetweight && flagscaley) SYScale = "_Yscale";	
		TString SFineBin="";
		if(!WIDEBIN) SFineBin = "_FineBin";
		TString SNoFake="";
		if(NoFakeOpt) SNoFake="NoFake";
		TString SNoLoss="";
		if(NoLossOpt) SNoLoss="NoLoss";
		TString STpcSys="";
		if(OPT_TPCSYS) {
			if(OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrPlusAbs0.05";
			else if(OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErrPlus0.05";
			else if(!OPT_TPCSYSPM && OPT_TPCSYSABS) STpcSys = "_TpcErrAbs0.05";
			else if(!OPT_TPCSYSPM && !OPT_TPCSYSABS) STpcSys = "_TpcErr0.05";
		}
		TString Spt="Pt02";
		if(OPT_RC02MC05) {
			Spt="PtRC02MC05";
		}
		TFile *fprior = new TFile(Form("Unfolding%s_%s%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data()));
		if(!fprior->IsOpen()) {
			fprior = new TFile(Form("Hist%s_%s%s%s%s%s%s%s_Mc%s%s%s%s.root",inputname.Data(),YvariableName.Data(),TrigName.Data(),TranCharge.Data(),SExclude.Data(),SNFWeight.Data(),SYScale.Data(),SFineBin.Data(),SNoFake.Data(),SNoLoss.Data(),Spt.Data(),STpcSys.Data()));
		}
		if(!fprior->IsOpen()) {
			cout<<"WARNING!! void Unfold2D::Weight2ChargePrior(float jetpt): Can not find root file to read "<<endl;
			return 1;
		}
		cout<<"Read "<<fprior->GetName()<<" to Change Prior."<<endl;
		TH2F *phtraintrue = (TH2F*)fprior->Get("htraintrue");
		TH2F *phtrain = (TH2F*)fprior->Get("htrain");
		TH2F *phmeas = (TH2F*)fprior->Get("hmeas");
		TH1F *phtraintrue_px = (TH1F*)phtraintrue->ProjectionX("phtraintrue_px");
		TH1F *phtrain_px = (TH1F*)phtrain->ProjectionX("phtrain_px");
		TH1F *phmeas_px = (TH1F*)phmeas->ProjectionX("phmeas_px");
		hpriorratio = (TH1F*) phmeas_px->Clone("hpriorratio");
		//hpriorratio->Divide(phtraintrue_px);
		hpriorratio->Divide(phtrain_px);
	}
	pbin = hpriorratio->FindBin(jetpt);

	if(pbin<=0 || pbin>hpriorratio->GetNbinsX()) { cout<<"WARNING!! void Unfold2D::Weight2ChargePrior(float jetpt): "<<jetpt<<" not in hpriorratio x range. "<<endl; return 1;  }
	
	float ratio = hpriorratio->GetBinContent(pbin);

	return ratio;
}


//==============================================================================
// Apply event weight to change prior for unfolding
// use a user defined function to generate weigth
//==============================================================================

float Unfold2D::Weight2ChargePrior(TF1 *f1, float jetpt) 
{
	float ratio = 1; 
	if(!f1) {
		f1 = new TF1("f1","exp(-[0]*x)",Wptbins[0],Wptbins[WNbins]);
		f1->SetParameter(0,0.1);
		ratio = f1->Eval(jetpt);
		f1->Delete();
	}
	else {
		ratio = f1->Eval(jetpt);
	}

	return ratio;
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
	hXRcVsMc = hYRcVsMc = 0;

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
	const TH2* h2= static_cast<const TH2*>(h);
	TH1D* h1= h2->ProjectionX (name, 1, h->GetNbinsY(), opt);
	if (title) h1->SetTitle (title);
	return h1;
}

TH1D* Unfold2D::ProjectionY (const TH1* h, const char* name, const char* title, Option_t* opt)
{
	const TH2* h2= static_cast<const TH2*>(h);
	TH1D* h1= h2->ProjectionY (name, 1, h->GetNbinsX(), opt);
	if (title) h1->SetTitle (title);
	return h1;
}

TProfile* Unfold2D::ProfileX (const TH1* h, const char* name, const char* title, Option_t* opt)
{
	const TH2* h2= static_cast<const TH2*>(h);
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

	if(FlagDefaultCalled) SetDefaultParms();		// If SetXYname called, need to update default parameters using YvariableName if we are using default values
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

void Unfold2D::SetChangePrior(Int_t change) {
	CHANGEPRIOR = change;
	cout<<"Set CHANGEPRIOR = "<<CHANGEPRIOR<<endl;
}

void Unfold2D::SetTpcSys(int opt_tpcsys, int opt_tpcsyspm, int opt_tpcsysabs) {
	//cout<<"opt_tpcsys = "<<opt_tpcsys<<" opt_tpcsyspm = "<<opt_tpcsyspm<<" opt_tpcsysabs = "<<opt_tpcsysabs<<endl;
	if(opt_tpcsys>0) {
		OPT_TPCSYS = true;
	}
	else {
		OPT_TPCSYS = false;
	}
	if(!OPT_TPCSYS) return;
	cout<<"Set OPT_TPCSYS = "<<OPT_TPCSYS<<" to apply "<<endl;
	if(opt_tpcsyspm>0) {
		OPT_TPCSYSPM = true;
		cout<<" + ";
	}
	else {
		OPT_TPCSYSPM = false;
		cout<<" - ";
	}
	if(opt_tpcsysabs>0) {
		OPT_TPCSYSABS = true;
		cout<<" absolute ";
	}
	else {
		OPT_TPCSYSABS = false;
		cout<<" relative ";
	}
	cout<<" 5% TPC tracking sys. err"<<endl;

}

void Unfold2D::SetRc02Mc05(int opt_rc02mc05) {
	if(opt_rc02mc05>0) {
		OPT_RC02MC05 = true;
		cout<<"OPT_RC02MC05 = true. Only unfold to pT>0.5GeV/c."<<endl;
	}
	else OPT_RC02MC05 = false;
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


TH2F* Unfold2D::CorrelationHist (const TMatrixD& cov,
		const char* name, const char* title,
		Double_t lo, Double_t hi)
{
	Int_t nb= cov.GetNrows();
	TH2F* h= new TH2F (name, title, nb, lo, hi, nb, lo, hi);
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
TH2F* Unfold2D::RebinAs(TH2F *old, TH2F *model)		// rebin old as the same binining as model
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

TH2F* Unfold2D::RebinX2DHisto(TH2F *old)	
{
	//const int WNbins = 15;
        //double Wptbins[WNbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};

	int Ny = old->GetNbinsY();
	double *ybins = new double[Ny+1];
	for(int ibin = 0; ibin<Ny; ibin++) {
		*(ybins+ibin) = old->GetYaxis()->GetBinLowEdge(ibin+1);
	}
	*(ybins+Ny) = old->GetYaxis()->GetBinLowEdge(Ny) + old->GetYaxis()->GetBinWidth(Ny);

	TH2F *hnew = Rebin2DHisto(old, WNbins, Wptbins, Ny, ybins);
	return hnew;
}


TH2F* Unfold2D::Rebin2DHisto(TH2F *old, int Nx, double *xbins, int Ny, double *ybins)
{
	TH2F *h = new TH2F(Form("Rebin%s",old->GetName()),old->GetTitle(),Nx,xbins,Ny,ybins);	
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

TH2F* Unfold2D::InitWideX2DHisto(TString name, TString title, int Ny, double ylo, double yhi) 
{

	if(Ny<=0) {cout<<"Err: "<<__PRETTY_FUNCTION__<<", Ny = "<<Ny<<"<=0"<<endl;return NULL;}
	double *ybins = new double[Ny+1];
	for(int ibin = 0; ibin<Ny; ibin++) {
		*(ybins+ibin) = ylo + ibin*(yhi-ylo)/Ny;	
	}
	*(ybins+Ny) = yhi;

	TH2F* h = new TH2F(name, title, WNbins, Wptbins, Ny,ybins);
	return h;
}

TH2F* Unfold2D::InitWideX2DHisto(TString name, TString title, int Ny, double *ybins) 
{

	if(Ny<=0) {cout<<"Err: "<<__PRETTY_FUNCTION__<<", Ny = "<<Ny<<"<=0"<<endl;return NULL;}

	TH2F* h = new TH2F(name, title, WNbins, Wptbins, Ny,ybins);
	return h;
}

TH2F* Unfold2D::InitWideXY2DHisto(TString name, TString title) 
{
	TH2F *h;

	cout<<"Init WideBin for "<<name<<" "<<title<<" with ";
	if(YvariableName.Contains("Ntrk")) 
	{
		const int Nltrkbins = 18;
		double ltrkbins[Nltrkbins+1] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,17.5,21.5,39.5};
		h = InitWideX2DHisto(name, title, Nltrkbins, ltrkbins);
		cout<<Nltrkbins<<" bins"<<endl;
	}

	else if(YvariableName.Contains("PtAve")) 
	{
		//const int Nptavebins = 178;
		//double ptavebins[Nptavebins+1];
		//if MC unfolded to pT==0
		//ptavebins[0] = 0;
		//for(int i = 1; i<132; i++) {
		//	ptavebins[i] = 0.2+(i-1)*0.01;
		//}
		//for(int i = 132; i<142; i++) {
		//	ptavebins[i] = 1.5+(i-131)*0.05;
		//}
		//for(int i = 142; i<162; i++) {
		//	ptavebins[i] = 2+(i-141)*0.1;
		//}
		//for(int i = 162; i<172; i++) {
		//	ptavebins[i] = 4+(i-161)*0.2;
		//}
		//for(int i = 172; i<179; i++) {
		//	ptavebins[i] = 6+(i-171)*2;
		//}


		//// if Mc unfolded to pT->0.2
		//const int Nptavebins = 198;
		//double ptavebins[Nptavebins+1];
		//for(int i = 0; i<198; i++) {
		//	ptavebins[i] = 0.2+i*0.1;
		//}
		//ptavebins[198] = 20;


		// if Mc unfolded to pT->0.2
		const int Nptavebins = 50;
		double ptavebins[Nptavebins+1];
		for(int i = 0; i<26; i++) {
			ptavebins[i] = 0.2+i*0.04;
		}
		for(int i = 26; i<34; i++) {
			ptavebins[i] = 1.2+(i-25)*0.1;
		}
		for(int i = 34; i<44; i++) {
			ptavebins[i] = 2+(i-33)*0.2;
		}
		for(int i = 44; i<48; i++) {
			ptavebins[i] = 4+(i-43)*0.5;
		}
		for(int i = 48; i<50; i++) {
			ptavebins[i] = 6+(i-47)*2;
		}
		ptavebins[50] = 20;
	
		//// if Mc unfolded to pT->0.2
		//// AND if we want to profile to pT>0.5 later on
		//const int Nptavebins = 50;
		//double ptavebins[Nptavebins+1];
		//for(int i = 0; i<7; i++) {
		//	ptavebins[i] = 0.2+i*0.04;
		//}
		//ptavebins[7] = 0.5;
		//for(int i = 8; i<26; i++) {
		//	ptavebins[i] = 0.2+i*0.04;
		//}
		//for(int i = 26; i<34; i++) {
		//	ptavebins[i] = 1.2+(i-25)*0.1;
		//}
		//for(int i = 34; i<44; i++) {
		//	ptavebins[i] = 2+(i-33)*0.2;
		//}
		//for(int i = 44; i<48; i++) {
		//	ptavebins[i] = 4+(i-43)*0.5;
		//}
		//for(int i = 48; i<50; i++) {
		//	ptavebins[i] = 6+(i-47)*2;
		//}
		//ptavebins[50] = 20;
	
		/*
		const int Nptavebins = 16;
                double ptavebins[Nptavebins+1];
		for(int i = 0; i<10; i++) {
			ptavebins[i] = 2+i*0.2;
		}
		for(int i = 10; i<14; i++) {
			ptavebins[i] = 4+(i-9)*0.5;
		}
		for(int i = 14; i<16; i++) {
			ptavebins[i] = 6+(i-13)*2;
		}
		ptavebins[16] = 20;
		*/

		h = new TH2F(name, title, WNbins, Wptbins, Nptavebins, ptavebins);
		cout<<Nptavebins<<" bins"<<endl;
		//fortest fine bin in x, wide bin in y
		//double *xbins = new double[ntx+1];		
		//for(int ibin = 0; ibin<ntx; ibin++) {
		//	*(xbins+ibin) = xlo + ibin*(xhi-ylo)/ntx;	
		//}
		//*(xbins+ntx) = xhi;
		//h = new TH2F(name, title, ntx, xbins, Nptavebins, ptavebins);
	
	}
	else 
	{
		h = InitWideX2DHisto(name, title, nty, ylo, yhi);
		cout<<nty<<" bins"<<endl;
	}



	return h;
	
}


TH2F* Unfold2D::InitWideXX2DHisto(TString name, TString title) 
{
	double MOnePlusWptbins[WNbins+2] = {-2};		// first bin for fake/loss
	for(int i = 1; i<WNbins+2; i++) {
		MOnePlusWptbins[i] = Wptbins[i-1];
	}

	cout<<"Init XX WideBin for "<<name<<" "<<title<<" with "<<WNbins+1<<" bins"<<endl;

	TH2F* h = new TH2F(name, title, WNbins+1, MOnePlusWptbins, WNbins+1, MOnePlusWptbins);
	return h;
}


TH2F* Unfold2D::InitWideYY2DHisto(TString name, TString title) 
{
	TH2F *h;

	cout<<"Init YY WideBin for "<<name<<" "<<title<<" with ";
	if(YvariableName.Contains("Ntrk")) 
	{
		const int Nltrkbins = 19;
		double ltrkbins[Nltrkbins+1] = {-2, -0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,17.5,21.5,39.5};	// first bin for fake/loss
		h = new TH2F(name, title, Nltrkbins, ltrkbins, Nltrkbins, ltrkbins);
		cout<<Nltrkbins<<" bins"<<endl;
	}

	else if(YvariableName.Contains("PtAve")) 
	{
		//const int Nptavebins = 178;
		//double ptavebins[Nptavebins+1];
		//if MC unfolded to pT==0
		//ptavebins[0] = 0;
		//for(int i = 1; i<132; i++) {
		//	ptavebins[i] = 0.2+(i-1)*0.01;
		//}
		//for(int i = 132; i<142; i++) {
		//	ptavebins[i] = 1.5+(i-131)*0.05;
		//}
		//for(int i = 142; i<162; i++) {
		//	ptavebins[i] = 2+(i-141)*0.1;
		//}
		//for(int i = 162; i<172; i++) {
		//	ptavebins[i] = 4+(i-161)*0.2;
		//}
		//for(int i = 172; i<179; i++) {
		//	ptavebins[i] = 6+(i-171)*2;
		//}


		//// if Mc unfolded to pT->0.2
		//const int Nptavebins = 198;
		//double ptavebins[Nptavebins+1];
		//for(int i = 0; i<198; i++) {
		//	ptavebins[i] = 0.2+i*0.1;
		//}
		//ptavebins[198] = 20;


		
		// if Mc unfolded to pT->0.2
		const int Nptavebins = 51;
		double ptavebins[Nptavebins+1];
		ptavebins[0] = -2;
		for(int i = 1; i<27; i++) {
			ptavebins[i] = 0.2+i*0.04;
		}
		for(int i = 27; i<35; i++) {
			ptavebins[i] = 1.2+(i-25)*0.1;
		}
		for(int i = 35; i<45; i++) {
			ptavebins[i] = 2+(i-33)*0.2;
		}
		for(int i = 45; i<49; i++) {
			ptavebins[i] = 4+(i-43)*0.5;
		}
		for(int i = 49; i<51; i++) {
			ptavebins[i] = 6+(i-47)*2;
		}
		ptavebins[51] = 20;

		h = new TH2F(name, title, Nptavebins, ptavebins, Nptavebins, ptavebins);
		cout<<Nptavebins<<" bins"<<endl;
		//fortest fine bin in x, wide bin in y
		//double *xbins = new double[ntx+1];		
		//for(int ibin = 0; ibin<ntx; ibin++) {
		//	*(xbins+ibin) = xlo + ibin*(xhi-ylo)/ntx;	
		//}
		//*(xbins+ntx) = xhi;
		//h = new TH2F(name, title, ntx, xbins, Nptavebins, ptavebins);
	
	}
	else 
	{
		h = new TH2F(name, title, nty, ylo, yhi, nty, ylo, yhi);
		cout<<nty<<" bins"<<endl;
	}


	return h;
	
}



int Unfold2D::Find(std::vector<int> vec, int val) {

        for(int i = 0; i<vec.size(); i++) {
                if(val==vec.at(i)) return i;
        }

        return -1;
}


#endif
