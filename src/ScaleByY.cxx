//======================================================================================================================================================
//
//		2017.01.18	Li Yi
//		Modified from ReWegithByNF
//		In order to do syst. err. check, we scale mean of Y for each jet pt bin (TranTotNtrk, LeadAreaNtrk
//		SubAreaNtrk..) from JP without VPD with the weight from MB with VPD/JP with VPD
//
//======================================================================================================================================================


#include "ScaleByY.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <algorithm>    // std::max, std::max_element
#include <vector>


using namespace std;

ScaleByY::ScaleByY() {		
	ScaleByY("NeutralFrac");
}


ScaleByY::ScaleByY(TString yname): YVariableName(yname), NeutralFracMax(0.9), NeutralFracMin(0.0), JetPtCorrMax(11) {		// the Max and Min values of NF will not be included

	// output
	h1w = NULL;		// weight 1D histogram. Used by other variables than Neutral Fraction, such as TranTotNtrk

	// intermediate histograms
	htmprx = NULL;		

	pfjp = NULL;		// Y vs jet pt TProfile
	pfmb = NULL;		// Y vs jet pt TProfile


	// input
	fjp = NULL;
	treejp = NULL;
	fmb = NULL;
	treemb = NULL;

	runidjp = 0;
	jptjp = 0;
	jneutralfracjp = 0;
	Yjp = 0;
	
	runidmb = 0;
	jptmb = 0;
	jneutralfracmb = 0;
	Ymb = 0;

	// control parameter
	mbptcut = JetPtCorrMax;		// for jet pt > mbptcut,  One need to check to make sure this mbptcut is used in pt bin so that it will no be cut into two bins.
					// at pt = 11, VPD cut result for JP and MB are consistent

}

ScaleByY::~ScaleByY() {

	if(fjp->IsOpen()) fjp->Close();
}

bool ScaleByY::FillYRatios() {		// for each pT bin, get MB/JP Y ratio w_i vs Y_i


	if(!treejp) return false;
	if(!treemb) return false;

	for(int ievt = 0; ievt<treejp->GetEntries(); ievt++) {
		treejp->GetEntry(ievt);
		
		if(jneutralfracjp>NeutralFracMax||jneutralfracjp<NeutralFracMin) continue;

		if(jptjp<=0) continue;

		if(runidjp<13048000) continue;

		if(YVariableName.Contains("LeadAreaNtrk", TString::kIgnoreCase)) Yjp = leadntrkjp;
		else if(YVariableName.Contains("SubAreaNtrk", TString::kIgnoreCase)) Yjp = subntrkjp;
		else if(YVariableName.Contains("TranMaxNtrk", TString::kIgnoreCase)) Yjp = maxntrkjp;
		else if(YVariableName.Contains("TranMinNtrk", TString::kIgnoreCase)) Yjp = minntrkjp;
		else if(YVariableName.Contains("NeutralFrac", TString::kIgnoreCase)) Yjp = jneutralfracjp;
		else if(YVariableName.Contains("TranTotNtrk", TString::kIgnoreCase)) Yjp = maxntrkjp+minntrkjp;
		else { cout<<"Wrong YVariableName for class ScaleByY!!!!"<<endl; return false;}

		pfjp->Fill(jptjp,Yjp,1);
	}

	vector<float> bincenter;
	for(int i = htmprx->FindBin(mbptcut); htmprx->GetBinCenter(i)>mbptcut && i<=htmprx->GetNbinsX(); i++) {
		bincenter.push_back(htmprx->GetBinCenter(i));
	}
	
	for(int ievt = 0; ievt<treemb->GetEntries(); ievt++) {
		treemb->GetEntry(ievt);
		
		if(jneutralfracmb>=NeutralFracMax||jneutralfracmb<=NeutralFracMin) continue;

		if(jptjp<=0) continue;

		if(runidmb<13048000) continue;


		if(YVariableName.Contains("LeadAreaNtrk", TString::kIgnoreCase)) Ymb = leadntrkmb;
		else if(YVariableName.Contains("SubAreaNtrk", TString::kIgnoreCase)) Ymb = subntrkmb;
		else if(YVariableName.Contains("TranMaxNtrk", TString::kIgnoreCase)) Ymb = maxntrkmb;
		else if(YVariableName.Contains("TranMinNtrk", TString::kIgnoreCase)) Ymb = minntrkmb;
		else if(YVariableName.Contains("NeutralFrac", TString::kIgnoreCase)) Ymb = jneutralfracmb;
		else if(YVariableName.Contains("TranTotNtrk", TString::kIgnoreCase)) Ymb = maxntrkmb+minntrkmb;
		else { cout<<"Wrong YVariableName for class ScaleByY!!!!"<<endl; return false;}


		// for MB, the last few bins are all filled same as intergrated value for the last few pT bins
		if(jptmb>mbptcut) {
			for(vector<float>::iterator it = bincenter.begin(); it!=bincenter.end(); it++) {
				pfmb->Fill(*it,Ymb,1);
			}
		}
		else {
			pfmb->Fill(jptmb,Ymb,1);
		}
	}

	for(int i = 0; i<h1w->GetNbinsX(); i++) {
		float yjp = pfjp->GetBinContent(i+1);
		float ymb = pfmb->GetBinContent(i+1);
		h1w->SetBinContent(i+1,(yjp>0?ymb/yjp:1));
	}

	return true;

}

double ScaleByY::GetYScale(double pt, double nf) {
	return (double) GetYScale((float)pt, (float)nf);
}

float ScaleByY::GetYScale(float pt, float nf) {

	if(!h1w) { cout<<"ERR!!!  float ScaleByY::GetYScale(float pt, float nf): hw is empty!!"<<endl;return 0;}

	if(nf<=NeutralFracMin || nf>=NeutralFracMax) return 0;
	if(pt>JetPtCorrMax) return 1;		// no correction for pT> 11 GeV/c

	int ithbin = h1w->FindBin(pt);
	return h1w->GetBinContent(ithbin);

}

bool ScaleByY::Init4Read(TString filenamejp, TString filenamemb, TString treename) {

	fjp = new TFile(filenamejp);
	if(fjp->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading for JP "<<filenamejp<<endl;
	else {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot open TFile "<<filenamejp<<endl;
		return 0;
	}
	treejp = (TTree*)fjp->Get(treename);
	if(!treejp) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot read TTree "<<treename<<endl;
		return 0;
	}
	fmb = new TFile(filenamemb);
	if(fmb->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading for MB "<<filenamemb<<endl;
	else {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot open TFile "<<filenamemb<<endl;
		return 0;
	}
	treemb = (TTree*)fmb->Get(treename);
	if(!treemb) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot read TTree "<<treename<<endl;
		return 0;
	}

	treejp->SetBranchAddress("runid",&runidjp);
	treejp->SetBranchAddress("j1pt",&jptjp);
	treejp->SetBranchAddress("j1neutralfrac",&jneutralfracjp);
	treejp->SetBranchAddress("LeadAreaNtrk",&leadntrkjp);
	treejp->SetBranchAddress("SubAreaNtrk",&subntrkjp);
	treejp->SetBranchAddress("TranMaxNtrk",&maxntrkjp);
	treejp->SetBranchAddress("TranMinNtrk",&minntrkjp);

	treemb->SetBranchAddress("runid",&runidmb);
	treemb->SetBranchAddress("j1pt",&jptmb);
	treemb->SetBranchAddress("j1neutralfrac",&jneutralfracmb);
	treemb->SetBranchAddress("LeadAreaNtrk",&leadntrkmb);
	treemb->SetBranchAddress("SubAreaNtrk",&subntrkmb);
	treemb->SetBranchAddress("TranMaxNtrk",&maxntrkmb);
	treemb->SetBranchAddress("TranMinNtrk",&minntrkmb);

        //const int WJbins = 14;
        //float WJptbinning[WJbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,65,100};
        const int WJbins = 10;			// default
        float WJptbinning[WJbins+1] = {0,2,3,4,5,7,9,11,15,20,100};		// default
        //const int WJbins = 5;
        //float WJptbinning[WJbins+1] = {0,2,3,4,5,100};
	//const int WJbins = 100;
	//float WJptbinning[WJbins+1];
	//for(int i = 0; i<WJbins+1; i++) {
	//	WJptbinning[i] = 1.*i;
	//}
 
	//const int NFbins = 30;			// Neutral Fraction from 0 to 0.90
	//float NFbinning[NFbins+1];
	int Ybins;			// Neutral Fraction from 0 to 0.90
	float *Ybinning;
	if(YVariableName.Contains("Ntrk",TString::kIgnoreCase)) {
		Ybins = 40;
		Ybinning = new float[Ybins];
	}
	else if(YVariableName.Contains("NeutralFrac",TString::kIgnoreCase)) {
		Ybins = 30;
		Ybinning = new float[Ybins];
	}

	for(int i = 0; i<Ybins+1; i++) {
		if(YVariableName.Contains("Ntrk",TString::kIgnoreCase)) {
			Ybinning[i] = i;
		}
		else if(YVariableName.Contains("NeutralFrac",TString::kIgnoreCase)) {
			Ybinning[i] = 1.*i*(NeutralFracMax-NeutralFracMin)/Ybins+NeutralFracMin;
		}
	}

	htmprx = new TH1D("htmprx","tmp", WJbins, WJptbinning);

	pfjp = new TProfile("pfjp", "JPs Number of Events vs pT and "+YVariableName, WJbins, WJptbinning);
	pfmb = new TProfile("pfmb", "JPs Number of Events vs pT and "+YVariableName, WJbins, WJptbinning);

	h1w = new TH1F("h1w", "Weights for "+YVariableName+" Vs jet pt", WJbins, WJptbinning);

	return true;
}

bool ScaleByY::WriteY(TString outfilename, bool opt) {		// if opt==true, write all intermediate histograms

	TFile *fout = new TFile(outfilename,"RECREATE");
	h1w->Write();
	if(opt) {
		pfjp->Write();
		pfmb->Write();
	}
	fout->Close();

	return true;
}

