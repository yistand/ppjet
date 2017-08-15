#include "ReWeightByNF.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <algorithm>    // std::max, std::max_element
#include <vector>


using namespace std;

ReWeightByNF::ReWeightByNF():NeutralFracMax(0.9), NeutralFracMin(0.0), JetPtCorrMax(20) {		// the Max and Min values of NF will not be included
	// output
	hw = NULL;		// final weight 2D histogram of pT bins and Neutral Fraction

	// intermediate histograms
	hjp = NULL;		// # of events per pT bin and neutral fraction bin for JP
	hmb = NULL;		// # of events per pT bin and neutral fraction bin for MB

	htmprx = NULL;		

	// input
	fjp = NULL;
	treejp = NULL;
	fmb = NULL;
	treemb = NULL;

	runidjp = 0;
	jptjp = 0;
	jneutralfracjp = 0;
	
	runidmb = 0;
	jptmb = 0;
	jneutralfracmb = 0;

	// control parameter
	mbptcut = 5;		// for jet pt > mbptcut, fill all the bins in MB histogram to be the same (intergral over pt > mbptcut). The reason is we want to increase statistis and NF distribution is stable for MB triggers after 5 GeV/c. One need to check to make sure this mbptcut is used in pt bin so that it will no be cut into two bins.

}

ReWeightByNF::~ReWeightByNF() {

	if(fjp->IsOpen()) fjp->Close();
}

bool ReWeightByNF::FillNFRatios() {		// for each pT bin, get MB/JP neutral fraction ratio w_i vs NF_i

	if(!treejp) return false;
	if(!treemb) return false;

	for(int ievt = 0; ievt<treejp->GetEntries(); ievt++) {
		treejp->GetEntry(ievt);
		
		if(jneutralfracjp>NeutralFracMax||jneutralfracjp<NeutralFracMin) continue;

		if(jptjp<=0) continue;

		if(runidjp<13048000) continue;

		hjp->Fill(jptjp,jneutralfracjp,1);
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

		// for MB, the last few bins are all filled same as intergrated value for the last few pT bins
		if(jptmb>mbptcut) {
			for(vector<float>::iterator it = bincenter.begin(); it!=bincenter.end(); it++) {
				hmb->Fill(*it,jneutralfracmb,1);
			}
		}
		else {
			hmb->Fill(jptmb,jneutralfracmb,1);
		}
	}

	hw = (TH2F*)hmb->Clone("hw");
	hw->SetName("hw");
	hw->SetTitle("JPs ratio MB NeutralFrac");
	hw->Divide(hjp);

	return true;

}

double ReWeightByNF::GetNFWeight(double pt, double nf) {
	return (double) GetNFWeight((float)pt, (float)nf);
}

float ReWeightByNF::GetNFWeight(float pt, float nf) {

	if(!hw) { cout<<"ERR!!!  float ReWeightByNF::GetNFWeight(float pt, float nf): hw is empty!!"<<endl;return 0;}

	if(nf<=NeutralFracMin || nf>=NeutralFracMax) return 0;
	if(pt>JetPtCorrMax) return 1;		// no correction for pT> 20 GeV/c
	int ithbin = hw->FindBin(pt,nf);
	return hw->GetBinContent(ithbin);
}

bool ReWeightByNF::AdjustNevts(TH2D *h2, TH1D *hx)	{		// adjust h2 x-axis has same distribution as hx. hx should have same binning as h2 x-axis
	TH1D *h1 = (TH1D*)h2->ProjectionX(Form("%sx",h2->GetName()));
	for(int i = 0; i<hx->GetNbinsX(); i++) {
		double modelNevts = hx->GetBinContent(i+1);
		double Aratio = 0;
		if(modelNevts>0) {
			Aratio = modelNevts/h1->GetBinContent(i+1);
		}
		for(int j = 0; j<h2->GetNbinsY(); j++) {
			h2->SetBinContent(i+1,j+1,h2->GetBinContent(i+1,j+1)*Aratio);
		}
	}
	delete h1;
	return true;
}


bool ReWeightByNF::AdjustNevts(TH2F *h2, TH1D *hx)	{		// adjust h2 x-axis has same distribution as hx. hx should have same binning as h2 x-axis
	TH1D *h1 = (TH1D*)h2->ProjectionX(Form("%sx",h2->GetName()));
	for(int i = 0; i<hx->GetNbinsX(); i++) {
		double modelNevts = hx->GetBinContent(i+1);
		double Aratio = 0;
		if(modelNevts>0) {
			Aratio = modelNevts/h1->GetBinContent(i+1);
		}
		for(int j = 0; j<h2->GetNbinsY(); j++) {
			h2->SetBinContent(i+1,j+1,h2->GetBinContent(i+1,j+1)*Aratio);
		}
	}
	delete h1;
	return true;
}

bool ReWeightByNF::Init4Read(TString filenamejp, TString filenamemb, TString treename) {

	fjp = new TFile(filenamejp);
	if(fjp->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading "<<filenamejp<<endl;
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
	if(fmb->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading "<<filenamejp<<endl;
	else {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot open TFile "<<filenamejp<<endl;
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

	treemb->SetBranchAddress("runid",&runidmb);
	treemb->SetBranchAddress("j1pt",&jptmb);
	treemb->SetBranchAddress("j1neutralfrac",&jneutralfracmb);

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
 
	const int NFbins = 30;			// Neutral Fraction from 0 to 0.90
	float NFbinning[NFbins+1];
	for(int i = 0; i<NFbins+1; i++) {
		NFbinning[i] = 1.*i*(NeutralFracMax-NeutralFracMin)/NFbins+NeutralFracMin;
	}

	hjp = new TH2F("hjp", "JPs Number of Events vs pT and NeutralFrac", WJbins, WJptbinning, NFbins, NFbinning);
	hjp->Sumw2();
	hmb = new TH2F("hmb", "MB Number of Events vs pT and NeutralFrac", WJbins, WJptbinning, NFbins, NFbinning);
	hmb->Sumw2();

	htmprx = new TH1D("htmprx","tmp", WJbins, WJptbinning);

	return true;
}

bool ReWeightByNF::WriteNF(TString outfilename, bool opt) {		// if opt==true, write all intermediate histograms

	TFile *fout = new TFile(outfilename,"RECREATE");
	hw->Write();
	if(opt) {
		hjp->Write();
		hmb->Write();
	}
	fout->Close();

	return true;
}


bool ReWeightByNF::ReadNFfile(TString infilename) {
	TFile *fin  = new TFile(infilename);
	cout<<__PRETTY_FUNCTION__<<" read in "<<infilename<<endl;
	if(!fin->IsOpen())	{
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<" Cannot open NFfile " <<infilename<<endl;
		return false;
	}
	hw = (TH2F*)fin->Get("hw");
	if(!hw) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<" Cannot find hw in "<<infilename<<endl;
		return false;
	}

	return true;
}
