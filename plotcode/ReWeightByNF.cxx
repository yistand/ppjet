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
	hsumNF = NULL;		// 1D vs pT		for JP
	hsumrNF = NULL;		// 1D vs pT		for JP

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
}

bool ReWeightByNF::FillNFRatios() {		// for each pT bin, get MB/JP neutral fraction ratio w_i vs NF_i, sum(w_i*NF_i), sum(NF_i) 

	if(!treejp) return false;
	if(!treemb) return false;

	for(int ievt = 0; ievt<treejp->GetEntries(); ievt++) {
		treejp->GetEntry(ievt);
		
		if(jneutralfracjp>NeutralFracMax||jneutralfracjp<NeutralFracMin) continue;

		if(runidjp<13048000) continue;

		hjp->Fill(jptjp,jneutralfracjp,1);
	}

	vector<float> bincenter;
	for(int i = hsumrNF->FindBin(mbptcut); hsumrNF->GetBinCenter(i)>mbptcut && i<=hsumrNF->GetNbinsX(); i++) {
		bincenter.push_back(hsumrNF->GetBinCenter(i));
	}
	
	for(int ievt = 0; ievt<treemb->GetEntries(); ievt++) {
		treemb->GetEntry(ievt);
		
		if(jneutralfracmb>=NeutralFracMax||jneutralfracmb<=NeutralFracMin) continue;

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

	TH2F *hr = (TH2F*)hmb->Clone("hr");
	hr->SetName("hr");
	hr->SetTitle("JPs ratio MB NeutralFrac");
	hr->Divide(hjp);

	hsumNF = (TH1D*)hjp->ProjectionX("hsumNF");
	hsumNF->SetName("hsumNF");
	hsumNF->SetTitle("Sum(NeutralFraction) for JPs vs p_{T}");

	for(int i = 0; i<hr->GetNbinsX(); i++) {
		float isum = 0;
		for(int j = 0; j<hr->GetNbinsY(); j++) {
			isum+=hsumNF->GetBinContent(i+1)*hr->GetBinContent(i+1,j+1);
		}
		hsumrNF->SetBinContent(i+1,isum);
	}

	for(int i = 0; i<hw->GetNbinsX(); i++) {
		for(int j = 0; j<hw->GetNbinsY(); j++) {
			if(hsumrNF->GetBinContent(i+1)) {
				hw->SetBinContent(i+1,j+1,hr->GetBinContent(i+1,j+1)*hsumNF->GetBinContent(i+1)/hsumrNF->GetBinContent(i+1));
			}
			else {
				hw->SetBinContent(i+1,j+1,0);
			}
		}
	}

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

bool ReWeightByNF::Init4Read(TString filenamejp, TString filenamemb, TString treename) {

	fjp = new TFile(filenamejp);
	treejp = (TTree*)fjp->Get(treename);
	fmb = new TFile(filenamemb);
	treemb = (TTree*)fmb->Get(treename);

	treejp->SetBranchAddress("runid",&runidjp);
	treejp->SetBranchAddress("j1pt",&jptjp);
	treejp->SetBranchAddress("j1neutralfrac",&jneutralfracjp);

	treemb->SetBranchAddress("runid",&runidmb);
	treemb->SetBranchAddress("j1pt",&jptmb);
	treemb->SetBranchAddress("j1neutralfrac",&jneutralfracmb);

        //const int WJbins = 14;
        //float WJptbinning[WJbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,65,100};
        const int WJbins = 10;
        float WJptbinning[WJbins+1] = {0,2,3,4,5,7,9,11,15,20,100};
        //const int WJbins = 5;
        //float WJptbinning[WJbins+1] = {0,2,3,4,5,100};
 
	const int NFbins = 30;			// Neutral Fraction from 0 to 0.90
	float NFbinning[NFbins+1];
	for(int i = 0; i<NFbins+1; i++) {
		NFbinning[i] = 1.*i*(NeutralFracMax-NeutralFracMin)/NFbins+NeutralFracMin;
	}

	hw = new TH2F("hw", "JPs weight to MB NeutralFrac", WJbins, WJptbinning, NFbins, NFbinning);

	hjp = new TH2F("hjp", "JPs Number of Events vs pT and NeutralFrac", WJbins, WJptbinning, NFbins, NFbinning);
	hjp->Sumw2();
	hmb = new TH2F("hmb", "MB Number of Events vs pT and NeutralFrac", WJbins, WJptbinning, NFbins, NFbinning);
	hmb->Sumw2();

	//hsumNF = new TH1F("hsumNF","Sum(NeutralFraction) for JPs vs p_{T}", WJbins, WJptbinning);
	hsumrNF = new TH1F("hsumrNF","Sum(weight*NeutralFraction) for JPs vs p_{T}", WJbins, WJptbinning);

	return true;
}

bool ReWeightByNF::WriteNF(TString outfilename, bool opt) {		// if opt==true, write all intermediate histograms

	TFile *fout = new TFile(outfilename,"RECREATE");
	hw->Write();
	if(opt) {
		hjp->Write();
		hmb->Write();
		hsumNF->Write();
		hsumrNF->Write();
	}
	fout->Close();

	return true;
}

