#include "RcVzWeight.h"

//#include "/home/fas/caines/ly247/code/ppjet/include/CrossSectionPerpT.h"
#include "CrossSectionPerpT.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <algorithm>    // std::max, std::max_element
#include <vector>


using namespace std;

RcVzWeight::RcVzWeight():VzRange(30) {
	// output
	hw = NULL;		// final weight 1D histogram of Vz 
	fpol = NULL;		// polynomial fit to hw

	// intermediate histograms
	hrc = NULL;		// # of events per Vz for Rc embedding data
	hdata = NULL;		// # of events per Vz for real data 

	// input
	frc = NULL;		// read embedding data, will loop over files and fill with cross section weight
	treerc = NULL;

	fdata = NULL;
	treedata = NULL;

	runidrc = 0;
	jptrc = 0;
	vzrc = 0;

	xsecw = 0;

	runiddata = 0;
	jptdata = 0;
	vzdata = 0;

}

RcVzWeight::~RcVzWeight() {

	if(fdata->IsOpen()) fdata->Close();
	if(frc->IsOpen()) frc->Close();
}


double RcVzWeight::GetVzWeight(double vz) {
	return (double) GetVzWeight((float)vz);
}

float RcVzWeight::GetVzWeight(float vz) {

	if(fpol) {
		return fpol->Eval(vz);
	}
	else { 
		cout<<"WARNING!!!  "<<__PRETTY_FUNCTION__<<": fpol is empty!!"<<endl;

		if(hw) {
			int ithbin = hw->FindBin(vz);
			return hw->GetBinContent(ithbin);
		}
		else {
			cout<<"ERR!!!  "<<__PRETTY_FUNCTION__<<": hw is empty!!"<<endl;
			return 0;
		}
	}
}


bool RcVzWeight::InitHist(TString filetagrc) {
	hrc = new TH1D("hrc","hrc"+filetagrc,VzRange*2,-VzRange,VzRange);
	hdata = new TH1D("hdata","hdata"+filetagrc,VzRange*2,-VzRange,VzRange);
	return true;
}

bool RcVzWeight::ReadAndFill(TString filetagrc, TString filenamedata) {

	InitHist(filetagrc);

	TString treename = "ResultTree";
	TString dirrc = "/home/fas/caines/ly247/Scratch/embedPythia/" + filetagrc + "/";

	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename = dirrc+Form("pt%s_underMcVsEmbed_FullJetTransCharged_McPt02.root",PTBINS[i]);
		frc = new TFile(ifilename);
		if(frc->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading "<<ifilename<<endl;
		else {
			cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot open TFile "<<ifilename<<endl;
			return 0;
		}
		treerc = (TTree*)frc->Get(treename);
		if(!treerc) {
			cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot read TTree "<<treename<<" from "<<ifilename<<endl;
			return 0;
		}

		treerc->SetBranchAddress("runid",&runidrc);
		treerc->SetBranchAddress("Rcj1pt",&jptrc);
		treerc->SetBranchAddress("Rcvz",&vzrc);

		// Weight per pT bin
		double xsecweight = XSEC[i]/NUMBEROFEVENT[i];
		//{ cout<<__PRETTY_FUNCTION__<<": pT bin "<<PTBINS[i]<<" Use "<<NUMBEROFEVENT[i]<<" for NUMBEROFEVENT in weight = XSEC["<<i<<"]/NUMBEROFEVENT["<<i<<"]"<<endl;}
		xsecw = xsecweight;

		//cout<<"Total events: "<<treerc->GetEntries()<<endl;
		cout<<"xsec weight for "<<PTBINS[i]<<": "<<xsecw<<endl;

		for(int ievt = 0; ievt<treerc->GetEntries(); ievt++) {
			treerc->GetEntry(ievt);

			if(jptrc<=0) continue;

			if(runidrc<13048000) continue;

			hrc->Fill(vzrc,xsecw);
		}
		frc->Close();
	}




	if(!filenamedata.Contains(filetagrc)) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<":  Data file "<<filenamedata<<" and Rc embedding file request "<<filetagrc<<" not consistent"<<endl;
		return 0;
	}

	fdata = new TFile(filenamedata);
	if(fdata->IsOpen()) cout<<__PRETTY_FUNCTION__<<"Reading "<<filenamedata<<endl;
	else {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot open TFile "<<filenamedata<<endl;
		return 0;
	}
	treedata = (TTree*)fdata->Get(treename);
	if(!treedata) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<": cannot read TTree "<<treename<<endl;
		return 0;
	}

	treedata->SetBranchAddress("runid",&runiddata);
	treedata->SetBranchAddress("j1pt",&jptdata);
	treedata->SetBranchAddress("vz",&vzdata);



	for(int ievt = 0; ievt<treedata->GetEntries(); ievt++) {
		treedata->GetEntry(ievt);

		if(jptdata<=0) continue;

		if(runiddata<13048000) continue;

		hdata->Fill(vzdata,1);
	}
	fdata->Close();

	hw = (TH1D*)hdata->Clone("hw");
	hw->SetTitle("hw"+filetagrc);
	hw->Divide(hrc);

	double sum = hw->Integral();
	hw->Scale(VzRange*2/sum);	// average 1

	fpol = new TF1("fpol","[0]+[1]*x+[2]*x*x",-30,30);
	hw->Fit(fpol);

	return true;
}

bool RcVzWeight::WriteVzFile(TString outfilename) {		

	TFile *fout = new TFile(outfilename,"RECREATE");
	hw->Write();
	fpol->Write();
	hrc->Write();
	hdata->Write();
	fout->Close();

	return true;
}


bool RcVzWeight::ReadVzfile(TString infilename) {
	TFile *fin  = new TFile(infilename);
	cout<<__PRETTY_FUNCTION__<<" read in "<<infilename<<endl;
	if(!fin->IsOpen())	{
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<" Cannot open Vz weight file " <<infilename<<endl;
		return false;
	}
	hw = (TH1D*)fin->Get("hw");
	if(!hw) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<" Cannot find hw in "<<infilename<<endl;
		return false;
	}
	fpol = (TF1*)fin->Get("fpol");
	if(!fpol) {
		cout<<"ERR!!! "<<__PRETTY_FUNCTION__<<" Cannot find fpol in "<<infilename<<endl;
		return false;
	}

	return true;
}

bool RcVzWeight::SetVzRange(double range) {
	VzRange = range;
	return true;
}
