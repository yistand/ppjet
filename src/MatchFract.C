//==========================================================================================================================================================================
//
//		2016.07.29		Li YI
//		what is the ratio of Rc jet leading jet matched to JP2 trigger patch
//
//==========================================================================================================================================================================
#include "include/CrossSectionPerpT.h"

////HC run6 pp 200GeV root files:
//const int NUMBEROFPT = 11;
//const char *PTBINS[NUMBEROFPT]={"3_4","4_5","5_7","7_9","9_11","11_15","15_25","25_35","35_45","45_55","55_65"};


void MatchFract() {
	
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	int color[11] = {1,2,kGreen+2,4,kOrange,6,kCyan+1,8,9,12,41};		// check consistence with src/MergeByXsec.cxx


	TFile *f1, *f2; 
	TH1D *h1, *h2;
	TString histname = "RcLeadJetPt";
	TString histname2 = "McLeadJetPt";

	TCanvas *c = new TCanvas();

	for(int i = 0; i<NUMBEROFPT; i++) {		// NUMBEROFPT defined in include/CrossSectionPerpT.h
		f1 = new TFile(Form("~/Scratch/embedPythia/pt%s_JetMcVsEmbedMatchTrig_nobbc.root",PTBINS[i]));
		//f1 = new TFile(Form("~/Scratch/embedPythia/pt%s_JetMcVsEmbedMatchTrig.root",PTBINS[i]));
		//f1 = new TFile(Form("~/Scratch/embedPythia/HCpt%s_JetMcVsEmbedMatchTrig.root",PTBINS[i]));
		h1 = (TH1D*)f1->Get(histname);	

		f2 = new TFile(Form("~/Scratch/embedPythia/pt%s_JetMcVsEmbed_nobbc.root",PTBINS[i]));
		//f2 = new TFile(Form("~/Scratch/embedPythia/pt%s_JetMcVsEmbedJP2.root",PTBINS[i]));
		//f2 = new TFile(Form("~/Scratch/embedPythia/pt%s_JetMcVsEmbed.root",PTBINS[i]));
		//f2 = new TFile(Form("~/Scratch/embedPythia/HCpt%s_JetMcVsEmbed.root",PTBINS[i]));
		h2 = (TH1D*)f2->Get(histname);	
		//h2 = (TH1D*)f1->Get("McLeadJetPt");	

		h1->Divide(h2);	

		h1->SetMaximum(1);
		h1->SetMinimum(0);

		h1->GetXaxis()->SetTitle("Leading jet p_{T}");	
		h1->GetYaxis()->SetTitle("Matched to JetPatch2 Trigger");	
		h1->GetXaxis()->SetRangeUser(0,60);
		h1->SetLineColor(color[i]);
		h1->SetMarkerColor(color[i]);
		h1->SetMarkerStyle(8);//7);

		c->cd();
		if(i==0) h1->Draw("");
		else h1->Draw("same");
	} 

	//c->SaveAs(Form("HCMatchRatiopt%s.png","2_-1"));
	//c->SaveAs(Form("MatchRatioInJP2pt%s.png","2_-1"));
	//c->SaveAs(Form("MatchRatiopt%s.png","2_-1"));
	//c->SaveAs(Form("JP2Ratiopt%s.png","2_-1"));
}
