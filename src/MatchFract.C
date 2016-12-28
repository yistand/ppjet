//==========================================================================================================================================================================
//
//		2016.09.19		Li YI
//		Update to reflect:
//		1. trigger efficiency
//		2. match to trigger patch efficiency
//		3. good (inside eta acceptance) Rc (not istrigger or trigger matched) to any Mc matching ratio 
//		4. good (inside eta acceptance) Mc to any Rc (not istrigger or trigger matched) matching ratio
//		Not plot seperately for each pT bin. Instead merging together according to their weights
//
//		2016.07.29		Li YI
//		what is the ratio of Rc jet leading jet matched to JP2 trigger patch
//
//==========================================================================================================================================================================

#include "include/CrossSectionPerpT.h"


#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"


#include <iostream>


Int_t ReadbyXsec (TH1D *hrc, TH1D *hrc_istrig, TH1D *hrc_trigmatch, TH1D *hrc_match1, TH1D *hrc_match2, TH1D *hmc, TH1D *hmc_match1, TH1D *hmc_match2, TString trigname);

//==================================================   Main   ==================================================
void MatchFract(TString trigname="JP1") {
	
	gStyle->SetOptStat(0);
	//gStyle->SetOptTitle(0);

	TFile *f1, *f2; 
	const int Nbins = 15;
        double Xbins[Nbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};
	TH1D *hrc = new TH1D("hrc","",Nbins, Xbins);
 	TH1D *hrc_istrig = new TH1D("hrc_istrig","",Nbins, Xbins);
 	TH1D *hrc_trigmatch = new TH1D("hrc_trigmatch","",Nbins, Xbins);
 	TH1D *hrc_match1 = new TH1D("hrc_match1","",Nbins, Xbins);
 	TH1D *hrc_match2 = new TH1D("hrc_match2","",Nbins, Xbins);
 	TH1D *hmc = new TH1D("hmc","",Nbins, Xbins);
 	TH1D *hmc_match1 = new TH1D("hmc_match1","",Nbins, Xbins);
 	TH1D *hmc_match2 = new TH1D("hmc_match2","",Nbins, Xbins);


	ReadbyXsec (hrc, hrc_istrig, hrc_trigmatch, hrc_match1, hrc_match2, hmc, hmc_match1, hmc_match2, trigname);
	
	hrc->Sumw2();
	hrc_istrig->Sumw2();
	hrc_trigmatch->Sumw2();
	hrc_match1->Sumw2();
	hrc_match2->Sumw2();
	hmc->Sumw2();
	hmc_match1->Sumw2();
	hmc_match2->Sumw2();

	TH1D *r_trigfrac, *r_matchfrac, *r_matchonly, *r_rc_1, *r_rc_2, *r_mc_1, *r_mc_2;

	double xmax = 60;

	r_trigfrac = (TH1D*)hrc_istrig->Clone("r_trigfrac");
	r_trigfrac->Divide(hrc);
	r_trigfrac->SetTitle("Trigger Efficiency");
	r_trigfrac->GetXaxis()->SetTitle("Detector-level leading jet p_{T}");
	r_trigfrac->GetXaxis()->SetRangeUser(0,xmax);
	r_trigfrac->SetMinimum(0);

	r_matchfrac = (TH1D*)hrc_trigmatch->Clone("r_matchfrac");
	r_matchfrac->Divide(hrc);
	r_matchfrac->SetTitle("Trigger Efficiency & Patch Matching Ratio");
	r_matchfrac->GetXaxis()->SetTitle("Detector-level leading jet p_{T}");
	r_matchfrac->GetXaxis()->SetRangeUser(0,xmax);
	r_matchfrac->SetMinimum(0);

	r_matchonly = (TH1D*)hrc_trigmatch->Clone("r_matchonly");
	r_matchonly->Divide(hrc_istrig);
	r_matchonly->SetTitle("Trigger Patch Matching Ratio");
	r_matchonly->GetXaxis()->SetTitle("Detector-level leading jet p_{T}");
	r_matchonly->GetXaxis()->SetRangeUser(0,xmax);
	r_matchonly->SetMinimum(0);


	r_rc_1 = (TH1D*)hrc_match1->Clone("r_rc_1");
	r_rc_1->Divide(hrc);
	r_rc_1->SetTitle("Leading Detector-level |#eta|<0.4 Matched to Leading Particle-level Jet Fraction");
	r_rc_1->GetXaxis()->SetTitle("Detector-level leading jet p_{T}");
	r_rc_1->GetXaxis()->SetRangeUser(0,xmax);
	r_rc_1->SetMinimum(0);

	r_rc_2 = (TH1D*)hrc_match2->Clone("r_rc_2");
	r_rc_2->Divide(hrc);
	r_rc_2->SetTitle("Leading Detector-level |#eta|<0.4 Matched to SubLeading Particle-level Jet Fraction");
	r_rc_2->GetXaxis()->SetTitle("Detector-level leading jet p_{T}");
	r_rc_2->GetXaxis()->SetRangeUser(0,xmax);
	r_rc_2->SetMinimum(0);

	r_mc_1 = (TH1D*)hmc_match1->Clone("r_mc_1");
	r_mc_1->Divide(hmc);
	r_mc_1->SetTitle("Leading Particle-level |#eta|<0.4 Matched to Leading Detector-level Jet Fraction");
	r_mc_1->GetXaxis()->SetTitle("Particle-level leading jet p_{T}");
	r_mc_1->GetXaxis()->SetRangeUser(0,xmax);
	r_mc_1->SetMinimum(0);

	r_mc_2 = (TH1D*)hmc_match2->Clone("r_mc_2");
	r_mc_2->Divide(hmc);
	r_mc_2->SetTitle("Leading Particle-level |#eta|<0.4 Matched to SubLeading Detector-level Jet Fraction");
	r_mc_2->GetXaxis()->SetTitle("Particle-level leading jet p_{T}");
	r_mc_2->GetXaxis()->SetRangeUser(0,xmax);
	r_mc_2->SetMinimum(0);



	
	TCanvas *c[7];
	for(int i = 0; i<7; i++) {
		c[i] = new TCanvas();
	}

	TLine *line = new TLine(0,1,xmax,1);

	c[0]->cd();
	r_trigfrac->Draw();
	line->Draw();
	c[0]->SaveAs("fig/r_trigfrac_"+trigname+".png");


	c[1]->cd();
	r_matchfrac->Draw();
	line->Draw();
	c[1]->SaveAs("fig/r_matchfrac_"+trigname+".png");

	c[2]->cd();
	r_matchonly->Draw();
	line->Draw();
	c[2]->SaveAs("fig/r_matchonly_"+trigname+".png");

	c[3]->cd();
	r_rc_1->Draw();
	c[3]->SaveAs("fig/r_rc_1_"+trigname+".png");

	c[4]->cd();
	r_rc_2->Draw();
	c[4]->SaveAs("fig/r_rc_2_"+trigname+".png");

	c[5]->cd();
	r_mc_1->Draw();
	c[5]->SaveAs("fig/r_mc_1_"+trigname+".png");

	c[6]->cd();
	r_mc_2->Draw();
	c[6]->SaveAs("fig/r_mc_2_"+trigname+".png");



	TFile *fout = new TFile("MatchFrac"+trigname+".root","RECREATE");
	r_trigfrac->Write();
	r_matchfrac->Write();
	r_matchonly->Write();
	r_rc_1->Write();
	r_rc_2->Write();
	r_mc_1->Write();
	r_mc_2->Write();


}



//==============================================================================
Int_t ReadbyXsec (TH1D *hrc, TH1D *hrc_istrig, TH1D *hrc_trigmatch, TH1D *hrc_match1, TH1D *hrc_match2, TH1D *hmc, TH1D *hmc_match1, TH1D *hmc_match2, TString trigname)
{
	Bool_t flagIsTrigger;
	Bool_t trigmatch;
	Int_t MatchedNthMcj;
	Int_t MatchedNthRcj;
	
	Float_t Mcj1pt;
	Float_t Rcj1pt;
	
	for(int i = 0; i<NUMBEROFPT; i++) {
		TString ifilename = Form("/home/fas/caines/ly247/Scratch/embedPythia/%s/pt%s_underMcVsEmbed_FullJetTransCharged.root",trigname.Data(),PTBINS[i]);
		cout<<"Read in "<<ifilename<<endl;
		TFile *ftrain = new TFile(ifilename);
		TTree *tree = (TTree*)ftrain->Get("ResultTree");
		if(!tree) {cout<<"Cannot find tree from input in Unfold2D::ReadbyXsec()"<<endl;exit;}
		tree->SetBranchAddress("flagIsTrigger",&flagIsTrigger);
		tree->SetBranchAddress("trigmatch",&trigmatch);
		tree->SetBranchAddress("MatchedNthMcj",&MatchedNthMcj);
		tree->SetBranchAddress("MatchedNthRcj",&MatchedNthRcj);
		tree->SetBranchAddress("Mcj1pt",&Mcj1pt);
		tree->SetBranchAddress("Rcj1pt",&Rcj1pt);

		// Weight per pT bin
		double weight = XSEC[i]/NUMBEROFEVENT[i];

		for (Int_t j= 0; j<tree->GetEntries(); j++) {		// loop over entries
			tree->GetEntry(j);

			if(Rcj1pt>0) {
				hrc->Fill(Rcj1pt);
				if(flagIsTrigger) hrc_istrig->Fill(Rcj1pt);
				if(trigmatch) hrc_trigmatch->Fill(Rcj1pt);
				if(MatchedNthMcj==0) hrc_match1->Fill(Rcj1pt);
				if(MatchedNthMcj==1) hrc_match2->Fill(Rcj1pt);
			}
			if(Mcj1pt>0) {
				hmc->Fill(Mcj1pt);
				if(MatchedNthRcj==0) hmc_match1->Fill(Mcj1pt);
				if(MatchedNthRcj==1) hmc_match2->Fill(Mcj1pt);
			}
		}
		ftrain->Close();
	}

	return 1;
}

