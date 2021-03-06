//=====================================================================================================
//
//	2018.08.13	Li YI
//	plot NFweighted JPs with noweight JP0, JP1, JP2 for Towards, Away & Transverse on the same canvas
//
//=====================================================================================================
//
//	2017.11.02	Li YI
//	compare the NFweighted JPs with no weight JP0, JP1, JP2 results
//
//=====================================================================================================
//
//	2016.02.26	Li YI
//
//	plot the various histograms obtained from plotpT_fromJetTree.C togehter on the same canvas
//
//
//=====================================================================================================

#include "TColor.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"


void plotcompareTrigUnfold_all() {

	int savefig = 0;

	//const int Nvar = 9;
	//char filenametag[Nvar][20] = {"TranTotNtrk","TranPtAve","TranTotPtSum","LeadAreaNtrk","LeadPtAve","LeadAreaPtSum","SubAreaNtrk","SubPtAve","SubAreaPtSum"};
	const int Nvar = 6;
	char filenametag[Nvar][20] = {"LeadAreaNtrk","SubAreaNtrk","TranTotNtrk","LeadPtAve","SubPtAve","TranPtAve"};

	const int Nunf = 4;
	char unfoldtag[Nunf][20] = {"2Charged","1Charged","0Charged","Charged_NFweight"};

	char filenameformat[200] = {"Unfolding_%sJP%s_BT170928_12JetBinv2_McPt02_embed%s_Baye5.root"};
	char filenameformat2[200] = {"Unfolding_%sJP%s_BT170928_RcVzW_12JetBinv2_McPt02_embed%s_Baye5.root"};

	char filetag[Nunf][100] = {"JP2","JP1","JP0","NFweighted"};
	char filetag1[Nunf][100] = {"JP2","JP1","JP0","MB"};

	double drawXrange[2] = {4,60}; //0,60};
	
	int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue,kMagenta,kPink-7};
	//int markerstyle[7] = {8,8,21,21,34,33,29};
	int markerstyle[7] = {4,8,34,21,34,33,29};
	float markersize[7] = {2,2,2,1,1,1,1};

	char histname[200] = {"hreco"};//"hmeas"};

	char histy[2][200] = {"<N/d#etad#phi>","<p_{T}>"};
	//char histy[Nvar][200] = {"Transverse Multiplicity","Transverse <p_{T}>","Transverse Sum p_{T}","Towards Multiplicity","Towards Track <p_{T}>","Towards Track Sum p_{T}","Away Multiplicity","Away Track <p_{T}>","Away Track Sum p_{T}"};
	char histx[200] = {"leading jet p_{T}"};

	double xmax = 45;  //60; // 40;//35;			// maximum x-axis range to draw

	TFile *fin[Nvar][Nunf];
	TProfile *profile[Nvar][Nunf];
	TH2D *h2d[Nvar][Nunf];
	for(int i = 0; i<Nvar; i++) {
		for(int j = 0; j<Nunf; j++) {
			if(j==(Nunf-1)) {
				fin[i][j] = new TFile(Form(filenameformat2,filenametag[i],unfoldtag[j],filetag1[j]));		// using RcVzW for NFweight one, other does not have it..
			}
			else {
				fin[i][j] = new TFile(Form(filenameformat,filenametag[i],unfoldtag[j],filetag1[j]));
			}
			h2d[i][j] = (TH2D*)fin[i][j]->Get(histname);
			h2d[i][j]->SetName(Form("h2d_%sJP%s",filenametag[i],unfoldtag[j]));
			profile[i][j] = (TProfile*)h2d[i][j]->ProfileX(Form("profile%sJP%s",filenametag[i],unfoldtag[j]));
			profile[i][j]->SetMarkerColor(filecolor[j]);
			profile[i][j]->SetLineColor(filecolor[j]);
			profile[i][j]->GetXaxis()->SetRangeUser(drawXrange[0],drawXrange[1]);
			profile[i][j]->GetXaxis()->SetTitle(Form("Unfolded %s",histx));
			profile[i][j]->GetYaxis()->SetTitle(Form("Unfolded %s",histy[i]));
			profile[i][j]->SetMarkerStyle(markerstyle[i%3]);
			profile[i][j]->SetMarkerSize(markersize[i%3]);

		}
	}

	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);    // TranTot: 2 area; eta 2; phi pi/3
	for(int i = 0; i<Nvar; i++) {
		for(int j = 0; j<Nunf; j++) {
			if(!(TString(filenametag[i]).Contains("LeadJetNtrk",TString::kIgnoreCase)||TString(filenametag[i]).Contains("PtAve",TString::kIgnoreCase))) {
				profile[i][j]->Scale(DeDpNorma);
			}
		}
	}
	
	TCanvas *c[2];
	for(int i = 0; i<2; i++) {
		c[i] = new TCanvas();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TLegend *l[2];
	for(int i = 0; i<2; i++) {
		if(strstr(filenametag[i], "Tran")) {
			l[i] = new TLegend(0.55,0.2,0.85,0.45);
		}
		else {
			l[i] = new TLegend(0.15,0.6,0.45,0.85);
		}
		l[i]->SetTextFont(42);
		l[i]->SetFillColor(0);
		l[i]->SetBorderSize(0);
		for(int j = 0; j<Nunf; j++) {
			l[i]->AddEntry(profile[i][j],filetag[j],"l");
		}
	}

	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV JP2 R=0.6 FullJet ");
	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	if(strstr(filenametag[0],"TransNeutral")) lat->SetText(0.1,0.94,"pp@200GeV R=0.6 FullJet TransNeutral NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


	// Scaling the tranphi30 one -- CAUTION!!! the error is not properly handled right now. Need to work on this. 
	//profile[Nfile-1][0]->Scale(60./30);

	for(int i = 0; i<2; i++) {
		c[i]->cd();
		TH1D *htmp = (TH1D*)profile[i][0]->Clone(Form("htmp%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,xmax);
		//htmp->SetMaximum(profile[i][0]->GetMaximum()*1.1);
		htmp->SetMaximum(2.8);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(Form("Unfolded %s",histx));
		htmp->GetYaxis()->SetTitle(Form("Unfolded %s", histy[i]));
		htmp->Draw("");
		for(int j = 0; j<Nunf; j++) {
			for(int k=0;k<3;k++) {
				profile[i*3+k][j]->Draw("eX0same");
			}
		}
		l[i]->Draw("same");
		lat->Draw("same");

		if(savefig) {
			TString outtag="";
			if(i==0) outtag = "Ntrk";
			else outtag = "PtAve";
			if(strstr(filenametag[0],"TransNeutral")) c[i]->SaveAs(Form("CompareUnfoldTrigs_all_TransNeutral_%s.png",outtag.Data()));
			else c[i]->SaveAs(Form("CompareUnfoldTrigs_all_%s.png",outtag.Data()));
		}
	}




}



