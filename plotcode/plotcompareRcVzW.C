//=====================================================================================================
//
//	2017.11.02	Li YI
//	compare RcVzWeight or not 
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

#include "plotcompare30Vs60.C"


void plotcompareRcVzW() {

	int savefig = 0;

	const int Nvar = 9;
	char filenametag[Nvar][20] = {"TranTotNtrk","TranPtAve","TranTotPtSum","LeadAreaNtrk","LeadPtAve","LeadAreaPtSum","SubAreaNtrk","SubPtAve","SubAreaPtSum"};

	const int Nrcv = 2;
	char unfoldtag[Nrcv][20] = {"_RcVzW",""};

	char unfoldnametag[Nrcv][20] = {"RcVzWeighted","No RcVzW"};

	char filenameformat[200] = {"Unfolding_%sJPCharged_NFWeight_BT170928%s_12JetBinv2_McPt02_embedMB_Baye5.root"};

	double drawXrange[2] = {4,60}; //0,60};
	
	int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue,kMagenta,kPink-7};
	int markerstyle[7] = {8,8,21,21,34,33,29};
	float markersize[7] = {1,1,1,1,1,1,1};

	char histname[200] = {"hreco"};//"hmeas"};

	char histy[Nvar][200] = {"Transverse Multiplicity","Transverse <p_{T}>","Transverse Sum p_{T}","Towards Multiplicity","Towards Track <p_{T}>","Towards Track Sum p_{T}","Away Multiplicity","Away Track <p_{T}>","Away Track Sum p_{T}"};
	char histx[200] = {"leading jet p_{T}"};

	double xmax = 60; // 40;//35;			// maximum x-axis range to draw
	double xmin = 4;

	TFile *fin[Nvar][Nrcv];
	TProfile *profile[Nvar][Nrcv];
	TH2D *h2d[Nvar][Nrcv];
	for(int i = 0; i<Nvar; i++) {
		for(int j = 0; j<Nrcv; j++) {
			fin[i][j] = new TFile(Form(filenameformat,filenametag[i],unfoldtag[j]));
			h2d[i][j] = (TH2D*)fin[i][j]->Get(histname);
			h2d[i][j]->SetName(Form("h2d_%sJP%s",filenametag[i],unfoldtag[j]));
			profile[i][j] = (TProfile*)h2d[i][j]->ProfileX(Form("profile%sJP%s",filenametag[i],unfoldtag[j]));
			profile[i][j]->SetMarkerColor(filecolor[j]);
			profile[i][j]->SetLineColor(filecolor[j]);
			profile[i][j]->GetXaxis()->SetRangeUser(drawXrange[0],drawXrange[1]);
			profile[i][j]->GetXaxis()->SetTitle(Form("Unfolded %s",histx));
			profile[i][j]->GetYaxis()->SetTitle(Form("Unfolded %s",histy[i]));
			profile[i][j]->SetMarkerStyle(markerstyle[j]);
			profile[i][j]->SetMarkerSize(markersize[j]);

		}
	}

	TH1D *hratio[Nvar];	
	for(int i = 0; i<Nvar; i++) {
		hratio[i] = Ratio2Profile(profile[i][0],profile[i][1]);
		hratio[i]->SetMarkerStyle(8);
	}
	
	TCanvas *c[Nvar];
	TCanvas *cratio[Nvar];
	for(int i = 0; i<Nvar; i++) {
		c[i] = new TCanvas();
	}
	for(int i = 0; i<Nvar; i++) {
		cratio[i] = new TCanvas();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TLegend *l[Nvar];
	for(int i = 0; i<Nvar; i++) {
		if(strstr(filenametag[i], "Tran")) {
			l[i] = new TLegend(0.55,0.2,0.85,0.45);
		}
		else {
			l[i] = new TLegend(0.15,0.6,0.45,0.85);
		}
		l[i]->SetTextFont(42);
		l[i]->SetFillColor(0);
		l[i]->SetBorderSize(0);
		for(int j = 0; j<Nrcv; j++) {
			l[i]->AddEntry(profile[i][j],unfoldnametag[j],"pl");
		}
	}

	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV JP2 R=0.6 FullJet ");
	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	if(strstr(filenametag[0],"TransNeutral")) lat->SetText(0.1,0.94,"pp@200GeV R=0.6 FullJet TransNeutral NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


	// Scaling the tranphi30 one -- CAUTION!!! the error is not properly handled right now. Need to work on this. 
	//profile[Nfile-1][0]->Scale(60./30);

	for(int i = 0; i<Nvar; i++) {
		c[i]->cd();
		TH1D *htmp = (TH1D*)profile[i][0]->Clone(Form("htmp%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,xmax);
		htmp->SetMaximum(profile[i][0]->GetMaximum()*1.1);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(Form("Unfolded %s",histx));
		htmp->GetYaxis()->SetTitle(Form("Unfolded %s", histy[i]));
		htmp->Draw();
		for(int j = 0; j<Nrcv; j++) {
			profile[i][j]->Draw("same");
		}
		l[i]->Draw("same");
		lat->Draw("same");

		if(savefig) {
			if(strstr(filenametag[0],"TransNeutral")) c[i]->SaveAs(Form("CompareRcVzw_TransNeutral_%s.png",filenametag[i]));
			else c[i]->SaveAs(Form("CompareRcVzw%s.png",filenametag[i]));
		}
	}


	TLine *line = new TLine(0,1,xmax+hratio[0]->GetBinWidth(hratio[0]->GetNbinsX())/2.,1);

	for(int i = 0; i<Nvar; i++) {
		cratio[i]->cd();
		TH1D *htmp2 = (TH1D*)profile[i][0]->Clone(Form("htmp2%d",i));
		htmp2->Reset();
		htmp2->GetXaxis()->SetRangeUser(0,xmax);
		htmp2->SetMaximum(1.5);
		htmp2->SetMinimum(0.5);
		htmp2->GetXaxis()->SetTitle(Form("Unfolded %s",histx));
		htmp2->GetYaxis()->SetTitle(Form("#frac{%s}{%s} Ratio for Unfolded %s", unfoldnametag[0],unfoldnametag[1],histy[i]));
		htmp2->Draw();
		line->Draw("same");
		hratio[i]->GetXaxis()->SetRangeUser(xmin,xmax);
		hratio[i]->Draw("same");
		lat->Draw("same");

		if(savefig) {
			if(strstr(filenametag[0],"TransNeutral")) cratio[i]->SaveAs(Form("RatioRcVzw_TransNeutral_%s.png",filenametag[i]));
			else cratio[i]->SaveAs(Form("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/Ratio_%s_RcVzw.pdf",filenametag[i]));
		}
	}




}



