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


void plotcompareMaxTrackMBvsJP() {

	int savefig = 0;

	const int Nfile = 2;
	//char filename[Nfile][200] = {"maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.root","maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.root"};
	//char filename[Nfile][200] = {"maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppJP_160811P12id_R06_HadrCorr_170418_WideBin_CDFcutPT05eta08.root","maxtrackpthist4MaxTrackNoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_170418_WideBin_CDFcutPT05eta08.root"};
	//char filename[Nfile][200] = {"maxtrackpthist4FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_WideBin_PT05eta1.root","maxtrackpthist4FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_170418_WideBin_PT05eta1.root"};
	char filename[Nfile][200] = {"maxtrackpthist4FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_WideBin.root","maxtrackpthist4FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_170418_WideBin.root"};
	char filetag[Nfile][100] = {"JP","MB"};
	char filetag4fig[Nfile][100] = {"JP","MB"};
	//int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue,kMagenta,kPink-7};
	int filecolor[7] = {kBlack,kBlue,kOrange+3,kGreen+1,kBlue,kMagenta,kPink-7};
	//int markerstyle[7] = {24,8,21,21,34,33,29};
	int markerstyle[7] = {24,25,30,8,21,29,33};
	float markersize[7] = {1.5,1,1,1,1,1,1};

	double drawXrange[Nfile][2] = {{0,60},{0,60}};

	const int Nhist = 9;		
	char histname[Nhist][200] = {"leadjetareantrkvsmaxtrackpt","subjetareantrkvsmaxtrackpt","trantotntrkvsmaxtrackpt","leadjetareaptavevsmaxtrackpt","subjetareaptavevsmaxtrackpt","tranptavevsmaxtrackpt","leadjetareaptsumvsmaxtrackpt","subjetareaptsumvsmaxtrackpt","trantotptsumvsmaxtrackpt"};
	char histtag[Nhist][200] = {"Towards","Away","Transverse","Towards","Away","Transverse","Towards","Away","Transverse"};
	char histy[Nhist][200] = {"Towards Region Multiplicity","Away Region Multiplicity","Transverse Region Multiplicity","Towards Region Track <p_{T}>","Away Region Track <p_{T}>","Transverse Region <p_{T}>","Towards Region Track Sum p_{T}","Away Region Track Sum p_{T}","Transverse Region Sum p_{T}"};
	//char histy[Nhist][200] = {"Towards Tower Multiplicity","Away Tower Multiplicity","Transverse Tower Multiplicity","Towards Region Track <E_{T}>","Away Region Track <E_{T}>","Transverse Region <E_{T}>","Towards Region Track Sum E_{T}","Away Region Track Sum E_{T}","Transverse Region Sum E_{T}"};
	char histx[200] = {"Detector-level Leading Charged Particle p_{T}"};
	
	double xmax = 40; // 40;//35;			// maximum x-axis range to draw

	const int Nfig=3;
	//char figy[Nfig][200] = {"Detector-level <dN_{ch}/d#etad#phi>","Detector-level <p_{T}^{ch}>","Detector-level #sump_{T}^{ch}"};
	char figy[Nfig][200] = {"Efficiency Corrected Detector-level <dN_{ch}/d#etad#phi>","Efficiency Corrected Detector-level <p_{T}^{ch}>","Efficiency Corrected Detector-level #sump_{T}^{ch}"};


	TFile *fin[Nfile];
	TProfile *profile[Nfile][Nhist];
	double DeDpNorma[Nhist] = {1./(2.*2.*TMath::Pi()/3.),1,1./(2.*2.*TMath::Pi()/3.)};
	//TH1D *profile[Nfile][Nhist];
	for(int i = 0; i<Nfile; i++) {
		fin[i] = new TFile(filename[i]);
		for(int j = 0; j<Nhist; j++) {
			///TH2D *htmp = (TH2D*)fin[i]->Get(histname[j]);
			///profile[i][j]  = htmp->ProfileX(Form("profile%d%d",i,j));
			profile[i][j] = (TProfile*)fin[i]->Get(histname[j]);
			//profile[i][j]->SetMarkerColor(filecolor[i]);
			//profile[i][j]->SetLineColor(filecolor[i]);
			profile[i][j]->GetXaxis()->SetRangeUser(0,xmax);
			//profile[i][j]->GetXaxis()->SetRangeUser(drawXrange[i][0],drawXrange[i][1]);
			profile[i][j]->GetXaxis()->SetTitle(histx);
			profile[i][j]->GetYaxis()->SetTitle(histy[j]);
			//profile[i][j]->SetMarkerStyle(markerstyle[i]);
			//profile[i][j]->SetMarkerSize(markersize[i]);
			profile[i][j]->Scale(DeDpNorma[j/3]);
			//profile[i][j]->Scale(prescale[i]);	// for leading jet pt distribution TH1D
			//if(i==0) {
			//	profile[i][j]->SetMarkerSize(1.5);
			//	profile[i][j]->SetMarkerStyle(4);
			//}
		}
	}
	
	TCanvas *c[Nfig];
	for(int i = 0; i<Nfig; i++) {
		c[i] = new TCanvas();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TLegend *l[Nfig][Nfile];
	for(int i = 0; i<Nfig; i++) {
		for(int j = 0; j<Nfile; j++) {
			l[i][j] = new TLegend(0.5+j*0.2,0.6,0.7+j*0.2,0.85);
			l[i][j]->SetTextFont(42);
			l[i][j]->SetFillColor(0);
			l[i][j]->SetBorderSize(0);
			l[i][j]->SetHeader(filetag[j]);
			for(int k=0; k<3; k++) {
				l[i][j]->AddEntry(profile[j][i*3+k],histtag[k],"pl");
			}
		}
	}

	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV TofMatch |#eta|<1 p_{T} > 0.2 GeV/#it{c}");
	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	if(strstr(filename[0],"TransNeutral")) lat->SetText(0.1,0.94,"pp@200GeV R=0.6 FullJet TransNeutral NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


	// Scaling the tranphi30 one -- CAUTION!!! the error is not properly handled right now. Need to work on this. 
	//profile[Nfile-1][0]->Scale(60./30);

	for(int i = 0; i<Nfig; i++) {
		c[i]->cd();
		TH1D *htmp = (TH1D*)profile[0][i]->Clone(Form("htmp%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,xmax);
		htmp->SetMaximum(profile[0][i*3]->GetMaximum()*1.1);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(histx);
		htmp->GetYaxis()->SetTitle(figy[i]);
		htmp->Draw();
		for(int j = 0; j<Nfile; j++) {
			for(int k=0; k<3; k++) {
				profile[j][i*3+k]->SetMarkerStyle(markerstyle[j*3+k]);
				profile[j][i*3+k]->SetMarkerSize(1.5);
				profile[j][i*3+k]->SetMarkerColor(filecolor[k]);
				profile[j][i*3+k]->SetLineColor(filecolor[k]);
				profile[j][i*3+k]->Draw("same");
			}
			l[i][j]->Draw("same");
		}
		lat->Draw("same");

		if(savefig) {
			c[i]->SaveAs(Form("Compare_%s.png",histname[i]));
		}
	}


	TCanvas *ct[Nfig];
	for(int i = 0; i<Nfig; i++) {
		ct[i] = new TCanvas();
	}
	
	TLegend *lt[Nfig];
	for(int i = 0; i<Nfig; i++) {
		lt[i] = new TLegend(0.6,0.6,0.87,0.87);
		lt[i]->SetHeader(histtag[2]);
		for(int j = 0; j<Nfile; j++) {
			lt[i]->SetTextFont(42);
			lt[i]->SetFillColor(0);
			lt[i]->SetBorderSize(0);
			lt[i]->AddEntry(profile[j][i*3+2],filetag[j],"pl");
		}
	}
	for(int i = 0; i<Nfig; i++) {
		ct[i]->cd();
		TH1D *htmp = (TH1D*)profile[0][i]->Clone(Form("htmptran%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,20);
		htmp->SetMaximum(profile[0][i*3+2]->GetMaximum()*1.5);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(histx);
		htmp->GetYaxis()->SetTitle(figy[i]);
		htmp->Draw();
		profile[1][i*3+2]->GetXaxis()->SetRangeUser(0,10);
		for(int j = 0; j<Nfile; j++) {
			profile[j][i*3+2]->Draw("eX0same");
		}
		lt[i]->Draw("same");
		lat->Draw("same");
	}



}



