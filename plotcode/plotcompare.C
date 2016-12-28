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


void plotcompare() {

	int savefig = 0;

	const int Nfile = 3;
	char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin_reweighted.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root"};
	char filetag[Nfile][100] = {"JP reweighted w/ VPD","JP w/ VPD","MB w/ VPD"};
	char filetag4fig[Nfile][100] = {"reJP","JP","MB"};
	double drawXrange[Nfile][2] = {{0,60},{0,60},{0,60}};

	//const int Nfile = 3;
	//char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_161209_NoEffCorr_WideBin_reweighted.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_161209_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_161209_NoEffCorr_WideBin.root"};
	//char filetag[Nfile][100] = {"JP reweighted w/o VPD","JP w/o VPD","VPDMB"};
	//char filetag4fig[Nfile][100] = {"reJP","JP","MB"};
	//double drawXrange[Nfile][2] = {{0,60},{0,60},{0,60}};


	int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue,kMagenta,kPink-7};
	int markerstyle[7] = {24,8,21,21,34,33,29};
	float markersize[7] = {1.5,1,1,1,1,1,1};


	//const int Nfile = 6;
	//char filename[Nfile][200] = {"leadjetpthist4FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","leadjetpthist4kT_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","leadjetpthist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","leadjetpthist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224_AsJetGt5.root","leadjetpthist4FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224_wR1.root", "leadjetpthist4TranPhi30_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root"};
	//char filename[Nfile][200] = {"multiplicityhist4FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","multiplicityhist4kT_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","multiplicityhist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","multiplicityhist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224_AsJetGt5.root","multiplicityhist4TranPhi30_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root"};
	//char filename[Nfile][200] = {"transntrkhist4FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","transntrkhist4kT_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","transntrkhist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root","transntrkhist4Dijet_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224_AsJetGt5.root","transntrkhist4FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224_wR1.root","transntrkhist4TranPhi30_FullJet_TransCharged_MatchTrig_ppJP2_R06_HadrCorr_160224.root"};
	//char filetag[Nfile][100] = {"default","fastjet k_{T}","Dijet angle","Dijet angle w/ recoil jet>5GeV/c","jet R=1","Transverse size < 30"};
	//char filetag4fig[Nfile][100] = {"default","kT","Dijet","DijetAsGt5","jetR1","TranPhi30"};
	//int filecolor[7] = {kBlack,kRed,kOrange-2,kGreen+1,kBlue,kMagenta,kPink-7};


	const int Nhist = 9;		
	char histname[Nhist][200] = {"leadjetareantrkvsleadjetpt","subjetareantrkvsleadjetpt","tranntrkvsleadjetpt","leadjetareaptavevsleadjetpt","subjetareaptavevsleadjetpt","tranptavevsleadjetpt","leadjetareaptsumvsleadjetpt","subjetareaptsumvsleadjetpt","tranptsumvsleadjetpt"};
	char histtag[Nhist][200] = {"Towards","Away","Transverse","Towards","Away","Transverse","Towards","Away","Transverse"};
	char histy[Nhist][200] = {"Towards Region Multiplicity","Away Region Multiplicity","Transverse Region Multiplicity","Towards Region Track <p_{T}>","Away Region Track <p_{T}>","Transverse Region <p_{T}>","Towards Region Track Sum p_{T}","Away Region Track Sum p_{T}","Transverse Region Sum p_{T}"};
	char histx[200] = {"leading jet p_{T}"};
	
	//const int Nhist = 1;		
	//char histname[Nhist][200] = {"leadjetpt"};
	//char histtag[Nhist][200] = {"NumberOfEvents"};
	//char histy[Nhist][200] = {"Number of Events"};
	//char histx[200] = {"leading jet p_{T}"};
	
	//const int Nhist = 2;		
	//char histname[Nhist][200] = {"tranntrkvsleadjetpt","tranptavevsleadjetpt"};
	///char histname[Nhist][200] = {"htranntrkvsleadjetpt","htranptavevsleadjetpt"};
	//char histx[200] = {"leading jet p_{T}"};

	//char histname[Nhist][200] = {"tranntrkvsmultiplicity","tranptavevsmultiplicity"};
	//char histx[200] = {"Total Uncorrected Charge Multiplicity"};

	//char histname[Nhist][200] = {"tranntrkvstransntrk","tranptavevstransntrk"};
	//char histx[200] = {"Transverse Uncorrected Charge Multiplicity"};
	
	//char histy[Nhist][200] = {"Transverse Region Multiplicity","Transverse Region <p_{T}>"};
	//char histtag[Nhist][200] = {"Transverse","Transverse"};

	double xmax = 60; // 40;//35;			// maximum x-axis range to draw

	TFile *fin[Nfile];
	TProfile *profile[Nfile][Nhist];
	//TH1D *profile[Nfile][Nhist];
	//double prescale[Nfile] = {1,2.5,141.35,110.52};
	for(int i = 0; i<Nfile; i++) {
		fin[i] = new TFile(filename[i]);
		for(int j = 0; j<Nhist; j++) {
			///TH2D *htmp = (TH2D*)fin[i]->Get(histname[j]);
			///profile[i][j]  = htmp->ProfileX(Form("profile%d%d",i,j));
			profile[i][j] = (TProfile*)fin[i]->Get(histname[j]);
			profile[i][j]->SetMarkerColor(filecolor[i]);
			profile[i][j]->SetLineColor(filecolor[i]);
			profile[i][j]->GetXaxis()->SetRangeUser(0,xmax);
			//profile[i][j]->GetXaxis()->SetRangeUser(drawXrange[i][0],drawXrange[i][1]);
			profile[i][j]->GetXaxis()->SetTitle(histx);
			profile[i][j]->GetYaxis()->SetTitle(histy[j]);
			profile[i][j]->SetMarkerStyle(markerstyle[i]);
			profile[i][j]->SetMarkerSize(markersize[i]);
			//profile[i][j]->Scale(prescale[i]);	// for leading jet pt distribution TH1D
			//if(i==0) {
			//	profile[i][j]->SetMarkerSize(1.5);
			//	profile[i][j]->SetMarkerStyle(4);
			//}
		}
	}
	
	TCanvas *c[Nhist];
	for(int i = 0; i<Nhist; i++) {
		c[i] = new TCanvas();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TLegend *l[Nhist];
	for(int i = 0; i<Nhist; i++) {
		l[i] = new TLegend(0.15,0.6,0.45,0.85);
		l[i]->SetTextFont(42);
		l[i]->SetFillColor(0);
		l[i]->SetBorderSize(0);
		for(int j = 0; j<Nfile; j++) {
			l[i]->AddEntry(profile[j][i],filetag[j],"pl");
		}
	}

	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV JP2 R=0.6 FullJet ");
	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


	// Scaling the tranphi30 one -- CAUTION!!! the error is not properly handled right now. Need to work on this. 
	//profile[Nfile-1][0]->Scale(60./30);

	for(int i = 0; i<Nhist; i++) {
		c[i]->cd();
		TH1D *htmp = (TH1D*)profile[0][i]->Clone(Form("htmp%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,xmax);
		htmp->SetMaximum(profile[0][i]->GetMaximum()*1.1);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(histx);
		htmp->GetYaxis()->SetTitle(histy[i]);
		htmp->Draw();
		for(int j = 0; j<Nfile; j++) {
			profile[j][i]->Draw("same");
		}
		l[i]->Draw("same");
		lat->Draw("same");

		if(savefig) c[i]->SaveAs(Form("Compare_%s.png",histname[i]));
	}




}



