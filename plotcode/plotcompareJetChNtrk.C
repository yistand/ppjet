//=====================================================================================================
//
//	2017.07.18	Li YI
//	modified from plotcompare.C
//	to check MB vs JP (after NF reweighted) jet charged multiplicity and compare with embedding (MB, no VPD)
//
//=====================================================================================================
//
//	2016.02.26	Li YI
//
//	plot the various histograms obtained from plotpT_fromJetTree.C togehter on the same canvas
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


void plotcompareJetChNtrk() {

	int savefig = 0;

	const int Nfile = 3;
	char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin_reweighted.root","Unfolding_LeadJetChargeNtrkJPCharged_NFweight_McPt02_embedMB_Baye5.root"};
	//char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransNeutral_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin_reweighted.root","leadjetpthist4NoTofMatch_FullJet_TransNeutral_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransNeutral_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root"};
	char filetag[Nfile][100] = {"VPDMB","JP reweighted","Embedding (MB)"};
	char filetag4fig[Nfile][100] = {"MB","reJP","Embed"};
	int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue,kMagenta,kPink-7};
	int markerstyle[7] = {24,8,21,21,34,33,29};
	float markersize[7] = {1.5,1,1,1,1,1,1};

	double drawXrange[Nfile][2] = {{0,60},{0,60},{0,60}};

	const int Nhist = 1;		
	char histname[Nhist][200] = {"jconstchargentrkvsleadjetpt"};
	char alterhistname[Nhist][200] = {"htrain"};
	char histtag[Nhist][200] = {"JetChargeMultiplicity"};
	char histy[Nhist][200] = {"Jet Charge Multiplicity"};
	//char histy[Nhist][200] = {"Towards Tower Multiplicity","Away Tower Multiplicity","Transverse Tower Multiplicity","Towards Region Track <E_{T}>","Away Region Track <E_{T}>","Transverse Region <E_{T}>","Towards Region Track Sum E_{T}","Away Region Track Sum E_{T}","Transverse Region Sum E_{T}"};
	char histx[200] = {"leading jet p_{T}"};
	
	double xmax = 45; // 40;//35;			// maximum x-axis range to draw

	TFile *fin[Nfile];
	TProfile *profile[Nfile][Nhist];
	//TH1D *profile[Nfile][Nhist];
	//double prescale[Nfile] = {1,2.5,141.35,110.52};
	for(int i = 0; i<Nfile; i++) {
		fin[i] = new TFile(filename[i]);
		for(int j = 0; j<Nhist; j++) {
			///TH2D *htmp = (TH2D*)fin[i]->Get(histname[j]);
			///profile[i][j]  = htmp->ProfileX(Form("profile%d%d",i,j));
			cout<<filename[i]<<endl;
			profile[i][j] = (TProfile*)fin[i]->Get(histname[j]);
			if(!profile[i][j]) {
				cout<<alterhistname[j]<<endl;
				TH2D *htmp = (TH2D*)fin[i]->Get(alterhistname[j]);
				profile[i][j]  = htmp->ProfileX(Form("profile%d%d",i,j));
			}
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
		l[i] = new TLegend(0.55,0.2,0.85,0.45);
		l[i]->SetTextFont(42);
		l[i]->SetFillColor(0);
		l[i]->SetBorderSize(0);
		for(int j = 0; j<Nfile-1; j++) {
			l[i]->AddEntry(profile[j][i],filetag[j],"pl");
		}
		l[i]->AddEntry(profile[Nfile-1][i],filetag[Nfile-1],"l");
	}

	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV JP2 R=0.6 FullJet ");
	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	if(strstr(filename[0],"TransNeutral")) lat->SetText(0.1,0.94,"pp@200GeV R=0.6 FullJet TransNeutral NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


	// Scaling the tranphi30 one -- CAUTION!!! the error is not properly handled right now. Need to work on this. 
	//profile[Nfile-1][0]->Scale(60./30);

	for(int i = 0; i<Nhist; i++) {
		c[i]->cd();
		TH1D *htmp = (TH1D*)profile[0][i]->Clone(Form("htmp%d",i));
		htmp->Reset();
		htmp->GetXaxis()->SetRangeUser(0,xmax);
		htmp->SetMaximum(profile[1][i]->GetMaximum()*1.1);
		htmp->SetMinimum(0);
		htmp->GetXaxis()->SetTitle(histx);
		htmp->GetYaxis()->SetTitle(histy[i]);
		htmp->Draw();
		for(int j = 0; j<Nfile-1; j++) {
			profile[j][i]->Draw("same");
		}
		profile[Nfile-1][i]->Draw("HISTsame");
		l[i]->Draw("same");
		lat->Draw("same");

		if(savefig) {
			c[i]->SaveAs(Form("CompareJetChargedMult_170418_ppMBvsJPvsEmbedMB.pdf"));
		}
	}




}



