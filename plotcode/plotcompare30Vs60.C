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

TH1D *Ratio2Profile(TProfile *p1, TProfile *p2) {		// p1/p2
	if(p1->GetNbinsX()!=p2->GetNbinsX()) {
		cout<<"WRONG!!! Ratio2Profile() in plotcompare30Vs60.C: p1 and p2 must have same bins"<<endl;
		return NULL;
	}
	int bins = p1->GetNbinsX();
	const double  *xbins = p1->GetXaxis()->GetXbins()->GetArray();
	TH1D *pr = new TH1D(Form("ratio_%s",p1->GetName()),"ratio",bins,xbins);
	for(int i = 1; i<p1->GetNbinsX()+1; i++) {
		float y1 = p1->GetBinContent(i);
		float e1 = p1->GetBinError(i);
		float y2 = p2->GetBinContent(i);
		float e2 = p2->GetBinError(i);

		float y = 0, e = 0;
		if(y1>0 || y2>0 ) {
			y = (y2>0?y1/y2:0);
			e = sqrt(pow(e1/y1,2)+pow(e2/y2,2))*y;
		}

		pr->SetBinContent(i,y);
		pr->SetBinError(i,e);
	}
	return pr;
}


void plotcompare30Vs60() {

	int savefig = 0;

	const int Nfile = 2;
	char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.root","leadjetpthist4TranPhi30NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.root"};
	char filetag[Nfile][100] = {"#Delta#phi=60","#Delta#phi=30"};
	char filetag4fig[Nfile][100] = {"60","30"};
	double drawXrange[Nfile][2] = {{0,60},{0,60}};

	//const int Nfile = 3;
	//char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin_reweighted.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root"};
	////char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransNeutral_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin_reweighted.root","leadjetpthist4NoTofMatch_FullJet_TransNeutral_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root","leadjetpthist4NoTofMatch_FullJet_TransNeutral_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209_NoEffCorr_WideBin.root"};
	//char filetag[Nfile][100] = {"JP reweighted w/ VPD","JP w/ VPD","MB w/ VPD"};
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
	//char histname[Nhist][200] = {"leadjetareantrkvsleadjetpt","subjetareantrkvsleadjetpt","tranntrkvsleadjetpt","leadjetareaptavevsleadjetpt","subjetareaptavevsleadjetpt","tranptavevsleadjetpt","leadjetareaptsumvsleadjetpt","subjetareaptsumvsleadjetpt","tranptsumvsleadjetpt"};
	char histname[Nhist][200] = {"leadjetareantrkvsleadjetpt","subjetareantrkvsleadjetpt","trantotntrkvsleadjetpt","leadjetareaptavevsleadjetpt","subjetareaptavevsleadjetpt","tranptavevsleadjetpt","leadjetareaptsumvsleadjetpt","subjetareaptsumvsleadjetpt","trantotptsumvsleadjetpt"};
	char histtag[Nhist][200] = {"Towards","Away","Transverse","Towards","Away","Transverse","Towards","Away","Transverse"};
	char histy[Nhist][200] = {"Towards Region Multiplicity","Away Region Multiplicity","Transverse Region Multiplicity","Towards Region Track <p_{T}>","Away Region Track <p_{T}>","Transverse Region <p_{T}>","Towards Region Track Sum p_{T}","Away Region Track Sum p_{T}","Transverse Region Sum p_{T}"};
	//char histy[Nhist][200] = {"Towards Tower Multiplicity","Away Tower Multiplicity","Transverse Tower Multiplicity","Towards Region Track <E_{T}>","Away Region Track <E_{T}>","Transverse Region <E_{T}>","Towards Region Track Sum E_{T}","Away Region Track Sum E_{T}","Transverse Region Sum E_{T}"};
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
	TH1D *rp[Nhist];
	//TH1D *profile[Nfile][Nhist];
	//double prescale[Nfile] = {1,2.5,141.35,110.52};
	for(int i = 0; i<Nfile; i++) {
		fin[i] = new TFile(filename[i]);
		for(int j = 0; j<Nhist; j++) {
			///TH2D *htmp = (TH2D*)fin[i]->Get(histname[j]);
			///profile[i][j]  = htmp->ProfileX(Form("profile%d%d",i,j));
			profile[i][j] = (TProfile*)fin[i]->Get(histname[j]);
			if(TString(filename[i]).Contains("TranPhi30",TString::kIgnoreCase)&&(TString(histname[j]).Contains("ntrk",TString::kIgnoreCase)||TString(histname[j]).Contains("ptsum",TString::kIgnoreCase)))	{
				if(TString(histname[j]).Contains("tran",TString::kIgnoreCase)) {
					profile[i][j]->Scale(2);	// tran phi 30 -> 60
				}
				else {
					profile[i][j]->Scale(60./75);	// toward, away 75->60
				}
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
	for(int j = 0; j<Nhist; j++) { 
		rp[j] = Ratio2Profile(profile[0][j],profile[1][j]);
		rp[j]->SetMarkerStyle(8);
		rp[j]->SetMarkerColor(1);
		rp[j]->SetLineColor(1);
		rp[j]->GetXaxis()->SetTitle(histx);
		rp[j]->GetYaxis()->SetTitle(Form("Ratio 60/30 degree %s", histy[j]));
	}
	
	TCanvas *c[Nhist];
	for(int i = 0; i<Nhist; i++) {
		c[i] = new TCanvas();
	}
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	TLegend *l[Nhist];
	for(int i = 0; i<Nhist; i++) {
		if(strstr(histname[i], "tran")) {
			l[i] = new TLegend(0.55,0.2,0.85,0.45);
		}
		else {
			l[i] = new TLegend(0.15,0.6,0.45,0.85);
		}
		l[i]->SetTextFont(42);
		l[i]->SetFillColor(0);
		l[i]->SetBorderSize(0);
		for(int j = 0; j<Nfile; j++) {
			l[i]->AddEntry(profile[j][i],filetag[j],"pl");
		}
	}

	//TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV JP2 R=0.6 FullJet ");
	TLatex *lat = new TLatex(0.2,0.94,"pp@200GeV R=0.6 FullJet NoTofMatch");
	if(strstr(filename[0],"TransNeutral")) lat->SetText(0.1,0.94,"pp@200GeV R=0.6 FullJet TransNeutral NoTofMatch");
	lat->SetNDC();
	lat->SetTextFont(42);


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

		if(savefig) {
			if(strstr(filename[0],"TransNeutral")) c[i]->SaveAs(Form("Compare_TransNeutral_%s.png",histname[i]));
			else c[i]->SaveAs(Form("Compare_%s_TranPhi30Vs60_NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.pdf",histname[i]));
		}
	}


	TCanvas *cr[Nhist];
        for(int i = 0; i<Nhist; i++) {
                cr[i] = new TCanvas();
        }

	TLine *line = new TLine(0,1,xmax,1);

	for(int i = 0; i<Nhist; i++) {
		cr[i]->cd();
		rp[i]->GetXaxis()->SetRangeUser(0,xmax);
		rp[i]->GetYaxis()->SetRangeUser(0.5,1.5);
		rp[i]->Draw();
		lat->Draw("same");
		line->Draw();
		if(savefig) {
			if(strstr(filename[0],"TransNeutral")) c[i]->SaveAs(Form("Compare_TransNeutral_%s.png",histname[i]));
			else cr[i]->SaveAs(Form("Ratio_%s_TranPhi30Vs60_NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418_NoEffCorr_WideBin.pdf",histname[i]));
		}
	}



}



