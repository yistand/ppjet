//=====================================================================================================
//
//
//	2016.07.13	Li YI
//	modified from plotcompare.C
//	compare two TProfile and do the ratio
//
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


void plotCompareTwo() {

	int savefig = 0;

	const int Nfile = 2;

	//char filename[Nfile][200] = {"leadjetpthist4FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160314_NoEffCorr.root","leadjetpthist4BemcOrTofMatch_FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_NoEffCorr_160712.root"};
	//char filetag[Nfile][100] = {"TOF Match","TOF||BEMC Match"};
	//char filetag4fig[Nfile][100] = {"TofOnly","BemcOrTof"};
	char filename[Nfile][200] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160711.root","leadjetpthist4FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160712.root"};
	char filetag[Nfile][100] = {"No Match","TOF Match"};
	char filetag4fig[Nfile][100] = {"NoTofMatch","TofMatch"};
	int filecolor[7] = {kBlack,kRed,kOrange-2,kGreen+1,kBlue,kMagenta,kPink-7};


	const int Nhist = 2;			// 6;
	
	char histname[Nhist][200] = {"tranntrkvsleadjetpt","tranptavevsleadjetpt"};
	char histx[200] = {"leading jet p_{T}"};

	//char histname[Nhist][200] = {"tranntrkvsmultiplicity","tranptavevsmultiplicity"};
	//char histx[200] = {"Total Uncorrected Charge Multiplicity"};

	//char histname[Nhist][200] = {"tranntrkvstransntrk","tranptavevstransntrk"};
	//char histx[200] = {"Transverse Uncorrected Charge Multiplicity"};
	
	char histy[Nhist][200] = {"Transverse Region Multiplicity","Transverse Region <p_{T}>"};
	char histtag[Nhist][200] = {"Transverse","Transverse"};

	double xmax = 40;//35;			// maximum x-axis range to draw

	TFile *fin[Nfile];
	TProfile *profile[Nfile][Nhist];
	for(int i = 0; i<Nfile; i++) {
		fin[i] = new TFile(filename[i]);
		for(int j = 0; j<Nhist; j++) {
			profile[i][j] = (TProfile*)fin[i]->Get(histname[j]);
			profile[i][j]->SetMarkerColor(filecolor[i]);
			profile[i][j]->SetLineColor(filecolor[i]);
			profile[i][j]->GetXaxis()->SetRangeUser(0,xmax);
			profile[i][j]->GetXaxis()->SetTitle(histx);
			profile[i][j]->GetYaxis()->SetTitle(histy[j]);
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

	TLatex *lat = new TLatex(0.05,0.94,"pp@200GeV JP2 R=0.6 FullJet, w/ Tracking/Matching eff. corr.");
	lat->SetNDC();
	lat->SetTextFont(42);



	for(int i = 0; i<Nhist; i++) {
		c[i]->cd();
		profile[0][i]->Draw();
		for(int j = 1; j<Nfile; j++) {
			profile[j][i]->Draw("same");
		}
		l[i]->Draw("same");
		lat->Draw("same");

		if(savefig) c[i]->SaveAs(Form("Compare%sVs%s_%s.png",filetag[0],filetag[1],histname[i]));
	}


	// do the ratio
	TH1D *hratio[Nhist];
	for(int i = 0; i<Nhist; i++) {
		hratio[i] = new TH1D(Form("hratio%s",histname[i]),Form("%s ratio of %s and %s",histname[i],filetag[0],filetag[1]),profile[0][i]->GetNbinsX(),profile[0][i]->GetXaxis()->GetBinLowEdge(1),profile[0][i]->GetXaxis()->GetBinLowEdge(profile[0][i]->GetNbinsX())+profile[0][i]->GetXaxis()->GetBinWidth(profile[0][i]->GetNbinsX()));
		hratio[i]->Reset();
		for(int j = 0; j<hratio[i]->GetNbinsX(); j++) {
			double y0 = profile[0][i]->GetBinContent(j+1);
			double y1 = profile[1][i]->GetBinContent(j+1);
			double e0 = profile[0][i]->GetBinError(j+1);
			double e1 = profile[1][i]->GetBinError(j+1);


			if(y1>0&&y0>0) {
				hratio[i]->SetBinContent(j+1,y0/y1);
				hratio[i]->SetBinError(j+1,sqrt(pow(e0/y0,2)+pow(e1/y1,2))*y0/y1);
			}
			else {
				hratio[i]->SetBinContent(j+1,0);
				hratio[i]->SetBinError(j+1,0);
			}

		}
		hratio[i]->SetLineColor(i+1);
		hratio[i]->GetYaxis()->SetTitle(Form("%s #frac{%s}{%s}",profile[0][i]->GetYaxis()->GetTitle(),filetag[0],filetag[1]));
		hratio[i]->GetXaxis()->SetRangeUser(0,xmax);
	}

	TCanvas *cr[Nhist]; 
	TLine *line = new TLine(0,1,xmax,1);
	for(int i = 0; i<Nhist; i++) {
		cr[i] = new TCanvas();

		hratio[i]->Draw("h");
		lat->Draw("same");
		line->Draw("same");

		if(savefig) cr[i]->SaveAs(Form("Ratio%sR%s_%s.png",filetag[0],filetag[1],histname[i]));
		
	}



}



