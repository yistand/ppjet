#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"


float foldpi(float phi) {
	double pi = TMath::Pi();
	while(phi>pi*(1+60./180.) || phi<-pi*(2*60./180.)) {
		if(phi>pi) phi=phi-2*pi;
		else phi = phi+2*pi;
	}
	return phi;
}


void plotDphi(float ptmin = 10, float ptmax = 50, float ymax = 6) {
	TFile *f = new TFile("NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170127.root");
	TTree *t = (TTree*)f->Get("ResultTree");
	float jphi, jpt;
	int nlead, nsub, nmax, nmin;
	const int MAX = 1000;
	float philead[MAX], phisub[MAX], phimax[MAX], phimin[MAX];
	t->SetBranchAddress("j1pt",&jpt);
	t->SetBranchAddress("j1phi",&jphi);
	t->SetBranchAddress("LeadAreaNtrk",&nlead);
	t->SetBranchAddress("SubAreaNtrk",&nsub);
	t->SetBranchAddress("TranMaxNtrk",&nmax);
	t->SetBranchAddress("TranMinNtrk",&nmin);
	t->SetBranchAddress("TrkLeadAreaPhi",philead);
	t->SetBranchAddress("TrkSubAreaPhi",phisub);
	t->SetBranchAddress("TrkTranMaxPhi",phimax);
	t->SetBranchAddress("TrkTranMinPhi",phimin);

	TH1D *h = new TH1D("h","#Delta#phi = #phi-#phi_{leading jet}",60,-TMath::Pi()*2./3., TMath::Pi()*4./3.);

	int count = 0;
	for(int i = 0 ;i<t->GetEntries(); i++) {
		t->GetEntry(i);
		if(jpt>ptmax||jpt<ptmin) continue;
		for(int j = 0; j<nlead; j++) {
			h->Fill(foldpi(philead[j]-jphi));
		}
		for(int j = 0; j<nsub; j++) {
			h->Fill(foldpi(phisub[j]-jphi));
		}
		for(int j = 0; j<nmax; j++) {
			h->Fill(foldpi(phimax[j]-jphi));
		}
		for(int j = 0; j<nmin; j++) {
			h->Fill(foldpi(phimin[j]-jphi));
		}
		count++;
	}

	float norm=0;
	float dphi = h->GetBinWidth(1);
	float deta = 2.;
       	if(count>0) {
		norm = 1./(count*dphi*deta);
	}
	h->Scale(norm);

	float ymin = 0 ;
	h->SetMaximum(ymax);
	h->SetMinimum(ymin);

	TCanvas *c1 = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	h->GetXaxis()->SetTitle("#Delta#phi");
	h->GetYaxis()->SetTitle("#frac{N_{particle}}{N_{evt}#Delta#phi#Delta#eta}");

	h->SetLineColor(1);
	h->SetLineWidth(2);
	h->Draw("HIST");

	TLatex *lat = new TLatex(0.6,0.8,Form("%g < p_{T} < %g GeV/#it{c}",ptmin, ptmax));
	lat->SetTextFont(42);
	lat->SetNDC();
	lat->DrawClone("same");

	TLine *l = new TLine();
	l->SetLineStyle(7);
	l->SetLineColor(4);
	l->DrawLine(TMath::Pi()*60./180.,ymin,TMath::Pi()*60./180.,ymax);
	l->DrawLine(-TMath::Pi()*60./180.,ymin,-TMath::Pi()*60./180.,ymax);
	l->DrawLine(TMath::Pi()*120./180.,ymin,TMath::Pi()*120./180.,ymax);
	l->DrawLine(-TMath::Pi()*120./180.,ymin,-TMath::Pi()*120./180.,ymax);

	lat->SetTextAngle(15);
	lat->SetTextColor(4);
	lat->SetNDC(0);
	lat->DrawLatex(-2,ymax,"Transverse");
	lat->DrawLatex(1.2,ymax,"Transverse");
	lat->DrawLatex(-0.5,ymax,"Towards");
	lat->DrawLatex(2.5,ymax,"Away");

}


