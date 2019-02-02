#include "plotUnfold.C"		// void SetHistStyle(TH1D *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) ;


void TranDiff() {
	const int Nfile = 2;
	int energy[Nfile] = {200,500};
	const char *filename[Nfile] = {"leadjetpthist4NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_161209_NoEffCorr_WideBin.root","GrantWebb500.root"};
	TFile *f[Nfile];
	const char *maxname[Nfile] = {"tranmaxntrkvsleadjetpt","tran500max"};
	const char *minname[Nfile] = {"tranminntrkvsleadjetpt","tran500min"};
	TProfile *max[Nfile];
	TProfile *min[Nfile];
	TProfile *diff[Nfile];

	for(int i = 0; i<Nfile; i++) {
		f[i] = new TFile(filename[i]);
		max[i] = (TProfile*)f[i]->Get(maxname[i]);
		min[i] = (TProfile*)f[i]->Get(minname[i]);
		diff[i] = (TProfile*)max[i]->Clone(Form("diff%d",energy[i]));
		diff[i]->Add(min[i],-1);
	}

	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	diff[0]->Scale(DeDpNorma);

	TH1D *htmp = new TH1D("htmp","",60,5,65);
	htmp->SetMaximum(0.6);
	htmp->SetMinimum(0);
	htmp->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T} (GeV/#it{c})");
	htmp->GetYaxis()->SetTitle("Detector-level #LTN_{ch}/#delta#eta#delta#phi#GT");
	htmp->GetYaxis()->SetNdivisions(505);
	htmp->GetXaxis()->SetTitleSize(0.05);
	htmp->GetYaxis()->SetTitleSize(0.05);
	htmp->GetXaxis()->SetLabelSize(0.05);
	htmp->GetYaxis()->SetLabelSize(0.05);

	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	int mcolor[Nfile] = {1,2};
	int mstyle[Nfile] = {24,8};
	htmp->Draw();
	diff[0]->GetXaxis()->SetRangeUser(10,65);
	diff[1]->GetXaxis()->SetRangeUser(14,62);
	for(int i = 0; i<Nfile; i++) {
		diff[i]->Draw("same");
		diff[i]->SetMarkerSize(2);
		diff[i]->SetMarkerColor(mcolor[i]);
		diff[i]->SetMarkerStyle(mstyle[i]);
		diff[i]->SetLineColor(mcolor[i]);
	}

	TLegend *leg = new TLegend(0.17,0.57,0.55,0.88);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetHeader("TranMax - TranMin");
	for(int i = Nfile-1; i>=0; i--) {
		leg->AddEntry(diff[i],Form("p+p@%dGeV",energy[i]),"pl");
	}
	leg->Draw();

	TLatex *latex = new TLatex();
	//latex->DrawLatex(40,0.5,"TranMax - TranMin");


}
