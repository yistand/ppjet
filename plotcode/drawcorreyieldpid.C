//====================================================================================================
//
//	2018.01.21	Li YI
//	quickly draw comparison between pid calculated and published one
//
//====================================================================================================


void drawcorreyieldpid(TString filename="corryieldpid.root") {
	const int NSpecies = 3;

	TH1D *h[NSpecies];
	TGraphErrors *gr[NSpecies];
	TString grname[NSpecies] = {"pion","kaon","proton"};
	int color[NSpecies] = {kRed,kBlue,kGreen+3};
	int style[NSpecies] = {8,21,47};

	TFile *fin = new TFile(filename);

	for(int i = 0; i<NSpecies; i++) {
		h[i] = (TH1D*)fin->Get(Form("hCorrYieldspid%d",i));
		gr[i] = (TGraphErrors*)fin->Get("gsum_"+grname[i]);

		h[i]->GetXaxis()->SetTitle("p_{T}");
		h[i]->GetYaxis()->SetTitle("#frac{d^{2}N}{2#pip_{T}dp_{T}d#eta}");

		h[i]->SetLineColor(color[i]);
		h[i]->SetMarkerColor(color[i]);
		h[i]->SetMarkerStyle(style[i]);
		gr[i]->SetLineColor(color[i]);
	}

	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	c->SetLogy();

	h[0]->Draw();
	for(int i = 0; i<NSpecies; i++) {
		h[i]->Draw("same");
		gr[i]->Draw("same");
	}

	TLegend *leg = new TLegend(0.7,0.6,0.85,0.85);
	leg->SetBorderSize(0);
	for(int i = 0; i<NSpecies; i++) {
		leg->AddEntry(h[i],grname[i],"p");
	}
	leg->Draw();
}
