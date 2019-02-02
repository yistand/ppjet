void comparepid() {
	const int Nfile = 2;
	char filename[Nfile][30] = {"corryieldpid.root","corryieldpid_1Dx.root"};
	int mstyle[Nfile] = {8,24};
	char filetag[Nfile][20] = {"TPC+TOF","TPC only"};

	const int Nhist = 3;
	int color[Nhist] = {kRed,kBlue,kGreen+3};
	char histtag[Nhist][20] = {"pion","kaon","proton"};

	TFile *f[Nfile];
	TH1D *h[Nfile][Nhist];

	for(int i = 0; i<Nfile; i++) {
		f[i] = new TFile(filename[i]);
		cout<<f[i]<<endl;
		for(int j = 0; j<Nhist; j++) {
			h[i][j] = (TH1D*)f[i]->Get(Form("hCorrYieldspid%d",j));
			cout<<h[i][j]<<endl;
			h[i][j]->SetLineColor(color[j]);
			h[i][j]->SetMarkerColor(color[j]);
			h[i][j]->SetMarkerStyle(mstyle[i]);
		}
	}

	TCanvas *c = new TCanvas();
	c->SetLogy();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	h[0][0]->GetXaxis()->SetTitle("p_{T}");
	h[0][0]->GetYaxis()->SetTitle("#frac{d^{2}N}{2#pip_{T}dp_{T}d#eta}");
	h[0][0]->Draw();
	for(int i = 0; i<Nfile; i++) {
		for(int j = 0; j<Nhist; j++) {
			h[i][j]->Draw("same");
		}
	}

	TLegend *l = new TLegend(0.7,0.85,0.7,0.85);
	for(int i = 0; i<Nfile; i++) {
		for(int j = 0; j<Nhist; j++) {
			l->AddEntry(h[i][j],Form("%s %s",filetag[i],histtag[j]),"pl");
		}
	}
	l->Draw();


}
