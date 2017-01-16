
void SetHistStyle(TH1D *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) {
	h->SetMarkerColor(mcolor);
	h->SetLineColor(lcolor);
	h->SetMarkerStyle(mstyle);
	h->SetLineStyle(lstyle);
	h->SetMarkerSize(msize);
	h->SetLineWidth(lwidth);

	h->GetXaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetNdivisions(510);
	h->GetYaxis()->SetNdivisions(505);
}



TH2D *GetBbB(TH2D* htt, TH2D* ht, TH2D* hm) {	// TrainTrue Train Measurement
	TH2D *hr = (TH2D*)hm->Clone("hrecobbb");
	hr->SetTitle("Bin-By-Bin corrected");

	TH2D *r = (TH2D*)htt->Clone("r");
	r->Divide(ht);

	for(int i = 0; i<hr->GetNbinsX(); i++) {
		for(int j=0; j<hr->GetNbinsY(); j++) {
			hr->SetBinContent(i+1,j+1,hr->GetBinContent(i+1,j+1)*r->GetBinContent(i+1,j+1));
		}
	}

	return hr;
}


void compare2BinByBin(TString filename) {
	
	TFile *f = new TFile(filename);
	if(!f->IsOpen()) {
		cout<<"cannot find "<<filename<<endl;
		return;
	}

	TH2D *htraint = (TH2D*)f->Get("htraintrue");
	TH2D *htrainm = (TH2D*)f->Get("htrain");
	TH2D *hmeas = (TH2D*)f->Get("hmeas");
	TH2D *hreco = (TH2D*)f->Get("hreco");

	TH2D *hbbb = (TH2D*)GetBbB(htraint,htrainm,hmeas);

	TH1D* hmeaspx = (TH1D*)hmeas->ProfileX("hmeaspx");
	TH1D* hrecopx = (TH1D*)hreco->ProfileX("hrecopx");

	TH1D* htrainmpx = (TH1D*)htrainm->ProfileX("htrainmpx");
	TH1D* htraintpx = (TH1D*)htraint->ProfileX("htraintpx");

	TH1D* hbbbpx = (TH1D*)hbbb->ProfileX("hbbbpx");

	SetHistStyle(hmeaspx, 2, 2, 25, 1, 2, 1);
	SetHistStyle(hrecopx, 1, 1, 20, 1, 2, 1);

	SetHistStyle(htrainmpx, 2, 2, 25, 1, 2, 3);
	SetHistStyle(htraintpx, 1, 1, 20, 1, 2, 3);

	SetHistStyle(hbbbpx, kAzure+7, kAzure+7, 20, 1, 2, 3);

	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	if(!filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) {
		hmeaspx->Scale(DeDpNorma);
		hrecopx->Scale(DeDpNorma);
		htrainmpx->Scale(DeDpNorma);
		htraintpx->Scale(DeDpNorma);
		hbbbpx->Scale(DeDpNorma);
	}



	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TString ftitle = f->GetName();
	TString sregion = "Transverse";
	if(ftitle.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(ftitle.Contains("Away",TString::kIgnoreCase) || ftitle.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) svariable = " Constituents Multiplicity";
	hmeaspx->GetYaxis()->SetTitle(sregion+svariable);

	hmeaspx->GetXaxis()->SetRangeUser(0,55);
	hmeaspx->SetMinimum(0);
	hmeaspx->SetMaximum(1.5);
	if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hmeaspx->SetMaximum(2.5);//Lead Ntrk density
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) hmeaspx->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)
	
	hmeaspx->Draw("p");
	hrecopx->Draw("psame");
	htrainmpx->Draw("HISTsame");
	htraintpx->Draw("HISTsame");
	hbbbpx->Draw("HISTsame");


	TLegend *leg;
	if(filename.Contains("LeadJetNtrk")) {
		leg = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else {
       		leg = new TLegend(0.5,0.6,0.88,0.88);	
	}
	leg->SetFillColor(0);
	leg->AddEntry(htrainmpx,"MC Detector-level","l");
	leg->AddEntry(htraintpx,"MC Particle-level","l");
	leg->AddEntry(hmeaspx,"Measured data","p");
	leg->AddEntry(hrecopx,"Bayesian Unfolded data","p");
	leg->AddEntry(hbbbpx,"Bin-by-Bin Unfolded data","l");
	leg->Draw();

	TString outtag = sregion;
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) outtag="LeadJetChargedAndNeutralNtrk";	// Lead Jet Ntrk (Charged+Neutral)
	c->SaveAs("compared2bbb_"+outtag+"TrkJPCharged_NFweight_McPt02_embedMB.png");
}
