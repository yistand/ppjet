int verbose = 0;

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

TH2D* GetTH2D(TFile *f, TString s) {
	if(!f->IsOpen()) {
		cout<<"Cannot open "<<f->GetName()<<endl;
	        return NULL;
	}
	if(verbose) cout<<"read "<<s<<endl;
	TH2D *h = (TH2D*)f->Get(s);
	if(!h) {
		cout<<"Cannot find "<<s<<" in "<<f->GetName()<<endl;
		return NULL;
	}
	return h;
}

bool ReadHist(TFile *f, TH1D* &hmeaspx, TH1D* &hrecopx, TH1D* &htrainmpx, TH1D* &htraintpx, TString tag="") {

	if(!f->IsOpen()) {
		cout<<"Cannot open "<<f->GetName()<<endl;
		return false;
	}

	cout<<"Read "<<f->GetName()<<endl;

	TH2D *hmeas = (TH2D*)GetTH2D(f,"hmeas");

	TH2D *hreco = (TH2D*)GetTH2D(f,"hreco");

	hmeaspx = (TH1D*)hmeas->ProfileX("hmeaspx"+tag);
	hrecopx = (TH1D*)hreco->ProfileX("hrecopx"+tag);

	TH2D *htrainm = (TH2D*)GetTH2D(f,"htrain");
	TH2D *htraint = (TH2D*)GetTH2D(f,"htraintrue");

	htrainmpx = (TH1D*)htrainm->ProfileX("htrainmpx"+tag);
	htraintpx = (TH1D*)htraint->ProfileX("htraintpx"+tag);

	SetHistStyle(hmeaspx, 2, 2, 25, 1, 2, 1);
	SetHistStyle(hrecopx, 1, 1, 20, 1, 2, 1);

	SetHistStyle(htrainmpx, 2, 2, 25, 1, 2, 3);
	SetHistStyle(htraintpx, 1, 1, 20, 1, 2, 3);

	TString ftitle = f->GetName();

	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	if(!ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) {
		hmeaspx->Scale(DeDpNorma);
		hrecopx->Scale(DeDpNorma);
		htrainmpx->Scale(DeDpNorma);
		htraintpx->Scale(DeDpNorma);
	}

	TString sregion = "Transverse";
	if(ftitle.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(ftitle.Contains("Away",TString::kIgnoreCase) || ftitle.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(ftitle.Contains("PtAve",TString::kIgnoreCase) || ftitle.Contains("AvePt",TString::kIgnoreCase)) svariable = "#LTp_{T}^{ch}#GT";
	if(ftitle.Contains("PtSum",TString::kIgnoreCase) || ftitle.Contains("SumPt",TString::kIgnoreCase)) svariable = "#sump_{T}^{ch}";
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) svariable = " Constituents Multiplicity";
	hmeaspx->GetYaxis()->SetTitle(sregion+svariable);
	hmeaspx->GetXaxis()->SetRangeUser(0,55);
	hmeaspx->SetMinimum(0);
	if(sregion.EqualTo("Transverse")) hmeaspx->SetMaximum(1.5);//Tran Ntrk density
	if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hmeaspx->SetMaximum(2.5);//Lead Ntrk density
	if(ftitle.Contains("PtSum",TString::kIgnoreCase) || ftitle.Contains("SumPt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hmeaspx->SetMaximum(0.9);//Tran PtSum
		if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hmeaspx->SetMaximum(10);//Lead PtSum
	}
	if(ftitle.Contains("PtAve",TString::kIgnoreCase) || ftitle.Contains("AvePt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hmeaspx->SetMaximum(0.68);//Tran PtAve
		if(sregion.EqualTo("Twoard")||sregion.EqualTo("Away"))  hmeaspx->SetMaximum(4);//Lead PtAve
	}
	if(ftitle.Contains("LeadJetNtrk",TString::kIgnoreCase)) hmeaspx->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)

	return true;

}


void plotUnfold(TString filename, TString filename2="") {		// filename is _NFweight, filename2 is not weight

	TFile *f = new TFile(filename);
	TH1D *hmeaspx, *hrecopx, *htrainmpx, *htraintpx;
	if(!ReadHist(f,hmeaspx, hrecopx, htrainmpx, htraintpx, "_NFweight")) {
		return;
	}


	TFile *f2 = new TFile(filename2);
	TH1D *hmeaspx2, *hrecopx2, *htrainmpx2, *htraintpx2;
	bool flagf2 = ReadHist(f2,hmeaspx2, hrecopx2, htrainmpx2, htraintpx2, "");
	if(flagf2) {
		SetHistStyle(hmeaspx2, 4, 4, 21, 1, 1.5, 1);
		SetHistStyle(hrecopx2, 4, 4, 30, 1, 2, 1);
		SetHistStyle(htrainmpx2, 2, 2, 25, 2, 2, 3);
		SetHistStyle(htraintpx2, 1, 1, 20, 2, 2, 3);
	}


	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	hmeaspx->Draw("p");
	hrecopx->Draw("psame");
	htrainmpx->Draw("HISTsame");
	htraintpx->Draw("HISTsame");
	if(flagf2) {
		hmeaspx2->Draw("eX0same");
		hrecopx2->Draw("eX0same");
		htrainmpx2->Draw("HISTsame");
		htraintpx2->Draw("HISTsame");
	}

	TLegend *leg;
	if(filename.Contains("LeadJetNtrk")) {
		leg = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else if(filename.Contains("Lead",TString::kIgnoreCase)||filename.Contains("Away",TString::kIgnoreCase)) {
		leg = new TLegend(0.5,0.2,0.88,0.48);		//Lead
	}
	else {
       		leg = new TLegend(0.5,0.6,0.88,0.88);	//Tran
	}
	leg->SetFillColor(0);
	leg->AddEntry(htrainmpx,"MC Detector-level","l");
	leg->AddEntry(htraintpx,"MC Particle-level","l");
	if(filename.Contains("_NFweight")) {
		leg->AddEntry(hmeaspx,"Measured data w/ NF weight","p");
		leg->AddEntry(hrecopx,"Unfolded data w/ NF weight","p");
	}
	else {
		leg->AddEntry(hmeaspx,"Measured data w/o NF weight","p");
		leg->AddEntry(hrecopx,"Unfolded data w/o NF weight","p");
	}
	if(flagf2) {
		if(filename2.Contains("FineBin")) {
			leg->AddEntry(hmeaspx2,"Measured data Fine Bin w/o NF weight","p");
			leg->AddEntry(hrecopx2,"Unfolded data Fine Bin w/o NF weight","p");
		}
		else {
			leg->AddEntry(htrainmpx2,"MC JP0 Detector-level","l");
			leg->AddEntry(htraintpx2,"MC JP0 Particle-level","l");
			leg->AddEntry(hmeaspx2,"Measured data w/o NF weight","p");
			leg->AddEntry(hrecopx2,"Unfolded data w/o NF weight","p");
		}
	}
	leg->Draw();

	TString sregion = "Tran";
	if(filename.Contains("Lead",TString::kIgnoreCase)) sregion = "Lead";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase)) sregion = "Sub";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion="LeadJetChargedAndNeutralNtrk";	// Lead Jet Ntrk (Charged+Neutral)
	c->SaveAs("Unfold"+sregion+"NtrkVsJetPt_JPs_NFweightedembedMBVsNoWeightembedJP0.png");


}
