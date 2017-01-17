
void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) {
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

void CompareBayesTimes(TString filename) {

	int savefig = 0;

	int DefaultTh = 4;
	const int Ntimes = 5;
	int OtherTh[Ntimes] = {2,3,5,6,7};

	TFile *f[Ntimes+1];
	f[0] = new TFile(filename);
	TH2D *hreco[Ntimes+1];
	hreco[0] = (TH2D*)f[0]->Get("hreco");
	hreco[0]->SetName(Form("hreco%d",DefaultTh));
	for(int i = 1; i<Ntimes+1; i++) {
		Ssiz_t toinsert = filename.Index(".root");
		TString ifilename = filename;
		ifilename.Insert(toinsert,Form("_Baye%d",OtherTh[i-1]));
		cout<<ifilename<<endl;
		f[i] = new TFile(ifilename);
		hreco[i] = (TH2D*)f[i]->Get("hreco");
		hreco[i]->SetName(Form("hreco%d",OtherTh[i-1]));
	}
	TProfile *hrecopfy[Ntimes+1];
	hrecopfy[0] = (TProfile*)hreco[0]->ProfileX(Form("hreco_pfx%d",DefaultTh));
	for(int i = 1; i<Ntimes+1; i++) {
		hrecopfy[i] = (TProfile*)hreco[i]->ProfileX(Form("hreco_pfx%d",OtherTh[i-1]));
	}

	//int filecolor[7] = {kBlack,kRed,kBlue,kGreen+1,kBlue+3,kMagenta,kPink-7};
	int filecolor[6] = {kBlack,kOrange+1, kRed,kBlue,kAzure+10, kTeal+3};

	SetHistStyle(hrecopfy[0], filecolor[0], filecolor[0], 8, 1, 2, 1);
	for(int i = 1; i<Ntimes+1; i++) {
		SetHistStyle(hrecopfy[i], filecolor[i], filecolor[i], 25, 1, 2, 1);
	}


	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	if(!filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) {
		for(int i = 0; i<Ntimes+1; i++) {
			hrecopfy[i]->Scale(DeDpNorma);
		}
	}


	TString sregion = "Transverse";
	if(filename.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) svariable = "#LTp_{T}^{ch}#GT";
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) svariable = "#sump_{T}^{ch}";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) svariable = " Constituents Multiplicity";
	hrecopfy[0]->GetYaxis()->SetTitle(sregion+svariable);
	hrecopfy[0]->GetXaxis()->SetRangeUser(0,55);
	hrecopfy[0]->SetMinimum(0);
	if(sregion.EqualTo("Transverse")) {
		hrecopfy[0]->SetMaximum(0.9);//Tran Ntrk density
		hrecopfy[0]->SetMinimum(0.3);
	}
	if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hrecopfy[0]->SetMaximum(2.5);//Lead Ntrk density
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hrecopfy[0]->SetMaximum(0.9);//Tran PtSum
		if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hrecopfy[0]->SetMaximum(10);//Lead PtSum
	}
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hrecopfy[0]->SetMaximum(0.68);//Tran PtAve
		if(sregion.EqualTo("Twoard")||sregion.EqualTo("Away"))  hrecopfy[0]->SetMaximum(4);//Lead PtAve
	}
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) hrecopfy[0]->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)





	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	hrecopfy[0]->Draw("p");
	for(int i = 1; i<Ntimes+1; i++) {
		hrecopfy[i]->Draw("psame");
	}


	TLegend *leg;
	if(filename.Contains("LeadJetNtrk")) {
		leg = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else if(filename.Contains("Lead",TString::kIgnoreCase)||filename.Contains("Away",TString::kIgnoreCase)||filename.Contains("Sub",TString::kIgnoreCase)) {
		leg = new TLegend(0.5,0.2,0.88,0.48);		//Lead
	}
	else {
       		leg = new TLegend(0.5,0.6,0.88,0.88);	//Tran
	}
	leg->SetFillColor(0);
	leg->AddEntry(hrecopfy[0],Form("iteration=%d",DefaultTh));
	for(int i = 1; i<Ntimes+1; i++) {
		leg->AddEntry(hrecopfy[i],Form("iteration=%d",OtherTh[i-1]));
	}
	leg->Draw();

	TString outtag = "TranTot";
	if(filename.Contains("Lead",TString::kIgnoreCase)) outtag="LeadArea";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase))  outtag="SubArea";
	//c->SaveAs("BayesIter_"+outtag+"NtrkJPCharged_NFweight_McPt02_embedMB.png");
	//c->SaveAs("BayesIter_"+outtag+"NtrkJPCharged_McPt02_embedJP0.png");

}
