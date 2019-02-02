
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

void CompareBayesTimes4Test(TString filename="TrainTest_TranTotNtrkMBChargedmPrior_BT170928_McPt02.root") {

	int savefig = 1;

	int DefaultTh = 5;
	const int Ntimes = 6;
	int OtherTh[Ntimes] = {2,3,4,6,7,8};

	TFile *f[Ntimes+1];

	Ssiz_t toremove = filename.Index("_Baye");
	TString Deffilename = filename;
	if(toremove!=-1) {
		Deffilename.Remove(toremove,6);	// assuming Bayes time no more than 9, other wise we will need 7 instead of 6
	}

	TString simfilename = Deffilename;
	Ssiz_t toremove2 = simfilename.Index(".root");
	if(toremove2!=-1) {
		simfilename.Remove(toremove2,5); // remove ".root"
	}

	//Ssiz_t toinsert = Deffilename.Index(".root");
	Ssiz_t toinsert; 
	if(filename.Contains("mPrior")) {toinsert = Deffilename.Index("mPrior_BT170928");}
	else {toinsert = Deffilename.Index("_BT170928");};

	if(DefaultTh!=4) {
		Deffilename.Insert(toinsert,Form("_Baye%d",DefaultTh));
	}
	f[0] = new TFile(Deffilename);
	cout<<Deffilename<<endl;
	
	TH2D *htraintrue;
	htraintrue = (TH2D*)f[0]->Get("htraintrue");
	TProfile *htraintruepfy;
	htraintruepfy = (TProfile*)htraintrue->ProfileX(Form("htraintrue_pfx"));

	TH2D *htrain;
	htrain = (TH2D*)f[0]->Get("htrain");
	TProfile *htrainpfy;
	htrainpfy = (TProfile*)htrain->ProfileX(Form("htrain_pfx"));

	TH2D *htrue;
	htrue = (TH2D*)f[0]->Get("htrue");
	TProfile *htruepfy;
	htruepfy = (TProfile*)htrue->ProfileX(Form("htrue_pfx"));

	TH2D *hreco[Ntimes+1];
	hreco[0] = (TH2D*)f[0]->Get("hreco");
	hreco[0]->SetName(Form("hreco%d",DefaultTh));
	for(int i = 1; i<Ntimes+1; i++) {
		TString ifilename = filename;
		if(toremove!=-1) {
			ifilename.Remove(toremove,6);	// assuming Bayes time no more than 9, other wise we will need 7 instead of 6
		}
		if(OtherTh[i-1]!=4) {
			ifilename.Insert(toinsert,Form("_Baye%d",OtherTh[i-1]));
		}
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

	// As xsec is used to merge different pT bins. we scale 

	// Take ratio between reco and true
	TH1D *hratio_recoRtrue[Ntimes+1+1];
        const TArrayD* xa = hrecopfy[0]->GetXaxis()->GetXbins();
        const Double_t* xbin = xa->GetArray();
	hratio_recoRtrue[0] = new TH1D(Form("hratio_recoRtrue%d",DefaultTh),Form("reco/true Bayes %d",DefaultTh),hrecopfy[0]->GetNbinsX(),xbin);
	for(int i = 0; i<Ntimes+1; i++) {
		hratio_recoRtrue[i] = new TH1D(Form("hratio_recoRtrue%d",OtherTh[i-1]),Form("reco/true Bayes %d",OtherTh[i-1]),hrecopfy[0]->GetNbinsX(),xbin);
	}
	hratio_recoRtrue[Ntimes+1] = new TH1D(Form("hratio_recoRtrue_train"),Form("traintrue/true Bayes"),hrecopfy[0]->GetNbinsX(),xbin);
	for(int j = 0; j<htrue->GetNbinsX(); j++) {
		double ytrue = htruepfy->GetBinContent(j+1);
		double eytrue = htruepfy->GetBinError(j+1);
		if(fabs(ytrue)<1e-13) continue;
		for(int i = 0; i<Ntimes+1; i++) {
			double yreco = hrecopfy[i]->GetBinContent(j+1);
			double eyreco = hrecopfy[i]->GetBinError(j+1);
			//double yr = (yreco-ytrue)/ytrue;
			//double eyr = (eyreco*ytrue-yreco*eytrue)/(ytrue*ytrue); 
			double yr = (yreco)/ytrue;
			double eyr = (eyreco*ytrue-yreco*eytrue)/(ytrue*ytrue); 
			hratio_recoRtrue[i]->SetBinContent(j+1,yr);
			hratio_recoRtrue[i]->SetBinError(j+1,eyr);
			//cout<<"Bayes "<<i<<" pt "<<j<<" ratio = "<<yreco<<"-"<<ytrue<<"/"<<ytrue<<" = "<<yreco-ytrue<<"/"<<ytrue<<" = "<<yr<<" +- "<<eyr<<endl;
		}
		double yreco = htraintruepfy->GetBinContent(j+1);
		double eyreco = htraintruepfy->GetBinError(j+1);
		//double yr = (yreco-ytrue)/ytrue;
		//double eyr = (eyreco*ytrue-yreco*eytrue)/(ytrue*ytrue); 
		double yr = (yreco)/ytrue;
		double eyr = (eyreco*ytrue-yreco*eytrue)/(ytrue*ytrue); 
		hratio_recoRtrue[Ntimes+1]->SetBinContent(j+1,yr);
		hratio_recoRtrue[Ntimes+1]->SetBinError(j+1,eyr);
	}

	int filecolor[8] = {kBlack,kRed,kBlue,kGreen+1,kBlue+3,kMagenta,kPink-7,kOrange-2};
	//int filecolor[6] = {kBlack,kOrange+1, kRed,kBlue,kAzure+10, kTeal+3};

	SetHistStyle(htruepfy, filecolor[0], filecolor[0], 23, 1, 2, 1);
	SetHistStyle(hrecopfy[0], filecolor[0], filecolor[0], 8, 1, 2, 1);
	SetHistStyle(hratio_recoRtrue[0], filecolor[0], filecolor[0], 8, 1, 2, 1);
	for(int i = 1; i<Ntimes+1; i++) {
		SetHistStyle(hrecopfy[i], filecolor[i], filecolor[i], 25, 1, 2, 1);
	}
	for(int i = 1; i<Ntimes+2; i++) {
		SetHistStyle(hratio_recoRtrue[i], filecolor[i], filecolor[i], 25, 1, 2, 1);
	}


	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	cout<<filename<<endl;
	if(!filename.Contains("PtAve",TString::kIgnoreCase)&&!filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) {
		htruepfy->Scale(DeDpNorma);
		for(int i = 0; i<Ntimes+1; i++) {
			hrecopfy[i]->Scale(DeDpNorma);
		}
	}


	// Set Y Title
	TString sregion = "Transverse";
	if(filename.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) svariable = "#LTp_{T}^{ch}#GT";
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) svariable = "#sump_{T}^{ch}";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) svariable = " Constituents Multiplicity";
	hrecopfy[0]->GetYaxis()->SetTitle(sregion+svariable);


	//hratio_recoRtrue[0]->GetYaxis()->SetTitle(sregion+svariable+Form(" (reco-true)/true"));
	hratio_recoRtrue[0]->GetYaxis()->SetTitle(sregion+svariable+Form(" reco/true"));


	// Set Draw range and so 
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
		hrecopfy[0]->SetMinimum(0);
		if(sregion.Contains("Tran")) hrecopfy[0]->SetMaximum(1);//Tran PtAve
		if(sregion.EqualTo("Toward")) hrecopfy[0]->SetMaximum(5);//Lead PtAve
		if(sregion.EqualTo("Away"))  hrecopfy[0]->SetMaximum(5);//Away PtAve
	}
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) hrecopfy[0]->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)



	// Draw
	// all Bayes times together
	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetErrorX(0);
	
	hrecopfy[0]->Draw("X0ep");
	htruepfy->Draw("X0epsame");
	for(int i = 1; i<Ntimes+1; i++) {
		hrecopfy[i]->Draw("X0epsame");
	}
	htraintruepfy->Draw("X0epsame");


	TLegend *leg;
	if(filename.Contains("LeadJetNtrk")) {
		leg = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else if(filename.Contains("Lead",TString::kIgnoreCase)||filename.Contains("Away",TString::kIgnoreCase)||filename.Contains("Sub",TString::kIgnoreCase)) {
		leg = new TLegend(0.7,0.2,0.88,0.48);		//Lead
	}
	else {
       		leg = new TLegend(0.7,0.6,0.88,0.88);	//Tran
	}
	leg->SetFillColor(0);
	leg->AddEntry(htruepfy,Form("Truth"));
	leg->AddEntry(hrecopfy[0],Form("iteration=%d",DefaultTh));
	for(int i = 1; i<Ntimes+1; i++) {
		leg->AddEntry(hrecopfy[i],Form("iteration=%d",OtherTh[i-1]));
	}
	leg->AddEntry(htraintruepfy,Form("TrainTruth"));
	leg->Draw();

	TString outtag = "TranTot";
	if(filename.Contains("Lead",TString::kIgnoreCase)) outtag="LeadArea";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase))  outtag="SubArea";

	if(savefig && false) {
		//c->SaveAs("BayesIter_"+outtag+"NtrkJPCharged_NFweight_McPt02_embedMB.png");
		//c->SaveAs("BayesIter_"+outtag+"NtrkJPCharged_McPt02_embedJP0.png");
		c->SaveAs("GPC-2/BayesIter_"+simfilename+".pdf");
	}


	// Draw ratio between reco and true

	hratio_recoRtrue[0]->GetYaxis()->SetTitleOffset(1.5);
	hratio_recoRtrue[0]->GetXaxis()->SetRangeUser(0,45);
	//hratio_recoRtrue[0]->SetMaximum(0.1);
	//hratio_recoRtrue[0]->SetMinimum(-0.1);
	hratio_recoRtrue[0]->SetMaximum(1.2);
	hratio_recoRtrue[0]->SetMinimum(0.8);

	TCanvas *c2 = new TCanvas("c2","c2",1000,800);
	c2->SetLeftMargin(0.12);
	c2->SetBottomMargin(0.12);
	c2->SetLeftMargin(0.15);
	c2->SetRightMargin(0.05);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	//hratio_recoRtrue[0]->Draw("p");
	TH1D *htemp2 = (TH1D*)hratio_recoRtrue[0]->Clone("htemp2");
	htemp2->Reset();
	htemp2->Draw();
	double rptmin = 4;
	double rptmax = 45;
	hratio_recoRtrue[0]->GetXaxis()->SetRangeUser(rptmin,rptmax);
	hratio_recoRtrue[0]->Draw("HIST][same");
	for(int i = 1; i<Ntimes+2; i++) {
		hratio_recoRtrue[i]->GetXaxis()->SetRangeUser(rptmin,rptmax);
		hratio_recoRtrue[i]->Draw("HIST][same");
	}

	TLegend *leg2;
	if(filename.Contains("LeadJetNtrk")) {
		leg2 = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else if(filename.Contains("Lead",TString::kIgnoreCase)||filename.Contains("Away",TString::kIgnoreCase)||filename.Contains("Sub",TString::kIgnoreCase)) {
		leg2 = new TLegend(0.7,0.2,0.88,0.48);		//Lead
	}
	else {
       		leg2 = new TLegend(0.7,0.6,0.88,0.88);	//Tran
	}
	leg2->SetFillColor(0);
	leg2->AddEntry(hratio_recoRtrue[0],Form("iteration=%d",DefaultTh),"l");
	for(int i = 1; i<Ntimes+1; i++) {
		leg2->AddEntry(hratio_recoRtrue[i],Form("iteration=%d",OtherTh[i-1]),"l");
	}
	leg2->AddEntry(hratio_recoRtrue[Ntimes+1],Form("TrainTruth"),"l");
	leg2->Draw();

	TLine *l = new TLine(0,0,45,0);
	l->Draw();

	if(savefig) {
		//cout<<TString("GPC-2/BayesIter_"+simfilename+".png")<<endl;
		c2->SaveAs("GPC-2/BayesIter_ratio_"+simfilename+".pdf");
	}

}
