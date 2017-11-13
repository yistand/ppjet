//====================================================================================================
//		2017.10.24	Li YI
//		There is jumps in TranPtAve, let's try to smooth it to see what happens	
//====================================================================================================			

TH1D* SmoothProfileX(TH2D *h2, TString name, double xmin, double xmax) {

	//cout<<"xmin = "<<xmin<<" xmax = "<<xmax<<endl;
	const double *binning = h2->GetXaxis()->GetXbins()->GetArray();
	TH1D *hp = new TH1D(name,name,h2->GetXaxis()->GetNbins(),binning);

	for(int i = h2->GetXaxis()->FindBin(xmin); i<=h2->GetXaxis()->FindBin(xmax); i++) {
		TH1D* htmp = (TH1D*)h2->ProjectionY(Form("htmp_projy%d",i),i,i);
		TH1D* htmp_sc = (TH1D*)htmp->Clone(Form("htmpsc_projy%d",i));	// scale width
		TH1D* htmp_scsm = (TH1D*)htmp->Clone(Form("htmpscsm_projy%d",i)); // scale width + smooth

		htmp_sc->Scale(1,"width");
		htmp_scsm->Scale(1,"width");
		htmp_scsm->Smooth();	// default smooth once
		TH1D* htmp_r = (TH1D*)htmp_scsm->Clone(Form("htmpr_projy%d",i));	// ratio for correction 
		htmp_r->Divide(htmp_sc);

		htmp->Multiply(htmp_r);
		hp->SetBinContent(i,htmp->GetMean());	
		hp->SetBinError(i,htmp->GetMeanError());	
		//cout<<name<<i<<" = "<<htmp->GetMean()<<"+-"<<htmp->GetMeanError()<<endl;
	}
	return hp;
}


TH1D* SmoothProfileX(TH2D *h2, TString name, double xmin) {
	double xmax = h2->GetXaxis()->GetXmax();
	return SmoothProfileX(h2,name,xmin,xmax);
}


TH1D* SmoothProfileX(TH2D *h2, TString name) {
	double xmin = h2->GetXaxis()->GetXmin();
	double xmax = h2->GetXaxis()->GetXmax();
	return SmoothProfileX(h2,name,xmin,xmax);
}

