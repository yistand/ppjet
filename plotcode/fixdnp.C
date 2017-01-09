// 	Li Yi		2017.01.08
// 	Found out plots shown on DNP and MPI,  Charged Density and Neutral Energy Density for toward, away and transverse ones are scaled by 1/(2*3/pi) instead of 1/(4*pi/3)
// 	Delta eta = 2, Delta phi = 120 degrees = 120/180*pi = 2/3*pi, Delta eta * Delta phi = 4/3*pi
// 	The root files are still the raw values. Here we redo the scaling with the right number

void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, float msize, int ydivi = 505, int lstyle=1, int lwidth=1) {
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
        h->GetYaxis()->SetNdivisions(ydivi);

	h->GetXaxis()->SetRangeUser(9,60);

}





void fixdnp() {

	const int nfig = 2;
	TString filename[nfig] = {"ChargedNtrk.root","NeutralSumPt.root"};
	TString cavname[nfig] = {"c1","c1_n3"};
	TString yname[nfig] = {"Detector-level #LTN_{ch}/#delta#eta#delta#phi#GT", "Detector-level #LTE^{ne}_{T}#GT (GeV/#it{c})"};
	TString xname = "Detector-level Leading Jet p_{T} (GeV/#it{c})";
	TFile *f[nfig];
	TCanvas *c[nfig];
	float ymax[nfig] = {3,6.7};
	int ydivi[nfig] = {505,510};
	
	
	const int nhist = 3;
	char *histname[nfig][nhist] = {{"leadjetareantrkvsleadjetpt","subjetareantrkvsleadjetpt","tranntrkvsleadjetpt"},{"leadjetareaptsumvsleadjetpt","subjetareaptsumvsleadjetpt","tranptsumvsleadjetpt"}};
	int mcolor[nhist] = {1,4,2};
	int mstyle[nhist] = {8,25,25};
	TProfile *p[nfig][nhist];
	TCanvas *c1;
	TProfile *ptmp;
	for(int i=0 ; i<nfig; i++) {
		cout<<"Open "<<filename[i]<<endl;
		f[i] = new TFile(filename[i]);	
		for(int j=0; j<nhist; j++) {
			cout<<"Read "<<histname[i][j]<<endl;
			c1 = (TCanvas*)f[i]->Get(cavname[i]);
			p[i][j] = (TProfile*)c1->FindObject(histname[i][j]);
			if(!p[i][j]) { cout<<"read NULL for "<<histname[i][j]<<endl; return;}
			p[i][j]->Scale(1./(4.*TMath::Pi()/3.));
			SetHistStyle(p[i][j], mcolor[j], mcolor[j], mstyle[j], 1.5, ydivi[i]);
		}
		c[i] = new TCanvas(Form("canvas%d",i));
		c[i]->SetFrameLineWidth(2);
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		p[i][0]->GetXaxis()->SetTitle(xname);
		p[i][0]->GetYaxis()->SetTitle(yname[i]);
		p[i][0]->GetYaxis()->SetTitleOffset(0.8);
		ptmp = (TProfile*)p[i][0]->Clone(Form("ptmp%d",i));
		ptmp->Reset();
		ptmp->GetXaxis()->SetRangeUser(5,60);
		ptmp->GetYaxis()->SetRangeUser(0,ymax[i]);
		ptmp->Draw();
		for(int j = 0; j<nhist; j++) {
			p[i][j]->Draw("same");
		}
	}



}


