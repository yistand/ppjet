//====================================================================================================
//
//	2017.04.07	Li YI
//	Test 3 cases for TranMax TranMin Ntrk
//	1. from same two-Gaussian distributions
//	2. Max from 1 Gaussian, Min from another Gaussian
//	3. Part of event from 1 Gauss, the rest from another Gauss
//
//====================================================================================================


void SampleMaxMin1(TH1D *hmax, TH1D *hmin, TF1 *f1, int Nevent = 50000) {
	for(int i = 0; i<Nevent; i++) {
		double a = f1->GetRandom();
		double b = f1->GetRandom();
		if(a>=b) {hmax->Fill(a); hmin->Fill(b);}
		else {hmax->Fill(b); hmin->Fill(a);}
	}
}

void SampleMaxMin2(TH1D *hmax, TH1D *hmin, TF1 *f1, TF1 *f2, int Nevent = 50000) {
	for(int i = 0; i<Nevent; i++) {
		double a = f1->GetRandom();
		double b = f2->GetRandom();
		if(a>=b) {hmax->Fill(a); hmin->Fill(b);}
		else {hmax->Fill(b); hmin->Fill(a);}
	}
}

void SampleMaxMin3(TH1D *hmax, TH1D *hmin, TF1 *f1, TF1 *f2, int Nevent = 50000) {
	for(int i = 0; i<Nevent/2; i++) {
		double a = f1->GetRandom();
		double b = f1->GetRandom();
		if(a>=b) {hmax->Fill(a); hmin->Fill(b);}
		else {hmax->Fill(b); hmin->Fill(a);}
	}
	for(int i = 0; i<(Nevent-Nevent/2); i++) {
		double a = f2->GetRandom();
		double b = f2->GetRandom();
		if(a>=b) {hmax->Fill(a); hmin->Fill(b);}
		else {hmax->Fill(b); hmin->Fill(a);}
	}
}



void TestMaxMin(int Nevent = 50000) {

	TF1 *g1 = new TF1("g1","[2]*TMath::Gaus(x,[0],[1])",0,20);
	TF1 *g2 = new TF1("g2","[2]*TMath::Gaus(x,[0],[1])",0,20);
	g1->SetParameters(4,3,1);
	g2->SetParameters(10,3,1.5);
	TF1 *gsum = new TF1("gsum","g1+g2",0,20);
	for(int i = 0; i<3; i++) {
		cout<<"Set par["<<i<<"] = "<<g1->GetParameter(i)<<endl;
		gsum->SetParameter(i,g1->GetParameter(i));
	}
	for(int i = 0; i<3; i++) {
		cout<<"Set par["<<i+3<<"] = "<<g2->GetParameter(i)<<endl;
		gsum->SetParameter(i+3,g2->GetParameter(i));
	}

	TCanvas *cg = new TCanvas("c1","cg");
	gsum->GetXaxis()->SetTitle("Ntrk");
	gsum->GetYaxis()->SetTitle("dist.");
	gsum->Draw();
	g1->SetLineColor(kYellow);
	g2->SetLineColor(kGreen);
	g1->Draw("same");
	g2->Draw("same");
	TLatex *lcg = new TLatex();
	lcg->SetNDC();
	lcg->SetTextColor(2);
	lcg->DrawLatex(0.6,0.8,"Gaus1+Gaus2");
	lcg->SetTextColor(g1->GetLineColor());
	lcg->DrawLatex(0.2,0.3,"Gaus1");
	lcg->SetTextColor(g2->GetLineColor());
	lcg->DrawLatex(0.45,0.5,"Gaus2");


	const int Ncase = 3;
	TH1D *hmax[Ncase], *hmin[Ncase];

	for(int i = 0; i<Ncase; i++) {
		hmax[i] = new TH1D(Form("hmax%d",i+1),Form("Max Case%d",i+1),20,0,20);
		hmin[i] = new TH1D(Form("hmin%d",i+1),Form("Min Case%d",i+1),20,0,20);
	}

	SampleMaxMin1(hmax[0],hmin[0],gsum,Nevent);
	SampleMaxMin2(hmax[1],hmin[1],g1,g2,Nevent);
	SampleMaxMin3(hmax[2],hmin[2],g1,g2,Nevent);

	int color[Ncase] = {1,2,4};
	for(int i = 0; i<Ncase; i++) {
		hmax[i]->GetXaxis()->SetTitle("Ntrk");
		hmax[i]->GetYaxis()->SetTitle("Nevents");
		hmin[i]->GetXaxis()->SetTitle("Ntrk");
		hmin[i]->GetYaxis()->SetTitle("Nevents");
		hmax[i]->SetLineColor(color[i]);
		hmin[i]->SetLineColor(color[i]);
		hmax[i]->SetMarkerColor(color[i]);
		hmin[i]->SetMarkerColor(color[i]);
		hmax[i]->SetMarkerStyle(4);
		hmin[i]->SetMarkerStyle(8);
	}
	
	TCanvas *cmaxmin[Ncase];
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
       	for(int i = 0; i<Ncase; i++) {
		cmaxmin[i] = new TCanvas(Form("case%d",i+1),Form("cmaxmin case %d",i+1));
	}
	for(int i = 0; i<Ncase; i++) {
		cmaxmin[i]->cd();
		hmin[i]->SetMaximum(7000);
		hmin[i]->Draw("p");
		hmax[i]->Draw("psame");
		cout<<"Case "<<i+1<<" Max-Min = "<<hmax[i]->GetMean()<<" - "<<hmin[i]->GetMean()<<" = "<<hmax[i]->GetMean() - hmin[i]->GetMean()<<endl;
	}
	TLatex *lm = new TLatex();
	lm->SetNDC();
	for(int i = 0; i<Ncase; i++) {
		cmaxmin[i]->cd();
		lm->SetTextColor(hmax[i]->GetMarkerColor());
		lm->SetTextFont(62);
		lm->DrawLatex(0.3,0.92,Form("Case %d: <Max>-<Min>=%.1f",i+1,hmax[i]->GetMean() - hmin[i]->GetMean()));
		lm->SetTextFont(62);
		lm->DrawLatex(g1->GetParameter(0)/20.,0.5,"Min");
		lm->SetTextFont(12);
		lm->DrawLatex(g2->GetParameter(0)/20.,0.5,"Max");
	}

}
