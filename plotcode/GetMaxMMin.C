//====================================================================================================
//
//	2017.04.06	Li Yi
//	Assume TransMax TransMin from the same distribution
//	How much is expected mean for TransMax - TransMin
//	P(max) = 2f(x)F(x)
//	P(min) = 2f(x)(1-F(x))
//	E(max-min) = Integral(x*2*f(x)*(2*Integral(f(t),t, -int, x)-1),x, -int, int)
//
//====================================================================================================


float Emean(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function EMean() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMean() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	float Mean = 0;
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		Mean+=x*2*y*(2*Fdist-1);
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return Mean;
}

float EMax(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function EMax() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMax() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	float Mean = 0;
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		Mean+=x*2*y*Fdist;
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return Mean;
}


float EMin(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function EMin() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMin() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	float Mean = 0;
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		Mean+=x*2*y*(1-Fdist);
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return Mean;
}

TH1I *inthistMin(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function EMin() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMin() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	const TArrayD* xa = fdist->GetXaxis()->GetXbins();
	const Double_t* xbin = xa->GetArray();
	TH1I *hmin = new TH1I(Form("%s_min",fdist->GetName()),Form("%s_min",fdist->GetName()),fdist->GetNbinsX(),xbin);
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		//float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		//Mean+=x*2*y*(1-Fdist);
		hmin->SetBinContent(i,(int)(2*y*(1-Fdist)*sum));
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return hmin;
}


TH1I *inthistMax(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function hMax() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMax() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	const TArrayD* xa = fdist->GetXaxis()->GetXbins();
	const Double_t* xbin = xa->GetArray();
	TH1I *hmax = new TH1I(Form("%s_max",fdist->GetName()),Form("%s_max",fdist->GetName()),fdist->GetNbinsX(),xbin);
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		//float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		//Mean+=x*2*y*Fdist;
		hmax->SetBinContent(i,(int)(2*y*Fdist*sum));
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return hmax;
}


TH1F *histMin(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function EMin() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMin() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	TH1F *hmin = (TH1F*)fdist->Clone(Form("%s_min",fdist->GetName()));
	hmin->Reset();
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		//float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		//Mean+=x*2*y*(1-Fdist);
		hmin->SetBinContent(i,(2*y*(1-Fdist)*sum));
		//hmin->SetBinContent(i,(2*y*(1-Fdist)));
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return hmin;
}


TH1F *histMax(const TH1F *fdist) {

	if(fdist->Integral()<=0) {cout<<"Warning!! function hMax() input fdist "<<fdist->GetName()<<" is empty "<<endl; return 0;}
	float sum = fdist->Integral();
	if(sum<=0) {cout<<"Warning!! function EMax() input fdist is not positive "<<endl; return 0;}
	TH1F *fcopy = (TH1F*)fdist->Clone(Form("%s_copy",fdist->GetName()));
	//cout<<sum<<endl;
	fcopy->Scale(1./sum);

	TH1F *hmax = (TH1F*)fdist->Clone(Form("%s_max",fdist->GetName()));
	hmax->Reset();
	for(int i = 1; i<fcopy->GetNbinsX()+1; i++) {
		float Fdist = fcopy->Integral(1,i);
		//float x = fcopy->GetBinCenter(i);
		float y = fcopy->GetBinContent(i);
		//Mean+=x*2*y*Fdist;
		hmax->SetBinContent(i,(2*y*Fdist*sum));
		//hmax->SetBinContent(i,(2*y*Fdist));
		//cout<<"x = "<<x<<" y = "<<y<<" Fdist = "<<Fdist<<" Mean+= "<<x*2*y*(2*Fdist-1)<<endl;
	}
	return hmax;
}

TH1I *TH1F2TH1I(TH1F* h) {
	const TArrayD* xa = h->GetXaxis()->GetXbins();
	const Double_t* xbin = xa->GetArray();
	TH1I *hi = new TH1I(Form("int%s",h->GetName()),Form("int%s",h->GetName()),h->GetNbinsX(),xbin);
	for(int i = 1; i<h->GetNbinsX()+1; i++) {
		hi->SetBinContent(i,(int)h->GetBinContent(i+1));
		hi->SetBinError(i,(int)h->GetBinError(i+1));
	}
	return hi;
}

float MaxpValue(TH1F *hmeas, TH1F *hmax_meas) {		// inputs are Ntrk distributions for one jet pt bin
	TH1F *hmax_calc = histMax(hmeas);
	//hmax_meas->Scale(1000);
	cout<<endl<<endl;
	hmax_meas->Print("range");
	hmax_calc->Print("range");
	//return hmax_meas->Chi2Test(hmax_calc,"UU P CHI2/NDF");
	return hmax_meas->Chi2Test(hmax_calc,"UU P");
}

float MinpValue(TH1F *hmeas, TH1F *hmin_meas) {		// inputs are Ntrk distributions for one jet pt bin
	//TH1F *hmin_calc = histMin(hmeas);
	TH1F *hmin_calc = histMax(hmeas);
	//hmin_meas->Scale(1000);
	cout<<endl<<endl;
	hmin_meas->Print("range");
	hmin_calc->Print("range");
	//return hmin_meas->Chi2Test(hmin_calc,"UU P CHI2/NDF");
	return hmin_meas->Chi2Test(hmin_calc,"UU P");
}

void PrintHistWErr(TH1 *h) {
	cout<<h->GetName()<<" "<<h->GetTitle()<<":"<<endl;
	for(int i = 1; i<h->GetNbinsX()+1; i++) {
		cout<<h->GetBinLowEdge(i)<<"-"<<h->GetBinLowEdge(i)+h->GetBinWidth(i)<<", "
			<<h->GetBinContent(i)<<" +- "
			<<h->GetBinError(i)
			<<endl;

	}
}


void PrintHistWoErr(TH1 *h) {
	cout<<h->GetName()<<" "<<h->GetTitle()<<":"<<endl;
	for(int i = 1; i<h->GetNbinsX()+1; i++) {
		cout<<h->GetBinLowEdge(i)<<"-"<<h->GetBinLowEdge(i)+h->GetBinWidth(i)<<", "<<h->GetBinContent(i)<<endl;
	}
}

//================================================  Call This one to draw p-Value for each jet pT for TransMax and Trans Min ====================================================
// > .L GetMaxMMin.C
// > Draw4pValue()
void Draw4pValue(TString filename = "Unfolding_TranTotNtrkJPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root") {
	TFile *file = new  TFile(filename);
	TH2F *hmeas = (TH2F*)file->Get("hmeas");
	
        TString filenameMax = filename;
        filenameMax.ReplaceAll("Tot","Max");
        TFile *filemax = new TFile(filenameMax);
	TH2F *hmaxmeas = (TH2F*)filemax->Get("hmeas");

        TString filenameMin = filename;
        filenameMin.ReplaceAll("Tot","Min");
        TFile *filemin = new TFile(filenameMin);
	TH2F *hminmeas = (TH2F*)filemin->Get("hmeas");

	TH1F *hpvaluemax = (TH1F*)hmaxmeas->ProjectionX("hpvaluemax");
	hpvaluemax->Reset();
	hpvaluemax->GetYaxis()->SetTitle("p-value TransMax");
	TH1F *hpvaluemin = (TH1F*)hminmeas->ProjectionX("hpvaluemin");
	hpvaluemin->Reset();
	hpvaluemin->GetYaxis()->SetTitle("p-value TransMin");

	for(int j = 1; j<=hmeas->GetNbinsX(); j++) {
		TH1F *hmeas_py = (TH1F*)hmeas->ProjectionY(Form("hmeas_py%d",j),j,j);
		TH1F *hmaxmeas_py = (TH1F*)hmaxmeas->ProjectionY(Form("hmaxmeas_py%d",j),j,j);
		TH1F *hminmeas_py = (TH1F*)hminmeas->ProjectionY(Form("hminmeas_py%d",j),j,j);
		hpvaluemax->SetBinContent(j,MaxpValue(hmeas_py,hmaxmeas_py));
		hpvaluemin->SetBinContent(j,MinpValue(hmeas_py,hminmeas_py));
	}

	hpvaluemax->GetXaxis()->SetRangeUser(3,40);
	hpvaluemin->GetXaxis()->SetRangeUser(3,40);
	TLine *l = new TLine(3,1,45,1);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *cmax = new TCanvas(); 
	//hpvaluemax->SetFillColor(4);
	//hpvaluemax->SetBarWidth(0.4);
	//hpvaluemax->Draw("bar");
	hpvaluemax->Draw("HIST");
	l->Draw("same");
	//hpvaluemax->Print("range");
	PrintHistWoErr(hpvaluemax);

	TCanvas *cmin = new TCanvas(); 
	//hpvaluemin->SetFillColor(4);
	//hpvaluemin->SetBarWidth(0.4);
	//hpvaluemin->Draw("bar");
	hpvaluemin->Draw("HIST");
	l->Draw("same");
	//hpvaluemin->Print("range");
	PrintHistWoErr(hpvaluemin);

}


TH1F *Diff(TString filename) {
	TString filenameMax = filename;
	filenameMax.ReplaceAll("Tot","Max");
	TString filenameMin = filename;
	filenameMin.ReplaceAll("Tot","Min");

	TFile *filemax = new TFile(filenameMax);
	TFile *filemin = new TFile(filenameMin);

	TProfile *pymax = (TProfile*)filemax->Get("pfxmeas");
	pymax->SetName("pymax_copy");
	TProfile *pymin = (TProfile*)filemin->Get("pfxmeas");
	pymin->SetName("pymin_copy");

	TH2F *hmeas = (TH2F*)filemax->Get("hmeas");
	hmeas->SetName("hmeas_copy");
	TH1F *diff_maxmin = (TH1F*)hmeas->ProjectionX("diff_maxmin");
	diff_maxmin->Reset();
	diff_maxmin->SetTitle("<TransMax>-<TransMin>");

	for(int i = 1; i<=pymax->GetNbinsX() ; i++) {
		diff_maxmin->SetBinContent(i,pymax->GetBinContent(i)-pymin->GetBinContent(i));
		diff_maxmin->SetBinError(i,sqrt(pow(pymax->GetBinError(i),2)+pow(pymin->GetBinError(i),2)));
	}
	diff_maxmin->SetDirectory(0);
	filemax->Close();
	filemin->Close();
	
	return diff_maxmin;
}

TProfile *ReadMax(TString filename) {
        TString filenameMax = filename;
        filenameMax.ReplaceAll("Tot","Max");
        TFile *filemax = new TFile(filenameMax);
        TProfile *pymax = (TProfile*)filemax->Get("pfxmeas");
	pymax->SetName("pymax");
	return pymax;
}

TProfile *ReadMin(TString filename) {
        TString filenameMin = filename;
        filenameMin.ReplaceAll("Tot","Min");
        TFile *filemin = new TFile(filenameMin);
        TProfile *pymin = (TProfile*)filemin->Get("pfxmeas");
	pymin->SetName("pymin");
	return pymin;
}

TProfile *ReadDiff(TString filename) {
        TString filenameDiff = filename;
        filenameDiff.ReplaceAll("Tot","Diff");
        TFile *filediff = new TFile(filenameDiff);
        TProfile *pydiff = (TProfile*)filediff->Get("pfxmeas");
	pydiff->SetName("pydiff");
	return pydiff;
}

//TH1F* GetMaxMDiff(TString filename = "Unfolding_TranTotNtrkJPCharged_NFweight_Baye5_McPt02.root") {
TH1F* GetMaxMMin(TString filename = "Unfolding_TranTotNtrkJPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root") {

	TFile *file = new TFile(filename);
	TH2F *hmeas = (TH2F*)file->Get("hmeas");

	TH1F *hdiff = (TH1F*)hmeas->ProjectionX("hdiff");
	hdiff->Reset();
	hdiff->SetTitle("TransMax - TransMin");

	TH1F *hmax = (TH1F*)hdiff->Clone("hmax");
	hmax->SetTitle("TransMax");

	TH1F *hmin = (TH1F*)hdiff->Clone("hmin");
	hmin->SetTitle("TransMin");

	for(int j = 1; j<=hmeas->GetNbinsX(); j++) {
		TH1F *hmeas_py = (TH1F*)hmeas->ProjectionY(Form("hmeas_py%d",j),j,j);
		hdiff->SetBinContent(j, Emean(hmeas_py));
		hmax->SetBinContent(j, EMax(hmeas_py));
		hmin->SetBinContent(j, EMin(hmeas_py));
	}



	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	hdiff->Scale(DeDpNorma);
	hmax->Scale(DeDpNorma);
	hmin->Scale(DeDpNorma);

	TCanvas *c1 = new TCanvas("c1","c1",1000,800);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	c1->GetFrame()->SetLineWidth(2);
	c1->SetLeftMargin(0.12);
	c1->SetBottomMargin(0.12);

	//hdiff->GetXaxis()->SetTitle("Detector-level Leading jet p_{T}(GeV/#it{c})");
	hdiff->GetXaxis()->SetTitle("Raw leading jet p_{T}(GeV/#it{c})");
	hdiff->SetMaximum(1.4);
	hdiff->SetMinimum(0);
	hmax->SetMaximum(1.7);
	hmax->SetMinimum(0);
	//hmax->GetXaxis()->SetTitle("Detector-level Leading jet p_{T}(GeV/#it{c})");
	hmax->GetXaxis()->SetTitle("Raw leading jet p_{T}(GeV/#it{c})");
//	hmax->GetYaxis()->SetTitle("Detector-level #LTN_{ch}/#delta#eta#delta#phi#GT");
	//hmax->GetYaxis()->SetTitle("Detector-level #LTdN_{ch}/d#etad#phi#GT");
	hmax->GetYaxis()->SetTitle("Raw #LTdN_{ch}/d#etad#phi#GT");
	hmax->GetXaxis()->SetTitleSize(0.05);
	hmax->GetXaxis()->SetLabelSize(0.05);
	hmax->GetYaxis()->SetTitleSize(0.05);
	hmax->GetYaxis()->SetLabelSize(0.05);
	hmax->GetYaxis()->SetTitleOffset(1);
	
	TH1D *htmp = (TH1D*)hmax->Clone("htmp");
	htmp->Reset();
	htmp->Draw();
	htmp->GetXaxis()->SetRangeUser(0,45);
	htmp->GetXaxis()->SetNdivisions(205);
	htmp->GetYaxis()->SetNdivisions(505);
	double XaxisMin = 3;
	double XaxisMax = 45;

	//hdiff->Draw("HIST");
	
	hmax->SetLineColor(kBlue);
	hmax->SetLineWidth(3);
	hmax->SetLineStyle(7);
	hmax->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
	hmax->Draw("sameHIST][");
	hmin->SetLineWidth(2);
	hmin->SetLineColor(kSpring-5);
	hmin->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
	hmin->Draw("sameHIST][");

	TH1F *diff_maxmin = Diff(filename);
	diff_maxmin->SetLineColor(1);
	diff_maxmin->SetMarkerColor(1);
	diff_maxmin->SetMarkerStyle(8);
	//diff_maxmin->Draw("same");
	diff_maxmin->Scale(2*DeDpNorma);

	TH1F *ratio_diff = (TH1F*)diff_maxmin->Clone("ratio_diff");
	ratio_diff->Divide(hdiff);
	ratio_diff->SetTitle("<TransMax>-<TransMin> Measured/Calculated");
	ratio_diff->GetYaxis()->SetTitle("TransDiff Measured/Calculated");
	PrintHistWErr(ratio_diff);

	TProfile *pmax = ReadMax(filename);
	pmax->Scale(2*DeDpNorma);
	pmax->SetLineColor(hmax->GetLineColor());
	pmax->SetMarkerColor(hmax->GetLineColor());
	pmax->SetMarkerStyle(4);
	pmax->SetMarkerSize(2);
	pmax->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
	pmax->Draw("sameEX0");
	TProfile *pmin = ReadMin(filename);
	pmin->Scale(2*DeDpNorma);
	pmin->SetLineColor(hmin->GetLineColor());
	pmin->SetMarkerColor(hmin->GetLineColor());
	pmin->SetMarkerStyle(8);
	pmin->SetMarkerSize(2);
	pmin->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
	pmin->Draw("sameEX0");
	TProfile *pdiff = ReadDiff(filename);
	pdiff->Scale(2*DeDpNorma);
	pdiff->SetMarkerColor(kYellow);
	pdiff->SetLineColor(kYellow);
	pdiff->SetMarkerStyle(23);
	//pdiff->Draw("same");

	TLegend *leg = new TLegend(0.485,0.58,0.885,0.89);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.043);
	//leg->AddEntry(pmax,"Pythia TransMax");
	//leg->AddEntry(pmin,"Pythia TransMin");
	//leg->AddEntry(pdiff,"Pythia <TransMax-TransMin>");
	//leg->AddEntry(diff_maxmin,"Pythia <TransMax>-<TransMin>");
	leg->AddEntry(pmax,"Measured TransMax","p");
	leg->AddEntry(pmin,"Measured TransMin","p");
	//leg->AddEntry(pdiff,"Measured <TransMax-TransMin>");
	//leg->AddEntry(diff_maxmin,"Measured <TransMax>-<TransMin>");
	leg->AddEntry(hmax,"Calculated TransMax","l");
	leg->AddEntry(hmin,"Calculated TransMin","l");
	//leg->AddEntry(hdiff,"Calculated TransMax-TransMin");
	leg->Draw();


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextFont(42);
	lat->DrawLatex(0.16,0.8,"#splitline{p+p@200 GeV}{p_{T} > 0.2 GeV/#it{c}, |#eta|<1}");
	lat->SetTextColor(1);
	lat->SetTextFont(62);
	lat->DrawLatex(0.16,0.16,"STAR");


        time_t rawtime;
        struct tm *timeinfo;
        char buffer4time[8];
        time(&rawtime);
        timeinfo=localtime(&rawtime);
        strftime(buffer4time,8,"%y%m%d",timeinfo);
        puts(buffer4time);



	if(1) {
		//c1->SaveAs("/Users/liyi/Research/Underlying/PaperDraft180920/TranMaxVsMin.pdf");
		c1->SaveAs(Form("fig/TranMaxVsMin_%s.pdf",buffer4time));
		c1->SaveAs("/Users/liyi/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/TranMaxVsMin.pdf");
	}


	return hdiff;

}



