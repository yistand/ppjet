#include <fstream>
#include "plothrecoWunfolderr.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 
							// void graphSystBand(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
							//                    Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBand=5,
							//                    Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1) 
							// void graphSystBox(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
							// 		  Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBox=5,
							// 		  Double_t dxaxis=1, Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1, Int_t fillstyleBox=1001)
							// 

TGraphAsymmErrors *Graph4UnfoldXErr(TGraphAsymmErrors *gas, TProfile *pf) {

	TGraphAsymmErrors *gr = new TGraphAsymmErrors();
	gr->SetName(Form("%s_xerr",gas->GetName()));

	gr->SetPoint(0,pf->GetBinLowEdge(1),gas->GetY()[0]);
	gr->SetPointEYlow(0,gas->GetErrorYlow(0));
	gr->SetPointEYhigh(0,gas->GetErrorYhigh(0));

	int gasNbins = pf->GetNbinsX();

	for(int i = 1; i<gasNbins; i++) {
		double tmp1, tmp2, tmpe1l, tmpe2l, tmpe1h, tmpe2h;
		tmp1 = gas->GetY()[i-1];
		tmp2 = gas->GetY()[i];
		tmpe1l = gas->GetErrorYlow(i-1);
		tmpe2l = gas->GetErrorYlow(i);
		tmpe1h = gas->GetErrorYhigh(i-1);
		tmpe2h = gas->GetErrorYhigh(i);
		double tmpmid = (tmp1+tmp2)/2.;
		double tmperr = (fabs(tmp1-tmp2))/2.;
		double tmperrl = (tmp1>tmp2?tmperr+tmpe2l:tmperr+tmpe1l);
		double tmperrh = (tmp1>tmp2?tmperr+tmpe1h:tmperr+tmpe2h);
		gr->SetPoint(i,pf->GetBinLowEdge(i+1),tmpmid);
		gr->SetPointEYlow(i,tmperrl);
		gr->SetPointEYhigh(i,tmperrh);
	}

	gr->SetPoint(gasNbins,pf->GetBinLowEdge(gasNbins)+pf->GetBinWidth(gasNbins),gas->GetY()[gasNbins-1]);
	gr->SetPointEYlow(gasNbins,gas->GetErrorYlow(gasNbins-1));
	gr->SetPointEYhigh(gasNbins,gas->GetErrorYhigh(gasNbins-1));

	return gr;

}


//TGraphErrors *Graph4UnfoldXErr(TProfile *pf) {
//
//	TGraphErrors *gr = new TGraphErrors();
//	gr->SetName(Form("%s_xerr",pf->GetName()));
//
//	gr->SetPoint(0,pf->GetBinLowEdge(1),pf->GetBinContent(1));
//	gr->SetPointError(0,0,0);
//
//	int pfNbins = pf->GetNbinsX();
//
//	for(int i = 1; i<pfNbins; i++) {
//		double tmp1, tmp2;
//		tmp1 = pf->GetBinContent(i);
//		tmp2 = pf->GetBinContent(i+1);
//		double tmpmid = (tmp1+tmp2)/2.;
//		double tmperr = fabs(tmp1-tmp2)/2.;
//		gr->SetPoint(i,pf->GetBinLowEdge(i+1),tmpmid);
//		gr->SetPointError(i,0,tmperr);
//	}
//
//	gr->SetPoint(pfNbins,pf->GetBinLowEdge(pfNbins)+pf->GetBinWidth(pfNbins),pf->GetBinContent(pfNbins));
//	gr->SetPointError(pfNbins,0,0);
//
//	return gr;
//
//}


TGraphAsymmErrors *Graph4UnfoldXErr(TProfile *pf) {

	TGraphAsymmErrors *gr = new TGraphAsymmErrors();
	gr->SetName(Form("%s_xerr",pf->GetName()));

	gr->SetPoint(0,pf->GetBinLowEdge(1),pf->GetBinContent(1));
	gr->SetPointEYlow(0,0);
	gr->SetPointEYhigh(0,0);

	int pfNbins = pf->GetNbinsX();

	for(int i = 1; i<pfNbins; i++) {
		double tmp1, tmp2;
		tmp1 = pf->GetBinContent(i);
		tmp2 = pf->GetBinContent(i+1);
		double tmpmid = (tmp1+tmp2)/2.;
		double tmperr = fabs(tmp1-tmp2)/2.;
		gr->SetPoint(i,pf->GetBinLowEdge(i+1),tmpmid);
		gr->SetPointEYlow(i,tmperr);
		gr->SetPointEYhigh(i,tmperr);
	}

	gr->SetPoint(pfNbins,pf->GetBinLowEdge(pfNbins)+pf->GetBinWidth(pfNbins),pf->GetBinContent(pfNbins));
	gr->SetPointEYlow(pfNbins,0);
	gr->SetPointEYhigh(pfNbins,0);

	return gr;

}






//======================================================================================================================================================

void PlotAllwSysErr(const char *Variable = "Ntrk", const char *filetag="NFWeight_BT170928_RcVzW_12JetBinv2_McPt02") {		// Ntrk, PtAve, PtSum

	int savefig = 0;
	int savecsv = 0;	// save output numbers to csv file
	int detplot = 0; 	// plot detector-level (no save)
	int saveroot = 0;

	int opt_drawdefpythia6 = 1;	// if 1, draw default pythia6 Perugia 2012 by Kevin

	const int NRegion = 3;
	const char *legtag[NRegion] = {"Toward","Away","Transverse"};

	//const char *Variable = "Ntrk";			// Ntrk, PtAve
	//const char *Region[NRegion] = {"LeadArea","SubArea","TranTot"};	// they are actually normalized to density

	//const char *Variable = "PtAve";			// Ntrk, PtAve
	const char *Region[NRegion] = {"Lead","Sub","Tran"};	// they are actually normalized to density

	double XaxisMin = 3; //2;
	double XaxisMax = 45;
	double YaxisMin = 0;
	double YaxisMax = 3.2;

	if(!strcmp(Variable,"PtSum")) YaxisMax = 7.5;

	int DefaultTh = 5;

	// Measurement data + STAR PYTHIA 6 Perugia2012
	TFile *f[NRegion];
	TProfile *hrecopf[NRegion];
	TProfile *hmeaspf[NRegion];
	TProfile *httpf[NRegion];
	TProfile *htpf[NRegion];
	TGraphAsymmErrors *grreco[NRegion];
       	for(int i = 0; i<NRegion; i++) {
		bool openi = false;
		//f[i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_NFWeight_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//f[i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_NFWeight_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//f[i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_NFWeight_BT170928_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		if(strcmp(Region[i],"Tran")==0 && (strcmp(Variable,"Ntrk")==0||strcmp(Variable,"PtSum")==0) ) {
			f[i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
		}
		else if(strcmp(Region[i],"Tran")!=0 && strcmp(Variable,"PtAve")!=0) {
			f[i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
		}
		else {
			f[i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
		}
		if(f[i]->IsOpen()) {
			openi=true;
		}
		//else {
		//	//f[i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_NFWeight_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//	//f[i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_NFWeight_BT170928_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//	f[i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
		//	if(f[i]->IsOpen()) {
		//		openi=true;
		//	}
		//	else {
		//		//f[i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_NFWeight_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//		//f[i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_NFWeight_BT170928_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		//		f[i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
		//		if(f[i]->IsOpen()) {
		//			openi=true;
		//		}
		//	}
		//}
		if(openi) {
			cout<<"Read "<<f[i]->GetName()<<endl;
			hrecopf[i] = (TProfile*)f[i]->Get(Form("hreco_default%d_pfx",DefaultTh));
			hmeaspf[i] = (TProfile*)f[i]->Get(Form("hmeas_pfx"));
			httpf[i] = (TProfile*)f[i]->Get(Form("htraintrue_pfx"));
			htpf[i] = (TProfile*)f[i]->Get(Form("htrain_pfx"));
			grreco[i] = (TGraphAsymmErrors*)f[i]->Get(Form("gSysErr%s",legtag[i]));
		}
		else {		// If cannot find the sys. err. one, open the unfold result file without sys. err. 
			//f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
			//f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_BinByBin_McPt02_embedMB.root",Region[i],Variable));
			//f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_2HalfMcPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
			//f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
			//f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_BT170928_12JetBinv2_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
			f[i] = new TFile(Form("Unfolding_%s%sJPCharged_%s_embedMB_Baye%d.root",Region[i],Variable,filetag,DefaultTh));
			if(f[i]->IsOpen()) {
				cout<<"Instead Open "<<f[i]->GetName()<<endl;
				TH2F *hreco = (TH2F*)f[i]->Get("hreco");
				TH2F *hmeas = (TH2F*)f[i]->Get("hmeas");
				TH2F *htrain = (TH2F*)f[i]->Get("htrain");
				TH2F *htraintrue = (TH2F*)f[i]->Get("htraintrue");

				if(hreco)
				hrecopf[i]=(TProfile*)hreco->ProfileX(Form("hreco_default%d_pfx",DefaultTh));
				else cout<<"Cannot find hreco"<<endl;
				if(hmeas)
				hmeaspf[i]=(TProfile*)hmeas->ProfileX(Form("hmeas_pfx"));
				else cout<<"Cannot find hmeas"<<endl;
				if(htraintrue)
				httpf[i]=(TProfile*)htraintrue->ProfileX(Form("htraintrue_pfx"));
				else cout<<"Cannot find htraintrue"<<endl;
				if(htrain)
				htpf[i]=(TProfile*)htrain->ProfileX(Form("htrain_pfx"));
				else cout<<"Cannot find htrain"<<endl;

				grreco[i] = NULL;
			}
		}
	}

	// PYTHIA 6 default Perugia 2012 without STAR's tune
	TFile *fp6def = new TFile("Profile_12JetBinv2_FullJet_TransCharged_pythia6default_pp200hard_PionDecayOff_180324.root");
	TProfile *httpf6def[NRegion];
	for(int i = 0 ;i <NRegion; i++) {
		httpf6def[i] = (TProfile*)fp6def->Get(Form("%s%s",Region[i],Variable));
	}


	// PYTHIA 8
	TFile *fp8 = new TFile("Profile_12JetBinv2_FullJet_TransCharged_pythia8215_pp200hard_PionDecayOff_seed134123_170422.root");//ProfileFullJet_TransCharged_pythia8215_pp200hard_PionDecayOff_seed134123_170422.root");
	TProfile *httpf8[NRegion];
	for(int i = 0 ;i <NRegion; i++) {
		httpf8[i] = (TProfile*)fp8->Get(Form("%s%s",Region[i],Variable));
	}


	// Unfolding x-position uncertainty
	TGraphAsymmErrors *xgr[NRegion];
	for(int i = 0; i<NRegion; i++) {
		if(grreco[i])
		xgr[i] = (TGraphAsymmErrors*) Graph4UnfoldXErr(grreco[i],hrecopf[i]);
		else 
		xgr[i] = NULL;
	}

	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	int mstyle[NRegion] = {8,21,24};
	int mcolor[NRegion] = {1,kBlue+1,2};
	float msize[NRegion] = {2,2,2};
	int lstyle[NRegion] = {1,7,3};			// Line Style
	int bstyle[NRegion] = {3001,3144,1001};			// Box Style
	int bcolor[NRegion] = {kGray, kBlue-10, kMagenta-10};

	for(int i = 0; i<NRegion; i++) {
		SetHistStyle(hrecopf[i],mcolor[i],mcolor[i],mstyle[i],1,msize[i],1);
		SetHistStyle(httpf[i],mcolor[i],mcolor[i],mstyle[i],lstyle[i],msize[i],7);	// last option is line width	// last option is line width
		SetHistStyle(httpf6def[i],mcolor[i],mcolor[i],mstyle[i],lstyle[i],msize[i],3);
		SetHistStyle(httpf8[i],mcolor[i],mcolor[i],mstyle[i],lstyle[i],msize[i],1);
		if(!xgr[i]) continue;
		xgr[i]->SetFillColor(bcolor[i]);
		xgr[i]->SetLineColor(bcolor[i]);
		xgr[i]->SetMarkerColor(bcolor[i]);
	}

	hrecopf[0]->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/#it{c})");
	if(!strcmp(Variable,"PtAve")) {
		hrecopf[0]->GetYaxis()->SetTitle("#LTp_{T}^{ch}#GT");
	}
	else if(!strcmp(Variable,"PtSum")) {
		//hrecopf[0]->GetYaxis()->SetTitle("#LT#sump_{T}^{ch}/#delta#eta#delta#phi#GT");
		hrecopf[0]->GetYaxis()->SetTitle("#LTd#scale[0.5]{#sum}p_{T}/d#etad#phi#GT");
	}
	else {		// default: Ntrk
		//hrecopf[0]->GetYaxis()->SetTitle("#LTN_{ch}/#delta#eta#delta#phi#GT");
		hrecopf[0]->GetYaxis()->SetTitle("#LTdN_{ch}/d#etad#phi#GT");
	}
	hrecopf[0]->SetMinimum(YaxisMin);
	hrecopf[0]->SetMaximum(YaxisMax);

	TH1D *htmp = (TH1D*) hrecopf[0]->Clone("htmp");
	htmp->Reset();
	htmp->GetXaxis()->SetRangeUser(0,45);
	htmp->GetXaxis()->SetNdivisions(205);
	htmp->Draw();

	for(int i = 0; i<NRegion; i++) {
	//for(int i = 2; i<NRegion; i++) {
		if(!grreco[i]) continue;

		double *ex;
		if(hrecopf[i]) {
			ex = new double[hrecopf[i]->GetNbinsX()];
			for( int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
				ex[j] = hrecopf[i]->GetBinWidth(j+1)/2.;
			}
		}

		double *eylow;
		eylow = new double[grreco[i]->GetN()];
		for(int j = 0; j<grreco[i]->GetN(); j++) {
			eylow[j] = - grreco[i]->GetErrorYlow(j);
		}
		grreco[i]->SetFillColor(kGray+i);

		// find the x-axis to plot points
		double *tmpx=grreco[i]->GetX();
		int skip=0;
		int total=0;
		for(int jg = 0; jg<grreco[i]->GetN(); jg++) {
			if(tmpx[jg]<XaxisMin) {
				skip++;
			}
			else {
				total++;
			}
			if(tmpx[jg]>XaxisMax) {
				total--;
				break;
			}
		}
		//cout<<"skip = "<<skip<<" total = "<<total<<endl;
		//cout<<"start at "<<tmpx[skip]<<" end at "<<tmpx[total+skip]<<endl;
		//graphSystBand(total,grreco[i]->GetX()+skip,grreco[i]->GetY()+skip,0,0, eylow+skip,  grreco[i]->GetEYhigh()+skip,bcolor[i], 0,0);
		graphSystBox(total,grreco[i]->GetX()+skip,grreco[i]->GetY()+skip,ex+skip,0, eylow+skip,  grreco[i]->GetEYhigh()+skip,bcolor[i], 0,0,1,0,1,bstyle[i]);	
	}
	if(grreco[0]) 
	{
		double *eylow;
		eylow = new double[grreco[0]->GetN()];
		for(int j = 0; j<grreco[0]->GetN(); j++) {
			eylow[j] = - grreco[0]->GetErrorYlow(j);
		}
		grreco[0]->SetFillColor(kGray+0);
		//graphSystBand(grreco[0]->GetN()-1,grreco[0]->GetX()+1,grreco[0]->GetY()+1,0,0, eylow+1,  grreco[0]->GetEYhigh()+1,bcolor[0], 0,0);
	}


	for(int i = 0; i<NRegion; i++) {
		if(!xgr[i]) continue;
		double *eylow;
		eylow = new double[xgr[i]->GetN()];
		for(int j = 0; j<xgr[i]->GetN(); j++) {
			eylow[j] = - xgr[i]->GetErrorYlow(j);
		}
		xgr[i]->SetMarkerStyle(8);
		//xgr[i]->Draw("psame");
		//graphSystBand(xgr[i]->GetN()-1,xgr[i]->GetX()+1,xgr[i]->GetY()+1,0,0, eylow+1,  xgr[i]->GetEYhigh()+1,bcolor[i], 0,0);	// make sure use the same color as greco...
	}
	if(xgr[0])
	{
		double *eylow;
		eylow = new double[xgr[0]->GetN()];
		for(int j = 0; j<xgr[0]->GetN(); j++) {
			eylow[j] = - xgr[0]->GetErrorYlow(j);
		}
		xgr[0]->SetMarkerStyle(8);
		//xgr[0]->Draw("psame");
		//graphSystBand(xgr[0]->GetN()-1,xgr[0]->GetX()+1,xgr[0]->GetY()+1,0,0, eylow+1,  xgr[0]->GetEYhigh()+1,bcolor[0], 0,0);	// make sure use the same color as greco...
	}

	for(int i=0;i<NRegion; i++) {		// Data unfolded
		hrecopf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
		hrecopf[i]->Draw("peX0same");
	}

	
	//hrecopf[2]->Draw("peX0same");//TEST

	// Normalized ProfileX to density
	if(!strcmp(Variable,"Ntrk")||!strcmp(Variable,"PtSum")) {
		float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
		for(int i = 0; i<NRegion; i++) {
			hmeaspf[i]->Scale(DeDpNorma);
			httpf[i]->Scale(DeDpNorma);
			htpf[i]->Scale(DeDpNorma);
			httpf6def[i]->Scale(DeDpNorma);
			httpf8[i]->Scale(DeDpNorma);
		}
	}


	for(int i=0;i<NRegion; i++) {		// PYTHIA
		httpf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
		httpf6def[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
		httpf8[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);

		httpf[i]->Draw("histcsame");
		if(opt_drawdefpythia6)  {
			httpf6def[i]->SetName(Form("%s_p6def",httpf6def[i]->GetName()));
			httpf6def[i]->Draw("histcsame");
		}
		httpf8[i]->SetName(Form("%s_p8",httpf8[i]->GetName()));
		httpf8[i]->Draw("histcsame");
		//httpf[i]->Draw("HISTsame");
		//httpf8[i]->Draw("HISTsame");
	}
	
	//httpf[2]->Draw("histcsame");//TEST


	TLegend *leg;
	leg = new TLegend(0.15,0.62,0.42,0.88);	
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	//leg->AddEntry(hrecopf[2],legtag[2],"p");//TEST
	for(int i = 0; i<NRegion; i++) {
		leg->AddEntry(hrecopf[i],legtag[i],"p");
	}
	leg->AddEntry(httpf[0],"Perugia 2012 (STAR)","l");
	//leg->AddEntry(httpf[0],"STAR Tune","l");
	if(opt_drawdefpythia6)  leg->AddEntry(httpf6def[0],"Perugia 2012","l");
	//leg->AddEntry(httpf8[0],"pythia 8.215","l");
	leg->AddEntry(httpf8[0],"Monash 2013","l");
	leg->Draw();

	TLegend *leg2;
	leg2 = new TLegend(0.42,0.55,0.69,0.88);	
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->SetHeader("perugia 2012");
	for(int i = 0; i<NRegion; i++) {
		leg2->AddEntry(httpf[i],Form("MC %s",legtag[i]),"l");
	}
	//leg2->Draw();

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextFont(42);
	if(!strcmp(Variable,"PtAve")) {
		//lat->DrawLatex(0.56,0.19,"#splitline{p+p@200 GeV}{#splitline{p_{T} > 0.2 GeV/#it{c}, |#eta|<1}{R = 0.6, |#eta_{jet}|<0.4}}");
		lat->DrawLatex(0.37,0.19,"p+p@200 GeV");
		lat->DrawLatex(0.56,0.18,"#splitline{p_{T} > 0.2 GeV/#it{c}, |#eta|<1}{R = 0.6, |#eta_{jet}|<0.4}");
	}
	else {
		lat->DrawLatex(0.56,0.8,"#splitline{p+p@200 GeV}{#splitline{p_{T} > 0.2 GeV/#it{c}, |#eta|<1}{R = 0.6, |#eta_{jet}|<0.4}}");
	}
	lat->SetTextColor(1);
	lat->SetTextFont(62);
	if(!strcmp(Variable,"PtSum")) {
		lat->DrawLatex(0.75,0.2,"STAR");
	}
	else {
		lat->DrawLatex(0.16,0.16,"STAR");
	}

	if(savefig) {
		//c->SaveAs(Form("/Users/li/Research/Underlying/PaperDraft170405/%s_DataVsPythia6Vs8_Box.pdf",Variable));
		c->SaveAs(Form("/Users/liyi/Research/Underlying/PaperDraft180920/%s_DataVsPythia6Vs8_Box.pdf",Variable));
		c->SaveAs(Form("/Users/liyi/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/%s_DataVsPythia6Vs8_Box.pdf",Variable));
		//c->SaveAs(Form("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/NoRcVzW_%s_DataVsPythia6Vs8_Box.pdf",Variable));
	}
	if(saveroot) {
		c->SaveAs(Form("/Users/liyi/Research/Underlying/PaperDraft180920/%s_DataVsPythia6Vs8_Box.root",Variable));
	}

	if(detplot) {
		TCanvas *c2 = new TCanvas("c2","c2",1000,800);
		c2->SetLeftMargin(0.12);
		c2->SetBottomMargin(0.12);
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		
		hmeaspf[0]->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/#it{c})");
		if(!strcmp(Variable,"PtAve")) {
			hmeaspf[0]->GetYaxis()->SetTitle("#LTp_{T}^{ch}#GT");
		}
		else if(!strcmp(Variable,"PtSum")) {
			//hmeaspf[0]->GetYaxis()->SetTitle("#LT#sump_{T}/#delta#eta#delta#phi#GT");
			hmeaspf[0]->GetYaxis()->SetTitle("#LTd#scale[0.5]{#sum}p_{T}/d#etad#phi#GT");
		}
		else {		// default: Ntrk
			//hmeaspf[0]->GetYaxis()->SetTitle("#LTN_{ch}/#delta#eta#delta#phi#GT");
			hmeaspf[0]->GetYaxis()->SetTitle("#LTdN_{ch}/d#etad#phi#GT");
		}
		hmeaspf[0]->SetMinimum(YaxisMin);
		hmeaspf[0]->SetMaximum(YaxisMax);


		TH1D *htmp2 = (TH1D*) hmeaspf[0]->Clone("htmp2");
		htmp2->Reset();
		htmp2->GetXaxis()->SetRangeUser(0,45);
		htmp2->GetXaxis()->SetNdivisions(205);
		htmp2->Draw();


		for(int i=0;i<NRegion; i++) {
			SetHistStyle(hmeaspf[i],mcolor[i],mcolor[i],mstyle[i],1,msize[i],1);
			SetHistStyle(htpf[i],mcolor[i],mcolor[i],mstyle[i],lstyle[i],msize[i],2);
			hmeaspf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			hmeaspf[i]->Draw("peX0same");
			htpf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			htpf[i]->Draw("histcsame");
		}


		TLegend *leg_det;
		leg_det = new TLegend(0.15,0.62,0.42,0.88);	
		leg_det->SetFillColor(0);
		leg_det->SetBorderSize(0);
		for(int i = 0; i<NRegion; i++) {
			leg_det->AddEntry(hmeaspf[i],legtag[i],"p");
		}
		leg_det->AddEntry(httpf[0],"Perugia 2012","l");
		leg_det->Draw();

		leg_det->Draw();

		TLatex *laxmeas = new TLatex(0.15,0.92,"Detector-level");
		laxmeas->SetTextFont(42);
		laxmeas->SetNDC();
		laxmeas->Draw();
	}


	// Print out data point
	for(int i=0;i<NRegion; i++) {
		if(!grreco[i]) continue;
		cout<<legtag[i]<<": "<<endl;
		for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
			cout<<hrecopf[i]->GetBinCenter(j+1)<<" "<<hrecopf[i]->GetBinContent(j+1)<<" +/- "<<hrecopf[i]->GetBinError(j+1)<<" (stat);  "<<grreco[i]->GetY()[j]<<" +"<<grreco[i]->GetEYhigh()[j]<<" - "<<grreco[i]->GetEYlow()[j]<<" (sys)"<<"\t\t\t PYTHIA6(star) = "<<httpf[i]->GetBinContent(j+1);
			if(opt_drawdefpythia6) cout<<"\t\t\t PYTHIA6(def) = "<<httpf6def[i]->GetBinContent(j+1);
			cout<<"\t\t\t PYTHIA8 = "<<httpf8[i]->GetBinContent(j+1)<<endl;
		}
	}


	ofstream tout;
	TString Stout=Form("%s_JP_Charged_%s_embedMB_Baye%d.csv",Variable,filetag,DefaultTh);
	if(savecsv) {
		tout.open(Stout);
		if(tout.is_open()) {
			cout<<"Output to "<<Stout<<endl;
		}
		else {
			cout<<"Error open file "<<Stout<<" to write"<<endl;
			return;
		}

		tout<<"pt,";
		for(int j = 0; j<hrecopf[0]->GetNbinsX(); j++) {
			tout<<hrecopf[0]->GetBinLowEdge(j+1)<<"-"<<hrecopf[0]->GetBinLowEdge(j+1)+hrecopf[0]->GetBinWidth(j+1)<<",";
		}
		tout<<endl;
		for(int i=0;i<NRegion; i++) {
			if(!grreco[i]) continue;
			tout<<legtag[i]<<",";
			for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
				tout<<hrecopf[i]->GetBinContent(j+1)<<" +/- "<<hrecopf[i]->GetBinError(j+1)<<" (stat) + "<<grreco[i]->GetEYhigh()[j]<<" - "<<grreco[i]->GetEYlow()[j]<<" (sys)"<<",";
			}
			tout<<endl;
		}
		tout<<"PYTHIA6(star)"<<endl;
		for(int i=0;i<NRegion; i++) {
			if(!grreco[i]) continue;
			tout<<legtag[i]<<",";
			for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
				tout<<httpf[i]->GetBinContent(j+1)<<",";
			}
			tout<<endl;
		}
		if(opt_drawdefpythia6) {
			tout<<"PYTHIA6(def)"<<endl;
			for(int i=0;i<NRegion; i++) {
				if(!grreco[i]) continue;
				tout<<legtag[i]<<",";
				for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
					tout<<httpf6def[i]->GetBinContent(j+1)<<",";
				}
				tout<<endl;
			}
		}
		tout<<"PYTHIA8"<<endl;
		for(int i=0;i<NRegion; i++) {
			if(!grreco[i]) continue;
			tout<<legtag[i]<<",";
			for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
				tout<<httpf8[i]->GetBinContent(j+1)<<",";
			}
			tout<<endl;
		}

		tout.close();
	}


}


