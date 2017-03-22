#include "plothrecoWunfolderr.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 
							// void graphSystBand(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
							//                    Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBand=5,
							//                    Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1) 


void graphSystBox(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
		  Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBox=5,
		  Double_t dxaxis=1, Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
   int utilityPrint =0;
  if (utilityPrint) printf("\n__________graphSystBox__________symb=%d sSize=%f color=%d colorBox=%d__________\n",symb,sSize,color,colorBox);
  if (utilityPrint) {
    for (int i=0; i<n; ++i) {
      printf("x[%d]=%f; y[%d]=%f; ",i,x[i],i,y[i]);
      if (utilityPrint<0) {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]/x[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]/y[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]/y[i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]/y[i]);
      }
      else {
	if (ex) printf("ex[%d]=%f; ",i,ex[i]);
	if (ey) printf("ey[%d]=%f; ",i,ey[i]);
	if (syPlus ) printf("sy1[%d]=%f; ",i,syPlus [i]);
	if (syMinus) printf("sy2[%d]=%f; ",i,syMinus[i]);
      }
      printf("\n");
    }
  }
  if (syPlus) {
    float scale=0.015*sSize;
    if (scale<0.01) scale=0.01;
    float dx=dxaxis*scale;
    float xbox[5], ybox[5];
    for (int i=0; i<n; ++i) {
      float y1=y[i]+syPlus[i];
      float y2;
      if (syMinus) y2=y[i]+syMinus[i];
      else         y2=y[i]-syPlus[i];
      xbox[0]=x[i]-dx;ybox[0]=y1;
      xbox[1]=x[i]+dx;ybox[1]=y1;
      xbox[2]=x[i]+dx;ybox[2]=y2;
      xbox[3]=x[i]-dx;ybox[3]=y2;
      xbox[4]=x[i]-dx;ybox[4]=y1;
      TGraph *grafbox=new TGraph(5,xbox,ybox);
      if(colorBox>=0){
	grafbox->SetFillColor(colorBox);
	grafbox->Draw("f");
      }
      else{
	grafbox->SetLineWidth(-colorBox);
	grafbox->SetLineColor(color);
	grafbox->Draw("l");
      }
    }
  }
  if(!symb)return;
  TGraphErrors *graf = new TGraphErrors(n,x,y,ex,ey);
  graf->SetMarkerStyle(symb);
  graf->SetMarkerSize(sSize);
  graf->SetMarkerColor(color);
  graf->SetLineColor(color);
  graf->SetLineWidth(lwidth);
  graf->Draw("p");
}



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

void PlotAllwSysErr() {

	const int NRegion = 3;
	const char *legtag[NRegion] = {"Toward","Away","Transverse"};

	//const char *Variable = "Ntrk";			// Ntrk, PtAve
	//const char *Region[NRegion] = {"LeadArea","SubArea","TranTot"};	// they are actually normalized to density

	const char *Variable = "PtAve";			// Ntrk, PtAve
	const char *Region[NRegion] = {"Lead","Sub","Tran"};	// they are actually normalized to density

	const double XaxisMin = 2;
	const double XaxisMax = 45;

	int DefaultTh = 5;
	TFile *f[NRegion];
	TProfile *hrecopf[NRegion];
	TProfile *hmeaspf[NRegion];
	TProfile *httpf[NRegion];
	TProfile *htpf[NRegion];
	TGraphAsymmErrors *grreco[NRegion];
       	for(int i = 0; i<NRegion; i++) {
		f[i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_NFWeight_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
		if(f[i]->IsOpen()) {
			hrecopf[i] = (TProfile*)f[i]->Get(Form("hreco_default%d_pfx",DefaultTh));
			hmeaspf[i] = (TProfile*)f[i]->Get(Form("hmeas_pfx"));
			httpf[i] = (TProfile*)f[i]->Get(Form("htraintrue_pfx"));
			htpf[i] = (TProfile*)f[i]->Get(Form("htrain_pfx"));
			grreco[i] = (TGraphAsymmErrors*)f[i]->Get(Form("gSysErr%s",legtag[i]));
		}
		else {		// If cannot find the sys. err. one, open the unfold result file without sys. err. 
			f[i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh));
			if(f[i]->IsOpen()) {
				cout<<"Instead Open "<<Form("Unfolding_%s%sJPCharged_NFWeight_McPt02_embedMB_Baye%d.root",Region[i],Variable,DefaultTh)<<endl;
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
	
	int mstyle[NRegion] = {8,21,8};
	int mcolor[NRegion] = {1,kBlue+1,2};
	float msize[NRegion] = {2,2,2};
	int lstyle[NRegion] = {1,7,3};
	int bcolor[NRegion] = {kGray, kBlue-10, kMagenta-10};
	for(int i = 0; i<NRegion; i++) {
		SetHistStyle(hrecopf[i],mcolor[i],mcolor[i],mstyle[i],1,msize[i],1);
		SetHistStyle(httpf[i],mcolor[i],mcolor[i],mstyle[i],lstyle[i],msize[i],2);
		if(!xgr[i]) continue;
		xgr[i]->SetFillColor(bcolor[i]);
		xgr[i]->SetLineColor(bcolor[i]);
		xgr[i]->SetMarkerColor(bcolor[i]);
	}

	hrecopf[0]->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/#it{c})");
	if(!strcmp(Variable,"PtAve")) {
		hrecopf[0]->GetYaxis()->SetTitle("#LTp_{T}#GT");
	}
	else {		// default: Ntrk
		hrecopf[0]->GetYaxis()->SetTitle("#LTN_{ch}/#delta#eta#delta#phi#GT");
	}
	hrecopf[0]->SetMinimum(0);
	hrecopf[0]->SetMaximum(3.1);

	TH1D *htmp = (TH1D*) hrecopf[0]->Clone("htmp");
	htmp->Reset();
	htmp->GetXaxis()->SetRangeUser(0,45);
	htmp->GetXaxis()->SetNdivisions(5);
	htmp->Draw();

	for(int i = 0; i<NRegion; i++) {
		if(!grreco[i]) continue;
		double *eylow;
		eylow = new double[grreco[i]->GetN()];
		for(int j = 0; j<grreco[i]->GetN(); j++) {
			eylow[j] = - grreco[i]->GetErrorYlow(j);
		}
		grreco[i]->SetFillColor(kGray+i);
		graphSystBand(grreco[i]->GetN()-1,grreco[i]->GetX()+1,grreco[i]->GetY()+1,0,0, eylow+1,  grreco[i]->GetEYhigh()+1,bcolor[i], 0,0);
		//graphSystBox(grreco[i]->GetN(),grreco[i]->GetX(),grreco[i]->GetY(),0,0, eylow,  grreco[i]->GetEYhigh(),kGray+i, 50,0);
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

	for(int i=0;i<NRegion; i++) {
		hrecopf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
		hrecopf[i]->Draw("peX0same");
	}


	// Normalized ProfileX to density
	if(!strcmp(Variable,"Ntrk")) {
		float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
		for(int i = 0; i<NRegion; i++) {
			hmeaspf[i]->Scale(DeDpNorma);
			httpf[i]->Scale(DeDpNorma);
			htpf[i]->Scale(DeDpNorma);
		}
	}
	for(int i=0;i<NRegion; i++) {
		httpf[i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
		httpf[i]->Draw("histcsame");
	}


	TLegend *leg;
	leg = new TLegend(0.15,0.62,0.42,0.88);	
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	for(int i = 0; i<NRegion; i++) {
		leg->AddEntry(hrecopf[i],legtag[i],"p");
	}
	leg->AddEntry(httpf[0],"perugia 2012","l");
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



	// Print out data point
	for(int i=0;i<NRegion; i++) {
		if(!grreco[i]) continue;
		cout<<legtag[i]<<": "<<endl;
		for(int j = 0; j<hrecopf[i]->GetNbinsX(); j++) {
			cout<<hrecopf[i]->GetBinCenter(j+1)<<" "<<hrecopf[i]->GetBinContent(j+1)<<" +/- "<<hrecopf[i]->GetBinError(j+1)<<" (stat);  "<<grreco[i]->GetY()[j]<<" +"<<grreco[i]->GetEYhigh()[j]<<" - "<<grreco[i]->GetEYlow()[j]<<" (sys)"<<"\t\t\t PYTHIA = "<<httpf[i]->GetBinContent(j+1)<<endl;
		}
	}
}


