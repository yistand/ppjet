#include <time.h>       /* time_t, struct tm, time, localtime, strftime */
#include "plothrecoWunfolderr.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 
							// void graphSystBand(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
							//                    Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBand=5,
							//                    Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1) 
							// void graphSystBox(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
							// 		  Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBox=5,
							// 		  Double_t dxaxis=1, Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1, Int_t fillstyleBox=1001)
							// {

/*
void graphSystBox(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
		  Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBox=5,
		  Double_t dxaxis=1, Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1, Int_t fillstyleBox=1001)
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
      if(fabs(dxaxis)<1e-12) { 				// if dxaxis is 0, use bin width as dx
	      	dx = ex[i];
      }
      xbox[0]=x[i]-dx;ybox[0]=y1;
      xbox[1]=x[i]+dx;ybox[1]=y1;
      xbox[2]=x[i]+dx;ybox[2]=y2;
      xbox[3]=x[i]-dx;ybox[3]=y2;
      xbox[4]=x[i]-dx;ybox[4]=y1;
      TGraph *grafbox=new TGraph(5,xbox,ybox);
      if(colorBox>=0){
	grafbox->SetFillColor(colorBox);
	grafbox->SetFillStyle(fillstyleBox);
	grafbox->Draw("f");
      }
      else{
	grafbox->SetLineWidth(-colorBox);
	grafbox->SetLineColor(color);
	grafbox->SetFillStyle(fillstyleBox);
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
*/


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

void PlotTranMc02Vs05wSysErr(const char *Variable = "Ntrk", TString filetag="NFWeight_BT170928_RcVzW_12JetBinv2") {		// Ntrk, PtAve, PtSum

	int savefig = 1;
	int savecsv = 1;	// save output numbers to csv file
	int detplot = 0;
	int saveroot = 1;

	int opt_drawdefpythia6 = 1;	// if 1, draw default pythia6 Perugia 2012 by Kevin

	const int Npt = 2;
	//const char *pttag[Npt] = {"McPt02","McPtRC02MC05"};		// for measured data + default pythia6 
	const char *pttag[Npt] = {"McPt02","McPtRC05MC05"};		// for measured data + default pythia6 
	const char *p8pttag[Npt] = {"","_PT05"};				// for pythia8
	const char *ptname[Npt] = {"pt>0.2","pT>0.5"};		// for measured data + default pythia6 

	const int NRegion = 1;
	const char *legtag[NRegion] = {"Transverse"};
	//const char *legtagpt[Npt] = {"p_{T} > 0.2 GeV/#it{c}", "p_{T} > 0.5 GeV/#it{c}"};
	const char *legtagpt[Npt] = {"0.2 GeV/#it{c}", "0.5 GeV/#it{c}"};

	const char *Region[NRegion] = {"Tran"};	// they are actually normalized to density

	double XaxisMin = 3; //2;
	double XaxisMax = 45;
	double YaxisMin = 0;
	double YaxisMax = 1.7;// 3.1;

	if(!strcmp(Variable,"PtSum")) YaxisMax = 25;
	if(!strcmp(Variable,"PtAve")) YaxisMax = 2;

	int DefaultTh = 4;

	// Measurement data + STAR PYTHIA 6 Perugia2012
	TFile *f[Npt][NRegion];
	TProfile *hrecopf[Npt][NRegion];
	TProfile *hmeaspf[Npt][NRegion];
	TProfile *httpf[Npt][NRegion];
	TProfile *htpf[Npt][NRegion];
	TGraphAsymmErrors *grreco[Npt][NRegion];
       	for(int j = 0; j<Npt; j++) {
       		for(int i = 0; i<NRegion; i++) {
			bool openi = false;
			if(DefaultTh==4) {
				if(strcmp(Region[i],"Tran")==0 && (strcmp(Variable,"Ntrk")==0||strcmp(Variable,"PtSum")==0) ) {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_%s_%s_embedMB.root",Region[i],Variable,filetag.Data(),pttag[j]));
				}
				else if(strcmp(Region[i],"Tran")!=0 && strcmp(Variable,"PtAve")!=0) {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_%s_%s_embedMB.root",Region[i],Variable,filetag.Data(),pttag[j]));
				}
				else {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_%s_%s_embedMB.root",Region[i],Variable,filetag.Data(),pttag[j]));
				}
			}
			else {
				//f[j][i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_NFWeight_12JetBinv2_%s_embedMB_Baye%d.root",Region[i],Variable,pttag[j],DefaultTh));
				if(strcmp(Region[i],"Tran")==0 && (strcmp(Variable,"Ntrk")==0||strcmp(Variable,"PtSum")==0) ) {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
				}
				else if(strcmp(Region[i],"Tran")!=0 && strcmp(Variable,"PtAve")!=0) {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
				}
				else {
					f[j][i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
				}
			}
			if(f[j][i]->IsOpen()) {
				openi=true;
			}

			//f[j][i] = new TFile(Form("SysErr4Unfolding_%s%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
			//if(f[j][i]->IsOpen()) {
			//	openi=true;
			//}
			//else {
			//	//f[j][i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_NFWeight_12JetBinv2_%s_embedMB_Baye%d.root",Region[i],Variable,pttag[j],DefaultTh));
			//	f[j][i] = new TFile(Form("SysErr4Unfolding_%sArea%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
			//	if(f[j][i]->IsOpen()) {
			//		openi=true;
			//	}
			//	else {
			//		//f[j][i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_NFWeight_12JetBinv2_%s_embedMB_Baye%d.root",Region[i],Variable,pttag[j],DefaultTh));
			//		f[j][i] = new TFile(Form("SysErr4Unfolding_%sTot%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
			//		if(f[j][i]->IsOpen()) {
			//			openi=true;
			//		}
			//	}
			//}
			if(openi) {
				cout<<"Read "<<f[j][i]->GetName()<<endl;
				hrecopf[j][i] = (TProfile*)f[j][i]->Get(Form("hreco_default%d_pfx",DefaultTh));
				hmeaspf[j][i] = (TProfile*)f[j][i]->Get(Form("hmeas_pfx"));
				httpf[j][i] = (TProfile*)f[j][i]->Get(Form("htraintrue_pfx"));
				htpf[j][i] = (TProfile*)f[j][i]->Get(Form("htrain_pfx"));
				grreco[j][i] = (TGraphAsymmErrors*)f[j][i]->Get(Form("gSysErr%s",legtag[i]));
			}
			else {		// If cannot find the sys. err. one, open the unfold result file without sys. err. 
				if(DefaultTh==4) {
					f[j][i] = new TFile(Form("Unfolding_%s%sJPCharged_%s_%s_embedMB.root",Region[i],Variable,filetag.Data(),pttag[j]));
				}
				else {
					//f[j][i] = new TFile(Form("Unfolding_%s%sJPCharged_NFWeight_12JetBinv2_%s_embedMB_Baye%d.root",Region[i],Variable,pttag[j],DefaultTh));
					f[j][i] = new TFile(Form("Unfolding_%s%sJPCharged_%s_%s_embedMB_Baye%d.root",Region[i],Variable,filetag.Data(),pttag[j],DefaultTh));
				}
				if(f[j][i]->IsOpen()) {
					cout<<"Instead Open "<<f[j][i]->GetName()<<endl;
					TH2F *hreco = (TH2F*)f[j][i]->Get("hreco");
					TH2F *hmeas = (TH2F*)f[j][i]->Get("hmeas");
					TH2F *htrain = (TH2F*)f[j][i]->Get("htrain");
					TH2F *htraintrue = (TH2F*)f[j][i]->Get("htraintrue");

					if(hreco)
					hrecopf[j][i]=(TProfile*)hreco->ProfileX(Form("hreco_default%d_pfx",DefaultTh));
					else cout<<"Cannot find hreco"<<endl;
					if(hmeas)
					hmeaspf[j][i]=(TProfile*)hmeas->ProfileX(Form("hmeas_pfx"));
					else cout<<"Cannot find hmeas"<<endl;
					if(htraintrue)
					httpf[j][i]=(TProfile*)htraintrue->ProfileX(Form("htraintrue_pfx"));
					else cout<<"Cannot find htraintrue"<<endl;
					if(htrain)
					htpf[j][i]=(TProfile*)htrain->ProfileX(Form("htrain_pfx"));
					else cout<<"Cannot find htrain"<<endl;

					grreco[j][i] = NULL;
				}
			}
		}
	}


	// PYTHIA 6 default Perugia 2012 without STAR's tune
	TFile *fp6def[Npt];
	TProfile *httpf6def[Npt][NRegion];
	for(int j = 0; j<Npt; j++) {
		for(int i = 0 ;i <NRegion; i++) {
 			fp6def[j] = new TFile(Form("/Users/liyi/Research/Underlying/Profile_12JetBinv2_FullJet_TransCharged_pythia6default_pp200hard_PionDecayOff_180324%s.root",p8pttag[j]));	// same tag as p8 for "" and "_PT05"
			httpf6def[j][i] = (TProfile*)fp6def[j]->Get(Form("%s%s",Region[i],Variable));
		}
	}

	// PYTHIA 8
	TFile *fp8[Npt];
	TProfile *httpf8[Npt][NRegion];
	for(int j = 0; j<Npt; j++) {
		fp8[j] = new TFile(Form("/Users/liyi/Research/Underlying/Profile_12JetBinv2_FullJet_TransCharged_pythia8215_pp200hard_PionDecayOff_seed134123_170422%s.root",p8pttag[j]));//Profile_12JetBinv2_FullJet_TransCharged_pythia8215_pp200hard_PionDecayOff_seed134123_170422%s.root",p8pttag[j]));
		cout<<"Read from "<<fp8[j]->GetName()<<" for ";
		for(int i = 0 ;i <NRegion; i++) {
			httpf8[j][i] = (TProfile*)fp8[j]->Get(Form("%s%s",Region[i],Variable));
			cout<<httpf8[j][i]->GetName()<<endl;
		}
	}


	// Unfolding x-position uncertainty
	TGraphAsymmErrors *xgr[Npt][NRegion];
	for(int j = 0; j<Npt; j++) {
		for(int i = 0; i<NRegion; i++) {
			if(grreco[j][i])
			xgr[j][i] = (TGraphAsymmErrors*) Graph4UnfoldXErr(grreco[j][i],hrecopf[j][i]);
			else 
			xgr[j][i] = NULL;
		}
	}

	TCanvas *c = new TCanvas("c","c1",1000,800);
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	
	int mstyle[Npt][NRegion] = {{24},{8}};//{8,21,24};
	int mcolor[Npt][NRegion] = {{2},{kRed+2}};//{1,kBlue+1,2};
	float msize[NRegion] = {2};//{2,2,2};
	int lstyle[NRegion] = {3};//{1,7,3};			// Line Style
	int bstyle[NRegion] = {1001};//{3001,3144,1001};			// Box Style
	int bcolor[Npt][NRegion] = {{kMagenta-10},{kMagenta-8}};//{kGray, kBlue-10, kMagenta-10};

	for(int j = 0; j<Npt; j++) {
		for(int i = 0; i<NRegion; i++) {
			SetHistStyle(hrecopf[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],1,msize[i],1);
			SetHistStyle(httpf[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],lstyle[i],msize[i],7);	// last option is line width	// last option is line width
			SetHistStyle(httpf6def[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],lstyle[i],msize[i],3);
			SetHistStyle(httpf8[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],lstyle[i],msize[i],1);
			if(!xgr[j][i]) continue;
			xgr[j][i]->SetFillColor(bcolor[j][i]);
			xgr[j][i]->SetLineColor(bcolor[j][i]);
			xgr[j][i]->SetMarkerColor(bcolor[j][i]);
		}
	}

	hrecopf[0][0]->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/#it{c})");
	hrecopf[0][0]->GetYaxis()->SetTitleOffset(1.0);
	if(!strcmp(Variable,"PtAve")) {
		hrecopf[0][0]->GetYaxis()->SetTitle("#LTp_{T}^{ch}#GT");
	}
	else if(!strcmp(Variable,"PtSum")) {
		//hrecopf[0][0]->GetYaxis()->SetTitle("#LT#scale[0.5]{#sum}p_{T}^{ch}/#delta#eta#delta#phi#GT");
		hrecopf[0][0]->GetYaxis()->SetTitle("#LTd#scale[0.5]{#sum}p_{T}/d#etad#phi#GT");
	}
	else {		// default: Ntrk
		//hrecopf[0][0]->GetYaxis()->SetTitle("#LTN_{ch}/#delta#eta#delta#phi#GT");
		hrecopf[0][0]->GetYaxis()->SetTitle("#LTdN_{ch}/d#etad#phi#GT");
	}
	hrecopf[0][0]->SetMinimum(YaxisMin);
	hrecopf[0][0]->SetMaximum(YaxisMax);

	TH1D *htmp = (TH1D*) hrecopf[0][0]->Clone("htmp");
	htmp->Reset();
	htmp->GetXaxis()->SetRangeUser(0,45);
	htmp->GetXaxis()->SetNdivisions(205);
	htmp->Draw();

	for(int j = 0; j<Npt; j++) {
		for(int i = 0; i<NRegion; i++) {
		//for(int i = 2; i<NRegion; i++) {
			if(!grreco[j][i]) continue;

			double *ex;
			if(hrecopf[j][i]) {
				ex = new double[hrecopf[j][i]->GetNbinsX()];
				for( int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
					ex[k] = hrecopf[j][i]->GetBinWidth(k+1)/2.;
				}
			}

			double *eylow;
			eylow = new double[grreco[j][i]->GetN()];
			for(int k = 0; k<grreco[j][i]->GetN(); k++) {
				eylow[k] = - grreco[j][i]->GetErrorYlow(k);
			}
			grreco[j][i]->SetFillColor(kGray+i);

			// find the x-axis to plot points
			double *tmpx=grreco[j][i]->GetX();
			int skip=0;
			int total=0;
			for(int kg = 0; kg<grreco[j][i]->GetN(); kg++) {
				if(tmpx[kg]<XaxisMin) {
					skip++;
				}
				else {
					total++;
				}
				if(tmpx[kg]>XaxisMax) {
					total--;
					break;
				}
			}
			//cout<<"skip = "<<skip<<" total = "<<total<<endl;
			//cout<<"start at "<<tmpx[skip]<<" end at "<<tmpx[total+skip]<<endl;
			//graphSystBand(total,grreco[j][i]->GetX()+skip,grreco[j][i]->GetY()+skip,0,0, eylow+skip,  grreco[j][i]->GetEYhigh()+skip,bcolor[j][i], 0,0);
			graphSystBox(total,grreco[j][i]->GetX()+skip,grreco[j][i]->GetY()+skip,ex+skip,0, eylow+skip,  grreco[j][i]->GetEYhigh()+skip,bcolor[j][i], 0,0,1,0,1,bstyle[i]);	
		}
		if(grreco[j][0]) 
		{
			double *eylow;
			eylow = new double[grreco[j][0]->GetN()];
			for(int k = 0; k<grreco[j][0]->GetN(); k++) {
				eylow[k] = - grreco[j][0]->GetErrorYlow(k);
			}
			grreco[j][0]->SetFillColor(kGray+0);
			//graphSystBand(grreco[0]->GetN()-1,grreco[0]->GetX()+1,grreco[0]->GetY()+1,0,0, eylow+1,  grreco[0]->GetEYhigh()+1,bcolor[0], 0,0);
		}


		for(int i = 0; i<NRegion; i++) {
			if(!xgr[j][i]) continue;
			double *eylow;
			eylow = new double[xgr[j][i]->GetN()];
			for(int k = 0; k<xgr[j][i]->GetN(); k++) {
				eylow[k] = - xgr[j][i]->GetErrorYlow(k);
			}
			xgr[j][i]->SetMarkerStyle(8);
			//xgr[j][i]->Draw("psame");
			//graphSystBand(xgr[j][i]->GetN()-1,xgr[i]->GetX()+1,xgr[i]->GetY()+1,0,0, eylow+1,  xgr[i]->GetEYhigh()+1,bcolor[j][i], 0,0);	// make sure use the same color as greco...
		}
		if(xgr[j][0])
		{
			double *eylow;
			eylow = new double[xgr[j][0]->GetN()];
			for(int k = 0; k<xgr[j][0]->GetN(); k++) {
				eylow[k] = - xgr[j][0]->GetErrorYlow(k);
			}
			xgr[j][0]->SetMarkerStyle(8);
			//xgr[0]->Draw("psame");
			//graphSystBand(xgr[0]->GetN()-1,xgr[0]->GetX()+1,xgr[0]->GetY()+1,0,0, eylow+1,  xgr[0]->GetEYhigh()+1,bcolor[0], 0,0);	// make sure use the same color as greco...
		}

		for(int i=0;i<NRegion; i++) {		// Data unfolded
			hrecopf[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			hrecopf[j][i]->Draw("peX0same");
		}

		
		//hrecopf[j][2]->Draw("peX0same");//TEST

		// Normalized ProfileX to density
		if(!strcmp(Variable,"Ntrk")||!strcmp(Variable,"PtSum")) {
			float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
			for(int i = 0; i<NRegion; i++) {
				hmeaspf[j][i]->Scale(DeDpNorma);
				httpf[j][i]->Scale(DeDpNorma);
				htpf[j][i]->Scale(DeDpNorma);
				httpf6def[j][i]->Scale(DeDpNorma);
				httpf8[j][i]->Scale(DeDpNorma);
			}
		}


		for(int i=0;i<NRegion; i++) {		// PYTHIA
			httpf[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			httpf6def[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			httpf8[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
			httpf[j][i]->Draw("histcsame");
			if(opt_drawdefpythia6)  httpf6def[j][i]->Draw("histcsame");
			httpf8[j][i]->Draw("histcsame");
			//httpf[j][i]->Draw("HISTsame");
			//httpf8[j][i]->Draw("HISTsame");
		}
		
		//httpf[j][2]->Draw("histcsame");//TEST
	}

	TLatex *lat0 = new TLatex();
	lat0->SetNDC();
	lat0->SetTextFont(42);
	lat0->DrawLatex(0.14,0.84,"p_{T}>");

	TLegend *leg[Npt];
	for(int j = 0; j<Npt ; j++) {
		//leg[j] = new TLegend(0.15+j*0.25,0.6,0.38+j*0.25,0.88);	
		leg[j] = new TLegend(0.197+j*0.16,0.6,0.48+j*0.16,0.88);	
		leg[j]->SetFillColor(0);
		leg[j]->SetBorderSize(0);
		leg[j]->SetHeader(legtagpt[j]);
		leg[j]->SetTextSize(0.045);
		for(int i = 0; i<NRegion; i++) {
			if(j==0) {
				//leg[i]->AddEntry(hrecopf[j][i],Form("%s %s",legtag[i],legtagpt[j]),"p");
				leg[j]->AddEntry(hrecopf[j][i]," ","p");
			}
			else {
				leg[j]->AddEntry(hrecopf[j][i],legtag[i],"p");
			}
		}
		//leg->AddEntry(httpf[j][0],Form("Perugia 2012 %s",legtagpt[j]),"l");
		////leg->AddEntry(httpf8[j][0],Form("pythia 8.215 %s",legtagpt[j]),"l");
		//leg->AddEntry(httpf8[j][0],Form("Monash 2013 %s",legtagpt[j]),"l");
		//leg[j]->AddEntry(httpf[j][0],Form("Perugia 2012 (STAR)"),"l");
		if(j==0) {
			leg[j]->AddEntry(httpf[j][0],Form(" "),"l");
		}
		else{
			leg[j]->AddEntry(httpf[j][0],Form("PYTHIA 6 (STAR)"),"l");
		}
		//if(opt_drawdefpythia6)  leg[j]->AddEntry(httpf6def[j][0],"Perugia 2012 (default)","l");
		//if(opt_drawdefpythia6)  leg[j]->AddEntry(httpf6def[j][0],"Perugia 2012","l");
		if(j==0) {
			if(opt_drawdefpythia6)  leg[j]->AddEntry(httpf6def[j][0]," ","l");
		}
		else {
			if(opt_drawdefpythia6)  leg[j]->AddEntry(httpf6def[j][0],"PYTHIA 6","l");
		}
		if(j==0) {
			//leg[j]->AddEntry(httpf8[j][0],Form("Monash 2013"),"l");
			leg[j]->AddEntry(httpf8[j][0],Form(" "),"l");
		}
		else{
			leg[j]->AddEntry(httpf8[j][0],Form("PYTHIA 8"),"l");
		}
		leg[j]->Draw();
	}

	TLegend *leg2;
	leg2 = new TLegend(0.42,0.55,0.69,0.88);	
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->SetHeader("perugia 2012");
	for(int j = 0; j<Npt ; j++) {
		for(int i = 0; i<NRegion; i++) {
			leg2->AddEntry(httpf[j][i],Form("MC %s",legtag[i]),"l");
		}
	}
	//leg2->Draw();


	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextFont(42);
	//if(!strcmp(Variable,"PtAve")) {
	//	//lat->DrawLatex(0.65,0.19,"#splitline{p+p@200 GeV}{|#eta|<1}");
	//	//lat->DrawLatex(0.6,0.83,"p+p@200 GeV");
	//	lat->DrawLatex(0.6,0.88,"p+p@200 GeV");
	//	lat->DrawLatex(0.6,0.19,"#splitline{|#eta|<1}{R_{anti-k_{T}}=0.6, |#eta_{jet}|<0.4}");
	//}
	//else {
		//lat->DrawLatex(0.65,0.8,"#splitline{p+p@200 GeV}{|#eta|<1}");
		//lat->DrawLatex(0.65,0.8,"#splitline{p+p@200 GeV}{#splitline{|#eta|<1}{R = 0.6, |#eta_{jet}|<0.4}}");
		lat->DrawLatex(0.6,0.17,"p+p@200 GeV");
		lat->DrawLatex(0.72,0.81,"#splitline{|#eta|<1}{R_{anti-k_{T}}=0.6}");
		lat->DrawLatex(0.72,0.72,"|#eta_{jet}|<0.4");
	//}
	lat->SetTextColor(1);
	lat->SetTextFont(62);
	if(!strcmp(Variable,"PtSum")) {
		lat->DrawLatex(0.75,0.2,"STAR");
	}
	else {
		lat->DrawLatex(0.16,0.16,"STAR");
	}

        time_t rawtime;
        struct tm *timeinfo;
        char buffer4time[8];
        time(&rawtime);
        timeinfo=localtime(&rawtime);
        strftime(buffer4time,8,"%y%m%d",timeinfo);
        puts(buffer4time);


	if(savefig) {
		if(filetag.Contains("RcVzW",TString::kIgnoreCase)) {
			////c->SaveAs(Form("/Users/li/Research/Underlying/PaperDraft170405/Tran%s_MC02Vs05_Box.pdf",Variable));
			//c->SaveAs(Form("/Users/liyi/Research/Underlying/PaperDraft180920/Tran%s_MC02Vs05_Box.pdf",Variable));
			//c->SaveAs(Form("/Users/liyi/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/Tran%s_MCpT02Vs05.pdf",Variable));
			c->SaveAs(Form("fig/Tran%s_MCpT02Vs05_%s.pdf",Variable,buffer4time));
		}
		else {
			//c->SaveAs(Form("/Users/liyi/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/NoRcVzW_Tran%s_MCpT02Vs05.pdf",Variable));
			c->SaveAs(Form("fig/NoRcVzW_Tran%s_MCpT02Vs05_%s.pdf",Variable,buffer4time));
		}
	}
	if(saveroot) {
		if(filetag.Contains("RcVzW",TString::kIgnoreCase)) {
			c->SaveAs(Form("Tran%s_MCpT02Vs05_%s.root",Variable,buffer4time));
		}
		else {
			c->SaveAs(Form("NoRcVzW_Tran%s_MCpT02Vs05_%s.root",Variable,buffer4time));
		}
	}


	if(detplot) {
		TCanvas *c2 = new TCanvas("c2","c2",1000,800);
		c2->SetLeftMargin(0.12);
		c2->SetBottomMargin(0.12);
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		
		hmeaspf[0][0]->GetXaxis()->SetTitle("Leading jet p_{T} (GeV/#it{c})");
		if(!strcmp(Variable,"PtAve")) {
			hmeaspf[0][0]->GetYaxis()->SetTitle("#LTp_{T}#GT");
		}
		else {		// default: Ntrk
			//hmeaspf[0][0]->GetYaxis()->SetTitle("#LTN_{ch}/#delta#eta#delta#phi#GT");
			hmeaspf[0][0]->GetYaxis()->SetTitle("#LTdN_{ch}/d#etad#phi#GT");
		}
		hmeaspf[0][0]->SetMinimum(YaxisMin);
		hmeaspf[0][0]->SetMaximum(YaxisMax);


		TH1D *htmp2 = (TH1D*) hmeaspf[0][0]->Clone("htmp2");
		htmp2->Reset();
		htmp2->GetXaxis()->SetRangeUser(0,45);
		htmp2->GetXaxis()->SetNdivisions(205);
		htmp2->Draw();


		for(int j = 0; j<Npt; j++) {
			for(int i=0;i<NRegion; i++) {
				SetHistStyle(hmeaspf[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],1,msize[i],1);
				SetHistStyle(htpf[j][i],mcolor[j][i],mcolor[j][i],mstyle[j][i],lstyle[i],msize[i],2);
				hmeaspf[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
				hmeaspf[j][i]->Draw("peX0same");
				htpf[j][i]->GetXaxis()->SetRangeUser(XaxisMin, XaxisMax);
				htpf[j][i]->Draw("histcsame");
			}
			leg[j]->Draw("same");
		}


		TLatex *laxmeas = new TLatex(0.15,0.92,"Detector-level");
		laxmeas->SetTextFont(42);
		laxmeas->SetNDC();
		laxmeas->Draw();
	}


	// Print out data point
	for(int j=0;j<Npt; j++) {
		for(int i=0;i<NRegion; i++) {
			if(!grreco[j][i]) continue;
			cout<<legtag[i]<<": "<<endl;
			for(int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
				cout<<hrecopf[j][i]->GetBinCenter(k+1)<<" "<<hrecopf[j][i]->GetBinContent(k+1)<<" +/- "<<hrecopf[j][i]->GetBinError(k+1)<<" (stat);  "<<grreco[j][i]->GetY()[k]<<" +"<<grreco[j][i]->GetEYhigh()[k]<<" - "<<grreco[j][i]->GetEYlow()[k]<<" (sys)"<<"\t\t\t PYTHIA6 = "<<httpf[j][i]->GetBinContent(k+1)<<"\t\t\t PYTHIA8 = "<<httpf8[j][i]->GetBinContent(k+1)<<endl;
			}
		}
	}


	ofstream tout;
	TString Stout=Form("TranMc02Vs05_%s_JP_Charged_%s_embedMB_Baye%d.csv",Variable,filetag.Data(),DefaultTh);
	if(savecsv) {
		tout.open(Stout);
		if(tout.is_open()) {
			cout<<"Output to "<<Stout<<endl;
		}
		else {
			cout<<"Error open file "<<Stout<<" to write"<<endl;
			return;
		}

		tout<<"Variable pt,";
		for(int j = 0; j<hrecopf[0][0]->GetNbinsX(); j++) {
			tout<<hrecopf[0][0]->GetBinLowEdge(j+1)<<"-"<<hrecopf[0][0]->GetBinLowEdge(j+1)+hrecopf[0][0]->GetBinWidth(j+1)<<",";
		}
		tout<<endl;
		for(int j=0;j<Npt; j++) {
			for(int i=0;i<NRegion; i++) {
				if(!grreco[j][i]) continue;
				tout<<legtag[i]<<" "<<ptname[j]<<",";
				for(int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
					tout<<hrecopf[j][i]->GetBinContent(k+1)<<" +/- "<<hrecopf[j][i]->GetBinError(k+1)<<" (stat) + "<<grreco[j][i]->GetEYhigh()[k]<<" - "<<grreco[j][i]->GetEYlow()[k]<<" (sys)"<<",";
				}
				tout<<endl;
			}
		}
		tout<<"PYTHIA6 (STAR)"<<endl;
		for(int j=0;j<Npt; j++) {
			for(int i=0;i<NRegion; i++) {
				if(!grreco[j][i]) continue;
				tout<<legtag[i]<<" "<<ptname[j]<<",";
				for(int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
					tout<<httpf[j][i]->GetBinContent(k+1)<<",";
				}
				tout<<endl;
			}
		}
		tout<<"PYTHIA6 (default)"<<endl;
		for(int j=0;j<Npt; j++) {
			for(int i=0;i<NRegion; i++) {
				if(!grreco[j][i]) continue;
				tout<<legtag[i]<<" "<<ptname[j]<<",";
				for(int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
					tout<<httpf6def[j][i]->GetBinContent(k+1)<<",";
				}
				tout<<endl;
			}
		}
		tout<<"PYTHIA8"<<endl;
		for(int j=0;j<Npt; j++) {
			for(int i=0;i<NRegion; i++) {
				if(!grreco[j][i]) continue;
				tout<<legtag[i]<<" "<<ptname[j]<<",";
				for(int k = 0; k<hrecopf[j][i]->GetNbinsX(); k++) {
					tout<<httpf8[j][i]->GetBinContent(k+1)<<",";
				}
				tout<<endl;
			}
		}

		tout.close();
	}

}


