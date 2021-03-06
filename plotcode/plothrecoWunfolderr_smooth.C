#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include "/Users/li/commonmacro/ClassSysErr.C"
#include "compare2BinByBin.C"		// TH2D *GetBbB(TH2D* htt, TH2D* ht, TH2D* hm)
#include "CompareBayesTimes.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 
#include "smooth.C"


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





void graphSystBand(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey, 
		Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBand=5,
		Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1)
{
	int utilityPrint = 0;
	if (utilityPrint) printf("\n__________graphSysErr2band__________symb=%d sSize=%f color=%d colorBand=%d__________\n",symb,sSize,color,colorBand);
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
		float xband[1000], yband[1000];
		for (int i=0; i<n; ++i) {
			xband[i]=x[i];
			xband[2*n-1-i]=x[i];
			float y1=y[i]+syPlus[i];
			float y2;
			if (syMinus)  y2=y[i]+syMinus[i];
			else          y2=y[i]-syPlus[i];
			if (y1<y2) {float yy=y2; y2=y1; y1=yy;} 
			yband[i]=y1;
			yband[2*n-1-i]=y2;
			//printf("%d %d\n",i,2*n-1-i);
		}
		//for(int iii=0; iii<2*n;++iii)printf("%f %f\n",xband[iii],yband[iii]);
		TGraph *grafband = new TGraph(2*n,xband,yband);
		grafband->SetFillColor(colorBand);
		grafband->Draw("f");
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



//=======================================================================================

ClassSysErr *SumTwoSSysErr(TGraphErrors *gr, double *e1l, double *e1h, double *e2l, double *e2h) {

	const int MAXPOINT = 1024;
	if(MAXPOINT<gr->GetN()) { cout<<"increase MAXPOINT > "<<gr->GetN()<<" in function SysErrMultSingleS()"<<endl;}
	double x[MAXPOINT] = {0}, y[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
	double *xx, *yy; 
	xx = gr->GetX();
	yy = gr->GetY();
	for(int i = 0; i<gr->GetN(); i++) {
		x[i] = xx[i];
		y[i] = yy[i];
		eyl[i] = -sqrt( pow(e1l[i],2) + pow(e2l[i],2) );
		eyh[i] = sqrt( pow(e1h[i],2) + pow(e2h[i],2) );
		//cout<<i<<" e1h = "<<e1h[i]<<" e2h = "<<e2h[i]<<" eh = "<<eyh[i]<<endl;
	}

	ClassSysErr *syserr = new ClassSysErr(gr->GetN(), x, y, eyl, eyh);

	return syserr;

}



ClassSysErr *SumTwoSSysErr(TGraphAsymmErrors *gr, double *e1l, double *e1h, double *e2l, double *e2h) {

	const int MAXPOINT = 1024;
	if(MAXPOINT<gr->GetN()) { cout<<"increase MAXPOINT > "<<gr->GetN()<<" in function SysErrMultSingleS()"<<endl;}
	double x[MAXPOINT] = {0}, y[MAXPOINT] = {0}, eyl[MAXPOINT] = {0}, eyh[MAXPOINT] = {0};
	double *xx, *yy; 
	xx = gr->GetX();
	yy = gr->GetY();
	for(int i = 0; i<gr->GetN(); i++) {
		x[i] = xx[i];
		y[i] = yy[i];
		eyl[i] = -sqrt( pow(e1l[i],2) + pow(e2l[i],2) );
		eyh[i] = sqrt( pow(e1h[i],2) + pow(e2h[i],2) );
	}

	ClassSysErr *syserr = new ClassSysErr(gr->GetN(), x, y, eyl, eyh);

	return syserr;

}





//================================================== Main =================================================
void plothrecoWunfolderr_smooth(TString filename="Unfolding_TranPtAveJPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root", bool smooth = true, bool flagMC05=false) {
//bool flagMC05 = false;		// if we want to extract MC pt>0.5 from profile of MC pt>0.2 for TranPtAve. Use McPt02 for input
//2017.10.18 at one point in time, I tried to unfold 0.5 using 0.2, it gave large fluctation in TranPtAve. Then there are two options: unfold to 0.2 use 0.2, and later project to 0.5 on final TranPtAve; or unfold 0.5 use 0.5, which is probably the way to go, as for TranTotNtrk, we need to do that. So it is better to be consistent

	int DefaultTh = 5;
	const int BayesTimes = 5;
	int OtherTh[BayesTimes] = {2,3,4,6,7};
	if(filename.Contains("LeadAreaNtrk")&&OtherTh[2]==4) OtherTh[2] = 5;		// There is a jump when iteration==4 for LeadAreaNtrk

	const int TpcTimes = 4;
	const char *TpcName[TpcTimes] = {"_TpcErrPlus005","_TpcErrMinus005","_TpcErrPlusAbs005","_TpcErrMinusAbs005"};

	const int BemcTimes = 2;
	const char *BemcName[BemcTimes] = {"_BemcErrPlus004","_BemcErrMinus004"};

	const int Nfile = (BayesTimes+2)*2+1+TpcTimes+BemcTimes;						// (Bayes iterations + Bin-by-Bin) * (NFweight vs Not weight)	+ YScale  + (TPC 5% tracking plus+minus) + (BEMC 4% tower gain)
	TFile *f[Nfile]={NULL};	
	TH2D *hreco[Nfile]={NULL};


	//bool flagMC05 = false;		// if we want to extract MC pt>0.5 from profile of MC pt>0.2 for TranPtAve

	//Make sure input default filename is NFweighted with MB embedding
	//if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_HalfMcPt02_embedJP0","_NFweight_HalfMcPt02_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_12JetBinv2_McPt02_embedJP0","_NFweight_12JetBinv2_McPt02_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_12JetBinv2_McPtRC02MC05_embedJP0","_NFweight_12JetBinv2_McPtRC02MC05_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_12JetBinv2_McPtRC05MC05_embedJP0","_NFweight_12JetBinv2_McPtRC05MC05_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_12JetBinv2_McPt02_embedJP0","_NFweight_BT170928_12JetBinv2_McPt02_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_12JetBinv2_McPtRC02MC05_embedJP0","_NFweight_BT170928_12JetBinv2_McPtRC02MC05_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_12JetBinv2_McPtRC05MC05_embedJP0","_NFweight_BT170928_12JetBinv2_McPtRC05MC05_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_RcVzW_12JetBinv2_McPt02_embedJP0","_NFweight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_RcVzW_12JetBinv2_McPtRC02MC05_embedJP0","_NFweight_BT170928_RcVzW_12JetBinv2_McPtRC02MC05_embedMB");
	if(!filename.Contains("_NFweight",TString::kIgnoreCase)) filename.ReplaceAll("_BT170928_RcVzW_12JetBinv2_McPtRC05MC05_embedJP0","_NFweight_BT170928_RcVzW_12JetBinv2_McPtRC05MC05_embedMB");
	//Make sure input default filename is consistent with DefaultTh 
	if(! (filename.Contains(Form("_Baye%d",DefaultTh))||(DefaultTh==4)) ) {
		Ssiz_t toinsert = filename.Index(".root");
		Ssiz_t toremove = filename.Index("_Baye");
		if(toremove!=-1) filename.Remove(toremove,(toinsert-toremove));		// WARNINING!!  Assuming _Baye%d.root
		if(DefaultTh!=4) {
			toinsert = filename.Index(".root");
			filename.Insert(toinsert,Form("_Baye%d",DefaultTh));
		}
	}
	//Default one
	cout<<filename<<endl;
	f[0] = new TFile(filename);
	if(!f[0]->IsOpen()) {  
		if(filename.Contains("TranPtAve")&&filename.Contains("McPtRC02MC05")) {
			cout<<"WARNING!!!: if you want to read McPtRC02MC05 case for TranPtAve, I am going to change to read input files as MC02 and profile it as MC05."<<endl;
			//filename = "Unfolding_TranPtAveJPCharged_NFWeight_12JetBinv2_McPt02_embedMB_Baye5.root";
			filename = "Unfolding_TranPtAveJPCharged_NFWeight_BT170928_12JetBinv2_McPt02_embedMB_Baye5.root";
			cout<<"WARNING!!!: I am reading "<<filename<<"."<<endl;
			f[0] = new TFile(filename);
			flagMC05 = true;		// so we are doing profile to MC05 later
			if(!f[0]->IsOpen()) {
				cout<<"ERROR: failed again.. cannot open file "<<filename<<"!!"<<endl;
				return;
			}
		}
		else {
			cout<<"Error: Cannot Open File "<<filename<<"!!"<<endl; 
			return; 
		}
	}
	hreco[0] = (TH2D*)f[0]->Get("hreco");
	hreco[0]->SetName(Form("hreco_default%d", DefaultTh));
	TH2D *hmeas = (TH2D*)f[0]->Get("hmeas");
	TH2D *htt = (TH2D*)f[0]->Get("htraintrue");
	TH2D *ht = (TH2D*)f[0]->Get("htrain");

	//Different bayes iteration
	for(int i = 1; i<BayesTimes+1; i++) {
		TString ifilename = filename;
		Ssiz_t toinsert = ifilename.Index(".root");
		Ssiz_t toremove = ifilename.Index("_Baye");
		if(toremove!=-1) {		// find _Baye in filename
			ifilename.Remove(toremove,(toinsert-toremove));		// WARNINING!!  Assuming _Baye%d.root
			//cout<<" -- "<<ifilename<<endl;
			toinsert = ifilename.Index(".root");
		}
		if(OtherTh[i-1]!=4) {		// originally I assume 4 as the default without _Baye4 in filename
			ifilename.Insert(toinsert,Form("_Baye%d",OtherTh[i-1]));
		}
		cout<<ifilename<<endl;
		f[i] = new TFile(ifilename);
		if(!f[i]->IsOpen()) {  cout<<"Error: Cannot Open File "<<ifilename<<"!!"<<endl; }
		else {
			hreco[i] = (TH2D*)f[i]->Get("hreco");
			hreco[i]->SetName(Form("hreco%d",OtherTh[i-1]));
		}
	}

	//Bin-by-Bin
	hreco[BayesTimes+1] = (TH2D*)GetBbB(htt,ht,hmeas);
	hreco[BayesTimes+1]->SetName("hreco_BbB");


	//--- No NFweighted one ---
	// default bayes iteration WITHOUT NFweight
	TString nfilename = filename;
	TH2D *hmeasNFF;
	TH2D *httNFF;
	TH2D *htNFF;
	//if((nfilename.Contains("_NFweight_12JetBinv2_McPt02_embedMB")||nfilename.Contains("_NFWeight_12JetBinv2_McPt02_embedMB")))  {
	//	nfilename.ReplaceAll("_NFweight_12JetBinv2_McPt02_embedMB","_12JetBinv2_McPt02_embedJP0");
	//	nfilename.ReplaceAll("_NFWeight_12JetBinv2_McPt02_embedMB","_12JetBinv2_McPt02_embedJP0");
	if((nfilename.Contains("embedMB")))  {
		nfilename.ReplaceAll("embedMB","embedJP0");
		nfilename.ReplaceAll("_NFWeight","");
		cout<<"replace to --> "<<endl<< nfilename<<endl;
		f[BayesTimes+2] = new TFile(nfilename);
		if(!f[BayesTimes+2]->IsOpen()) {  cout<<"Error: Cannot Open File "<<nfilename<<"!!"<<endl; }
		else {
			hreco[BayesTimes+2] = (TH2D*)f[BayesTimes+2]->Get("hreco");
			hreco[BayesTimes+2]->SetName(Form("hreco_NoNFw_default"));
			hmeasNFF = (TH2D*)f[BayesTimes+2]->Get("hmeas");
			httNFF = (TH2D*)f[BayesTimes+2]->Get("htraintrue");
			htNFF = (TH2D*)f[BayesTimes+2]->Get("htrain");
		}

		// Different bayes iteration WITHOUT NFweight
		for(int i = BayesTimes+3; i<(BayesTimes+2)*2-1; i++) {
			TString ifilename = nfilename;
			Ssiz_t toinsert = ifilename.Index(".root");
			Ssiz_t toremove = ifilename.Index("_Baye");
			if(toremove!=-1) {		// find _Baye in filename
				ifilename.Remove(toremove,(toinsert-toremove));		// WARNINING!!  Assuming _Baye%d.root
				toinsert = ifilename.Index(".root");
			}
			if(OtherTh[i-(BayesTimes+3)]!=4) {		// originally I assume 4 as the default without _Baye4 in filename
				ifilename.Insert(toinsert,Form("_Baye%d",OtherTh[i-(BayesTimes+3)]));
			}
			cout<<ifilename<<endl;
			f[i] = new TFile(ifilename);
			if(!f[i]->IsOpen()) {  cout<<"Error: Cannot Open File "<<ifilename<<"!!"<<endl; }
			else {
				hreco[i] = (TH2D*)f[i]->Get("hreco");
				hreco[i]->SetName(Form("hrecoNNF%d",OtherTh[i-(BayesTimes+3)]));
			}
		}

		// Bin-by-Bin WITHOUT NFweight
		if(f[BayesTimes+2]->IsOpen()) {
			hreco[(BayesTimes+2)*2-1] = (TH2D*)GetBbB(httNFF,htNFF,hmeasNFF);
			hreco[(BayesTimes+2)*2-1]->SetName("hrecoNNF_BbB");
		}
	}



	//--- YScaled one ---
	TString ysfilename = filename;
	ysfilename.ReplaceAll("_NFweight","_YScale");
	ysfilename.ReplaceAll("_NFWeight","_YScale");
	cout<<"Read Yscaled file: "<<ysfilename<<endl;
	f[(BayesTimes+2)*2] = new TFile(ysfilename);
	if(!f[(BayesTimes+2)*2]->IsOpen()) {  cout<<"Error: Cannot Open File "<<ysfilename<<"!!"<<endl; }
	else {
		hreco[(BayesTimes+2)*2] = (TH2D*)f[(BayesTimes+2)*2]->Get("hreco");
		hreco[(BayesTimes+2)*2]->SetName(Form("hreco_Yscale"));
	}



	//--- TPC tracking 5%: relative + absolute ---
	cout<<"tpc ones: "<<endl;
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		TString ifilename = filename;
		Ssiz_t toinsert = ifilename.Index(".root");
		ifilename.Insert(toinsert,TpcName[i-((BayesTimes+2)*2+1)]);
		cout<<ifilename<<endl;
		f[i] = new TFile(ifilename);
		if(!f[i]->IsOpen()) {  cout<<"Error: Cannot Open File "<<ifilename<<"!!"<<endl; }
		else {
			hreco[i] = (TH2D*)f[i]->Get("hreco");
			hreco[i]->SetName(Form("hreco%s",TpcName[i-((BayesTimes+2)*2+1)]));
		}
	}


	//--- BEMC tracking 5%: relative + absolute ---
	cout<<"bemc ones: "<<endl;
	for(int i = (BayesTimes+2)*2+1+TpcTimes; i<(BayesTimes+2)*2+1+TpcTimes+BemcTimes; i++) {
		TString ifilename = filename;
		Ssiz_t toinsert = ifilename.Index(".root");
		ifilename.Insert(toinsert,BemcName[i-((BayesTimes+2)*2+1+TpcTimes)]);
		cout<<ifilename<<endl;
		f[i] = new TFile(ifilename);
		if(!f[i]->IsOpen()) {  cout<<"Error: Cannot Open File "<<ifilename<<"!!"<<endl; }
		else {
			hreco[i] = (TH2D*)f[i]->Get("hreco");
			hreco[i]->SetName(Form("hreco%s",BemcName[i-((BayesTimes+2)*2+1+TpcTimes)]));
		}
	}



	// Take ProfileX
	cout<<"Take ProfileX"<<endl;
	TProfile *hp[Nfile]={NULL};
	for(int i = 0; i<Nfile; i++) {
		if(hreco[i]) {
			if(flagMC05) {
				if(!filename.Contains("PtAve")) {cout<<"ERR!!! when flagMC05==1, we also need to have filename for PtAve, other filename variable type NOT acceptance!!!"<<endl; return;}
				TProfile *hpy = hreco[i]->ProfileY();
				if(smooth) {
					hp[i] = (TProfile*)SmoothProfileX(hreco[i],Form("%s_pfx",hreco[i]->GetName()),hpy->FindBin(0.5));
				}
				else {
					hp[i] = (TProfile*)hreco[i]->ProfileX(Form("%s_pfx",hreco[i]->GetName()),hpy->FindBin(0.5));	// profile starting from 0.5 instead of 0.2 GeV/c
				}
			}
			else {
				if(smooth) {
					hp[i] = (TProfile*)SmoothProfileX(hreco[i],Form("%s_pfx",hreco[i]->GetName()));
				}
				else {
					hp[i] = (TProfile*)hreco[i]->ProfileX(Form("%s_pfx",hreco[i]->GetName()));
				}
			}
			cout<<hp[i]->GetName()<<endl;
		}
	}

	// Normalized ProfileX to density
	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	if(!(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)||filename.Contains("PtAve",TString::kIgnoreCase))) {
		for(int i = 0; i<Nfile; i++) {
			if(hp[i]) {
				hp[i]->Scale(DeDpNorma);
			}
		}
	}


	// Convert ProfileX to TGraphErrors for ClassSysErr
	//const int Nbins = 13;
	int Nbins = hreco[0]->GetNbinsX();
	TGraphErrors *gr[Nfile]={NULL};
	for(int i = 0; i<Nfile; i++) {
		if(hreco[i]) {
			gr[i] = new TGraphErrors();
			gr[i]->SetName(Form("gr_%s",hreco[i]->GetName()));
			for(int j = 0; j<hp[i]->GetNbinsX(); j++) {
				gr[i]->SetPoint(j, hp[i]->GetBinCenter(j+1), hp[i]->GetBinContent(j+1));
				gr[i]->SetPointError(j, 0, hp[i]->GetBinError(j+1));
			}
			gr[i]->SetLineColor(kViolet-22+i*2);
			gr[i]->SetMarkerColor(kViolet-22+i*2);
			gr[i]->SetMarkerStyle(24);
			gr[i]->SetMarkerSize(2);
		}
	}
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		if(gr[i]) {
			gr[i]->SetMarkerStyle(23);
		}
	}
	for(int i = (BayesTimes+2)*2+1+TpcTimes; i<(BayesTimes+2)*2+1+TpcTimes+BemcTimes; i++) {
		if(gr[i]) {
			gr[i]->SetMarkerStyle(25);
			gr[i]->SetMarkerColor(kGreen+(i-((BayesTimes+2)*2+1+TpcTimes))*2);
		}
	}

	TGraphErrors *grNonNULL_unfold[Nfile] = {NULL};
	TGraphErrors *grNonNULL_tpc[Nfile] = {NULL};
	TGraphErrors *grNonNULL_bemc[Nfile] = {NULL};
	int NgrNonNull_unfold=0;
	int NgrNonNull_tpc=0;
	int NgrNonNull_bemc=0;
	for(int i = 1 ; i<(BayesTimes+2)*2+1; i++) {
		if(gr[i]&&gr[i]->GetN()!=0) {
			grNonNULL_unfold[NgrNonNull_unfold] = gr[i];
			NgrNonNull_unfold++;
		}
	}
	for(int i =(BayesTimes+2)*2+1 ; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		if(gr[i]&&gr[i]->GetN()!=0) {
			grNonNULL_tpc[NgrNonNull_tpc] = gr[i];
			NgrNonNull_tpc++;
		}
	}
	for(int i =(BayesTimes+2)*2+1+TpcTimes ; i<(BayesTimes+2)*2+1+TpcTimes+BemcTimes; i++) {
		if(gr[i]&&gr[i]->GetN()!=0) {
			grNonNULL_bemc[NgrNonNull_bemc] = gr[i];
			NgrNonNull_bemc++;
		}
	}
	TString opt_unfold = "s";			// will take the max of all as the final unfolding sys err
	//ClassSysErr *syserr_unfold = new ClassSysErr((BayesTimes+2)*2,gr[0], gr+1,opt_unfold);			// Bayes, Bin-by-Bin, NFweight or not, YScale
	ClassSysErr *syserr_unfold=NULL;
	if(NgrNonNull_unfold>0) {
		syserr_unfold = new ClassSysErr(NgrNonNull_unfold,gr[0], grNonNULL_unfold,opt_unfold);			// Bayes, Bin-by-Bin, NFweight or not, YScale
	}
	TString opt_tracking = "s";
	//ClassSysErr *syserr_tpc = new ClassSysErr(TpcTimes,gr[0], gr+(BayesTimes+2)*2+1,opt_tracking);
	ClassSysErr *syserr_tpc=NULL;
	if(NgrNonNull_tpc>0) {
		syserr_tpc = new ClassSysErr(TpcTimes,gr[0], grNonNULL_tpc, opt_tracking);
	}
	TString opt_towergain = "s";
	ClassSysErr *syserr_bemc=NULL;
	if(NgrNonNull_bemc>0) {
		syserr_bemc = new ClassSysErr(BemcTimes,gr[0], grNonNULL_bemc, opt_towergain);
	}

	ClassSysErr *syserr=NULL;
	if(syserr_unfold&&syserr_tpc) {
		syserr = SumTwoSSysErr(gr[0], syserr_unfold->GetEYlow(), syserr_unfold->GetEYhigh(), syserr_tpc->GetEYlow(), syserr_tpc->GetEYhigh());
		if(syserr_bemc) {
			ClassSysErr *syserr2 = syserr;
			syserr = SumTwoSSysErr(gr[0], syserr2->GetEYlow(), syserr2->GetEYhigh(), syserr_bemc->GetEYlow(), syserr_bemc->GetEYhigh());
		}
	}
	else if(syserr_unfold) {
		syserr = syserr_unfold;
	}
	else if(syserr_tpc) {
		syserr = syserr_tpc;
	}
	else if(syserr_bemc) {
		syserr = syserr_bemc;
	}

	TString sregion = "Transverse";
	if(filename.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) svariable = " #LTp_{T}^{ch}#GT";
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) svariable = " #LT#sump_{T}^{ch}/#delta#eta#delta#phi#GT";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) svariable = " Constituents Multiplicity";
	hp[0]->GetYaxis()->SetTitle(sregion+svariable);
	hp[0]->GetXaxis()->SetRangeUser(0,55);
	hp[0]->SetMinimum(0);
	if(sregion.EqualTo("Transverse")) {
		hp[0]->SetMaximum(0.9);//Tran Ntrk density
		hp[0]->SetMinimum(0.3);
	}
	if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hp[0]->SetMaximum(2.3);//Lead Ntrk density
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hp[0]->SetMaximum(3);//Tran PtSum
		if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hp[0]->SetMaximum(30);//Lead PtSum
	}
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hp[0]->SetMaximum(0.68);//Tran PtAve
		if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hp[0]->SetMaximum(4);//Lead PtAve
	}
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) hp[0]->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)
	if((filename.Contains("MC05")||flagMC05)&&filename.Contains("TranPtAve")) hp[0]->SetMaximum(1.8);	// Tran PtAve MC05



	TCanvas *c = new TCanvas();
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	c->SetFrameLineWidth(2);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);


	SetHistStyle(hp[0],1,1,8,1,2,1);
	hp[0]->GetXaxis()->SetRangeUser(0,55);
	hp[0]->Draw("epX0");

//#define SUMMARYPLOT

// draw sysband of sum and each individual data set
//#ifdef SUMMARYPLOT
	graphSystBand(Nbins,syserr->GetX(),syserr->GetY(),0,0,syserr->GetEYlow(),  syserr->GetEYhigh(),kGray);
	hp[0]->Draw("psame");
	for(int i = 1; i<Nfile; i++) {
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i>(BayesTimes+2)?i-(BayesTimes+3):i-1]==DefaultTh)&&gr[i])
			gr[i]->Draw("psame");
	}


	TLegend *leg;
	if(filename.Contains("LeadJetNtrk")) {
		leg = new TLegend(0.17,0.6,0.55,0.88);		//LeadJetNtrk
	}
	else if(filename.Contains("Lead",TString::kIgnoreCase)||filename.Contains("Away",TString::kIgnoreCase)||filename.Contains("Sub",TString::kIgnoreCase)) {
		leg = new TLegend(0.5,0.15,0.95,0.65);		//Lead
	}
	else {
		leg = new TLegend(0.5,0.6,0.88,0.88);	//Tran
	}
	leg->SetFillColor(0);
	leg->AddEntry(hp[0],Form("Default w/ NFweight iter=%d",DefaultTh),"p");
	for(int i = 1; i<BayesTimes+1; i++) {
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i-1]==DefaultTh) && gr[i])
			leg->AddEntry(gr[i],Form("w/ NFweight iter=%d",OtherTh[i-1]),"p");
	}
	leg->AddEntry(gr[BayesTimes+1],"w/ NFweight Bin-by-Bin","p");
	leg->AddEntry(gr[BayesTimes+2],Form("w/o NFweight iter=%d",DefaultTh),"p");
	for(int i = BayesTimes+3; i<(BayesTimes+2)*2-1; i++) {
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i-(BayesTimes+3)]==DefaultTh) && gr[i])
			leg->AddEntry(gr[i],Form("w/o NFweight iter=%d",OtherTh[i-(BayesTimes+3)]),"p");
	}
	if(gr[(BayesTimes+2)*2-1])
		leg->AddEntry(gr[(BayesTimes+2)*2-1],"w/o NFweight Bin-by-Bin","p");
	if(gr[(BayesTimes+2)*2])
		leg->AddEntry(gr[(BayesTimes+2)*2],"Scaled meas to match MB","p");
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		if(gr[i])
			leg->AddEntry(gr[i],TpcName[i-((BayesTimes+2)*2+1)],"p");
	}
	for(int i = (BayesTimes+2)*2+1+TpcTimes; i<(BayesTimes+2)*2+1+TpcTimes+BemcTimes; i++) {
		if(gr[i])
			leg->AddEntry(gr[i],BemcName[i-((BayesTimes+2)*2+1+TpcTimes)],"p");
	}
	leg->Draw();

//#endif

//#ifndef SUMMARYPLOT

	TCanvas *c2 = new TCanvas();
	c2->SetLeftMargin(0.12);
	c2->SetBottomMargin(0.12);
	c2->SetFrameLineWidth(2);

	hp[0]->Draw("epX0");

	double *ex;
	if(hp[0]) {
		ex = new double[hp[0]->GetNbinsX()];
		for( int j = 0; j<hp[0]->GetNbinsX(); j++) {
			ex[j] = hp[0]->GetBinWidth(j+1)/2.;
		}
	}
	int bandcolor[3] = {kRed, kBlue-10, kSpring};
	int bandstyle[3] = {3001, 3344, 3002};
	graphSystBox(Nbins,syserr_unfold->GetX(),syserr_unfold->GetY(),ex,0, syserr_unfold->GetEYlow(),  syserr_unfold->GetEYhigh(), bandcolor[0], 0,0,1,0,1,bandstyle[0]);		
	graphSystBox(Nbins,syserr_tpc->GetX(),syserr_tpc->GetY(),ex,0, syserr_tpc->GetEYlow(),  syserr_tpc->GetEYhigh(), bandcolor[1], 0,0,1,0,1,bandstyle[1]);	
	graphSystBox(Nbins,syserr_bemc->GetX(),syserr_bemc->GetY(),ex,0, syserr_bemc->GetEYlow(),  syserr_bemc->GetEYhigh(), bandcolor[2], 0,0,1,0,1,bandstyle[2]);	
	hp[0]->Draw("peX0same");

	TH1D *htmp_unfold = new TH1D("htmp_unfold","dummy",100,0,100);
	TH1D *htmp_tpc = new TH1D("htmp_tpc","dummy",100,0,10);
	TH1D *htmp_bemc = new TH1D("htmp_bemc","dummy",100,0,10);
	htmp_unfold->SetFillColor(bandcolor[0]);
	htmp_unfold->SetFillStyle(bandstyle[0]);
	htmp_tpc->SetFillColor(bandcolor[1]);
	htmp_tpc->SetFillStyle(bandstyle[1]);
	htmp_bemc->SetFillColor(bandcolor[2]);
	htmp_bemc->SetFillStyle(bandstyle[2]);

	TLegend *leg2 = new TLegend(0.55,0.2,0.85,0.5);
	leg2->AddEntry(hp[0],Form("Default w/ NFweight iter=%d",DefaultTh),"p");
	leg2->AddEntry(htmp_unfold,"Unfolding Sys. Err.","f");
	leg2->AddEntry(htmp_tpc,"TPC Tracking Sys. Err.","f");
	leg2->AddEntry(htmp_bemc,"BEMC Tower Sys. Err.","f");
	leg2->Draw("same");


	// print out sys. err. seperately for unfold, tpc, bemc
	for(int i = 0; i<Nbins; i++) {
		double *x_syserr = syserr->GetX();
		double *y_syserr = syserr->GetY();
		double *eyl_syserr = syserr->GetEYlow();
		double *eyh_syserr = syserr->GetEYhigh();
		double *eyl_syserr_unfold = syserr_unfold->GetEYlow();
		double *eyh_syserr_unfold = syserr_unfold->GetEYhigh();
		double *eyl_syserr_tpc = syserr_tpc->GetEYlow();
		double *eyh_syserr_tpc = syserr_tpc->GetEYhigh();
		double *eyl_syserr_bemc = syserr_bemc->GetEYlow();
		double *eyh_syserr_bemc = syserr_bemc->GetEYhigh();
		cout<<"At "<<x_syserr[i]<<": "<<y_syserr[i]<<" +"<<eyh_syserr[i]<<""<<eyl_syserr[i]<<"(+"<<eyh_syserr[i]/y_syserr[i]*100<<"%"<<eyl_syserr[i]/y_syserr[i]*100<<"%)"<<endl;
		cout<< std::setprecision(2)<<" -> +"<<eyh_syserr_unfold[i]<<""<<eyl_syserr_unfold[i]<<"(unfold +"<<eyh_syserr_unfold[i]/y_syserr[i]*100<<"%"<<eyl_syserr_unfold[i]/y_syserr[i]*100<<"%) +"<<eyh_syserr_tpc[i]<<""<<eyl_syserr_tpc[i]<<"(tpc +"<<eyh_syserr_tpc[i]/y_syserr[i]*100<<"%"<<eyl_syserr_tpc[i]/y_syserr[i]*100<<"%) +"<<eyh_syserr_bemc[i]<<""<<eyl_syserr_bemc[i]<<"(bemc +"<<eyh_syserr_bemc[i]/y_syserr[i]*100<<"%"<<eyl_syserr_bemc[i]/y_syserr[i]*100<<"%)"<<endl;
	}

//#endif


	if(1) {
		TString outname = "SysErr4"+filename;
		if(flagMC05) outname.ReplaceAll("McPt02","McPtRC02MC05");
		if(smooth) outname.ReplaceAll(".root","_smooth.root");
		TFile *fout = new TFile(outname,"RECREATE");
		c->Write();
		hp[0]->Write();
		double *eylow;
		eylow = new double[Nbins];
		for(int i = 0; i<Nbins; i++) {
			eylow[i] = fabs(syserr->GetEYlow()[i]);	//syserr->GetEYlow() return negative values
		}
		TGraphAsymmErrors *gsyserr = new TGraphAsymmErrors(Nbins, syserr->GetX(),syserr->GetY(),0,0,eylow,  syserr->GetEYhigh());
		gsyserr->SetName("gSysErr"+sregion);
		gsyserr->Write();
		if(!flagMC05) {
			hmeas->ProfileX()->Write();
			htt->ProfileX()->Write();
			ht->ProfileX()->Write();
		}
		else {
			TProfile *hpy = hmeas->ProfileY();
			hmeas->ProfileX(Form("%s_pfx",hmeas->GetName()),hpy->FindBin(0.5))->Write();
			htt->ProfileX(Form("%s_pfx",htt->GetName()),hpy->FindBin(0.5))->Write();
			ht->ProfileX(Form("%s_pfx",ht->GetName()),hpy->FindBin(0.5))->Write();
		}

		fout->Close();
	}

}
