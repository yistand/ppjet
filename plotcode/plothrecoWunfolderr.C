#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include "/Users/li/commonmacro/ClassSysErr.C"
#include "compare2BinByBin.C"		// TH2D *GetBbB(TH2D* htt, TH2D* ht, TH2D* hm)
#include "CompareBayesTimes.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 

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

void plothrecoWunfolderr(TString filename="Unfolding_TranTotNtrkJPCharged_NFweight_McPt02_embedMB_Baye5.root") {

	int DefaultTh = 5;
	const int BayesTimes = 5;
	int OtherTh[BayesTimes] = {2,3,4,6,7};
	if(filename.Contains("LeadAreaNtrk")&&OtherTh[2]==4) OtherTh[2] = 5;		// There is a jump when iteration==4 for LeadAreaNtrk

	const int TpcTimes = 2;
	const char *TpcName[TpcTimes] = {"_TpcErrPlus005","_TpcErrMinus005"};

	const int Nfile = (BayesTimes+2)*2+1+TpcTimes;						// (Bayes iterations + Bin-by-Bin) * (NFweight vs Not weight)	+ YScale  + (TPC 5% tracking plus+minus)
	TFile *f[Nfile];	
	TH2D *hreco[Nfile];


	//Make sure input default filename is NFweighted with MB embedding
	if(!filename.Contains("_NFweight")) filename.ReplaceAll("_McPt02_embedJP0","_NFweight_McPt02_embedMB");
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
		hreco[i] = (TH2D*)f[i]->Get("hreco");
		hreco[i]->SetName(Form("hreco%d",OtherTh[i-1]));
	}

	//Bin-by-Bin
	hreco[BayesTimes+1] = (TH2D*)GetBbB(htt,ht,hmeas);
	hreco[BayesTimes+1]->SetName("hreco_BbB");


	//--- No NFweighted one ---
	// default bayes iteration WITHOUT NFweight
	TString nfilename = filename;
	nfilename.ReplaceAll("_NFweight_McPt02_embedMB","_McPt02_embedJP0");
	cout<<"replace to --> "<<endl<< nfilename<<endl;
	f[BayesTimes+2] = new TFile(nfilename);
	hreco[BayesTimes+2] = (TH2D*)f[BayesTimes+2]->Get("hreco");
	hreco[BayesTimes+2]->SetName(Form("hreco_NoNFw_default"));
	TH2D *hmeasNFF = (TH2D*)f[BayesTimes+2]->Get("hmeas");
	TH2D *httNFF = (TH2D*)f[BayesTimes+2]->Get("htraintrue");
	TH2D *htNFF = (TH2D*)f[BayesTimes+2]->Get("htrain");

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
		hreco[i] = (TH2D*)f[i]->Get("hreco");
		hreco[i]->SetName(Form("hrecoNNF%d",OtherTh[i-(BayesTimes+3)]));
	}

	// Bin-by-Bin WITHOUT NFweight
	hreco[(BayesTimes+2)*2-1] = (TH2D*)GetBbB(httNFF,htNFF,hmeasNFF);
	hreco[(BayesTimes+2)*2-1]->SetName("hrecoNNF_BbB");



	//--- YScaled one ---
	TString ysfilename = filename;
	ysfilename.ReplaceAll("_NFweight","_YScale");
	cout<<"Read Yscaled file: "<<ysfilename<<endl;
	f[(BayesTimes+2)*2] = new TFile(ysfilename);
	hreco[(BayesTimes+2)*2] = (TH2D*)f[(BayesTimes+2)*2]->Get("hreco");
	hreco[(BayesTimes+2)*2]->SetName(Form("hreco_Yscale"));



	//--- TPC tracking 5% ---
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		TString ifilename = filename;
		Ssiz_t toinsert = ifilename.Index(".root");
		ifilename.Insert(toinsert,TpcName[i-((BayesTimes+2)*2+1)]);
		f[i] = new TFile(ifilename);
		hreco[i] = (TH2D*)f[i]->Get("hreco");
		hreco[i]->SetName(Form("hreco%s",TpcName[i-((BayesTimes+2)*2+1)]));
	}



	// Take ProfileX
	cout<<"Take ProfileX"<<endl;
	TProfile *hp[Nfile];
	for(int i = 0; i<Nfile; i++) {
		hp[i] = (TProfile*)hreco[i]->ProfileX(Form("%s_pfx",hreco[i]->GetName()));
		cout<<hp[i]->GetName()<<endl;
	}

	// Normalized ProfileX to density
	float DeDpNorma = 1./(2.*2.*TMath::Pi()/3.);	// TranTot: 2 area; eta 2; phi pi/3
	if(!filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) {
		for(int i = 0; i<Nfile; i++) {
			hp[i]->Scale(DeDpNorma);
		}
	}


	// Convert ProfileX to TGraphErrors for ClassSysErr
	const int Nbins = 13;
	TGraphErrors *gr[Nfile];
	for(int i = 0; i<Nfile; i++) {
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
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		gr[i]->SetMarkerStyle(23);
	}

	TString opt_unfold = "s";			// will take the max of all as the final unfolding sys err
	ClassSysErr *syserr_unfold = new ClassSysErr((BayesTimes+2)*2,gr[0], gr+1,opt_unfold);			// Bayes, Bin-by-Bin, NFweight or not, YScale
	TString opt_tracking = "s";
	ClassSysErr *syserr_tpc = new ClassSysErr(TpcTimes,gr[0], gr+(BayesTimes+2)*2+1,opt_tracking);
	// estimate 7% from 5% TPC tracking uncertainty 
	int N5 = syserr_tpc->GetN();
	double *x5, *y5, *eyl5, *eyh5;
	double *x7, *y7, *eyl7, *eyh7;
	x5 = syserr_tpc->GetX();
	y5 = syserr_tpc->GetY();
	eyl5 = syserr_tpc->GetEYlow();
	eyh5 = syserr_tpc->GetEYhigh();
	x7 = new double[N5];
	y7 = new double[N5];
	eyl7 = new double[N5];
	eyh7 = new double[N5];
	for(int i = 0; i<N5; i++) {
		x7[i] = x5[i];
		y7[i] = y5[i];
		eyl7[i] = eyl5[i]*1.4;
		eyh7[i] = eyh5[i]*1.4;
	}
	ClassSysErr *syserr_tpc7 = new ClassSysErr(N5, x7, y7, eyl7, eyh7);

	ClassSysErr *syserr = SumTwoSSysErr(gr[0], syserr_unfold->GetEYlow(), syserr_unfold->GetEYhigh(), syserr_tpc7->GetEYlow(), syserr_tpc7->GetEYhigh());

	//ClassSysErr *syserr = SumTwoSSysErr(gr[0], syserr_unfold->GetEYlow(), syserr_unfold->GetEYhigh(), syserr_tpc->GetEYlow(), syserr_tpc->GetEYhigh());

	TString sregion = "Transverse";
	if(filename.Contains("Lead",TString::kIgnoreCase)) sregion = "Toward";
	if(filename.Contains("Away",TString::kIgnoreCase) || filename.Contains("Sub",TString::kIgnoreCase)) sregion = "Away";
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) sregion = "Leading Jet";
	TString svariable = " #LTN_{ch}/#delta#eta#delta#phi#GT";		// Ntrk
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) svariable = "#LTp_{T}^{ch}#GT";
	if(filename.Contains("PtSum",TString::kIgnoreCase) || filename.Contains("SumPt",TString::kIgnoreCase)) svariable = "#sump_{T}^{ch}";
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
		if(sregion.Contains("Tran")) hp[0]->SetMaximum(0.9);//Tran PtSum
		if(sregion.EqualTo("Toward")||sregion.EqualTo("Away"))  hp[0]->SetMaximum(10);//Lead PtSum
	}
	if(filename.Contains("PtAve",TString::kIgnoreCase) || filename.Contains("AvePt",TString::kIgnoreCase)) {
		if(sregion.Contains("Tran")) hp[0]->SetMaximum(0.68);//Tran PtAve
		if(sregion.EqualTo("Twoard")||sregion.EqualTo("Away"))  hp[0]->SetMaximum(4);//Lead PtAve
	}
	if(filename.Contains("LeadJetNtrk",TString::kIgnoreCase)) hp[0]->SetMaximum(31);	// Lead Jet Ntrk (Charged+Neutral)



	TCanvas *c = new TCanvas();
	c->SetLeftMargin(0.12);
	c->SetBottomMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);


	SetHistStyle(hp[0],1,1,8,1,2,1);
	hp[0]->GetXaxis()->SetRangeUser(0,55);
	hp[0]->Draw("p");
	graphSystBand(Nbins,syserr->GetX(),syserr->GetY(),0,0,syserr->GetEYlow(),  syserr->GetEYhigh(),kGray);
	hp[0]->Draw("psame");
	for(int i = 1; i<Nfile; i++) {
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i>(BayesTimes+2)?i-(BayesTimes+3):i-1]==DefaultTh))
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
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i-1]==DefaultTh))
		leg->AddEntry(gr[i],Form("w/ NFweight iter=%d",OtherTh[i-1]),"p");
	}
	leg->AddEntry(gr[BayesTimes+1],"w/ NFweight Bin-by-Bin","p");
	leg->AddEntry(gr[BayesTimes+2],Form("w/o NFweight iter=%d",DefaultTh),"p");
	for(int i = BayesTimes+3; i<(BayesTimes+2)*2-1; i++) {
		if(!(filename.Contains("LeadAreaNtrk")&&OtherTh[i-(BayesTimes+3)]==DefaultTh))
		leg->AddEntry(gr[i],Form("w/o NFweight iter=%d",OtherTh[i-(BayesTimes+3)]),"p");
	}
	leg->AddEntry(gr[(BayesTimes+2)*2-1],"w/o NFweight Bin-by-Bin","p");
	leg->AddEntry(gr[(BayesTimes+2)*2],"Scaled meas to match MB","p");
	for(int i = (BayesTimes+2)*2+1; i<(BayesTimes+2)*2+1+TpcTimes; i++) {
		leg->AddEntry(gr[i],TpcName[i-((BayesTimes+2)*2+1)],"p");
	}
	leg->Draw();


	if(1) {
		TFile *fout = new TFile(("SysErr4"+filename),"RECREATE");
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
		hmeas->ProfileX()->Write();
		htt->ProfileX()->Write();
		ht->ProfileX()->Write();

		fout->Close();
	}

}
