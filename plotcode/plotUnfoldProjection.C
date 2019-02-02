//==========================================================================
//
//	2018.03.24	Li YI
//	plot unfolded result (Ntrk or PtAve distribution) for each pT bins.
//
//==========================================================================
#include "plothrecoWunfolderr.C"			// void SetHistStyle(TProfile *h, int mcolor, int lcolor, int mstyle, int lstyle, float msize, int lwidth) 
                                                        // void graphSystBand(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
                                                        //                    Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBand=5,
                                                        //                    Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1) 
                                                        // void graphSystBox(int n, Double_t *x, Double_t *y, Double_t *ex, Double_t *ey,
                                                        //                Double_t *syPlus=0, Double_t *syMinus=0, Int_t colorBox=5,
                                                        //                Double_t dxaxis=1, Int_t symb=20, Double_t sSize=1, Int_t color=1, Int_t lwidth=1, Int_t fillstyleBox=1001)
                                                        // 
							// include "compare2BinByBin.C"           // TH2D *GetBbB(TH2D* htt, TH2D* ht, TH2D* hm)
							// include /Users/li/commonmacro/ClassSysErr.C

const int BayesTimes = 6;
const int Nfunfold = BayesTimes*2+1;		// NFweighted vs NoNFweighted, + Yscaled
const int Nftpc = 4;
const int Nfbemc = 2;
const int Nfile = Nfunfold+Nftpc+Nfbemc;

const char *BayesName[BayesTimes] = {"_Baye2","_Baye3","","_Baye5","_Baye6","_Baye7"};		// 2,3,4,5,6,7
const char *TpcName[Nftpc] = {"_TpcErrPlus004","_TpcErrMinus004","_TpcErrPlusAbs004","_TpcErrMinusAbs004"};
const char *BemcName[Nfbemc] = {"_BemcErrPlus004","_BemcErrMinus004"};


TString GetFileName(TString filenamedef="Unfolding_%sJPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root", TString tag="TranPtAve", int ith=0) {

	TString outfilename = "";

	if(ith<0) {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<" ith ("<<ith<<") cannot be smaller than 0!!!"<<endl;
		return "";
	}

	TString basefilename = Form(filenamedef,tag.Data());
	if(basefilename.Contains("Baye",TString::kIgnoreCase)) {
		Ssiz_t toinsert = basefilename.Index(".root");
		Ssiz_t toremove = basefilename.Index("_Baye");
                if(toremove!=-1) basefilename.Remove(toremove,(toinsert-toremove));         // WARNINING!!  Assuming _Baye%d.root

	}

	outfilename = basefilename;

	Ssiz_t toinsert = basefilename.Index(".root");
	if(ith>=0 && ith<BayesTimes) {	//_NFWeight
		outfilename.Insert(toinsert,BayesName[ith]);
	}
	else if(ith<BayesTimes*2) {	//no NFWeight
		outfilename.Insert(toinsert,BayesName[ith-BayesTimes]);
		outfilename.ReplaceAll("embedMB","embedJP0");
		outfilename.ReplaceAll("_NFWeight","");
	}
	else if(ith<Nfunfold) {		// Yscale
		outfilename = Form(filenamedef,tag.Data());
		outfilename.ReplaceAll("_NFWeight","_YScale");
	}
	else if(ith<Nfunfold+Nftpc) {	// TPC
                outfilename.Insert(toinsert,TpcName[ith-Nfunfold]);
	}
	else if(ith<Nfile) {		// BEMC
                outfilename.Insert(toinsert,BemcName[ith-Nfunfold-Nftpc]);
	}
	else {
		cout<<"ERR!! "<<__PRETTY_FUNCTION__<<" ith ("<<ith<<") cannot be larger than Nfile!!!"<<endl;
		return "";
	}

	cout<<"to read "<<outfilename<<endl;
	return outfilename;
}


TGraphErrors *ConvertHist2Graph(TH1D *h) {

	TGraphErrors *gr = new TGraphErrors();
	gr->SetTitle(Form("gr_%s",h->GetTitle()));
	gr->SetName(Form("gr_%s",h->GetName()));
	int Nbins = h->GetNbinsX();
	for(int i = 0; i<Nbins; i++) {
		gr->SetPoint(i,h->GetBinCenter(i+1),h->GetBinContent(i+1));	
		gr->SetPointError(i,0,h->GetBinError(i+1));	
	}
	gr->SetMarkerStyle(h->GetMarkerStyle());
	gr->SetMarkerColor(h->GetMarkerColor());
	gr->SetLineStyle(h->GetLineStyle());
	gr->SetLineColor(h->GetLineColor());

	return gr;
}


void plotUnfoldProjection(TString tag = "TranTotNtrk", TString filename="Unfolding_%sJPCharged_NFWeight_BT170928_RcVzW_12JetBinv2_McPt02_embedMB_Baye5.root") {

	//Sys. Err. handling
	TString opt_unfold = "s";                       // will take the max of all as the final unfolding sys err
	TString opt_tracking = "s";
	TString opt_towergain = "s";


	TFile *fin[Nfile];

	for(int i = 0; i<Nfile; i++) {
       		fin[i] = new TFile(GetFileName(filename,tag,i));
	}

	int idefault = 3;	// _Bayes5
	if(!filename.Contains("Baye",TString::kIgnoreCase)) idefault = 2;	// Bayes == 4

	const int Nhist = Nfile+2; 	// + 2 BinByBin
	TH2D *hmeas[Nhist];		
	TH2D *hreco[Nhist];	
	TH2D *ht[Nhist];		
	TH2D *htt[Nhist];	
	for(int i = 0; i<Nfile; i++) {
       		hmeas[i] = (TH2D*)fin[i]->Get("hmeas");
		hreco[i] = (TH2D*)fin[i]->Get("hreco");
		ht[i] = (TH2D*)fin[i]->Get("htrain");
		htt[i] = (TH2D*)fin[i]->Get("htraintrue");
	}

	// BinByBin
	//TH2D *htt[2], *ht[2];		// NFWeighted and no weight
       	//htt[0] = (TH2D*)fin[idefault]->Get("htrain");
       	//ht[0] = (TH2D*)fin[idefault]->Get("htraintrue");
       	//htt[1] = (TH2D*)fin[BayesTimes+idefault]->Get("htrain");
       	//ht[1] = (TH2D*)fin[BayesTimes+idefault]->Get("htraintrue");
	ht[Nfile] = ht[idefault];
	ht[Nfile+1] = ht[BayesTimes+idefault];
	htt[Nfile] = htt[idefault];
	htt[Nfile+1] = htt[BayesTimes+idefault];
	hmeas[Nfile] = hmeas[idefault];
	hmeas[Nfile+1] = hmeas[BayesTimes+idefault];
	hreco[Nfile] = (TH2D*)GetBbB(htt[idefault],ht[idefault],hmeas[idefault]);
	hreco[Nfile+1] = (TH2D*)GetBbB(htt[BayesTimes+idefault],ht[BayesTimes+idefault],hmeas[BayesTimes+idefault]);

	TH1D *hmeaspy[Nhist];	
	TH1D *hrecopy[Nhist];	
	TH1D *htpy[Nhist];	
	TH1D *httpy[Nhist];	
	
	TList *listC = new TList();	// canvas list
	gStyle->SetOptStat(0);

	float pt1,pt2;
	TProfile *pfxmeas = (TProfile*)fin[idefault]->Get("pfxmeas");

	TString xtitle;
	if(tag.Contains("Tran",TString::kIgnoreCase)) xtitle="Transverse ";
	if(tag.Contains("PtAve",TString::kIgnoreCase)) xtitle+="#LTp_{T}#GT";
	if(tag.Contains("TotNtrk",TString::kIgnoreCase)) xtitle+="Total Multiplicity";

	TString ytitle;
	if(tag.Contains("PtAve",TString::kIgnoreCase)) ytitle = "Prob./#Deltap_{T}";
	if(tag.Contains("Ntrk",TString::kIgnoreCase)) ytitle = "Prob.";

	int fillcolor_meas = kGray, fillcolor_reco = kOrange;	// not really use for plot hist. It is for legend with band (for sys. err plot)
	int fillstyle = 3002;

	//for(int i = 0; i<1; i++) {
	for(int i = 0; i<hmeas[idefault]->GetNbinsX(); i++) {

		for(int j = 0; j<Nhist; j++) {
			//cout<<"i = "<<i<<" j = "<<j<<" hmeas = "<<hmeas[j]<<endl;
			hmeaspy[j] = (TH1D*)hmeas[j]->ProjectionY(Form("hmeaspy%d_at%d",j+1,i+1),i+1,i+1);
			hrecopy[j] = (TH1D*)hreco[j]->ProjectionY(Form("hrecopy%d_at%d",j+1,i+1),i+1,i+1);
			if(j<BayesTimes) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d %s NFWeighted",tag.Data(),i+1,BayesName[j]));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d %s NFWeighted",tag.Data(),i+1,BayesName[j]));
			}
			else if(j<2*BayesTimes) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d %s no NFWeighted",tag.Data(),i+1,BayesName[j-BayesTimes]));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d %s no NFWeighted",tag.Data(),i+1,BayesName[j-BayesTimes]));
			}
			else if(j<Nfunfold) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d YScaled",tag.Data(),i+1));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d YScaled",tag.Data(),i+1));
			}
			else if(j<Nfunfold+Nftpc) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d tpc%d NFWeighted",tag.Data(),i+1,j+1-Nfunfold));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d tpc%d NFWeighted",tag.Data(),i+1,j+1-Nfunfold));
			}
			else if(j<Nfunfold+Nftpc+Nfbemc) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d bemc%d NFWeighted",tag.Data(),i+1,j+1-Nfunfold-Nftpc));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d bemc%d NFWeighted",tag.Data(),i+1,j+1-Nfunfold-Nftpc));
			}
			else if(j<Nfunfold+Nftpc+Nfbemc+1) {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d BinbyBin NFWeighted",tag.Data(),i+1));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d BinbyBin NFWeighted",tag.Data(),i+1));
			}
			else {
				hmeaspy[j]->SetTitle(Form("Measured %s at pt%d BinbyBin no NFWeighted",tag.Data(),i+1));
				hrecopy[j]->SetTitle(Form("Unfolded %s at pt%d BinbyBin no NFWeighted",tag.Data(),i+1));
			}
			//cout<<"i = "<<i<<" j = "<<j<<" hmeaspy = "<<hmeaspy[j]<<endl;
			
			htpy[j] = (TH1D*)ht[j]->ProjectionY(Form("htpy%d_at%d",j+1,i+1),i+1,i+1);
			httpy[j] = (TH1D*)htt[j]->ProjectionY(Form("httpy%d_at%d",j+1,i+1),i+1,i+1);

			hmeaspy[j]->GetXaxis()->SetTitle(xtitle);
			hrecopy[j]->GetXaxis()->SetTitle(xtitle);
			hmeaspy[j]->GetYaxis()->SetTitle(ytitle);
			hrecopy[j]->GetYaxis()->SetTitle(ytitle);
			hmeaspy[j]->GetYaxis()->SetTitleOffset(0.8);
			hrecopy[j]->GetYaxis()->SetTitleOffset(0.8);

			//Normalization
			double sum_meas = hmeaspy[j]->Integral();
			double sum_reco = hrecopy[j]->Integral();
			if(sum_meas>0) hmeaspy[j]->Scale(1./sum_meas);
			if(sum_reco>0) hrecopy[j]->Scale(1./sum_reco);
			double sum_t = htpy[j]->Integral();
			double sum_tt = httpy[j]->Integral();
			if(sum_t>0) htpy[j]->Scale(1./sum_t);
			if(sum_tt>0) httpy[j]->Scale(1./sum_tt);

			// Set Style
			int reco_MarkerStyle = 4;
			int meas_MarkerStyle = 25;
			if(j == idefault) {
				reco_MarkerStyle = 8;
				meas_MarkerStyle = 21;
			}
			SetHistStyle(hmeaspy[j],1,1,meas_MarkerStyle,1,1,1);
			SetHistStyle(hrecopy[j],2,2,reco_MarkerStyle,1,1,1);
			hmeaspy[j]->SetFillColor(fillcolor_meas);
			hrecopy[j]->SetFillColor(fillcolor_reco);
			hmeaspy[j]->SetFillStyle(fillstyle);
			hrecopy[j]->SetFillStyle(fillstyle);
			SetHistStyle(htpy[j],1,1,meas_MarkerStyle,2,1,1);
			SetHistStyle(httpy[j],2,2,reco_MarkerStyle,2,1,1);

			// for pt dist, normalized by bin width
			if(tag.Contains("PtAve",TString::kIgnoreCase)) {
				hmeaspy[j]->Scale(1,"width"); 
				hrecopy[j]->Scale(1,"width");
				htpy[j]->Scale(1,"width"); 
				httpy[j]->Scale(1,"width");
			}

		}


		// Convert to TGraphErrors for syst. err. ClassSysErr
		TGraphErrors *gr_measpy[Nhist], *gr_recopy[Nhist];
		for(int j = 0; j<Nhist; j++) {
			gr_measpy[j] = (TGraphErrors*)ConvertHist2Graph(hmeaspy[j]);
			gr_recopy[j] = (TGraphErrors*)ConvertHist2Graph(hrecopy[j]);

			//cout<<"i = "<<i<<" j = "<<j<<" gr_measpy Npoints = "<<gr_measpy[j]->GetN()<<endl;
		}
		// reorder 
		TGraphErrors *gr_measpy_unfold[Nfunfold+2+2-1];	//  NFweighted vs NoNFweighted, + Yscaled, + BinByBin (for default NFweighted and no weight ones)
		//TGraphErrors *gr_measpy_tpc[Nftpc];
		//TGraphErrors *gr_measpy_bemc[Nfbemc];
		TGraphErrors *gr_recopy_unfold[Nfunfold+2+2-1];	//  NFweighted vs NoNFweighted, + Yscaled, + BinByBin (for default NFweighted and no weight ones)

		for(int k = 0, g = 0; k<Nfunfold; k++) {
			if(k!=idefault) {
				//cout<<"g = "<<g<<endl;
				gr_measpy_unfold[g] = gr_measpy[k];
				gr_recopy_unfold[g] = gr_recopy[k];
				g++;
			}
		}
		gr_measpy_unfold[Nfunfold-1] = gr_measpy[Nfile];
		gr_measpy_unfold[Nfunfold] = gr_measpy[Nfile+1];
		gr_recopy_unfold[Nfunfold-1] = gr_recopy[Nfile];
		gr_recopy_unfold[Nfunfold] = gr_recopy[Nfile+1];

		//for(int k = 0; k<Nfunfold+1; k++) {
		//	//cout<<k<<" "<<gr_measpy_unfold[k]->GetN()<<endl;
		//	gr_measpy_unfold[k]->Print();
		//	gr_recopy_unfold[k]->Print();
		//}

		// Syst. Err.
		ClassSysErr *syserr_measunfold = new ClassSysErr(Nfunfold+1, gr_measpy[idefault], gr_measpy_unfold, opt_unfold);
		ClassSysErr *syserr_recounfold = new ClassSysErr(Nfunfold+1, gr_recopy[idefault], gr_recopy_unfold, opt_unfold);
		//cout<<"syserr_unfold"<<endl;
		//syserr_unfold->Print();
		//cout<<endl;
		ClassSysErr *syserr_meastpc = new ClassSysErr(Nftpc, gr_measpy[idefault], &gr_measpy[Nfunfold], opt_tracking);
		ClassSysErr *syserr_recotpc = new ClassSysErr(Nftpc, gr_recopy[idefault], &gr_recopy[Nfunfold], opt_tracking);
		//cout<<"syserr_meastpc"<<endl;
		//syserr_meastpc->Print();
		//cout<<endl;
		ClassSysErr *syserr_measbemc = new ClassSysErr(Nfbemc, gr_measpy[idefault], &gr_measpy[Nfunfold+Nftpc], opt_towergain);
		ClassSysErr *syserr_recobemc = new ClassSysErr(Nfbemc, gr_recopy[idefault], &gr_recopy[Nfunfold+Nftpc], opt_towergain);
		//cout<<"syserr_bemc"<<endl;
		//syserr_bemc->Print();
		//cout<<endl;

        	ClassSysErr *syserr_meas=NULL;
        	if(syserr_measunfold&&syserr_meastpc) {
        	        syserr_meas = SumTwoSSysErr(gr_measpy[idefault], syserr_measunfold->GetEYlow(), syserr_measunfold->GetEYhigh(), syserr_meastpc->GetEYlow(), syserr_meastpc->GetEYhigh());
        	        if(syserr_measbemc) {
        	                ClassSysErr *syserr2 = syserr_meas;
        	                syserr_meas = SumTwoSSysErr(gr_measpy[idefault], syserr2->GetEYlow(), syserr2->GetEYhigh(), syserr_measbemc->GetEYlow(), syserr_measbemc->GetEYhigh());
        	        }
        	}
        	else if(syserr_measunfold) {
        	        syserr_meas = syserr_measunfold;
        	}
        	else if(syserr_meastpc) {
        	        syserr_meas = syserr_meastpc;
        	}
        	else if(syserr_measbemc) {
        	        syserr_meas = syserr_measbemc;
        	}

		//syserr_meas->Print();

        	ClassSysErr *syserr_reco=NULL;
        	if(syserr_recounfold&&syserr_recotpc) {
        	        syserr_reco = SumTwoSSysErr(gr_recopy[idefault], syserr_recounfold->GetEYlow(), syserr_recounfold->GetEYhigh(), syserr_recotpc->GetEYlow(), syserr_recotpc->GetEYhigh());
        	        if(syserr_recobemc) {
        	                ClassSysErr *syserr_reco2 = syserr_reco;
        	                syserr_reco = SumTwoSSysErr(gr_recopy[idefault], syserr_reco2->GetEYlow(), syserr_reco2->GetEYhigh(), syserr_recobemc->GetEYlow(), syserr_recobemc->GetEYhigh());
        	        }
        	}
        	else if(syserr_recounfold) {
        	        syserr_reco = syserr_recounfold;
        	}
        	else if(syserr_recotpc) {
        	        syserr_reco = syserr_recotpc;
        	}
        	else if(syserr_recobemc) {
        	        syserr_reco = syserr_recobemc;
        	}

		//syserr_reco->Print();

		pt1 = pfxmeas->GetBinLowEdge(i+1);
		pt2 = pt1+pfxmeas->GetBinWidth(i+1);
		cout<<"pt "<<pt1<<"-"<<pt2<<endl;



		// Draw
		TCanvas *ic = new TCanvas();
		gStyle->SetOptTitle(0);
		ic->Clear();
		if(tag.Contains("PtAve",TString::kIgnoreCase)) {ic->SetLogy(); ic->SetLogx();}			
		if(tag.Contains("Ntrk",TString::kIgnoreCase)) {ic->SetLogy(); }			
		listC->Add(ic);
		hmeaspy[idefault]->GetXaxis()->SetRangeUser(0,20);
		if(tag.Contains("PtAve",TString::kIgnoreCase)) hmeaspy[idefault]->GetYaxis()->SetRangeUser(5e-5,9);
		if(tag.Contains("Ntrk",TString::kIgnoreCase)) hmeaspy[idefault]->GetYaxis()->SetRangeUser(5e-5,0.9);
		hmeaspy[idefault]->Draw("eX0");
		hrecopy[idefault]->Draw("eX0same");
		//for(int ik = 0; ik<Nhist; ik++) {
		//	hmeaspy[ik]->Draw("same");
		//	hrecopy[ik]->Draw("same");
		//}

		//draw sys. err. box
	        double *ex = new double[hmeaspy[idefault]->GetNbinsX()];
        	for( int ix = 0; ix<hmeaspy[idefault]->GetNbinsX(); ix++) {
        	        ex[ix] = hmeaspy[idefault]->GetBinWidth(ix+1)/2.;
        	}
		//as we don't consider sys.err at measured level. All QA uncertainties are in unfolding procedure. graphSystBox(hmeaspy[idefault]->GetNbinsX(),syserr_meas->GetX(),syserr_meas->GetY(),ex,0, syserr_meas->GetEYlow(), syserr_meas->GetEYhigh(), fillcolor_meas, 0,0,1,0,1,fillstyle);
		graphSystBox(hrecopy[idefault]->GetNbinsX(),syserr_reco->GetX(),syserr_reco->GetY(),ex,0, syserr_reco->GetEYlow(), syserr_reco->GetEYhigh(), fillcolor_reco, 0,0,1,0,1,fillstyle);

		htpy[idefault]->Draw("HISTsame");
		httpy[idefault]->Draw("HISTsame");

		hmeaspy[idefault]->Draw("eX0same");
		hrecopy[idefault]->Draw("eX0same");

		TLegend *il = new TLegend(0.6,0.55,0.88,0.86);
		il->AddEntry(hmeaspy[idefault],"Measured","p");
		il->AddEntry(hrecopy[idefault],"Unfolded","pf");
		il->AddEntry(htpy[idefault],"MC Detector-level","l");
		il->AddEntry(httpy[idefault],"MC Particle-level","l");
		il->SetBorderSize(0);
		il->SetFillColor(0);
		il->SetHeader(Form("%g < p_{T} < %g GeV/#it{c}",pt1,pt2));
		il->Draw("same");


		//ic->SaveAs(Form("%sDistUnfolded_pt%g-%g.png",tag.Data(),pt1,pt2));
		//ic->SaveAs(Form("%sDistUnfolded_pt%g-%g.pdf",tag.Data(),pt1,pt2));
		//ic->SaveAs(Form("/Users/li/Documents/paperproposal/UnderlyingEvent/AnaNote/fig_ananote/%sDistUnfolded_pt%g-%g.pdf",tag.Data(),pt1,pt2));
	}


}
