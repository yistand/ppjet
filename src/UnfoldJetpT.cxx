//======================================================================================================================================================
//
//		2016.07.25		Li YI
//		use RooUnfold to unfold leading jet pT using Mc - Rc from embedding
//
//		To run this code:
//		root -l src/UnfoldJetpT.cxx
//
//
//		2016.08.01		Li YI
//		Reconstruct responsive matrix with fine binning, so that we can correct change prior (weight within pT bins).
//		Unfold with wide binning to reduce the bin-to-bin migration due to statistics.
//
//======================================================================================================================================================
#if !(defined(__CINT__) || defined(__CLING__))
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldBinByBin.h"
#endif

void UseAnotherPrior(TH1D *truth, TH2D *cov, TH1D *measured, TH1D *modtruth, TH2D *modcov, TH1D *modmeasured, float par=0.01) {
	
	TF1 *f1 = new TF1("f1","exp(-[0]*x)",0,100);
	f1->SetParameter(0,par);
	int nbins=truth->GetNbinsX();
	for(int i = 0; i<nbins; i++) {
		float iscale = f1->Eval(i);
		modtruth->SetBinContent(i+1,iscale*truth->GetBinContent(i+1)); 
		modtruth->SetBinError(i+1,iscale*truth->GetBinError(i+1));
		for(int j = 0; j<nbins; j++) {
			modcov->SetBinContent(j+1,i+1,iscale*cov->GetBinContent(j+1,i+1));
			modcov->SetBinError(j+1,i+1,iscale*cov->GetBinError(j+1,i+1));
		}
	}

	TH1D *modcov_px=(TH1D*)modcov->ProjectionX("modcov_px");
	TH1D *cov_px=(TH1D*)cov->ProjectionX("cov_px");
	for(int j=0; j<nbins; j++) {
		float jscale = cov_px->GetBinContent(j+1)>0?modcov_px->GetBinContent(j+1)/cov_px->GetBinContent(j+1):0;
		modmeasured->SetBinContent(j+1,jscale*measured->GetBinContent(j+1));
		modmeasured->SetBinError(j+1,jscale*measured->GetBinError(j+1));
	}
}

TH2D *Rebin2DHisto(TH2D *old, int Nx, double *xbins, int Ny, double *ybins){
	TH2D *h = new TH2D(Form("Rebin%s",Form("Rebin%s",old->GetName())),old->GetTitle(),Nx,xbins,Ny,ybins);	
	h->Sumw2();
	TAxis *xaxis = old->GetXaxis();
	TAxis *yaxis = old->GetYaxis();
	for (int j=1; j<=yaxis->GetNbins();j++) {
		for (int i=1; i<=xaxis->GetNbins();i++) {
			h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j));
		}
	}
	return h;
}



//========================================= MAIN ===========================================
void UnfoldJetpT(int par=4) {

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	cout << "==================================== TRAIN ====================================" << endl;
	TFile *fin = new TFile("UnfoldMatrxLead_Online_nobbc.root");
	TH1D *measured = (TH1D*)fin->Get("Rc");
	TH1D *truth = (TH1D*)fin->Get("Mc");
	TH2D *cov = (TH2D*)fin->Get("CovMatrix");
	measured->Sumw2();
	truth->Sumw2();
	cov->Sumw2();

	TH1D *measured_save = (TH1D*)measured->Clone("Rc");
	TH1D *truth_save = (TH1D*)truth->Clone("Mc");
	TH2D *cov_save = (TH2D*)cov->Clone("CovMatrix");

	RooUnfoldResponse response (measured, truth, cov);
	//RooUnfoldResponse response (measured, cov->ProjectionY(), cov); // test
	//RooUnfoldResponse response (cov->ProjectionX(),truth, cov); // test
	//RooUnfoldResponse response (cov->ProjectionX(), cov->ProjectionY(), cov); // test

	int iter = par;			// for RooUnfoldBayes		// default 4
	int kterm = par;			// for RooUnfoldSvd	// default 4

	cout << "==================================== self TEST =====================================" << endl;
	RooUnfoldBayes   unfold (&response, measured, iter);  
	//RooUnfoldSvd   unfold (&response, measured, kterm);  
	//RooUnfoldBinByBin unfold(&response, measured);
	TH1D* reco= (TH1D*) unfold.Hreco();

	//#define CHANGE_PRIOR

	// if change the prior
	#ifdef CHANGE_PRIOR
	cout << "INFO: Use the modified prior " << endl;
	reco->Reset();
	TH1D *modtruth = (TH1D*)truth->Clone("modtruth");
	modtruth->Reset();
	TH1D *modmeasured = (TH1D*)measured->Clone("modmeasured");
	modmeasured->Reset();
	TH2D *modcov = (TH2D*)cov->Clone("modcov");
	modcov->Reset();
	UseAnotherPrior(truth, cov, measured, modtruth, modcov, modmeasured, -0.1);
	RooUnfoldResponse modresponse (modmeasured, modtruth, modcov);
	RooUnfoldBayes   modunfold (&modresponse, measured, iter);  
	//RooUnfoldSvd   modunfold (&modresponse, measured, kterm);  
	//RooUnfoldBinByBin modunfold(&modresponse, measured);
	reco= (TH1D*) modunfold.Hreco();
	#endif

	reco->SetName("unfold_train");

	TCanvas *c1 = new TCanvas("ctest");
	c1->SetLogy();
	reco->SetLineColor(1);
	reco->GetXaxis()->SetTitle("Leading Jet p_{T}");
	reco->GetYaxis()->SetTitle("d(Number of events)/dp_{T}");
	reco->DrawClone("h");
	measured->SetLineColor(kOrange);
	measured->SetLineWidth(2);
	measured->DrawClone("psame");
	truth->SetLineColor(kRed);
	truth->DrawClone("psame");
	#ifdef CHANGE_PRIOR
	modtruth->Draw("psame");
	#endif
	TLegend *leg0 = new TLegend(0.6,0.7,0.85,0.86);
	leg0->SetFillColor(0);
	leg0->AddEntry(reco,"MC unfolded","pl");
	leg0->AddEntry(truth,"MC particle-level truth","pl");
	leg0->AddEntry(measured,"MC detector-level","pl");
	leg0->DrawClone("same");


	// plot ratio
	TH1D *hr = (TH1D*)reco->Clone("Ratio");
	hr->Divide(truth);

	TCanvas *cr1 = new TCanvas("cmcratio");
	hr->GetXaxis()->SetTitle("Leading Jet p_{T}");
	hr->GetYaxis()->SetTitle("MC unfolded/particle-level truth Ratio");
	hr->DrawClone();
	TLine *line = new TLine(0,1,hr->GetBinLowEdge(hr->GetNbinsX()+1),1);
	line->SetLineWidth(2);
	line->DrawClone();



	cout << "==================================== UNFOLDING  =====================================" << endl;
	TFile *freal = new TFile("/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP2_151030_P12id.root");
	TH1D *LeadJetPt_data = (TH1D*)freal->Get("LeadingJetPt");
	TH1D *LeadJetPt_data_save = (TH1D*)LeadJetPt_data->Clone("LeadJetPt_data");

	//#define CUT_HIGH_PT
	//#define CUT_LOW_PT
	
	double cutoffhigh=50;
	double cutofflow=7;

	#ifdef CUT_HIGH_PT
	cout << "INFO: CUT off high pT in Rc jet " << endl;
	for(int i = 0; i<LeadJetPt_data->GetNbinsX(); i++) {
		if(LeadJetPt_data->GetBinCenter(i+1)>cutoffhigh) {
			LeadJetPt_data->SetBinContent(i+1,0);	
			LeadJetPt_data->SetBinError(i+1,0);	
		}
	}
	#endif

	#ifdef CUT_LOW_PT
	cout << "INFO: CUT off low pT in Rc jet " << endl;
	for(int i = 0; i<LeadJetPt_data->GetNbinsX(); i++) {
		if(LeadJetPt_data->GetBinCenter(i+1)<cutofflow) {
			LeadJetPt_data->SetBinContent(i+1,0);	
			LeadJetPt_data->SetBinError(i+1,0);	
		}
	}
	#endif

	#if defined(CUT_HIGH_PT) && defined(CUT_LOW_PT) 
	for(int i = 0; i<measured->GetNbinsX(); i++) {
                if((measured->GetBinCenter(i+1)>cutoffhigh) || (measured->GetBinCenter(i+1)<cutofflow)) {
                        measured->SetBinContent(i+1,0);
                        measured->SetBinError(i+1,0);
                }
        }
	for(int i = 0; i<truth->GetNbinsX(); i++) {
                if((truth->GetBinCenter(i+1)>cutoffhigh) || (truth->GetBinCenter(i+1)<cutofflow)) {
                        truth->SetBinContent(i+1,0);
                        truth->SetBinError(i+1,0);
                }
        }
	for(int i = 0; i<cov->GetNbinsX(); i++) {
		for(int j = 0; j<cov->GetNbinsY(); j++) {
                	if((cov->GetXaxis()->GetBinCenter(i+1)>cutoffhigh) || (cov->GetXaxis()->GetBinCenter(i+1)<cutofflow) || (cov->GetYaxis()->GetBinCenter(j+1)>cutoffhigh) || (cov->GetYaxis()->GetBinCenter(j+1)<cutofflow)) {
                	        cov->SetBinContent(i+1,j+1,0);
                	        cov->SetBinError(i+1,j+1,0);
                	}
		}
        }
	#endif


	#define WIDEBIN
	#ifdef WIDEBIN
	// ====================== Use wide binning for actually unfolding ======================
        const int WNbins = 15;          // wide binning for unfolding procedure
        double Wptbins[WNbins+1] = {0,2,3,4,5,7,9,11,15,20,25,35,45,55,65,100};

	TH2D *Wcov = Rebin2DHisto(cov,WNbins,Wptbins,WNbins,Wptbins);
	TH1D *Wmeasured = (TH1D*)measured->Rebin(WNbins,"Wmeasured",Wptbins);
	TH1D *Wtruth = (TH1D*)truth->Rebin(WNbins,"Wtruth",Wptbins);
	RooUnfoldResponse Wresponse (Wmeasured, Wtruth, Wcov);
	TH1D *WLeadJetPt_data = (TH1D*)LeadJetPt_data->Rebin(WNbins,"WLeadJetPt_data",Wptbins);

	TH1D *Wtruth_save = (TH1D*)Wtruth->Clone("Wtruth");
	TH1D *Wmeasured_save = (TH1D*)Wmeasured->Clone("Wmeasured");
	TH2D *Wcov_save = (TH2D*)Wcov->Clone("Wcov");
	#endif


	#if !defined(CHANGE_PRIOR) && !defined(WIDEBIN)
	RooUnfoldBayes   unfold_data (&response, LeadJetPt_data , iter);
	//RooUnfoldSvd   unfold_data (&response, LeadJetPt_data , kterm);
	//RooUnfoldBinByBin unfold_data(&response, LeadJetPt_data);

	#elif !defined(CHANGE_PRIOR) && defined(WIDEBIN)
	RooUnfoldBayes   unfold_data (&Wresponse, WLeadJetPt_data , iter);
	//RooUnfoldSvd   unfold_data (&Wresponse, WLeadJetPt_data , kterm);
	//RooUnfoldBinByBin unfold_data(&Wresponse, WLeadJetPt_data);

	#elif defined(CHANGE_PRIOR) && !defined(WIDEBIN) 
	cout << "INFO: Use the modified prior " << endl;
	RooUnfoldBayes   unfold_data (&modresponse, LeadJetPt_data , iter);
	//RooUnfoldSvd   unfold_data (&modresponse, LeadJetPt_data , kterm);
	//RooUnfoldBinByBin unfold_data(&modresponse, LeadJetPt_data);

	#elif defined(CHANGE_PRIOR) && defined(WIDEBIN) 
	cout << "INFO: Use the modified prior " << endl;
	TH2D *Wmodcov = Rebin2DHisto(modcov,WNbins,Wptbins,WNbins,Wptbins);
	TH1D *Wmodmeasured = (TH1D*)modmeasured->Rebin(WNbins,"Wmodmeasured",Wptbins);
	TH1D *Wmodtruth = (TH1D*)modtruth->Rebin(WNbins,"Wmodtruth",Wptbins);
	RooUnfoldResponse Wmodresponse (Wmodmeasured, Wmodtruth, Wmodcov);
	RooUnfoldBayes   unfold_data (&Wmodresponse, WLeadJetPt_data , iter);
	//RooUnfoldSvd   unfold_data (&Wmodresponse, WLeadJetPt_data , kterm);
	//RooUnfoldBinByBin unfold_data(&Wmodresponse, WLeadJetPt_data);

	#endif

	TH1D *UnfoldedLeadJetPt_data = (TH1D*)unfold_data.Hreco();
	UnfoldedLeadJetPt_data->SetName("UnfoldedLeadJetPt_data");

	//unfold_data.Print();		// print detailed bin-to-bin information

	cout << "Get nFake = "<<unfold_data.GetMeasureFake()<< endl;
	
	TH1D *UnfoldedLeadJetPt_data_save  = (TH1D*)UnfoldedLeadJetPt_data->Clone("unfold");


	float ptmin = 15; // 12; 10;
	float ptmax = 40; // 40; 50;
	int ixmin=truth->FindBin(ptmin);
	int ixmax=truth->FindBin(ptmax);
	int Uxmin=UnfoldedLeadJetPt_data->FindBin(ptmin);
	int Uxmax=UnfoldedLeadJetPt_data->FindBin(ptmax);
	float truthscale = truth->Integral(ixmin,ixmax)>0?UnfoldedLeadJetPt_data->Integral(Uxmin,Uxmax)/truth->Integral(ixmin,ixmax):0;
	cout << " scale Mc by " << truthscale << endl;
	if(truthscale>0)truth->Scale(truthscale);
	#ifdef WIDEBIN
	float Wtruthscale = Wtruth->Integral(Uxmin,Uxmax)>0?UnfoldedLeadJetPt_data->Integral(Uxmin,Uxmax)/Wtruth->Integral(Uxmin,Uxmax):0;
	cout << " scale WMc by " << Wtruthscale << endl;
	if(Wtruthscale>0)Wtruth->Scale(Wtruthscale);
	#endif
	float measuredscale = measured->Integral(ixmin,ixmax)>0?LeadJetPt_data->Integral(ixmin,ixmax)/measured->Integral(ixmin,ixmax):0;
	cout << " scale Rc by " << measuredscale << endl;
	if(measuredscale>0)measured->Scale(measuredscale);

	// Normalized by bin width
	#ifdef WIDEBIN
	UnfoldedLeadJetPt_data->Scale(1,"width");			// wide bin, scale it by width (NEED TO CHECK, the fine bining uses 1 as bin width, (100,0,100)
	Wtruth->Scale(1,"width");
	#endif

	TCanvas *c2 = new TCanvas("cunfold");
	c2->SetLogy();
	UnfoldedLeadJetPt_data->SetLineColor(1);
	UnfoldedLeadJetPt_data->SetMarkerStyle(8);
	LeadJetPt_data->SetLineColor(kBlue);
	LeadJetPt_data->SetMarkerColor(kBlue);
	LeadJetPt_data->SetMarkerStyle(8);
	truth->SetLineColor(kRed);
	truth->SetLineWidth(2);
	measured->SetLineColor(kOrange);
	measured->SetLineWidth(2);
	#if !defined(WIDEBIN)
	truth->GetXaxis()->SetTitle("Leading Jet p_{T}");
	truth->GetYaxis()->SetTitle("d(Number of events)/dp_{T}");
	truth->DrawClone("h");
	#else 
	Wtruth->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Wtruth->GetYaxis()->SetTitle("d(Number of events)/dp_{T}");
	Wtruth->DrawClone("h");
	#endif
	measured->DrawClone("hsame");
	UnfoldedLeadJetPt_data->DrawClone("psame");
	LeadJetPt_data->DrawClone("psame");
	#if defined(CHANGE_PRIOR) && !defined(WIDEBIN)
	modtruth->Scale(truthscale);
	modtruth->SetLineColor(kGreen+1);
	modtruth->SetMarkerColor(kGreen+1);
	modtruth->Draw("csame");
	#elif defined(CHANGE_PRIOR) && defined(WIDEBIN)
	Wmodtruth->Scale(Wtruthscale);
	Wmodtruth->SetLineColor(kGreen+1);
	Wmodtruth->SetMarkerColor(kGreen+1);
	Wmodtruth->Draw("csame");
	#endif
	TLegend *leg1 = new TLegend(0.6,0.7,0.85,0.86);
	leg1->SetFillColor(0);
	leg1->AddEntry(LeadJetPt_data,"Measured data","pl");
	leg1->AddEntry(UnfoldedLeadJetPt_data,"Unfolded data","pl");
	leg1->AddEntry(measured,"MC detector-level (scaled)","pl");
	leg1->AddEntry(truth,"MC particle-level (scaled)","pl");
	#if defined(CHANGE_PRIOR) && !defined(WIDEBIN)
	leg1->AddEntry(modtruth,"prior for unfolding","l");
	#elif defined(CHANGE_PRIOR) && defined(WIDEBIN)
	leg1->AddEntry(Wmodtruth,"prior for unfolding","l");
	#endif
	leg1->DrawClone("same");

	// plot ratio
	TH1D *hmcr = (TH1D*)UnfoldedLeadJetPt_data->Clone("ParticleLev_Ratio");
	hmcr->Sumw2();
	#if !defined(WIDEBIN)
	hmcr->Divide(truth);
	#else 
	hmcr->Divide(Wtruth);
	#endif
	TH1D *hrcr = (TH1D*)LeadJetPt_data->Clone("DetectorLev_Ratio");
	hrcr->Sumw2();
	hrcr->Divide(measured);

	TCanvas *cr2 = new TCanvas("cratio");
	hrcr->GetXaxis()->SetTitle("Leading Jet p_{T}");
	hrcr->GetYaxis()->SetTitle("Ratio");
	hmcr->SetMarkerSize(1.5);
	hrcr->SetMaximum(2);
	hrcr->SetMinimum(0);
	hrcr->DrawClone();
	hmcr->DrawClone("same");
	TLine *line = new TLine(0,1,hrcr->GetBinLowEdge(hrcr->GetNbinsX()+1),1);
	line->SetLineWidth(2);
	line->DrawClone();
	TLegend *leg2 = new TLegend(0.6,0.7,0.85,0.86);
	leg2->SetFillColor(0);
	leg2->AddEntry(hrcr,"data/MC detector-level","pl");
	leg2->AddEntry(hmcr,"data/MC particle-level","pl");
	leg2->DrawClone("same");


	if(0) {
		TFile *fout = new TFile(Form("testBayes%d_W.root",iter),"RECREATE");
		//if(Wmodtruth) Wmodtruth->Write();
		reco->Write();

		hrcr->Write();
		UnfoldedLeadJetPt_data_save->Write();
		
		#ifdef WIDEBIN
		Wcov_save->Write();
		Wmeasured_save->Write();
		Wtruth_save->Write();
		WLeadJetPt_data->Write();
		#else
		cov_save->Write();
		measured_save->Write();
		truth_save->Write();
		LeadJetPt_data_save->Write();
		#endif

		fout->Close();
	}


}



#ifndef __CINT__
int main () { UnfoldJetpT(); return 0; }  // Main program when run stand-alone
#endif
