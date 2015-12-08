//==================================================================================================================================
//
//		2015.11.11	Li Yi
//		Read event info from ResultTree and Draw the Profile for Ntrk, Sum pT, Track pT
//		
//		This is to solve the issue that some runs have problematic behavior which ruins the original histograms  
//
//
//==================================================================================================================================


#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TF1.h"

#include <iostream>

using namespace std;

float inveff_pion(float pt) {
	return 1;			// test

	if(pt<0.1) return 0;

	float eff = 0;
	TF1* feff=new TF1("feff","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);
	feff->SetParameters(0.874739, 0.156624, 5.67316); 
	eff = feff->Eval(pt);

	if(pt>4.5) eff = feff->Eval(4.5);

	delete feff;		// prevent memory leak

	if(eff>0) return 1./eff;
	else return 0;
}

void plotTree2Histo(TString what2fill="multiplicity") {
// what2fill: refmult, leadjetpt, multiplicity	(no space)

	int savefig = 0;
	int saveroot = 0; 

	// reader
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root");
	TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_MB_R06_LeadJetAngle_MatchTrig_151206.root");
	if(!f) { cout<<"Cannot find input file"<<endl; return; }

	TTree *t = (TTree*)f->Get("ResultTree");
	if(!t) { cout<<"Cannot find Tree"<<endl; return; }

	float jpt, leadpt, subpt, tranmaxpt,tranminpt, tranpt;
	int leadntrk, subntrk, tranmaxntrk, tranminntrk, tranntrk;
	double refmult;
	int runid;

	const int MAXARRAY = 1000;
	float pt_min[MAXARRAY], pt_max[MAXARRAY], pt_jet[MAXARRAY], pt_sub[MAXARRAY];

	t->SetBranchAddress("runid",&runid);
	t->SetBranchAddress("refmult",&refmult);
	t->SetBranchAddress("j1pt",&jpt);
	t->SetBranchAddress("LeadAreaPtSum",&leadpt);
        t->SetBranchAddress("SubLeadAreaPtSum",&subpt);
        t->SetBranchAddress("TranPtSum",&tranpt);
        t->SetBranchAddress("TranMaxPtSum",&tranmaxpt);
        t->SetBranchAddress("TranMinPtSum",&tranminpt);
        t->SetBranchAddress("LeadAreaNtrk",&leadntrk);
        t->SetBranchAddress("SubAreaNtrk",&subntrk);
        t->SetBranchAddress("TranNtrk",&tranntrk);
        t->SetBranchAddress("TranMaxNtrk",&tranmaxntrk);
        t->SetBranchAddress("TranMinNtrk",&tranminntrk);
        t->SetBranchAddress("TrkLeadAreaPt",pt_jet);
        t->SetBranchAddress("TrkSubAreaPt",pt_sub);
        t->SetBranchAddress("TrkTranMaxPt",pt_max);
        t->SetBranchAddress("TrkTranMinPt",pt_min);

	
	// define histograms
	double maxpt = 50; // for refmult range
	int nbinning = 50; // for refmult range
	if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
		maxpt = 100;
		nbinning = 100;
	}
	TString xvariablename = "refmult";
	if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
		xvariablename = "Leading Jet Pt";
	}
	if(what2fill.Contains("multiplicity",TString::kIgnoreCase)) {
		xvariablename = "Total multiplicity";
	}

	TH1D *leadjetpt = new TH1D("leadjetpt",xvariablename,nbinning,0,maxpt);

        TProfile *leadjetntrkvsleadjetpt = new TProfile("leadjetareantrkvsleadjetpt","Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetntrkvsleadjetpt = new TProfile("subjetareantrkvsleadjetpt","SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxntrkvsleadjetpt = new TProfile("tranmaxntrkvsleadjetpt","Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminntrkvsleadjetpt = new TProfile("tranminntrkvsleadjetpt","Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranntrkvsleadjetpt = new TProfile("tranntrkvsleadjetpt","Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptsumvsleadjetpt = new TProfile("leadjetareaptsumvsleadjetpt","Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptsumvsleadjetpt = new TProfile("subjetareaptsumvsleadjetpt","SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptsumvsleadjetpt = new TProfile("tranmaxptsumvsleadjetpt","Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptsumvsleadjetpt = new TProfile("tranminptsumvsleadjetpt","Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptsumvsleadjetpt = new TProfile("tranptsumvsleadjetpt","Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptavevsleadjetpt = new TProfile("leadjetareaptavevsleadjetpt","Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptavevsleadjetpt = new TProfile("subjetareaptavevsleadjetpt","SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptavevsleadjetpt = new TProfile("tranmaxptavevsleadjetpt","Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptavevsleadjetpt = new TProfile("tranminptavevsleadjetpt","Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptavevsleadjetpt = new TProfile("tranptavevsleadjetpt","Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt);

	
	int nbinning_tpt = 1000;
	double max_tpt = 50;
        TH2D *hleadjetntrkvsleadjetpt = new TH2D("hleadjetareantrkvsleadjetpt","Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetntrkvsleadjetpt = new TH2D("hsubjetareantrkvsleadjetpt","SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxntrkvsleadjetpt = new TH2D("htranmaxntrkvsleadjetpt","Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminntrkvsleadjetpt = new TH2D("htranminntrkvsleadjetpt","Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranntrkvsleadjetpt = new TH2D("htranntrkvsleadjetpt","Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptsumvsleadjetpt = new TH2D("hleadjetareaptsumvsleadjetpt","Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptsumvsleadjetpt = new TH2D("hsubjetareaptsumvsleadjetpt","SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptsumvsleadjetpt = new TH2D("htranmaxptsumvsleadjetpt","Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptsumvsleadjetpt = new TH2D("htranminptsumvsleadjetpt","Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptsumvsleadjetpt = new TH2D("htranptsumvsleadjetpt","Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptavevsleadjetpt = new TH2D("hleadjetareaptavevsleadjetpt","Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptavevsleadjetpt = new TH2D("hsubjetareaptavevsleadjetpt","SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptavevsleadjetpt = new TH2D("htranmaxptavevsleadjetpt","Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptavevsleadjetpt = new TH2D("htranminptavevsleadjetpt","Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptavevsleadjetpt = new TH2D("htranptavevsleadjetpt","Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

	// loop over events
	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		t->GetEntry(ievt);

		if(runid>13052000&& runid<13060000) continue;		// problematic runs, need future investigation
		if(((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase)))&&(jpt<10)) continue;		// when studying multiplicity dependence, exclude jet pT<10GeV/c to ensure the real jet found

		if(ievt%100000==0) cout<<"event "<<ievt<<endl;

		double xvariable = 0;
		if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
			xvariable = jpt;
		}
		else if(what2fill.Contains("multiplicity",TString::kIgnoreCase)) {	// total mulitplicity
			xvariable = leadntrk+subntrk+tranntrk;
		}
		else {
			xvariable = refmult;
		}

		leadjetpt->Fill(xvariable);	
		
        	leadjetntrkvsleadjetpt->Fill(xvariable,leadntrk);
        	subjetntrkvsleadjetpt->Fill(xvariable,subntrk);
        	tranmaxntrkvsleadjetpt->Fill(xvariable,tranmaxntrk);
        	tranminntrkvsleadjetpt->Fill(xvariable,tranminntrk);
        	tranntrkvsleadjetpt->Fill(xvariable,tranntrk);

        	hleadjetntrkvsleadjetpt->Fill(xvariable,leadntrk);
        	hsubjetntrkvsleadjetpt->Fill(xvariable,subntrk);
        	htranmaxntrkvsleadjetpt->Fill(xvariable,tranmaxntrk);
        	htranminntrkvsleadjetpt->Fill(xvariable,tranminntrk);
        	htranntrkvsleadjetpt->Fill(xvariable,tranntrk);

        	//leadjetptsumvsleadjetpt->Fill(xvariable,leadpt);
        	//subjetptsumvsleadjetpt->Fill(xvariable,subpt);
        	//tranmaxptsumvsleadjetpt->Fill(xvariable,tranmaxpt);
        	//tranminptsumvsleadjetpt->Fill(xvariable,tranminpt);
        	//tranptsumvsleadjetpt->Fill(xvariable,tranpt);

		float sumleadpt = 0, sumsubpt = 0, sumtranmaxpt = 0, sumtranminpt = 0, sumtranpt = 0;
		//cout<<"tranmaxntrk = "<<tranmaxntrk<<"\ttranminntrk = "<<tranminntrk<<"\tleadntrk = "<<leadntrk<<"\tsubntrk = "<<subntrk<<endl;
		float w;
		for(int it = 0; it<tranmaxntrk; it++) {
			//cout<<"pt_max["<<it<<"] = "<<pt_max[it];
			w = inveff_pion(pt_max[it]);
			//cout<<" w = "<<w<<endl;
			tranmaxptavevsleadjetpt->Fill(xvariable,pt_max[it],w);
			tranptavevsleadjetpt->Fill(xvariable,pt_max[it],w);

			htranmaxptavevsleadjetpt->Fill(xvariable,pt_max[it],w);
			htranptavevsleadjetpt->Fill(xvariable,pt_max[it],w);

			sumtranmaxpt+=pt_max[it]*w;		// will be w tracks effectively ..
			sumtranpt+=pt_max[it]*w;		// will be w tracks effectively ..
		}
		for(int it = 0; it<tranminntrk; it++) {
			//cout<<"pt_min["<<it<<"] = "<<pt_min[it];
			w = inveff_pion(pt_min[it]);
			//cout<<" w = "<<w<<endl;
			tranminptavevsleadjetpt->Fill(xvariable,pt_min[it],w);
			tranptavevsleadjetpt->Fill(xvariable,pt_min[it],w);

			htranminptavevsleadjetpt->Fill(xvariable,pt_min[it],w);
			htranptavevsleadjetpt->Fill(xvariable,pt_min[it],w);

			sumtranminpt+=pt_min[it]*w;		// will be w tracks effectively ..
			sumtranpt+=pt_min[it]*w;		// will be w tracks effectively ..
		}
		for(int it = 0; it<leadntrk; it++) {
			//cout<<"pt_jet["<<it<<"] = "<<pt_jet[it];
			w = inveff_pion(pt_jet[it]);
			//cout<<" w = "<<w<<endl;
			leadjetptavevsleadjetpt->Fill(xvariable,pt_jet[it],w);

			hleadjetptavevsleadjetpt->Fill(xvariable,pt_jet[it],w);

			sumleadpt+=pt_jet[it]*w;		// will be w tracks effectively ..
		}
		for(int it = 0; it<subntrk; it++) {
			//cout<<"pt_sub["<<it<<"] = "<<pt_sub[it];
			w = inveff_pion(pt_sub[it]);
			//cout<<" w = "<<w<<endl;
			subjetptavevsleadjetpt->Fill(xvariable,pt_sub[it],w);

			hsubjetptavevsleadjetpt->Fill(xvariable,pt_sub[it],w);

			sumsubpt+=pt_sub[it]*w;		// will be w tracks effectively ..
		}
        	leadjetptsumvsleadjetpt->Fill(xvariable,sumleadpt);
        	subjetptsumvsleadjetpt->Fill(xvariable,sumsubpt);
        	tranmaxptsumvsleadjetpt->Fill(xvariable,sumtranmaxpt);
        	tranminptsumvsleadjetpt->Fill(xvariable,sumtranminpt);
        	tranptsumvsleadjetpt->Fill(xvariable,sumtranpt);

        	hleadjetptsumvsleadjetpt->Fill(xvariable,sumleadpt);
        	hsubjetptsumvsleadjetpt->Fill(xvariable,sumsubpt);
        	htranmaxptsumvsleadjetpt->Fill(xvariable,sumtranmaxpt);
        	htranminptsumvsleadjetpt->Fill(xvariable,sumtranminpt);
        	htranptsumvsleadjetpt->Fill(xvariable,sumtranpt);

        	//leadjetptavevsleadjetpt->Fill(xvariable,((leadntrk>0)?leadpt/leadntrk:0),leadntrk);
        	//subjetptavevsleadjetpt->Fill(xvariable,((subntrk>0)?subpt/subntrk:0),subntrk);
        	//tranmaxptavevsleadjetpt->Fill(xvariable,((tranmaxntrk>0)?tranmaxpt/tranmaxntrk:0),tranmaxntrk);
        	//tranminptavevsleadjetpt->Fill(xvariable,((tranminntrk>0)?tranminpt/tranminntrk:0),tranminntrk);
        	//tranptavevsleadjetpt->Fill(xvariable,((tranntrk>0)?tranpt/tranntrk:0),tranntrk);

	}

	// set histogram draw properties
	leadjetpt->GetXaxis()->SetTitle(xvariablename);
	leadjetpt->GetYaxis()->SetTitle("events");
	leadjetpt->SetLineColor(1);
	//leadjetpt->SetMarkerStyle(8);
	//leadjetpt->SetMarkerColor(1);
	
	int clead = 1, csub = 9, cmax = 2, cmin = 8, ctran = 28;
	int slead = 20, ssub = 25, smax = 24, smin = 20, stran = 21;

        leadjetntrkvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	leadjetntrkvsleadjetpt->GetYaxis()->SetTitle("Leading Jet Area Multiplicity");
	leadjetntrkvsleadjetpt->SetLineColor(clead);
	leadjetntrkvsleadjetpt->SetMarkerStyle(slead);
	leadjetntrkvsleadjetpt->SetMarkerColor(clead);

        subjetntrkvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	subjetntrkvsleadjetpt->GetYaxis()->SetTitle("Away Side Area Multiplicity");
	subjetntrkvsleadjetpt->SetLineColor(csub);
	subjetntrkvsleadjetpt->SetMarkerStyle(ssub);
	subjetntrkvsleadjetpt->SetMarkerColor(csub);

        tranmaxntrkvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranmaxntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Multiplicity");
	tranmaxntrkvsleadjetpt->SetLineColor(cmax);
	tranmaxntrkvsleadjetpt->SetMarkerStyle(smax);
	tranmaxntrkvsleadjetpt->SetMarkerColor(cmax);

        tranminntrkvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranminntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Multiplicity");
	tranminntrkvsleadjetpt->SetLineColor(cmin);
	tranminntrkvsleadjetpt->SetMarkerStyle(smin);
	tranminntrkvsleadjetpt->SetMarkerColor(cmin);

        tranntrkvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Area Multiplicity");
	tranntrkvsleadjetpt->SetLineColor(ctran);
	tranntrkvsleadjetpt->SetMarkerStyle(stran);
	tranntrkvsleadjetpt->SetMarkerColor(ctran);


        leadjetptsumvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	leadjetptsumvsleadjetpt->GetYaxis()->SetTitle("Lead Jet Area Sum Pt");
	leadjetptsumvsleadjetpt->SetLineColor(clead);
	leadjetptsumvsleadjetpt->SetMarkerStyle(slead);
	leadjetptsumvsleadjetpt->SetMarkerColor(clead);

        subjetptsumvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	subjetptsumvsleadjetpt->GetYaxis()->SetTitle("Away Side Are Sum Pt");
	subjetptsumvsleadjetpt->SetLineColor(csub);
	subjetptsumvsleadjetpt->SetMarkerStyle(ssub);
	subjetptsumvsleadjetpt->SetMarkerColor(csub);

        tranmaxptsumvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranmaxptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Sum Pt");
	tranmaxptsumvsleadjetpt->SetLineColor(cmax);
	tranmaxptsumvsleadjetpt->SetMarkerStyle(smax);
	tranmaxptsumvsleadjetpt->SetMarkerColor(cmax);

        tranminptsumvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranminptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Sum Pt");
	tranminptsumvsleadjetpt->SetLineColor(cmin);
	tranminptsumvsleadjetpt->SetMarkerStyle(smin);
	tranminptsumvsleadjetpt->SetMarkerColor(cmin);

        tranptsumvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Area Sum Pt");
	tranptsumvsleadjetpt->SetLineColor(ctran);
	tranptsumvsleadjetpt->SetMarkerStyle(stran);
	tranptsumvsleadjetpt->SetMarkerColor(ctran);


        leadjetptavevsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	leadjetptavevsleadjetpt->GetYaxis()->SetTitle("Leading Area Average Track Pt");
	leadjetptavevsleadjetpt->SetLineColor(clead);
	leadjetptavevsleadjetpt->SetMarkerStyle(slead);
	leadjetptavevsleadjetpt->SetMarkerColor(clead);

        subjetptavevsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	subjetptavevsleadjetpt->GetYaxis()->SetTitle("Away side Area Average Track Pt");
	subjetptavevsleadjetpt->SetLineColor(csub);
	subjetptavevsleadjetpt->SetMarkerStyle(ssub);
	subjetptavevsleadjetpt->SetMarkerColor(csub);

        tranmaxptavevsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranmaxptavevsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Average Track Pt");
	tranmaxptavevsleadjetpt->SetLineColor(cmax);
	tranmaxptavevsleadjetpt->SetMarkerStyle(smax);
	tranmaxptavevsleadjetpt->SetMarkerColor(cmax);

        tranminptavevsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranminptavevsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Average Track Pt");
	tranminptavevsleadjetpt->SetLineColor(cmin);
	tranminptavevsleadjetpt->SetMarkerStyle(smin);
	tranminptavevsleadjetpt->SetMarkerColor(cmin);

        tranptavevsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	tranptavevsleadjetpt->GetYaxis()->SetTitle("Transverse Area Average Track Pt");
	tranptavevsleadjetpt->SetLineColor(ctran);
	tranptavevsleadjetpt->SetMarkerStyle(stran);
	tranptavevsleadjetpt->SetMarkerColor(ctran);


	// Drawing
	const int ncanv = 4;
	TCanvas *c[ncanv];
	TH2D *htmp[ncanv];
	gStyle->SetOptStat(0);

	for(int i = 0; i<ncanv; i++) {
 		c[i] = new TCanvas();
	}
	
	c[0]->cd();
	c[0]->SetLogy(1);
	leadjetpt->Draw();
	if(savefig) c[0]->SaveAs(Form("%sVs%s.png",leadjetpt->GetName(),what2fill.Data()));
	
	TLegend *leg = new TLegend(0.16,0.6,0.35,0.86);
	leg->AddEntry(leadjetntrkvsleadjetpt,"Toward","pl");
	leg->AddEntry(subjetntrkvsleadjetpt,"Away","pl");
	leg->AddEntry(tranmaxntrkvsleadjetpt,"TransMax","pl");
	leg->AddEntry(tranminntrkvsleadjetpt,"TransMin","pl");
	leg->AddEntry(tranntrkvsleadjetpt,"Trans","pl");
	leg->SetFillColor(10);
	leg->SetBorderSize(0);
	leg->SetTextFont(42);

	TLatex *lat = new TLatex();
	lat->SetNDC(1);
	lat->SetTextFont(42);
	lat->SetTextAlign(21);	// center alignment
	
	c[1]->cd();
	float drawxmax = 50;
	htmp[1] = new TH2D("htmp1","",1000,0,drawxmax,1000,0,10);//5);
	htmp[1]->GetXaxis()->SetTitle(xvariablename);
	htmp[1]->GetYaxis()->SetTitle("Event Multiplicity #LTNtrk#GT");
	htmp[1]->Draw();
	lat->DrawLatex(0.5,0.94,"Event Multiplicity #LTNtrk#GT");
	leadjetntrkvsleadjetpt->Draw("same");
	subjetntrkvsleadjetpt->Draw("same");
	tranmaxntrkvsleadjetpt->Draw("same");
	tranminntrkvsleadjetpt->Draw("same");	
	tranntrkvsleadjetpt->Draw("same");	
	leg->Draw("same");
	if(savefig) c[1]->SaveAs(Form("%sVs%s.png","Ntrk",what2fill.Data()));

	c[2]->cd();
	htmp[2] = new TH2D("htmp2","",1000,0,drawxmax,1000,0,20);
	htmp[2]->GetXaxis()->SetTitle(xvariablename);
	htmp[2]->GetYaxis()->SetTitle("Event #LT#sum p_{T}#GT");
	htmp[2]->Draw();
	lat->DrawLatex(0.5,0.94,"Event #LT#sum p_{T}#GT");
	leadjetptsumvsleadjetpt->Draw("same");
	subjetptsumvsleadjetpt->Draw("same");
	tranmaxptsumvsleadjetpt->Draw("same");
	tranminptsumvsleadjetpt->Draw("same");	
	tranptsumvsleadjetpt->Draw("same");	
	leg->Draw("same");
	if(savefig) c[2]->SaveAs(Form("%sVs%s.png","PtSum",what2fill.Data()));
	
	c[3]->cd();
	htmp[3] = new TH2D("htmp3","",1000,0,drawxmax,1000,0,7);
	htmp[3]->GetXaxis()->SetTitle(xvariablename);
	htmp[3]->GetYaxis()->SetTitle("Track #LTp_{T}#GT");
	htmp[3]->Draw();
	lat->DrawLatex(0.5,0.94,"Track #LTp_{T}#GT");
	leadjetptavevsleadjetpt->Draw("same");
	subjetptavevsleadjetpt->Draw("same");
	tranmaxptavevsleadjetpt->Draw("same");
	tranminptavevsleadjetpt->Draw("same");	
	tranptavevsleadjetpt->Draw("same");	
	leg->Draw("same");
	if(savefig) c[3]->SaveAs(Form("%sVs%s.png","PtAve",what2fill.Data()));


	if(saveroot) {
		TFile *fout = new TFile("TwoHisto4underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root","RECREATE");
		leadjetpt->Write();
		
        	leadjetntrkvsleadjetpt->Write();
        	subjetntrkvsleadjetpt->Write();
        	tranmaxntrkvsleadjetpt->Write();
        	tranminntrkvsleadjetpt->Write();
        	tranntrkvsleadjetpt->Write();

        	leadjetptsumvsleadjetpt->Write();
        	subjetptsumvsleadjetpt->Write();
        	tranmaxptsumvsleadjetpt->Write();
        	tranminptsumvsleadjetpt->Write();
        	tranptsumvsleadjetpt->Write();

        	leadjetptavevsleadjetpt->Write();
        	subjetptavevsleadjetpt->Write();
        	tranmaxptavevsleadjetpt->Write();
        	tranminptavevsleadjetpt->Write();
        	tranptavevsleadjetpt->Write();

	
        	hleadjetntrkvsleadjetpt->Write();
        	hsubjetntrkvsleadjetpt->Write();
        	htranmaxntrkvsleadjetpt->Write();
        	htranminntrkvsleadjetpt->Write();
        	htranntrkvsleadjetpt->Write();

        	hleadjetptsumvsleadjetpt->Write();
        	hsubjetptsumvsleadjetpt->Write();
        	htranmaxptsumvsleadjetpt->Write();
        	htranminptsumvsleadjetpt->Write();
        	htranptsumvsleadjetpt->Write();

        	hleadjetptavevsleadjetpt->Write();
        	hsubjetptavevsleadjetpt->Write();
        	htranmaxptavevsleadjetpt->Write();
        	htranminptavevsleadjetpt->Write();
        	htranptavevsleadjetpt->Write();
			
		fout->Close();
	}

}


