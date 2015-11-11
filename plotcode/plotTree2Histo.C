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

#include <iostream>

using namespace std;

void plotTree2Histo() {

	// reader
	TFile *f = new TFile("underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root");
	if(!f) { cout<<"Cannot find input file"<<endl; return; }

	TTree *t = (TTree*)f->Get("ResultTree");
	if(!t) { cout<<"Cannot find Tree"<<endl; return; }

	float jpt, leadpt, subpt, tranmaxpt,tranminpt, tranpt;
	int leadntrk, subntrk, tranmaxntrk, tranminntrk, tranntrk;
	int runid;

	t->SetBranchAddress("runid",&runid);
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

	
	// define histograms
	TH1D *leadjetpt = new TH1D("leadjetpt","Leading Jet Pt",100,0,100);

        TProfile *leadjetntrkvsleadjetpt = new TProfile("leadjetareantrkvsleadjetpt","Leading Jet Area Ntrk vs Leading Jet Pt",100,0,100);
        TProfile *subjetntrkvsleadjetpt = new TProfile("subjetareantrkvsleadjetpt","SubLeading Jet Area Ntrk vs Leading Jet Pt",100,0,100);
        TProfile *tranmaxntrkvsleadjetpt = new TProfile("tranmaxntrkvsleadjetpt","Transverse Max Ntrk vs Leading Jet Pt",100,0,100);
        TProfile *tranminntrkvsleadjetpt = new TProfile("tranminntrkvsleadjetpt","Transverse Min Ntrk vs Leading Jet Pt",100,0,100);
        TProfile *tranntrkvsleadjetpt = new TProfile("tranntrkvsleadjetpt","Transverse Ntrk vs Leading Jet Pt",100,0,100);

        TProfile *leadjetptsumvsleadjetpt = new TProfile("leadjetareaptsumvsleadjetpt","Leading Jet Area Sum Pt vs Leading Jet Pt",100,0,100);
        TProfile *subjetptsumvsleadjetpt = new TProfile("subjetareaptsumvsleadjetpt","SubLeading Jet Area Sum Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranmaxptsumvsleadjetpt = new TProfile("tranmaxptsumvsleadjetpt","Transverse Max Sum Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranminptsumvsleadjetpt = new TProfile("tranminptsumvsleadjetpt","Transverse Min Sum Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranptsumvsleadjetpt = new TProfile("tranptsumvsleadjetpt","Transverse Sum Pt vs Leading Jet Pt",100,0,100);

        TProfile *leadjetptavevsleadjetpt = new TProfile("leadjetareaptavevsleadjetpt","Leading Jet Area Average Pt vs Leading Jet Pt",100,0,100);
        TProfile *subjetptavevsleadjetpt = new TProfile("subjetareaptavevsleadjetpt","SubLeading Jet Area Average Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranmaxptavevsleadjetpt = new TProfile("tranmaxptavevsleadjetpt","Transverse Max Average Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranminptavevsleadjetpt = new TProfile("tranminptavevsleadjetpt","Transverse Min Average Pt vs Leading Jet Pt",100,0,100);
        TProfile *tranptavevsleadjetpt = new TProfile("tranptavevsleadjetpt","Transverse Average Pt vs Leading Jet Pt",100,0,100);

	// loop over events
	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		if(runid>13052000&& runid<13060000) continue;		// problematic runs, need future investigation
		if(ievt%100000==0) cout<<"event "<<ievt<<endl;

		t->GetEntry(ievt);

		leadjetpt->Fill(jpt);	
		
        	leadjetntrkvsleadjetpt->Fill(jpt,leadntrk);
        	subjetntrkvsleadjetpt->Fill(jpt,subntrk);
        	tranmaxntrkvsleadjetpt->Fill(jpt,tranmaxntrk);
        	tranminntrkvsleadjetpt->Fill(jpt,tranminntrk);
        	tranntrkvsleadjetpt->Fill(jpt,tranntrk);

        	leadjetptsumvsleadjetpt->Fill(jpt,leadpt);
        	subjetptsumvsleadjetpt->Fill(jpt,subpt);
        	tranmaxptsumvsleadjetpt->Fill(jpt,tranmaxpt);
        	tranminptsumvsleadjetpt->Fill(jpt,tranminpt);
        	tranptsumvsleadjetpt->Fill(jpt,tranpt);

        	leadjetptavevsleadjetpt->Fill(jpt,((leadntrk>0)?leadpt/leadntrk:0),leadntrk);
        	subjetptavevsleadjetpt->Fill(jpt,((subntrk>0)?subpt/subntrk:0),subntrk);
        	tranmaxptavevsleadjetpt->Fill(jpt,((tranmaxntrk>0)?tranmaxpt/tranmaxntrk:0),tranmaxntrk);
        	tranminptavevsleadjetpt->Fill(jpt,((tranminntrk>0)?tranminpt/tranminntrk:0),tranminntrk);
        	tranptavevsleadjetpt->Fill(jpt,((tranntrk>0)?tranpt/tranntrk:0),tranntrk);

	}

	// set histogram draw properties
	leadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	leadjetpt->GetYaxis()->SetTitle("events");
	leadjetpt->SetLineColor(1);
	//leadjetpt->SetMarkerStyle(8);
	//leadjetpt->SetMarkerColor(1);
	
	int clead = 1, csub = 9, cmax = 2, cmin = 8, ctran = 28;
	int slead = 20, ssub = 25, smax = 24, smin = 20, stran = 21;

        leadjetntrkvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	leadjetntrkvsleadjetpt->GetYaxis()->SetTitle("Leading Jet Area Multiplicity");
	leadjetntrkvsleadjetpt->SetLineColor(clead);
	leadjetntrkvsleadjetpt->SetMarkerStyle(slead);
	leadjetntrkvsleadjetpt->SetMarkerColor(clead);

        subjetntrkvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	subjetntrkvsleadjetpt->GetYaxis()->SetTitle("Away Side Area Multiplicity");
	subjetntrkvsleadjetpt->SetLineColor(csub);
	subjetntrkvsleadjetpt->SetMarkerStyle(ssub);
	subjetntrkvsleadjetpt->SetMarkerColor(csub);

        tranmaxntrkvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranmaxntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Multiplicity");
	tranmaxntrkvsleadjetpt->SetLineColor(cmax);
	tranmaxntrkvsleadjetpt->SetMarkerStyle(smax);
	tranmaxntrkvsleadjetpt->SetMarkerColor(cmax);

        tranminntrkvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranminntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Multiplicity");
	tranminntrkvsleadjetpt->SetLineColor(cmin);
	tranminntrkvsleadjetpt->SetMarkerStyle(smin);
	tranminntrkvsleadjetpt->SetMarkerColor(cmin);

        tranntrkvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranntrkvsleadjetpt->GetYaxis()->SetTitle("Transverse Area Multiplicity");
	tranntrkvsleadjetpt->SetLineColor(ctran);
	tranntrkvsleadjetpt->SetMarkerStyle(stran);
	tranntrkvsleadjetpt->SetMarkerColor(ctran);


        leadjetptsumvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	leadjetptsumvsleadjetpt->GetYaxis()->SetTitle("Lead Jet Area Sum Pt");
	leadjetptsumvsleadjetpt->SetLineColor(clead);
	leadjetptsumvsleadjetpt->SetMarkerStyle(slead);
	leadjetptsumvsleadjetpt->SetMarkerColor(clead);

        subjetptsumvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	subjetptsumvsleadjetpt->GetYaxis()->SetTitle("Away Side Are Sum Pt");
	subjetptsumvsleadjetpt->SetLineColor(csub);
	subjetptsumvsleadjetpt->SetMarkerStyle(ssub);
	subjetptsumvsleadjetpt->SetMarkerColor(csub);

        tranmaxptsumvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranmaxptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Sum Pt");
	tranmaxptsumvsleadjetpt->SetLineColor(cmax);
	tranmaxptsumvsleadjetpt->SetMarkerStyle(smax);
	tranmaxptsumvsleadjetpt->SetMarkerColor(cmax);

        tranminptsumvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranminptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Sum Pt");
	tranminptsumvsleadjetpt->SetLineColor(cmin);
	tranminptsumvsleadjetpt->SetMarkerStyle(smin);
	tranminptsumvsleadjetpt->SetMarkerColor(cmin);

        tranptsumvsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranptsumvsleadjetpt->GetYaxis()->SetTitle("Transverse Area Sum Pt");
	tranptsumvsleadjetpt->SetLineColor(ctran);
	tranptsumvsleadjetpt->SetMarkerStyle(stran);
	tranptsumvsleadjetpt->SetMarkerColor(ctran);


        leadjetptavevsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	leadjetptavevsleadjetpt->GetYaxis()->SetTitle("Leading Area Average Track Pt");
	leadjetptavevsleadjetpt->SetLineColor(clead);
	leadjetptavevsleadjetpt->SetMarkerStyle(slead);
	leadjetptavevsleadjetpt->SetMarkerColor(clead);

        subjetptavevsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	subjetptavevsleadjetpt->GetYaxis()->SetTitle("Away side Area Average Track Pt");
	subjetptavevsleadjetpt->SetLineColor(csub);
	subjetptavevsleadjetpt->SetMarkerStyle(ssub);
	subjetptavevsleadjetpt->SetMarkerColor(csub);

        tranmaxptavevsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranmaxptavevsleadjetpt->GetYaxis()->SetTitle("Transverse Max Area Average Track Pt");
	tranmaxptavevsleadjetpt->SetLineColor(cmax);
	tranmaxptavevsleadjetpt->SetMarkerStyle(smax);
	tranmaxptavevsleadjetpt->SetMarkerColor(cmax);

        tranminptavevsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
	tranminptavevsleadjetpt->GetYaxis()->SetTitle("Transverse Min Area Average Track Pt");
	tranminptavevsleadjetpt->SetLineColor(cmin);
	tranminptavevsleadjetpt->SetMarkerStyle(smin);
	tranminptavevsleadjetpt->SetMarkerColor(cmin);

        tranptavevsleadjetpt->GetXaxis()->SetTitle("Leading Jet Pt");
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
	c[0]->SaveAs(Form("%s.png",leadjetpt->GetName()));
	
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
	htmp[1] = new TH2D("htmp1","",1000,0,50,1000,0,5);
	htmp[1]->GetXaxis()->SetTitle("Leading Jet p_{T} (GeV/c)");
	htmp[1]->GetYaxis()->SetTitle("Event Multiplicity #LTNtrk#GT");
	htmp[1]->Draw();
	lat->DrawLatex(0.5,0.94,"Event Multiplicity #LTNtrk#GT");
	leadjetntrkvsleadjetpt->Draw("same");
	subjetntrkvsleadjetpt->Draw("same");
	tranmaxntrkvsleadjetpt->Draw("same");
	tranminntrkvsleadjetpt->Draw("same");	
	tranntrkvsleadjetpt->Draw("same");	
	leg->Draw("same");
	c[1]->SaveAs(Form("%s.png","Ntrk"));

	c[2]->cd();
	htmp[2] = new TH2D("htmp2","",1000,0,50,1000,0,20);
	htmp[2]->GetXaxis()->SetTitle("Leading Jet p_{T} (GeV/c)");
	htmp[2]->GetYaxis()->SetTitle("Event #LT#sum p_{T}#GT");
	htmp[2]->Draw();
	lat->DrawLatex(0.5,0.94,"Event #LT#sum p_{T}#GT");
	leadjetptsumvsleadjetpt->Draw("same");
	subjetptsumvsleadjetpt->Draw("same");
	tranmaxptsumvsleadjetpt->Draw("same");
	tranminptsumvsleadjetpt->Draw("same");	
	tranptsumvsleadjetpt->Draw("same");	
	leg->Draw("same");
	c[2]->SaveAs(Form("%s.png","PtSum"));
	
	c[3]->cd();
	htmp[3] = new TH2D("htmp3","",1000,0,50,1000,0,7);
	htmp[3]->GetXaxis()->SetTitle("Leading Jet p_{T} (GeV/c)");
	htmp[3]->GetYaxis()->SetTitle("Track #LTp_{T}#GT");
	htmp[3]->Draw();
	lat->DrawLatex(0.5,0.94,"Track #LTp_{T}#GT");
	leadjetptavevsleadjetpt->Draw("same");
	subjetptavevsleadjetpt->Draw("same");
	tranmaxptavevsleadjetpt->Draw("same");
	tranminptavevsleadjetpt->Draw("same");	
	tranptavevsleadjetpt->Draw("same");	
	leg->Draw("same");
	c[3]->SaveAs(Form("%s.png","PtAve"));

}


