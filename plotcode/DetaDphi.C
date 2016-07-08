//==================================================================================================================================
//
//
//		2016.03.01	Li Yi
//		plot Delta eta - Delta phi
//		for jet and underlying event in the same eta region or not
//		
///==================================================================================================================================


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
#include "TMath.h"

#include <iostream>
#include <algorithm>    // std::max, std::max_element

#define MINPTCUT	0.2				// CHECK!!!! need to cut consistence with production code in: 
//UnderlyingAna::UnderlyingAna ( double R,
//                double max_const_rap, //double PtConsLo, double PtConsHi,
//                double min_const_pt,				<------------ This one
//                double dPhiCut,
//                TString name
//                )

using namespace std;

float foldphi(float phi) {

	int npi = floor(fabs(phi)/(2*TMath::Pi()));

	if(phi>0) phi = phi-npi*2*TMath::Pi();
	if(phi<0) phi = phi+(npi+1)*2*TMath::Pi();

	if(phi>TMath::Pi()) phi = 2*TMath::Pi()-phi;

	return phi;
}


void DetaDphi(TString dir="~/Scratch/pp200Y12_jetunderlying/", TString filetag = "underlyingevent_MB_R06_LeadJetAngle_FullJetFraclt90_160116", double jetptmin = 10, double jetptmax= 200, int ExclusiveEta = 0, int DijetSameSide = 0, double subjetptmin = 5, double subjetptmax = 200, int DoDijet = 0) {
// ExclusiveEta: 
// 		==1: jet and underlying event in seperate eta region. for example: jet in [-0.6, 0), then undelrying event [0.6,1)
// 		==0: no requirement for underlying event and jet eta
// DijetSameSide:
// 		==1: Dijet in the same eta side
// 		==0: Dijet in the opposite eta side
// DoDijet:
// 		==1: Fill Dijet no matter "Dijet" in the file name or not (When "Dijet" is in the file name, we will do Dijet cuts with subjetptmin cuts)
// 		==0: Only do Dijet if there is "Dijet" in file name

	int MCflag=0;
	if(filetag.Contains("pythia",TString::kIgnoreCase)) 	{
		MCflag = 1;
		cout<<endl<<"INFO: process as MC data: not efficiency correction will be applied"<<endl;
	}


	int chargeflag = 1;
	if(filetag.Contains("Charge0",TString::kIgnoreCase)||filetag.Contains("TransNeutral",TString::kIgnoreCase)) chargeflag = 0;
	cout<<endl<<"INFO: underlying event chargeflag == "<<chargeflag<<endl;
	if(chargeflag==1&&MCflag==0) cout<<"Apply TPC tracking & TOF matching efficiency"<<endl;
	cout<<endl;

	int savefig = 0;
	int saveroot = 1; 

	// reader
	TString filepath = dir+filetag+".root";
	TFile *f = new TFile(filepath);
	if(!f) { cout<<"Cannot find input file"<<endl; return; }

	TTree *t = (TTree*)f->Get("ResultTree");
	if(!t) { cout<<"Cannot find Tree"<<endl; return; }

	float jpt, jeta, jphi, jaspt, jaseta, jasphi, leadpt, subpt, tranmaxpt,tranminpt, tranpt;
	int leadntrk, subntrk, tranmaxntrk, tranminntrk, tranntrk;
	double refmult;
	int runid;

	float j1neutralfrac;

	const int MAXARRAY = 1000;
	float pt_min[MAXARRAY], pt_max[MAXARRAY], pt_jet[MAXARRAY], pt_sub[MAXARRAY];
	float eta_min[MAXARRAY], eta_max[MAXARRAY], eta_jet[MAXARRAY], eta_sub[MAXARRAY];
	float phi_min[MAXARRAY], phi_max[MAXARRAY], phi_jet[MAXARRAY], phi_sub[MAXARRAY];

	t->SetBranchAddress("runid",&runid);
	t->SetBranchAddress("refmult",&refmult);
	t->SetBranchAddress("j1pt",&jpt);
	//t->SetBranchAddress("j1r1pt",&jpt);		// same jet with R=1
	t->SetBranchAddress("j1eta",&jeta);
	t->SetBranchAddress("j1phi",&jphi);
	t->SetBranchAddress("jaspt",&jaspt);
	t->SetBranchAddress("jasphi",&jasphi);
	t->SetBranchAddress("jaseta",&jaseta);
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
        t->SetBranchAddress("TrkLeadAreaEta",eta_jet);
        t->SetBranchAddress("TrkSubAreaEta",eta_sub);
        t->SetBranchAddress("TrkTranMaxEta",eta_max);
        t->SetBranchAddress("TrkTranMinEta",eta_min);
        t->SetBranchAddress("TrkLeadAreaPhi",phi_jet);
        t->SetBranchAddress("TrkSubAreaPhi",phi_sub);
        t->SetBranchAddress("TrkTranMaxPhi",phi_max);
        t->SetBranchAddress("TrkTranMinPhi",phi_min);

	t->SetBranchAddress("j1neutralfrac",&j1neutralfrac);
	
	// define histograms
	TH2D *hetaphi = new TH2D("hetaphi","#Delta#eta-#Delta#phi",100,0,TMath::Pi(),100,-2,2);
	hetaphi->Sumw2();
	TH1D *hjet = new TH1D("hjet","Number of leading jet",1,0,1);
	TH2D *hetaphi2 = new TH2D("hetaphi2","#Delta#eta-#Delta#phi for subleading jet * associated particle",100,0,TMath::Pi(),100,-2,2);
	hetaphi2->Sumw2();
	TH1D *hjet2 = new TH1D("hjet2","Number of subleading jet",1,0,1);
	
	// loop over events
	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
	int processedevent = 0;


	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		t->GetEntry(ievt);

		if(jpt>jetptmax || jpt<jetptmin) continue;

		//if(fabs(jasptmax)>1e-6 && fabs(fasptmin)>1e-6) {
		//	if(jaspt>jasptmax || jaspt<jasptmin) continue;
		//	if(DijetSameSide && jeta*jaseta<0) continue;
		//}

		if(j1neutralfrac>0.9) continue;			// jet neutral pT fraction < 90%

		//if(fabs(jeta)>0.4) continue;			// test for R=0.2
		
		if(filetag.Contains("Dijet")||DoDijet) {
			if(fabs(jaspt)<subjetptmin) continue;			
			//if(fabs(jpt-jaspt)>0.05*jpt) continue;		// dijet balanced
			//if(fabs(jpt-jaspt)>0.3*jpt) continue;		// dijet balanced
			//if(fabs(jpt-jaspt)<=0.3*jpt) continue;		// dijet unbalanced

			if(DijetSameSide && jeta*jaseta<0) continue;
		}

		hjet->Fill(0);		// this is count++ for jet
		if(jaspt>0) hjet2->Fill(0);

		if(ievt%1000000==0) cout<<"event "<<ievt<<endl;
		
		int etaflag = 0;
		if( ExclusiveEta>0 && jeta>0 )	etaflag = -1;
		if( ExclusiveEta>0 && jeta<0 )	etaflag = 1;

		for(int it = 0; it<tranmaxntrk; it++) {
			if(etaflag==1&&eta_max[it]<0.6) continue;
			if(etaflag==-1&&eta_max[it]>-0.6) continue;

			hetaphi->Fill(foldphi(phi_max[it]-jphi),eta_max[it]-jeta);
			if(jaspt>0) hetaphi2->Fill(foldphi(phi_max[it]-jasphi),eta_max[it]-jaseta);

		}
		for(int it = 0; it<tranminntrk; it++) {
			if(etaflag==1&&eta_min[it]<0.6) continue;
			if(etaflag==-1&&eta_min[it]>-0.6) continue;

			hetaphi->Fill(foldphi(phi_min[it]-jphi),eta_min[it]-jeta);
			if(jaspt>0) hetaphi2->Fill(foldphi(phi_min[it]-jasphi),eta_min[it]-jaseta);
		}
		for(int it = 0; it<leadntrk; it++) {
			if(etaflag==1&&eta_jet[it]<0.6) continue;
			if(etaflag==-1&&eta_jet[it]>-0.6) continue;

			hetaphi->Fill(foldphi(phi_jet[it]-jphi),eta_jet[it]-jeta);
			if(jaspt>0) hetaphi2->Fill(foldphi(phi_jet[it]-jasphi),eta_jet[it]-jaseta);
		}
		for(int it = 0; it<subntrk; it++) {
			if(etaflag==1&&eta_sub[it]<0.6) continue;
			if(etaflag==-1&&eta_sub[it]>-0.6) continue;

			hetaphi->Fill(foldphi(phi_sub[it]-jphi),eta_sub[it]-jeta);
			if(jaspt>0) hetaphi2->Fill(foldphi(phi_sub[it]-jasphi),eta_sub[it]-jaseta);
		}

      		processedevent++;
	}
	cout<<"Total processed # of event "<<processedevent<<endl;

	// set histogram draw properties
	int clead = 1, csub = 9, cmax = 2, cmin = 8, ctran = 28;
	int slead = 20, ssub = 25, smax = 24, smin = 20, stran = 21;

	hetaphi->GetXaxis()->SetTitle("#Delta#phi");
	hetaphi->GetYaxis()->SetTitle("#Delta#eta");
	hetaphi->SetLineColor(clead);
	hetaphi->SetMarkerStyle(slead);
	hetaphi->SetMarkerColor(clead);

      	if(saveroot) {
		TString jettag=Form("_jet%g-%g",jetptmin,jetptmax);	
		TString dijettag="";
		if(DoDijet&&(!filetag.Contains("Dijet"))) dijettag = "DoDijet";

		TString outfilepath = dir+"Deta"+"hist4"+dijettag+filetag+jettag+".root";
		if(ExclusiveEta>0) {
			outfilepath = dir+"Deta"+"hist4"+dijettag+filetag+jettag+"_EtaExcl"+".root";
			if(DijetSameSide) {
				outfilepath = dir+"Deta"+"hist4"+dijettag+filetag+jettag+"_EtaExcl_SamesideDijet"+".root";
			}
		}
		TFile *fout = new TFile(outfilepath,"RECREATE");

		hetaphi->Write();
		hjet->Write();
		hetaphi2->Write();
		hjet2->Write();
			
		fout->Close();
	}

}


void MaxOrMin(float &max, float &min) {		// switch max or min
	if(max>=min) return;
	float tmp = max;
	max = min;
	min = tmp;
}

void MaxOrMin(int &max, int &min) {		// switch max or min
	if(max>=min) return;
	int tmp = max;
	max = min;
	min = tmp;
}



float getweight(float pt, int charge=1, int MC=0) {		// tpc and tof efficiency only apply to charged particle. set 'charge' to 0 if bemc neutral particles are used.
	
	if(MC>0) return 1;					// pythia data, no correction needed

	if(charge==0) return 1;

	if(pt<MINPTCUT) return 0;				// set minimum pt as 0.2 GeV !!!!! CHECK

	float eff = 0;
	TF1* feff=new TF1("feff","[0]*(exp(-pow([1]/x,[2])-pow([3]/x,[4])))",0.1,4.5);
	feff->SetParameters(6.372e-01,1.52059e-01, 5.28031,0.156624, 5.67316);		// get from the below inveff_tof & inveff_tpc function parameters
	eff = feff->Eval(pt);

	if(pt>4.5) eff = feff->Eval(4.5);

	delete feff;		// prevent memory leak

	if(eff>0) return 1./eff;
	else return 0;
	
}

float inveff_tof(float pt, int charge=1, int MC=0) {		// pion embedding
	//return 1;			// test

	if(MC>0) return 1;					// pythia data, no correction needed

	if(charge==0) return 1;

	if(pt<MINPTCUT) return 0;

	float eff = 0;
	TF1* feff=new TF1("feff","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);
	feff->SetParameters(7.28477e-01, 1.52059e-01, 5.28031); 
	eff = feff->Eval(pt);

	if(pt>4.5) eff = feff->Eval(4.5);

	delete feff;		// prevent memory leak

	if(eff>0) return 1./eff;
	else return 0;
}


float inveff_tpc(float pt, int charge=1, int MC=0) {		// pion embedding
	//return 1;			// test

	if(MC>0) return 1;					// pythia data, no correction needed

	if(charge==0) return 1;

	if(pt<MINPTCUT) return 0;

	float eff = 0;
	TF1* feff=new TF1("feff","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);
	feff->SetParameters(0.874739, 0.156624, 5.67316); 
	eff = feff->Eval(pt);

	if(pt>4.5) eff = feff->Eval(4.5);

	delete feff;		// prevent memory leak

	if(eff>0) return 1./eff;
	else return 0;
}


