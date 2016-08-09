//==================================================================================================================================
//
//		2015.11.11	Li Yi
//		Read event info from ResultTree and Draw the Profile for Ntrk, Sum pT, Track pT
//		
//		This is to solve the issue that some runs have problematic behavior which ruins the original histograms  
//
//
//		2016.01.23	Li Yi
//		Previously TransMax or TransMin is defined as sum pt.
//		TransMax or TransMin definition depends on each variable.
//		If variable is Ntrk, TransMax or TransMin is determined by whose Ntrk is larger/smaller
//		If variable is sum pt, TransMax or TransMin is determined by whose sum pt is larger/smaller
//		However, currently <pT> TransMax or TransMin is still determined by whose sum pt is larger/smaller. as right now, it is track average value, not event average value
//
//		2016.02.01	Li Yi
//		TransPt and TransNtrk now are the average of TransMax and TranMin, Previously it was the sum of these two
//
//
//
//		2016.07.11	Li Yi
//		change sumleadntrk, sumsubntrk, sumtranmaxntrk, sumtranminntrk to float from int:
//		Previously, TOF matching or TPC efficiency is only corrected for pT, Sum pT, but not for multiplicity. After changing the Ntrk to float type, we can also apply correction for those. 
//		Note: for jet pT, we didn't apply any correction. Currently, because jet pT dependence is not large, we neglect it.. But in principle, we can unfolding those together in later steps.

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


float foldphi(float phi) {

	int npi = floor(fabs(phi)/(2*TMath::Pi()));

	if(phi>0) phi = phi-npi*2*TMath::Pi();
	if(phi<0) phi = phi+(npi+1)*2*TMath::Pi();

	if(phi>TMath::Pi()) phi = 2*TMath::Pi()-phi;

	return phi;
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

float inveff_tof(float pt, int charge=1, int MC=0) {		// real data 2000<zdc<3000 
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


double etafun_tpc(double *x, double *par) {

	double eta = x[0];
	double y = 0; 
		
	double turnpoint = 0.5;
	double maxeta = 5.3;
	double p0 = par[0];
	double p1 = par[1];
	double p2 = par[2];
	double p3 = par[3];
	double p4 = par[4];
	double p5 = par[5];
	if(eta<-turnpoint&&eta>=-1) {
		y = p0*exp(-pow(p1/(eta+maxeta),p2));
	}
	if(fabs(eta)<=turnpoint) {
		double k = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))-(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/(fabs(turnpoint)*2.);
		double b = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))+(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/2.;
		y = k*eta+b;
	}
	if(eta>turnpoint&&eta<=1) {
		y = p3*exp(-pow(p4/(maxeta-eta),p5));
	}

	return y;

}

double eta_weight_tpc(double eta) {	// efficiency eta dependence (integral as 1, since its value has been included in the pT eff one)

	TF1 *f = new TF1("f",etafun_tpc,-1,1,6); 
	f->SetParameters(1.00723e+00,4.18242e+00,8.10867e+01,1.00394e+00,4.20742e+00,9.46997e+01);

	return f->Eval(eta);			 
}

void plotTree2Histo(TString what2fill="multiplicity", TString dir="~/Scratch/pp200Y12_jetunderlying/", TString filetag = "underlyingevent_MB_R06_LeadJetAngle_FullJetFraclt90_160116", double jetptmin = 10, double jetptmax= 200, int ExclusiveEta = 0) {
// what2fill: refmult, leadjetpt, multiplicity, transntrk	(no space)
// ExclusiveEta: 
// 		==1: jet and underlying event in seperate eta region. for example: jet in [-0.6, 0), then undelrying event [0.6,1)
// 		==0: no requirement for underlying event and jet eta

	if( (!what2fill.EqualTo("refmult",TString::kIgnoreCase)) && (!what2fill.EqualTo("leadjetpt",TString::kIgnoreCase)) && (!what2fill.EqualTo("multiplicity",TString::kIgnoreCase)) && (!(what2fill.EqualTo("transntrk",TString::kIgnoreCase)||(what2fill.EqualTo("tranntrk",TString::kIgnoreCase)))) ) {
  	  cout<<"ERR!! call plotTree2Histo(TString what2fill, TString filepath): what2fill should be \"refmult\", \"leadjetpt\", \"multiplicity\", \"transntrk\" or \"tranntrk\"."<<endl;
  	  return;
  	}

	int MCflag=1;			// if you don't want to correct for TOF/TPC efficency at this point, for example, do it as part of unfolding procedure, you can just set MCflag=1. Then it will be treated as pure MC (not geant) and no correction will be applied. 
	// if you are working on embedding (MC + geant), you can apply corrections. 
	if(filetag.Contains("pythia",TString::kIgnoreCase)) 	{
		MCflag = 1;
		cout<<endl<<"INFO: process as MC data: no efficiency correction will be applied"<<endl;
	}

	int chargeflag = 1;
	if(filetag.Contains("Charge0",TString::kIgnoreCase)||filetag.Contains("TransNeutral",TString::kIgnoreCase)) chargeflag = 0;
	cout<<endl<<"INFO: underlying event chargeflag == "<<chargeflag<<endl;
	if(chargeflag==1&&MCflag==0) cout<<"Apply TPC tracking efficiency"<<endl;

	// TOFMATCH:
	// 		==1: tracks matched to TOF and the matching efficiency will be corrected
	// 		==0: No TOF matching efficiency correction
	int TOFMATCH = 1;		// default one: do the TOF Matching and correction needed.
	if(filetag.Contains("NoTofMatch",TString::kIgnoreCase)) TOFMATCH=0;
	if(chargeflag==1&&MCflag==0&&TOFMATCH==1) cout<<"Apply TOF matching efficiency"<<endl;
	cout<<endl;

	int savefig = 0;
	int saveroot = 1; 

	// reader
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root");
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_MB_R06_LeadJetAngle_MatchTrig_151206.root");
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/Charge0underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_JetCharge0Fraclt90_151208.root");
	TString filepath = dir+filetag+".root";
	TFile *f = new TFile(filepath);
	if(!f) { cout<<"Cannot find input file"<<endl; return; }

	TTree *t = (TTree*)f->Get("ResultTree");
	if(!t) { cout<<"Cannot find Tree"<<endl; return; }

	float jpt, jeta, jphi, jaspt, jaseta, leadpt, subpt, tranmaxpt,tranminpt, tranpt;
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
	if(what2fill.Contains("transntrk",TString::kIgnoreCase)||what2fill.Contains("tranntrk",TString::kIgnoreCase)) {
		xvariablename = "Transverse multiplicity";
	}

	TH1D *xaxis = new TH1D(what2fill,xvariablename,nbinning,0,maxpt);
	xaxis->Sumw2();

	TH2D *hetaphi = new TH2D("hetaphi","#Delta#eta-#Delta#phi",100,0,TMath::Pi(),100,-2,2);
	hetaphi->Sumw2();

	TProfile *j1ptvsxaxis = new TProfile("j1ptvs"+what2fill,"Leading Jet p_{T} vs "+xvariablename,nbinning,0,maxpt);
	TProfile *jasptvsxaxis = new TProfile("jasptvs"+what2fill,"Leading Jet p_{T} vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetntrkvsxaxis = new TProfile("leadjetareantrkvs"+what2fill,"Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetntrkvsxaxis = new TProfile("subjetareantrkvs"+what2fill,"SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxntrkvsxaxis = new TProfile("tranmaxntrkvs"+what2fill,"Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminntrkvsxaxis = new TProfile("tranminntrkvs"+what2fill,"Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranntrkvsxaxis = new TProfile("tranntrkvs"+what2fill,"Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptsumvsxaxis = new TProfile("leadjetareaptsumvs"+what2fill,"Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptsumvsxaxis = new TProfile("subjetareaptsumvs"+what2fill,"SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptsumvsxaxis = new TProfile("tranmaxptsumvs"+what2fill,"Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptsumvsxaxis = new TProfile("tranminptsumvs"+what2fill,"Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptsumvsxaxis = new TProfile("tranptsumvs"+what2fill,"Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptavevsxaxis = new TProfile("leadjetareaptavevs"+what2fill,"Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptavevsxaxis = new TProfile("subjetareaptavevs"+what2fill,"SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptavevsxaxis = new TProfile("tranmaxptavevs"+what2fill,"Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptavevsxaxis = new TProfile("tranminptavevs"+what2fill,"Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptavevsxaxis = new TProfile("tranptavevs"+what2fill,"Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt);

	TProfile *maxtranptvsxaxis = new TProfile("maxtranptvs"+what2fill,"Maximum Transverse Pt vs "+xvariablename,nbinning,0,maxpt);
	
	int nbinning_tpt = 1000;
	double max_tpt = 50;
	TH2D *hj1ptvsxaxis = new TH2D("hj1ptvs"+what2fill,"Leading Jet p_{T} vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
	TH2D *hjasptvsxaxis = new TH2D("hjasptvs"+what2fill,"Leading Jet p_{T} vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
	
        TH2D *hleadjetntrkvsxaxis = new TH2D("hleadjetareantrkvs"+what2fill,"Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetntrkvsxaxis = new TH2D("hsubjetareantrkvs"+what2fill,"SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxntrkvsxaxis = new TH2D("htranmaxntrkvs"+what2fill,"Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminntrkvsxaxis = new TH2D("htranminntrkvs"+what2fill,"Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranntrkvsxaxis = new TH2D("htranntrkvs"+what2fill,"Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptsumvsxaxis = new TH2D("hleadjetareaptsumvs"+what2fill,"Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptsumvsxaxis = new TH2D("hsubjetareaptsumvs"+what2fill,"SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptsumvsxaxis = new TH2D("htranmaxptsumvs"+what2fill,"Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptsumvsxaxis = new TH2D("htranminptsumvs"+what2fill,"Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptsumvsxaxis = new TH2D("htranptsumvs"+what2fill,"Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptavevsxaxis = new TH2D("hleadjetareaptavevs"+what2fill,"Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptavevsxaxis = new TH2D("hsubjetareaptavevs"+what2fill,"SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptavevsxaxis = new TH2D("htranmaxptavevs"+what2fill,"Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptavevsxaxis = new TH2D("htranminptavevs"+what2fill,"Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptavevsxaxis = new TH2D("htranptavevs"+what2fill,"Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hmaxtranptvsxaxis = new TH2D("hmaxtranptvs"+what2fill,"Maximum Transverse Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

	hj1ptvsxaxis->Sumw2();
	hjasptvsxaxis->Sumw2();
	
        hleadjetntrkvsxaxis->Sumw2();
        hsubjetntrkvsxaxis->Sumw2();
        htranmaxntrkvsxaxis->Sumw2();
        htranminntrkvsxaxis->Sumw2();
        htranntrkvsxaxis->Sumw2();

        hleadjetptsumvsxaxis->Sumw2();
        hsubjetptsumvsxaxis->Sumw2();
        htranmaxptsumvsxaxis->Sumw2();
        htranminptsumvsxaxis->Sumw2();
        htranptsumvsxaxis->Sumw2();

        hleadjetptavevsxaxis->Sumw2();
        hsubjetptavevsxaxis->Sumw2();
        htranmaxptavevsxaxis->Sumw2();
        htranminptavevsxaxis->Sumw2();
        htranptavevsxaxis->Sumw2();

	hmaxtranptvsxaxis->Sumw2();


	// problematic runs, need future investigation		--> already exclude in ttree production code
	//const int NoBadRun = 186;
	//int badrun[NoBadRun] = {13044118, 13044123, 13044124, 13044125, 13045001, 13045003, 13045005, 13045006, 13045007, 13045012, 13045029, 13046002, 13046008, 13046010, 13046029, 13046118, 13046119, 13046120, 13047004, 13047014, 13047018, 13047036, 13047037, 13047039, 13047040, 13047041, 13047042, 13047043, 13047044, 13047045, 13047046, 13047047, 13047048, 13047049, 13047050, 13047051, 13047052, 13047053, 13047054, 13047055, 13048007, 13048022, 13048046, 13049004, 13049005, 13049050, 13049052, 13049075, 13049086, 13049087, 13049088, 13049089, 13050007, 13050025, 13050026, 13050027, 13050033, 13050039, 13050043, 13050044, 13050046, 13050047, 13050049, 13050050, 13051068, 13051080, 13051088, 13051095, 13051102, 13052021, 13052022, 13052054, 13052063, 13052068, 13053010, 13053021, 13054004, 13054005, 13054006, 13054007, 13054008, 13054009, 13054011, 13054012, 13054013, 13054014, 13054015, 13054016, 13054017, 13054018, 13054019, 13054020, 13054022, 13054042, 13054045, 13054046, 13054057, 13055015, 13055072, 13055081, 13055082, 13055086, 13055087, 13055088, 13055089, 13055090, 13056011, 13056012, 13056034, 13056035, 13056037, 13056038, 13056039, 13057038, 13057039, 13058019, 13058030, 13058047, 13058048, 13059003, 13059004, 13059005, 13059006, 13059007, 13059008, 13059009, 13059010, 13059019, 13059035, 13059082, 13059083, 13059084, 13059085, 13059086, 13059087, 13060001, 13060002, 13060003, 13060009, 13060012, 13061026, 13063033, 13064030, 13064057, 13064059, 13064074, 13066035, 13066036, 13066101, 13066102, 13066104, 13066109, 13066110, 13067001, 13067002, 13067003, 13067004, 13067005, 13067006, 13067007, 13067008, 13067009, 13067010, 13067011, 13067012, 13067013, 13067014, 13067015, 13067017, 13068017, 13068022, 13068027, 13068029, 13068034, 13068036, 13068037, 13069006, 13069009, 13069029, 13070030, 13070056, 13071034, 13071037, 13071038, 13071040};

	// loop over events
	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
	int processedevent = 0;
	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		t->GetEntry(ievt);

		if(j1neutralfrac>0.9) continue;			// jet neutral pT fraction < 90%
		//if(filetag.Contains("MB") && j1neutralfrac>0.75) continue;			// jet neutral pT fraction < 75% for MB data (there is a hot area for MB > 75%) -- UPDATE: It is due to the hot stripe in dataset at the begining of the run, shall mark as bad run
		if(runid<=13046029) continue;		// hot BEMC stripe in MB dataset
		
		if(fabs(jeta)>0.6) continue;			// I miss this cut in jet finding selector

		//if(fabs(jaspt)<5) continue;			
		//if(fabs(jpt-jaspt)>0.05*jpt) continue;		// dijet balanced

		//if(runid>=13058000&& runid<13061000) continue;		// a dip in TPC primary tracks. problematic runs, need future investigation		--> already exclude in ttree production code
		//test if(runid>13052000&& runid<13060000) continue;		// problematic runs, need future investigation
		//for(int i = 0; i<NoBadRun; i++) {
		//	if(runid==badrun[i]) continue;
		//}

		//if(((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase)))&&((jpt<10)||(jpt>40))) continue;		// when studying multiplicity dependence, exclude jet pT<10GeV/c to ensure the real jet found + jet pT<40GeV/c there is some unphysics structure seen from Ntrk, SumPt, pT vs leadjetpt distribution(solved -> due to hot tower missed in QA, now excluded)
		//if(((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase)))&&(jpt<10)) continue;		// when studying multiplicity dependence, exclude jet pT<10GeV/c to ensure the real jet found 
		if(((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase))||(what2fill.Contains("transntrk",TString::kIgnoreCase))||(what2fill.Contains("tranntrk",TString::kIgnoreCase)))&&(jpt<jetptmin||jpt>=jetptmax)) continue;	

		if(ievt%1000000==0) cout<<"event "<<ievt<<endl;

		double xvariable = 0;
		if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
			xvariable = jpt;
		}
		else if(what2fill.Contains("multiplicity",TString::kIgnoreCase)) {	// total mulitplicity
			xvariable = leadntrk+subntrk+tranmaxntrk+tranminntrk;
		}
		else if(what2fill.Contains("transntrk",TString::kIgnoreCase)||what2fill.Contains("tranntrk",TString::kIgnoreCase)) {	// transverse mulitplicity
			xvariable = tranmaxntrk+tranminntrk;
		}
		else {
			xvariable = refmult;
		}

		xaxis->Fill(xvariable);	

		j1ptvsxaxis->Fill(xvariable,jpt);
		jasptvsxaxis->Fill(xvariable,jaspt);
		
		hj1ptvsxaxis->Fill(xvariable,jpt);
		hjasptvsxaxis->Fill(xvariable,jaspt);

		float maxtranspt = max(*max_element(pt_min,pt_min+tranminntrk),*max_element(pt_max,pt_max+tranmaxntrk));	
		maxtranptvsxaxis->Fill(xvariable,maxtranspt);

		hmaxtranptvsxaxis->Fill(xvariable,maxtranspt);


		int etaflag = 0;
		if( ExclusiveEta>0 && jeta>0 )	etaflag = -1;
		if( ExclusiveEta>0 && jeta<0 )	etaflag = 1;

		//if( jeta>0 && jaseta<0 ) 	continue;		// requirement dijet on the same side	testly
		//if( jeta<=0 && jaseta>=0 ) 	continue;		// requirement dijet on the same side	testly


		float sumleadpt = 0, sumsubpt = 0, sumtranmaxpt = 0, sumtranminpt = 0, sumtranpt = 0;
		float sumleadntrk = 0, sumsubntrk = 0, sumtranmaxntrk = 0, sumtranminntrk = 0;
		//cout<<"tranmaxntrk = "<<tranmaxntrk<<"\ttranminntrk = "<<tranminntrk<<"\tleadntrk = "<<leadntrk<<"\tsubntrk = "<<subntrk<<endl;
		float w;
		for(int it = 0; it<tranmaxntrk; it++) {
			//cout<<"pt_max["<<it<<"] = "<<pt_max[it];
			if(etaflag==1&&eta_max[it]<0.6) continue;
			if(etaflag==-1&&eta_max[it]>-0.6) continue;

			if(TOFMATCH) {
				w = getweight(pt_max[it],chargeflag,MCflag);		// If TOF matched
			}
			else {
				w = inveff_tpc(pt_max[it],chargeflag,MCflag);		// if NO TOF match
			}
			//cout<<" w = "<<w<<endl;
			tranmaxptavevsxaxis->Fill(xvariable,pt_max[it],w);
			tranptavevsxaxis->Fill(xvariable,pt_max[it],w);

			htranmaxptavevsxaxis->Fill(xvariable,pt_max[it],w);
			htranptavevsxaxis->Fill(xvariable,pt_max[it],w);

			hetaphi->Fill(foldphi(phi_max[it]-jphi),eta_max[it]-jeta);

			sumtranmaxpt+=pt_max[it]*w;		// will be w tracks effectively ..
			sumtranpt+=pt_max[it]*w;		// will be w tracks effectively ..

			sumtranmaxntrk+=w;
		}
		for(int it = 0; it<tranminntrk; it++) {
			//cout<<"pt_min["<<it<<"] = "<<pt_min[it];
			if(etaflag==1&&eta_min[it]<0.6) continue;
			if(etaflag==-1&&eta_min[it]>-0.6) continue;

			if(TOFMATCH) {
				w = getweight(pt_min[it],chargeflag,MCflag);
			}
			else {
				w = inveff_tpc(pt_min[it],chargeflag,MCflag);
			}
			//cout<<" w = "<<w<<endl;
			tranminptavevsxaxis->Fill(xvariable,pt_min[it],w);
			tranptavevsxaxis->Fill(xvariable,pt_min[it],w);

			htranminptavevsxaxis->Fill(xvariable,pt_min[it],w);
			htranptavevsxaxis->Fill(xvariable,pt_min[it],w);

			hetaphi->Fill(foldphi(phi_min[it]-jphi),eta_min[it]-jeta);

			sumtranminpt+=pt_min[it]*w;		// will be w tracks effectively ..
			sumtranpt+=pt_min[it]*w;		// will be w tracks effectively ..

			sumtranminntrk+=w;
		}
		for(int it = 0; it<leadntrk; it++) {
			//cout<<"pt_jet["<<it<<"] = "<<pt_jet[it];
			if(etaflag==1&&eta_jet[it]<0.6) continue;
			if(etaflag==-1&&eta_jet[it]>-0.6) continue;

			if(TOFMATCH) {
				w = getweight(pt_jet[it],chargeflag,MCflag);
			}
			else {
				w = inveff_tpc(pt_jet[it],chargeflag,MCflag);
			}
			//cout<<" w = "<<w<<endl;
			leadjetptavevsxaxis->Fill(xvariable,pt_jet[it],w);

			hleadjetptavevsxaxis->Fill(xvariable,pt_jet[it],w);
			
			hetaphi->Fill(foldphi(phi_jet[it]-jphi),eta_jet[it]-jeta);

			sumleadpt+=pt_jet[it]*w;		// will be w tracks effectively ..

			sumleadntrk+=w;
		}
		for(int it = 0; it<subntrk; it++) {
			//cout<<"pt_sub["<<it<<"] = "<<pt_sub[it];
			if(etaflag==1&&eta_sub[it]<0.6) continue;
			if(etaflag==-1&&eta_sub[it]>-0.6) continue;

			if(TOFMATCH) {
				w = getweight(pt_sub[it],chargeflag,MCflag);
			}
			else {
				w = inveff_tpc(pt_sub[it],chargeflag,MCflag);
			}
			//cout<<" w = "<<w<<endl;
			subjetptavevsxaxis->Fill(xvariable,pt_sub[it],w);

			hsubjetptavevsxaxis->Fill(xvariable,pt_sub[it],w);

			hetaphi->Fill(foldphi(phi_sub[it]-jphi),eta_sub[it]-jeta);

			sumsubpt+=pt_sub[it]*w;		// will be w tracks effectively ..

			sumsubntrk+=w;
		}

        	leadjetntrkvsxaxis->Fill(xvariable,sumleadntrk);
        	hleadjetntrkvsxaxis->Fill(xvariable,sumleadntrk);
        	subjetntrkvsxaxis->Fill(xvariable,sumsubntrk);
        	hsubjetntrkvsxaxis->Fill(xvariable,sumsubntrk);
        	tranntrkvsxaxis->Fill(xvariable,(sumtranmaxntrk+sumtranminntrk)/2.);		// do not use tranntrk directly because it is int type and it is the sum divided by 2 
        	htranntrkvsxaxis->Fill(xvariable,(sumtranmaxntrk+sumtranminntrk)/2.);

		MaxOrMin(sumtranmaxntrk,sumtranminntrk);
        	tranmaxntrkvsxaxis->Fill(xvariable,sumtranmaxntrk);
        	htranmaxntrkvsxaxis->Fill(xvariable,sumtranmaxntrk);
        	tranminntrkvsxaxis->Fill(xvariable,sumtranminntrk);
        	htranminntrkvsxaxis->Fill(xvariable,sumtranminntrk);


		sumtranpt=sumtranpt/2.;		// take the average 
        	leadjetptsumvsxaxis->Fill(xvariable,sumleadpt);
        	subjetptsumvsxaxis->Fill(xvariable,sumsubpt);
        	tranptsumvsxaxis->Fill(xvariable,sumtranpt);
        	hleadjetptsumvsxaxis->Fill(xvariable,sumleadpt);
        	hsubjetptsumvsxaxis->Fill(xvariable,sumsubpt);
        	htranptsumvsxaxis->Fill(xvariable,sumtranpt);

		MaxOrMin(sumtranmaxpt,sumtranminpt);
        	tranmaxptsumvsxaxis->Fill(xvariable,sumtranmaxpt);
        	tranminptsumvsxaxis->Fill(xvariable,sumtranminpt);
        	htranmaxptsumvsxaxis->Fill(xvariable,sumtranmaxpt);
        	htranminptsumvsxaxis->Fill(xvariable,sumtranminpt);

		processedevent++;
	}
	cout<<"Total processed # of event "<<processedevent<<endl;

	// set histogram draw properties
	xaxis->GetXaxis()->SetTitle(xvariablename);
	xaxis->GetYaxis()->SetTitle("events");
	xaxis->SetLineColor(1);
	//xaxis->SetMarkerStyle(8);
	//xaxis->SetMarkerColor(1);
	
	int clead = 1, csub = 9, cmax = 2, cmin = 8, ctran = 28;
	int slead = 20, ssub = 25, smax = 24, smin = 20, stran = 21;

        j1ptvsxaxis->GetXaxis()->SetTitle(xvariablename);
	j1ptvsxaxis->GetYaxis()->SetTitle("Leading Jet p_{T}");
	j1ptvsxaxis->SetLineColor(clead);
	j1ptvsxaxis->SetMarkerStyle(slead);
	j1ptvsxaxis->SetMarkerColor(clead);

        jasptvsxaxis->GetXaxis()->SetTitle(xvariablename);
	jasptvsxaxis->GetYaxis()->SetTitle("Recoil Jet p_{T}");
	jasptvsxaxis->SetLineColor(clead);
	jasptvsxaxis->SetMarkerStyle(slead);
	jasptvsxaxis->SetMarkerColor(clead);

	hetaphi->GetXaxis()->SetTitle("#Delta#phi");
	hetaphi->GetYaxis()->SetTitle("#Delta#eta");
	hetaphi->SetLineColor(clead);
	hetaphi->SetMarkerStyle(slead);
	hetaphi->SetMarkerColor(clead);

        leadjetntrkvsxaxis->GetXaxis()->SetTitle(xvariablename);
	leadjetntrkvsxaxis->GetYaxis()->SetTitle("Leading Jet Area Multiplicity");
	leadjetntrkvsxaxis->SetLineColor(clead);
	leadjetntrkvsxaxis->SetMarkerStyle(slead);
	leadjetntrkvsxaxis->SetMarkerColor(clead);

        subjetntrkvsxaxis->GetXaxis()->SetTitle(xvariablename);
	subjetntrkvsxaxis->GetYaxis()->SetTitle("Away Side Area Multiplicity");
	subjetntrkvsxaxis->SetLineColor(csub);
	subjetntrkvsxaxis->SetMarkerStyle(ssub);
	subjetntrkvsxaxis->SetMarkerColor(csub);

        tranmaxntrkvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranmaxntrkvsxaxis->GetYaxis()->SetTitle("Transverse Max Area Multiplicity");
	tranmaxntrkvsxaxis->SetLineColor(cmax);
	tranmaxntrkvsxaxis->SetMarkerStyle(smax);
	tranmaxntrkvsxaxis->SetMarkerColor(cmax);

        tranminntrkvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranminntrkvsxaxis->GetYaxis()->SetTitle("Transverse Min Area Multiplicity");
	tranminntrkvsxaxis->SetLineColor(cmin);
	tranminntrkvsxaxis->SetMarkerStyle(smin);
	tranminntrkvsxaxis->SetMarkerColor(cmin);

        tranntrkvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranntrkvsxaxis->GetYaxis()->SetTitle("Transverse Area Multiplicity");
	tranntrkvsxaxis->SetLineColor(ctran);
	tranntrkvsxaxis->SetMarkerStyle(stran);
	tranntrkvsxaxis->SetMarkerColor(ctran);


        leadjetptsumvsxaxis->GetXaxis()->SetTitle(xvariablename);
	leadjetptsumvsxaxis->GetYaxis()->SetTitle("Lead Jet Area Sum Pt");
	leadjetptsumvsxaxis->SetLineColor(clead);
	leadjetptsumvsxaxis->SetMarkerStyle(slead);
	leadjetptsumvsxaxis->SetMarkerColor(clead);

        subjetptsumvsxaxis->GetXaxis()->SetTitle(xvariablename);
	subjetptsumvsxaxis->GetYaxis()->SetTitle("Away Side Are Sum Pt");
	subjetptsumvsxaxis->SetLineColor(csub);
	subjetptsumvsxaxis->SetMarkerStyle(ssub);
	subjetptsumvsxaxis->SetMarkerColor(csub);

        tranmaxptsumvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranmaxptsumvsxaxis->GetYaxis()->SetTitle("Transverse Max Area Sum Pt");
	tranmaxptsumvsxaxis->SetLineColor(cmax);
	tranmaxptsumvsxaxis->SetMarkerStyle(smax);
	tranmaxptsumvsxaxis->SetMarkerColor(cmax);

        tranminptsumvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranminptsumvsxaxis->GetYaxis()->SetTitle("Transverse Min Area Sum Pt");
	tranminptsumvsxaxis->SetLineColor(cmin);
	tranminptsumvsxaxis->SetMarkerStyle(smin);
	tranminptsumvsxaxis->SetMarkerColor(cmin);

        tranptsumvsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranptsumvsxaxis->GetYaxis()->SetTitle("Transverse Area Sum Pt");
	tranptsumvsxaxis->SetLineColor(ctran);
	tranptsumvsxaxis->SetMarkerStyle(stran);
	tranptsumvsxaxis->SetMarkerColor(ctran);


        leadjetptavevsxaxis->GetXaxis()->SetTitle(xvariablename);
	leadjetptavevsxaxis->GetYaxis()->SetTitle("Leading Area Average Track Pt");
	leadjetptavevsxaxis->SetLineColor(clead);
	leadjetptavevsxaxis->SetMarkerStyle(slead);
	leadjetptavevsxaxis->SetMarkerColor(clead);

        subjetptavevsxaxis->GetXaxis()->SetTitle(xvariablename);
	subjetptavevsxaxis->GetYaxis()->SetTitle("Away side Area Average Track Pt");
	subjetptavevsxaxis->SetLineColor(csub);
	subjetptavevsxaxis->SetMarkerStyle(ssub);
	subjetptavevsxaxis->SetMarkerColor(csub);

        tranmaxptavevsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranmaxptavevsxaxis->GetYaxis()->SetTitle("Transverse Max Area Average Track Pt");
	tranmaxptavevsxaxis->SetLineColor(cmax);
	tranmaxptavevsxaxis->SetMarkerStyle(smax);
	tranmaxptavevsxaxis->SetMarkerColor(cmax);

        tranminptavevsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranminptavevsxaxis->GetYaxis()->SetTitle("Transverse Min Area Average Track Pt");
	tranminptavevsxaxis->SetLineColor(cmin);
	tranminptavevsxaxis->SetMarkerStyle(smin);
	tranminptavevsxaxis->SetMarkerColor(cmin);

        tranptavevsxaxis->GetXaxis()->SetTitle(xvariablename);
	tranptavevsxaxis->GetYaxis()->SetTitle("Transverse Area Average Track Pt");
	tranptavevsxaxis->SetLineColor(ctran);
	tranptavevsxaxis->SetMarkerStyle(stran);
	tranptavevsxaxis->SetMarkerColor(ctran);

        maxtranptvsxaxis->GetXaxis()->SetTitle(xvariablename);
	maxtranptvsxaxis->GetYaxis()->SetTitle("Maximum Transverse Pt");
	maxtranptvsxaxis->SetLineColor(ctran);
	maxtranptvsxaxis->SetMarkerStyle(stran);
	maxtranptvsxaxis->SetMarkerColor(ctran);


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
	xaxis->Draw();
	if(savefig) c[0]->SaveAs(Form("figs/%sVs%s_%s.png",xaxis->GetName(),what2fill.Data(),filetag.Data()));
	
	TLegend *leg = new TLegend(0.16,0.6,0.35,0.86);
	leg->AddEntry(leadjetntrkvsxaxis,"Toward","pl");
	leg->AddEntry(subjetntrkvsxaxis,"Away","pl");
	leg->AddEntry(tranmaxntrkvsxaxis,"TransMax","pl");
	leg->AddEntry(tranminntrkvsxaxis,"TransMin","pl");
	leg->AddEntry(tranntrkvsxaxis,"Trans","pl");
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
	leadjetntrkvsxaxis->Draw("same");
	subjetntrkvsxaxis->Draw("same");
	tranmaxntrkvsxaxis->Draw("same");
	tranminntrkvsxaxis->Draw("same");	
	tranntrkvsxaxis->Draw("same");	
	leg->Draw("same");
	if(savefig) c[1]->SaveAs(Form("figs/%sVs%s_%s.png","Ntrk",what2fill.Data(),filetag.Data()));

	c[2]->cd();
	htmp[2] = new TH2D("htmp2","",1000,0,drawxmax,1000,0,20);
	htmp[2]->GetXaxis()->SetTitle(xvariablename);
	htmp[2]->GetYaxis()->SetTitle("Event #LT#sum p_{T}#GT");
	htmp[2]->Draw();
	lat->DrawLatex(0.5,0.94,"Event #LT#sum p_{T}#GT");
	leadjetptsumvsxaxis->Draw("same");
	subjetptsumvsxaxis->Draw("same");
	tranmaxptsumvsxaxis->Draw("same");
	tranminptsumvsxaxis->Draw("same");	
	tranptsumvsxaxis->Draw("same");	
	leg->Draw("same");
	if(savefig) c[2]->SaveAs(Form("figs/%sVs%s_%s.png","PtSum",what2fill.Data(),filetag.Data()));
	
	c[3]->cd();
	htmp[3] = new TH2D("htmp3","",1000,0,drawxmax,1000,0,7);
	htmp[3]->GetXaxis()->SetTitle(xvariablename);
	htmp[3]->GetYaxis()->SetTitle("Track #LTp_{T}#GT");
	htmp[3]->Draw();
	lat->DrawLatex(0.5,0.94,"Track #LTp_{T}#GT");
	leadjetptavevsxaxis->Draw("same");
	subjetptavevsxaxis->Draw("same");
	tranmaxptavevsxaxis->Draw("same");
	tranminptavevsxaxis->Draw("same");	
	tranptavevsxaxis->Draw("same");	
	leg->Draw("same");
	if(savefig) c[3]->SaveAs(Form("figs/%sVs%s_%s.png","PtAve",what2fill.Data(),filetag.Data()));


	if(saveroot) {
		//TFile *fout = new TFile("TwoHisto4underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root","RECREATE");
		//TString outfilepath = dir+what2fill+"hist4"+filetag+".root";
		TString jettag="";
		if( ((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase))||(what2fill.Contains("transntrk",TString::kIgnoreCase))||(what2fill.Contains("tranntrk",TString::kIgnoreCase))) && (fabs(jetptmin-10)>1e-6||fabs(jetptmax-200)>1e-6)) {
			jettag = Form("_jet%g-%g",jetptmin,jetptmax);	
		} 
		if(MCflag==1) {
			jettag+="_NoEffCorr";
		}
		TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_160712.root";
		//TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_NeutralFrac75.root";
		//TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_wR1"+".root";
		//TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_AsJetGt5"+".root";
		//TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_balance005"+".root";
		if(ExclusiveEta>0) {
			outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_EtaExcl"+".root";
			//outfilepath = dir+what2fill+"hist4"+filetag+jettag+"_EtaExcl_SamesideDijet"+".root";
		}
		TFile *fout = new TFile(outfilepath,"RECREATE");
		xaxis->Write();
		
		j1ptvsxaxis->Write();
		jasptvsxaxis->Write();

		hetaphi->Write();

        	leadjetntrkvsxaxis->Write();
        	subjetntrkvsxaxis->Write();
        	tranmaxntrkvsxaxis->Write();
        	tranminntrkvsxaxis->Write();
        	tranntrkvsxaxis->Write();

        	leadjetptsumvsxaxis->Write();
        	subjetptsumvsxaxis->Write();
        	tranmaxptsumvsxaxis->Write();
        	tranminptsumvsxaxis->Write();
        	tranptsumvsxaxis->Write();

        	leadjetptavevsxaxis->Write();
        	subjetptavevsxaxis->Write();
        	tranmaxptavevsxaxis->Write();
        	tranminptavevsxaxis->Write();
        	tranptavevsxaxis->Write();

		maxtranptvsxaxis->Write();


		hj1ptvsxaxis->Write();
		hjasptvsxaxis->Write();
	
        	hleadjetntrkvsxaxis->Write();
        	hsubjetntrkvsxaxis->Write();
        	htranmaxntrkvsxaxis->Write();
        	htranminntrkvsxaxis->Write();
        	htranntrkvsxaxis->Write();

        	hleadjetptsumvsxaxis->Write();
        	hsubjetptsumvsxaxis->Write();
        	htranmaxptsumvsxaxis->Write();
        	htranminptsumvsxaxis->Write();
        	htranptsumvsxaxis->Write();

        	hleadjetptavevsxaxis->Write();
        	hsubjetptavevsxaxis->Write();
        	htranmaxptavevsxaxis->Write();
        	htranminptavevsxaxis->Write();
        	htranptavevsxaxis->Write();

		hmaxtranptvsxaxis->Write();
			
		fout->Close();
	}

}


