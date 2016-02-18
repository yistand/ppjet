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

void plotTree2Histo(TString what2fill="multiplicity", TString dir="~/Scratch/pp200Y12_jetunderlying/", TString filetag = "underlyingevent_MB_R06_LeadJetAngle_FullJetFraclt90_160116", double jetptmin = 10, double jetptmax= 200) {
// what2fill: refmult, leadjetpt, multiplicity, transntrk	(no space)
	if( (!what2fill.EqualTo("refmult",TString::kIgnoreCase)) && (!what2fill.EqualTo("leadjetpt",TString::kIgnoreCase)) && (!what2fill.EqualTo("multiplicity",TString::kIgnoreCase)) && (!(what2fill.EqualTo("transntrk",TString::kIgnoreCase)||(what2fill.EqualTo("tranntrk",TString::kIgnoreCase)))) ) {
  	  cout<<"ERR!! call plotTree2Histo(TString what2fill, TString filepath): what2fill should be \"refmult\", \"leadjetpt\", \"multiplicity\", \"transntrk\" or \"tranntrk\"."<<endl;
  	  return;
  	}

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
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root");
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/underlyingevent_MB_R06_LeadJetAngle_MatchTrig_151206.root");
	//TFile *f = new TFile("~/Scratch/pp200Y12_jetunderlying/Charge0underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_JetCharge0Fraclt90_151208.root");
	TString filepath = dir+filetag+".root";
	TFile *f = new TFile(filepath);
	if(!f) { cout<<"Cannot find input file"<<endl; return; }

	TTree *t = (TTree*)f->Get("ResultTree");
	if(!t) { cout<<"Cannot find Tree"<<endl; return; }

	float jpt, jeta, leadpt, subpt, tranmaxpt,tranminpt, tranpt;
	int leadntrk, subntrk, tranmaxntrk, tranminntrk, tranntrk;
	double refmult;
	int runid;

	float j1neutralfrac;

	const int MAXARRAY = 1000;
	float pt_min[MAXARRAY], pt_max[MAXARRAY], pt_jet[MAXARRAY], pt_sub[MAXARRAY];

	t->SetBranchAddress("runid",&runid);
	t->SetBranchAddress("refmult",&refmult);
	t->SetBranchAddress("j1pt",&jpt);
	t->SetBranchAddress("j1eta",&jeta);
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

	TH1D *leadjetpt = new TH1D(what2fill,xvariablename,nbinning,0,maxpt);
	leadjetpt->Sumw2();

        TProfile *leadjetntrkvsleadjetpt = new TProfile("leadjetareantrkvs"+what2fill,"Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetntrkvsleadjetpt = new TProfile("subjetareantrkvs"+what2fill,"SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxntrkvsleadjetpt = new TProfile("tranmaxntrkvs"+what2fill,"Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminntrkvsleadjetpt = new TProfile("tranminntrkvs"+what2fill,"Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranntrkvsleadjetpt = new TProfile("tranntrkvs"+what2fill,"Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptsumvsleadjetpt = new TProfile("leadjetareaptsumvs"+what2fill,"Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptsumvsleadjetpt = new TProfile("subjetareaptsumvs"+what2fill,"SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptsumvsleadjetpt = new TProfile("tranmaxptsumvs"+what2fill,"Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptsumvsleadjetpt = new TProfile("tranminptsumvs"+what2fill,"Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptsumvsleadjetpt = new TProfile("tranptsumvs"+what2fill,"Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt);

        TProfile *leadjetptavevsleadjetpt = new TProfile("leadjetareaptavevs"+what2fill,"Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *subjetptavevsleadjetpt = new TProfile("subjetareaptavevs"+what2fill,"SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranmaxptavevsleadjetpt = new TProfile("tranmaxptavevs"+what2fill,"Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranminptavevsleadjetpt = new TProfile("tranminptavevs"+what2fill,"Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt);
        TProfile *tranptavevsleadjetpt = new TProfile("tranptavevs"+what2fill,"Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt);

	TProfile *maxtranptvsleadjetpt = new TProfile("maxtranptvs"+what2fill,"Maximum Transverse Pt vs "+xvariablename,nbinning,0,maxpt);
	
	int nbinning_tpt = 1000;
	double max_tpt = 50;
        TH2D *hleadjetntrkvsleadjetpt = new TH2D("hleadjetareantrkvs"+what2fill,"Leading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetntrkvsleadjetpt = new TH2D("hsubjetareantrkvs"+what2fill,"SubLeading Jet Area Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxntrkvsleadjetpt = new TH2D("htranmaxntrkvs"+what2fill,"Transverse Max Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminntrkvsleadjetpt = new TH2D("htranminntrkvs"+what2fill,"Transverse Min Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranntrkvsleadjetpt = new TH2D("htranntrkvs"+what2fill,"Transverse Ntrk vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptsumvsleadjetpt = new TH2D("hleadjetareaptsumvs"+what2fill,"Leading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptsumvsleadjetpt = new TH2D("hsubjetareaptsumvs"+what2fill,"SubLeading Jet Area Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptsumvsleadjetpt = new TH2D("htranmaxptsumvs"+what2fill,"Transverse Max Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptsumvsleadjetpt = new TH2D("htranminptsumvs"+what2fill,"Transverse Min Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptsumvsleadjetpt = new TH2D("htranptsumvs"+what2fill,"Transverse Sum Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hleadjetptavevsleadjetpt = new TH2D("hleadjetareaptavevs"+what2fill,"Leading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *hsubjetptavevsleadjetpt = new TH2D("hsubjetareaptavevs"+what2fill,"SubLeading Jet Area Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranmaxptavevsleadjetpt = new TH2D("htranmaxptavevs"+what2fill,"Transverse Max Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranminptavevsleadjetpt = new TH2D("htranminptavevs"+what2fill,"Transverse Min Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
        TH2D *htranptavevsleadjetpt = new TH2D("htranptavevs"+what2fill,"Transverse Average Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        TH2D *hmaxtranptvsleadjetpt = new TH2D("hmaxtranptvs"+what2fill,"Maximum Transverse Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

        hleadjetntrkvsleadjetpt->Sumw2();
        hsubjetntrkvsleadjetpt->Sumw2();
        htranmaxntrkvsleadjetpt->Sumw2();
        htranminntrkvsleadjetpt->Sumw2();
        htranntrkvsleadjetpt->Sumw2();

        hleadjetptsumvsleadjetpt->Sumw2();
        hsubjetptsumvsleadjetpt->Sumw2();
        htranmaxptsumvsleadjetpt->Sumw2();
        htranminptsumvsleadjetpt->Sumw2();
        htranptsumvsleadjetpt->Sumw2();

        hleadjetptavevsleadjetpt->Sumw2();
        hsubjetptavevsleadjetpt->Sumw2();
        htranmaxptavevsleadjetpt->Sumw2();
        htranminptavevsleadjetpt->Sumw2();
        htranptavevsleadjetpt->Sumw2();

	hmaxtranptvsleadjetpt->Sumw2();


	// problematic runs, need future investigation		--> already exclude in ttree production code
	//const int NoBadRun = 186;
	//int badrun[NoBadRun] = {13044118, 13044123, 13044124, 13044125, 13045001, 13045003, 13045005, 13045006, 13045007, 13045012, 13045029, 13046002, 13046008, 13046010, 13046029, 13046118, 13046119, 13046120, 13047004, 13047014, 13047018, 13047036, 13047037, 13047039, 13047040, 13047041, 13047042, 13047043, 13047044, 13047045, 13047046, 13047047, 13047048, 13047049, 13047050, 13047051, 13047052, 13047053, 13047054, 13047055, 13048007, 13048022, 13048046, 13049004, 13049005, 13049050, 13049052, 13049075, 13049086, 13049087, 13049088, 13049089, 13050007, 13050025, 13050026, 13050027, 13050033, 13050039, 13050043, 13050044, 13050046, 13050047, 13050049, 13050050, 13051068, 13051080, 13051088, 13051095, 13051102, 13052021, 13052022, 13052054, 13052063, 13052068, 13053010, 13053021, 13054004, 13054005, 13054006, 13054007, 13054008, 13054009, 13054011, 13054012, 13054013, 13054014, 13054015, 13054016, 13054017, 13054018, 13054019, 13054020, 13054022, 13054042, 13054045, 13054046, 13054057, 13055015, 13055072, 13055081, 13055082, 13055086, 13055087, 13055088, 13055089, 13055090, 13056011, 13056012, 13056034, 13056035, 13056037, 13056038, 13056039, 13057038, 13057039, 13058019, 13058030, 13058047, 13058048, 13059003, 13059004, 13059005, 13059006, 13059007, 13059008, 13059009, 13059010, 13059019, 13059035, 13059082, 13059083, 13059084, 13059085, 13059086, 13059087, 13060001, 13060002, 13060003, 13060009, 13060012, 13061026, 13063033, 13064030, 13064057, 13064059, 13064074, 13066035, 13066036, 13066101, 13066102, 13066104, 13066109, 13066110, 13067001, 13067002, 13067003, 13067004, 13067005, 13067006, 13067007, 13067008, 13067009, 13067010, 13067011, 13067012, 13067013, 13067014, 13067015, 13067017, 13068017, 13068022, 13068027, 13068029, 13068034, 13068036, 13068037, 13069006, 13069009, 13069029, 13070030, 13070056, 13071034, 13071037, 13071038, 13071040};

	// loop over events
	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
	int processedevent = 0;
	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		t->GetEntry(ievt);

		if(j1neutralfrac>0.9) continue;			// jet neutral pT fraction < 90%
		
		if(fabs(jeta)>0.6) continue;			// I miss this cut in jet finding selector

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

		leadjetpt->Fill(xvariable);	
		
        	leadjetntrkvsleadjetpt->Fill(xvariable,leadntrk);
        	subjetntrkvsleadjetpt->Fill(xvariable,subntrk);
        	tranntrkvsleadjetpt->Fill(xvariable,(tranmaxntrk+tranminntrk)/2.);		// do not use tranntrk directly because it is int type and it is the sum divided by 2 
		MaxOrMin(tranmaxntrk,tranminntrk);
        	tranmaxntrkvsleadjetpt->Fill(xvariable,tranmaxntrk);
        	tranminntrkvsleadjetpt->Fill(xvariable,tranminntrk);

		float maxtranspt = max(*max_element(pt_min,pt_min+tranminntrk),*max_element(pt_max,pt_max+tranmaxntrk));	
		maxtranptvsleadjetpt->Fill(xvariable,maxtranspt);

        	hleadjetntrkvsleadjetpt->Fill(xvariable,leadntrk);
        	hsubjetntrkvsleadjetpt->Fill(xvariable,subntrk);
        	htranntrkvsleadjetpt->Fill(xvariable,(tranmaxntrk+tranminntrk)/2.);
		MaxOrMin(tranmaxntrk,tranminntrk);
        	htranmaxntrkvsleadjetpt->Fill(xvariable,tranmaxntrk);
        	htranminntrkvsleadjetpt->Fill(xvariable,tranminntrk);

		hmaxtranptvsleadjetpt->Fill(xvariable,maxtranspt);

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
			w = getweight(pt_max[it],chargeflag,MCflag);
			//w = inveff_tpc(pt_max[it],chargeflag,MCflag);
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
			w = getweight(pt_min[it],chargeflag,MCflag);
			//w = inveff_tpc(pt_min[it],chargeflag,MCflag);
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
			w = getweight(pt_jet[it],chargeflag,MCflag);
			//w = inveff_tpc(pt_jet[it],chargeflag,MCflag);
			//cout<<" w = "<<w<<endl;
			leadjetptavevsleadjetpt->Fill(xvariable,pt_jet[it],w);

			hleadjetptavevsleadjetpt->Fill(xvariable,pt_jet[it],w);

			sumleadpt+=pt_jet[it]*w;		// will be w tracks effectively ..
		}
		for(int it = 0; it<subntrk; it++) {
			//cout<<"pt_sub["<<it<<"] = "<<pt_sub[it];
			w = getweight(pt_sub[it],chargeflag,MCflag);
			//w = inveff_tpc(pt_sub[it],chargeflag,MCflag);
			//cout<<" w = "<<w<<endl;
			subjetptavevsleadjetpt->Fill(xvariable,pt_sub[it],w);

			hsubjetptavevsleadjetpt->Fill(xvariable,pt_sub[it],w);

			sumsubpt+=pt_sub[it]*w;		// will be w tracks effectively ..
		}

		sumtranpt=sumtranpt/2.;		// take the average 
        	leadjetptsumvsleadjetpt->Fill(xvariable,sumleadpt);
        	subjetptsumvsleadjetpt->Fill(xvariable,sumsubpt);
        	tranptsumvsleadjetpt->Fill(xvariable,sumtranpt);
		MaxOrMin(sumtranmaxpt,sumtranminpt);
        	tranmaxptsumvsleadjetpt->Fill(xvariable,sumtranmaxpt);
        	tranminptsumvsleadjetpt->Fill(xvariable,sumtranminpt);

        	hleadjetptsumvsleadjetpt->Fill(xvariable,sumleadpt);
        	hsubjetptsumvsleadjetpt->Fill(xvariable,sumsubpt);
        	htranptsumvsleadjetpt->Fill(xvariable,sumtranpt);
		MaxOrMin(sumtranmaxpt,sumtranminpt);
        	htranmaxptsumvsleadjetpt->Fill(xvariable,sumtranmaxpt);
        	htranminptsumvsleadjetpt->Fill(xvariable,sumtranminpt);

		processedevent++;
	}
	cout<<"Total processed # of event "<<processedevent<<endl;

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

        maxtranptvsleadjetpt->GetXaxis()->SetTitle(xvariablename);
	maxtranptvsleadjetpt->GetYaxis()->SetTitle("Maximum Transverse Pt");
	maxtranptvsleadjetpt->SetLineColor(ctran);
	maxtranptvsleadjetpt->SetMarkerStyle(stran);
	maxtranptvsleadjetpt->SetMarkerColor(ctran);


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
	if(savefig) c[0]->SaveAs(Form("figs/%sVs%s_%s.png",leadjetpt->GetName(),what2fill.Data(),filetag.Data()));
	
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
	if(savefig) c[1]->SaveAs(Form("figs/%sVs%s_%s.png","Ntrk",what2fill.Data(),filetag.Data()));

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
	if(savefig) c[2]->SaveAs(Form("figs/%sVs%s_%s.png","PtSum",what2fill.Data(),filetag.Data()));
	
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
	if(savefig) c[3]->SaveAs(Form("figs/%sVs%s_%s.png","PtAve",what2fill.Data(),filetag.Data()));


	if(saveroot) {
		//TFile *fout = new TFile("TwoHisto4underlyingevent_JP2_R06_LeadJetAngle_MatchTrig_151110.root","RECREATE");
		//TString outfilepath = dir+what2fill+"hist4"+filetag+".root";
		TString jettag="";
		if( ((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase))||(what2fill.Contains("transntrk",TString::kIgnoreCase))||(what2fill.Contains("tranntrk",TString::kIgnoreCase))) && (fabs(jetptmin-10)>1e-6||fabs(jetptmax-200)>1e-6)) {
			jettag = Form("_jet%g-%g",jetptmin,jetptmax);	
		} 
		TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+".root";
		TFile *fout = new TFile(outfilepath,"RECREATE");
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

		maxtranptvsleadjetpt->Write();

	
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

		hmaxtranptvsleadjetpt->Write();
			
		fout->Close();
	}

}


