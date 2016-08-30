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
//
//		2016.08.30	Li Yi
//		Add option to use max pT track's phi as region reference instead of leading jet phi
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
#define MAXARRAY	1000
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


float foldphi(float phi) {			// fold to [0, pi]

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

int whichregion(float phi, float transize, float iphi) {		// 0: leading, 1: away, 2: tran
// phi: the direction of reference (leading jet or Max pT track)
// transize: transverse region phi size, usually it is 60 degree
// iphi: the phi direction the track for investigation its region
	
	if(foldphi(iphi-phi)<((180.-transize)/2.)/180.*TMath::Pi()) {		// leading (track) direction
		return 0;
	}
	else if(foldphi(iphi-phi)>((180.+transize)/2.)/180.*TMath::Pi()) {		// away direction
		return 1;
	}
	else return 2;
}

void printAll(int ntrk, const float *pt, const float *phi, const float *eta) {
	cout<<"Total "<<ntrk<<" particles:"<<endl;
	for(int i = 0; i<ntrk; i++) {
		cout<<"phi = "<<phi[i]<<" pt = "<<pt[i]<<" eta = "<<eta[i]<<endl;
	}
}

void RedistriRegions(float phi, float transize, int &leadntrk, int &subntrk, int &tranmaxntrk, int &tranminntrk, float *pt_jet, float *pt_sub, float *pt_max, float *pt_min, float *phi_jet, float *phi_sub, float *phi_max, float *phi_min, float *eta_jet, float *eta_sub, float *eta_max, float *eta_min) {		// redistribution jet, sub, tranmax, tranmin according to the input phi direction, this function will modify the array and size
// phi: is the axis to define jet. 
// transize: is the transverse region phi size, usually it is 60 degree. 

	int tmpleadntrk=0,  tmpsubntrk=0,  tmptranmaxntrk=0,  tmptranminntrk=0; 
	float tmpleadsumpt=0,  tmpsubsumpt=0,  tmptranmaxsumpt=0,  tmptranminsumpt=0; 
	float tmppt_min[MAXARRAY]={0}, tmppt_max[MAXARRAY]={0}, tmppt_jet[MAXARRAY]={0}, tmppt_sub[MAXARRAY]={0};
	float tmpeta_min[MAXARRAY]={0}, tmpeta_max[MAXARRAY]={0}, tmpeta_jet[MAXARRAY]={0}, tmpeta_sub[MAXARRAY]={0};
	float tmpphi_min[MAXARRAY]={0}, tmpphi_max[MAXARRAY]={0}, tmpphi_jet[MAXARRAY]={0}, tmpphi_sub[MAXARRAY]={0};

	for(int itr = 0; itr<leadntrk; itr++) {
		switch (whichregion(phi, transize, phi_jet[itr])) {
			case 0: 
				tmppt_jet[tmpleadntrk]=pt_jet[itr];
				tmpleadsumpt+=pt_jet[itr];
				tmpphi_jet[tmpleadntrk]=phi_jet[itr];
				tmpeta_jet[tmpleadntrk]=eta_jet[itr];
				tmpleadntrk++;
				break;
			case 1: 
				tmppt_sub[tmpsubntrk]=pt_jet[itr];
				tmpsubsumpt+=pt_jet[itr];
				tmpphi_sub[tmpsubntrk]=phi_jet[itr];
				tmpeta_sub[tmpsubntrk]=eta_jet[itr];
				tmpsubntrk++;
				break;
			case 2: 
				if(phi_jet[itr]-phi>0) {
					tmppt_max[tmptranmaxntrk]=pt_jet[itr];
					tmptranmaxsumpt+=pt_jet[itr];
					tmpphi_max[tmptranmaxntrk]=phi_jet[itr];
					tmpeta_max[tmptranmaxntrk]=eta_jet[itr];
					tmptranmaxntrk++;
					break;
				}
				else {
					tmppt_min[tmptranminntrk]=pt_jet[itr];
					tmptranminsumpt+=pt_jet[itr];
					tmpphi_min[tmptranminntrk]=phi_jet[itr];
					tmpeta_min[tmptranminntrk]=eta_jet[itr];
					tmptranminntrk++;
					break;
				}
		}
	}

	for(int itr = 0; itr<subntrk; itr++) {
		switch (whichregion(phi, transize, phi_sub[itr])) {
			case 0: 
				tmppt_jet[tmpleadntrk]=pt_sub[itr];
				tmpleadsumpt+=pt_sub[itr];
				tmpphi_jet[tmpleadntrk]=phi_sub[itr];
				tmpeta_jet[tmpleadntrk]=eta_sub[itr];
				tmpleadntrk++;
				break;
			case 1: 
				tmppt_sub[tmpsubntrk]=pt_sub[itr];
				tmpsubsumpt+=pt_sub[itr];
				tmpphi_sub[tmpsubntrk]=phi_sub[itr];
				tmpeta_sub[tmpsubntrk]=eta_sub[itr];
				tmpsubntrk++;
				break;
			case 2: 
				if(phi_sub[itr]-phi>0) {
					tmppt_max[tmptranmaxntrk]=pt_sub[itr];
					tmptranmaxsumpt+=pt_sub[itr];
					tmpphi_max[tmptranmaxntrk]=phi_sub[itr];
					tmpeta_max[tmptranmaxntrk]=eta_sub[itr];
					tmptranmaxntrk++;
					break;
				}
				else {
					tmppt_min[tmptranminntrk]=pt_sub[itr];
					tmptranminsumpt+=pt_sub[itr];
					tmpphi_min[tmptranminntrk]=phi_sub[itr];
					tmpeta_min[tmptranminntrk]=eta_sub[itr];
					tmptranminntrk++;
					break;
				}
		}
	}

	for(int itr = 0; itr<tranmaxntrk; itr++) {
		switch (whichregion(phi, transize, phi_max[itr])) {
			case 0: 
				tmppt_jet[tmpleadntrk]=pt_max[itr];
				tmpleadsumpt+=pt_max[itr];
				tmpphi_jet[tmpleadntrk]=phi_max[itr];
				tmpeta_jet[tmpleadntrk]=eta_max[itr];
				tmpleadntrk++;
				break;
			case 1: 
				tmppt_sub[tmpsubntrk]=pt_max[itr];
				tmpsubsumpt+=pt_max[itr];
				tmpphi_sub[tmpsubntrk]=phi_max[itr];
				tmpeta_sub[tmpsubntrk]=eta_max[itr];
				tmpsubntrk++;
				break;
			case 2: 
				if(phi_max[itr]-phi>0) {
					tmppt_max[tmptranmaxntrk]=pt_max[itr];
					tmptranmaxsumpt+=pt_max[itr];
					tmpphi_max[tmptranmaxntrk]=phi_max[itr];
					tmpeta_max[tmptranmaxntrk]=eta_max[itr];
					tmptranmaxntrk++;
					break;
				}
				else {
					tmppt_min[tmptranminntrk]=pt_max[itr];
					tmptranminsumpt+=pt_max[itr];
					tmpphi_min[tmptranminntrk]=phi_max[itr];
					tmpeta_min[tmptranminntrk]=eta_max[itr];
					tmptranminntrk++;
					break;
				}
		}
	}

	for(int itr = 0; itr<tranminntrk; itr++) {
		switch (whichregion(phi, transize, phi_min[itr])) {
			case 0: 
				tmppt_jet[tmpleadntrk]=pt_min[itr];
				tmpleadsumpt+=pt_min[itr];
				tmpphi_jet[tmpleadntrk]=phi_min[itr];
				tmpeta_jet[tmpleadntrk]=eta_min[itr];
				tmpleadntrk++;
				break;
			case 1: 
				tmppt_sub[tmpsubntrk]=pt_min[itr];
				tmpsubsumpt+=pt_min[itr];
				tmpphi_sub[tmpsubntrk]=phi_min[itr];
				tmpeta_sub[tmpsubntrk]=eta_min[itr];
				tmpsubntrk++;
				break;
			case 2: 
				if(phi_min[itr]-phi>0) {
					tmppt_max[tmptranmaxntrk]=pt_min[itr];
					tmptranmaxsumpt+=pt_min[itr];
					tmpphi_max[tmptranmaxntrk]=phi_min[itr];
					tmpeta_max[tmptranmaxntrk]=eta_min[itr];
					tmptranmaxntrk++;
					break;
				}
				else {
					tmppt_min[tmptranminntrk]=pt_min[itr];
					tmptranminsumpt+=pt_min[itr];
					tmpphi_min[tmptranminntrk]=phi_min[itr];
					tmpeta_min[tmptranminntrk]=eta_min[itr];
					tmptranminntrk++;
					break;
				}
		}
	}

	leadntrk = tmpleadntrk;
	std::copy(tmppt_jet, tmppt_jet+tmpleadntrk, pt_jet);		
	std::copy(tmpphi_jet, tmpphi_jet+tmpleadntrk, phi_jet);		
	std::copy(tmpeta_jet, tmpeta_jet+tmpleadntrk, eta_jet);		

	subntrk = tmpsubntrk;
	std::copy(tmppt_sub,tmppt_sub+tmpsubntrk,pt_sub);		
	std::copy(tmpphi_sub,tmpphi_sub+tmpsubntrk,phi_sub);		
	std::copy(tmpeta_sub,tmpeta_sub+tmpsubntrk,eta_sub);		

	//cout<<"transumpt "<<tmptranmaxsumpt<<" "<< tmptranminsumpt<<endl;
	if(tmptranmaxsumpt>=tmptranminsumpt) {
		tranmaxntrk = tmptranmaxntrk;
		std::copy(tmppt_max,tmppt_max+tmptranmaxntrk,pt_max);		
		std::copy(tmpphi_max,tmpphi_max+tmptranmaxntrk,phi_max);		
		std::copy(tmpeta_max,tmpeta_max+tmptranmaxntrk,eta_max);		

		tranminntrk = tmptranminntrk;
		std::copy(tmppt_min,tmppt_min+tmptranminntrk,pt_min);		
		std::copy(tmpphi_min,tmpphi_min+tmptranminntrk,phi_min);		
		std::copy(tmpeta_min,tmpeta_min+tmptranminntrk,eta_min);		
	}
	else {	
		tranmaxntrk = tmptranminntrk;
		std::copy(tmppt_min,tmppt_min+tmptranminntrk,pt_max);		
		std::copy(tmpphi_min,tmpphi_min+tmptranminntrk,phi_max);		
		std::copy(tmpeta_min,tmpeta_min+tmptranminntrk,eta_max);		

		tranminntrk = tmptranmaxntrk;
		std::copy(tmppt_max,tmppt_max+tmptranmaxntrk,pt_min);		
		std::copy(tmpphi_max,tmpphi_max+tmptranmaxntrk,phi_min);		
		std::copy(tmpeta_max,tmpeta_max+tmptranmaxntrk,eta_min);		
	}

	//cout<<"Reference phi = "<<phi<<" transize = "<<transize<<endl;
	//cout<<"Lead --- "<<endl;
	//printAll(leadntrk, pt_jet, phi_jet, eta_jet);
	//cout<<"Away --- "<<endl;
	//printAll(subntrk, pt_sub, phi_sub, eta_sub);
	//cout<<"TranMax --- "<<endl;
	//printAll(tranmaxntrk, pt_max, phi_max, eta_max);
	//cout<<"TranMin --- "<<endl;
	//printAll(tranminntrk, pt_min, phi_min, eta_min);

}



void plotTree2Histo(TString what2fill="multiplicity", TString dir="~/Scratch/pp200Y12_jetunderlying/", TString filetag = "underlyingevent_MB_R06_LeadJetAngle_FullJetFraclt90_160116", double jetptmin = 10, double jetptmax= 200, int ExclusiveEta = 0) {
// what2fill: refmult, leadjetpt, maxtrackpt, multiplicity, transntrk	(no space)
// ExclusiveEta: 
// 		==1: jet and underlying event in seperate eta region. for example: jet in [-0.6, 0), then undelrying event [0.6,1)
// 		==0: no requirement for underlying event and jet eta

	if( (!what2fill.EqualTo("refmult",TString::kIgnoreCase)) && (!what2fill.EqualTo("leadjetpt",TString::kIgnoreCase)) && (!what2fill.EqualTo("maxtrackpt",TString::kIgnoreCase)) && (!what2fill.EqualTo("multiplicity",TString::kIgnoreCase)) && (!(what2fill.EqualTo("transntrk",TString::kIgnoreCase)||(what2fill.EqualTo("tranntrk",TString::kIgnoreCase)))) ) {
  	  cout<<"ERR!! call plotTree2Histo(TString what2fill, TString filepath): what2fill should be \"refmult\", \"leadjetpt\", \"maxtrackpt\", \"multiplicity\", \"transntrk\" or \"tranntrk\"."<<endl;
  	  return;
  	}

	int MCflag=1;			// if you don't want to correct for TOF/TPC efficency at this point, for example, do it as part of unfolding procedure, you can just set MCflag=1. Then it will be treated as pure MC (not geant) and no correction will be applied. 
	// if you are working on embedding (MC + geant), you can apply corrections. 
	if(filetag.Contains("pythia",TString::kIgnoreCase)) 	{
		MCflag = 1;
	}
	if(MCflag==1) {
		cout<<endl<<"INFO: process as MC data: No efficiency correction will be applied"<<endl;
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
	int saveroot = 0; 

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

	float pt_min[MAXARRAY], pt_max[MAXARRAY], pt_jet[MAXARRAY], pt_sub[MAXARRAY];
	float eta_min[MAXARRAY], eta_max[MAXARRAY], eta_jet[MAXARRAY], eta_sub[MAXARRAY];
	float phi_min[MAXARRAY], phi_max[MAXARRAY], phi_jet[MAXARRAY], phi_sub[MAXARRAY];

	float maxtrackpt, maxtrackphi, maxtracketa;		// if Max pt track is used 

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
	
	if(what2fill.Contains("trackpt",TString::kIgnoreCase)||what2fill.Contains("maxtrack",TString::kIgnoreCase)) {
		t->SetBranchAddress("maxpt",&maxtrackpt);
		t->SetBranchAddress("maxphi",&maxtrackphi);
		t->SetBranchAddress("maxeta",&maxtracketa);
	}



	// define histograms
	double maxpt = 50; // for refmult range
	int nbinning = 50; // for refmult range
	if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
		maxpt = 100;
		nbinning = 100;
	}  
	else if(what2fill.Contains("trackpt",TString::kIgnoreCase)||what2fill.Contains("maxtrack",TString::kIgnoreCase)) {
		maxpt = 25;			// max track pt cut is 20 
		nbinning = 25;
	}
	TString xvariablename = "refmult";
	if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
		xvariablename = "Leading Jet p_{T}";
	}
	else if(what2fill.Contains("trackpt",TString::kIgnoreCase)||what2fill.Contains("maxtrack",TString::kIgnoreCase)) {
		xvariablename = "Max Track p_{T}";
	}
	else if(what2fill.Contains("multiplicity",TString::kIgnoreCase)) {
		xvariablename = "Total multiplicity";
	}
	else if(what2fill.Contains("transntrk",TString::kIgnoreCase)||what2fill.Contains("tranntrk",TString::kIgnoreCase)) {
		xvariablename = "Transverse multiplicity";
	}

	TH1D *xaxis = new TH1D(what2fill,xvariablename,nbinning,0,maxpt);
	xaxis->Sumw2();

	TH2D *hetaphi = new TH2D("hetaphi","#Delta#eta-#Delta#phi",100,0,TMath::Pi(),100,-2,2);
	hetaphi->Sumw2();

	TProfile *j1ptvsxaxis = new TProfile("j1ptvs"+what2fill,"Leading Jet p_{T} vs "+xvariablename,nbinning,0,maxpt);
	TProfile *jasptvsxaxis = new TProfile("jasptvs"+what2fill,"Recoild Jet p_{T} vs "+xvariablename,nbinning,0,maxpt);

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
	TH2D *hjasptvsxaxis = new TH2D("hjasptvs"+what2fill,"Recoil Jet p_{T} vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);
	
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

        TH2D *hmaxtranptvsxaxis = new TH2D("hmaxtranptvs"+what2fill,"Maximum Transverse Track/Tower Pt vs "+xvariablename,nbinning,0,maxpt,nbinning_tpt,0,max_tpt);

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

	TH1D *hphidiff = new TH1D("hphidiff","|Max pt Track phi - Jet phi|",100,0,3.15); // for xvariable is maxtrackpt case. 
	TH2D *hphidiffVsJetpT = new TH2D("hphidiffVsJetpT","|Max pt Track phi - Jet phi| vs jet pt",100,0,100,100,0,3.15); // for xvariable is maxtrackpt case. 
	TH2D *hphidiffVsMaxTrackpT = new TH2D("hphidiffVsMaxTrackpT","|Max pt Track phi - Jet phi| vs Max Track pt",20,0,20,100,0,3.15); // for xvariable is maxtrackpt case. 


	// loop over events
	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
	int processedevent = 0;
	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {	
		t->GetEntry(ievt);
		//if(ievt>1000) break;		// test

		if( !(what2fill.Contains("trackpt",TString::kIgnoreCase)&&filetag.Contains("MaxTrack",TString::kIgnoreCase)) ) {	//  apply the following cut if we are not using the max track pt files
			if(j1neutralfrac>0.9) continue;			// jet neutral pT fraction < 90%
			if(runid<=13046029) continue;		// hot BEMC stripe in MB dataset
			
			if(fabs(jeta)>0.6) continue;			// I miss this cut in jet finding selector in old files
			
			if(jpt<=0) continue;				// no jet. 2016.08.28 I didn't require non-ghost in jet finding selector: there were some events with jet pt 0, I suspected that is the reason. but anyway remove zero jet events for futher calculation

			//if(fabs(jaspt)<5) continue;			
			//if(fabs(jpt-jaspt)>0.05*jpt) continue;		// dijet balanced

		}

		if(((what2fill.Contains("refmult",TString::kIgnoreCase))||(what2fill.Contains("multiplicity",TString::kIgnoreCase))||(what2fill.Contains("transntrk",TString::kIgnoreCase))||(what2fill.Contains("tranntrk",TString::kIgnoreCase)))&&(jpt<jetptmin||jpt>=jetptmax)) continue;	// when studying multiplicity dependence, exclude low pT jet to ensure the real jet found 

		if(ievt%1000000==0) cout<<"event "<<ievt<<endl;

		double xvariable = 0;
		float highesttrackphi = -999;	// for xvariable is maxtrackpt case and the file uses jetpt need to reassign the regions
		if(what2fill.Contains("jetpt",TString::kIgnoreCase)) {
			xvariable = jpt;
		}
		else if(what2fill.Contains("trackpt",TString::kIgnoreCase)) {	// max track pt
			if(filetag.Contains("MaxTrack",TString::kIgnoreCase)) {
				xvariable = maxtrackpt;
				highesttrackphi = maxtrackphi;
			}
			else {	// which means we are try to use track phi with max pt as reference, but we are reading from a file with tree using jet phi as reference. so some work needed to redistribute the regions
				// locate the max pt track 
				int LocMaxjettrackpt = TMath::LocMax(leadntrk,pt_jet);
				int LocMaxsubtrackpt = TMath::LocMax(subntrk,pt_sub);
				int LocMaxtranmaxtrackpt = TMath::LocMax(tranmaxntrk,pt_max);
				int LocMaxtranmintrackpt = TMath::LocMax(tranminntrk,pt_min);
				//cout<<"ntrk: "<<leadntrk<<" "<<subntrk<<" "<<tranmaxntrk<<" "<<tranminntrk<<endl;
				//cout<<"LocMax: "<<LocMaxjettrackpt<<" "<<LocMaxsubtrackpt<<" "<<LocMaxtranmaxtrackpt<<" "<<LocMaxtranmintrackpt<<endl;
				//cout<<"LocMaxPt: "<<(LocMaxjettrackpt>=0?pt_jet[LocMaxjettrackpt]:0)<<" "<< (LocMaxsubtrackpt>=0?pt_sub[LocMaxsubtrackpt]:0)<<" "<< (LocMaxtranmaxtrackpt>=0?pt_max[LocMaxtranmaxtrackpt]:0)<<" "<< (LocMaxtranmintrackpt>=0?pt_min[LocMaxtranmintrackpt]:0)<<endl;
				float tmp[4] = {LocMaxjettrackpt>=0?pt_jet[LocMaxjettrackpt]:0, LocMaxsubtrackpt>=0?pt_sub[LocMaxsubtrackpt]:0, LocMaxtranmaxtrackpt>=0?pt_max[LocMaxtranmaxtrackpt]:0, LocMaxtranmintrackpt>=0?pt_min[LocMaxtranmintrackpt]:0}; 
				int LocMaxalltrack = TMath::LocMax(4,tmp);
				if(LocMaxalltrack==0) {//if(LocMaxjettrackpt==-1) {//cout<<"jet pt="<<jpt<<" phi="<<jphi<<" eta="<<jeta<<endl;printAll(leadntrk,pt_jet,phi_jet,eta_jet); printAll(subntrk,pt_sub,phi_sub,eta_sub);printAll(tranmaxntrk, pt_max, phi_max, eta_max);printAll(tranminntrk, pt_min, phi_min, eta_min);} 
					xvariable = pt_jet[LocMaxjettrackpt]; highesttrackphi = phi_jet[LocMaxjettrackpt]; 
				}
				else if(LocMaxalltrack==1) {xvariable = pt_sub[LocMaxsubtrackpt]; highesttrackphi = phi_sub[LocMaxsubtrackpt]; }
				else if(LocMaxalltrack==2) {
					xvariable = pt_max[LocMaxtranmaxtrackpt]; highesttrackphi = phi_max[LocMaxtranmaxtrackpt]; //cout<<"Caution!!! Highest Transverse track/tower pT in TransMax "<<"pt="<<xvariable<<" phi = "<<phi_max[LocMaxtranmaxtrackpt]<<endl;cout<<"jetpt = "<<jpt<<" jetphi = "<<jphi<<endl;printAll(leadntrk,pt_jet,phi_jet,eta_jet);
				}
				else if(LocMaxalltrack==3) {
					xvariable = pt_min[LocMaxtranmintrackpt]; highesttrackphi = phi_min[LocMaxtranmintrackpt]; //cout<<"Highest Transverse track/tower pT in TransMin "<<"pt="<<xvariable<<endl;
				}
				else {
					cout<<"Something is not right here to find the max pt tracks"<<endl;
					break;
				}

				//cout<<"Max track pt = "<<xvariable<<" at phi = "<<highesttrackphi<<endl;
				// redistribution jet, sub, tranmax, tranmin according to new phi direction defined by max track pt 
				// Caution, this will rewrite the old region difintion. Don't expect to go back to the old one/call one array after this.
				float transize = 60; 
				if(filetag.Contains("TranPhi30_",TString::kIgnoreCase)) transize = 30;
				RedistriRegions(highesttrackphi, transize, leadntrk, subntrk, tranmaxntrk, tranminntrk, pt_jet, pt_sub, pt_max, pt_min, phi_jet, phi_sub, phi_max, phi_min, eta_jet, eta_sub, eta_max, eta_min);

			}
				
			float diffphi = foldphi(highesttrackphi-jphi);
			if(jphi>-999) {			// leading jet exists
				hphidiff->Fill(diffphi);
				hphidiffVsJetpT->Fill(jpt, diffphi);
				hphidiffVsMaxTrackpT->Fill(xvariable, diffphi);
			}
			else {
				hphidiff->Fill(-999.);			// put an entry in overflow bin
				hphidiffVsJetpT->Fill(0.,-999.);		// put an entry in overflow bin
				hphidiffVsMaxTrackpT->Fill(xvariable,-999.);	// put an entry in overflow bin
			}
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


		//cout<<"event "<<ievt<<" x="<<xvariable<<" ntrklead = "<<sumleadntrk<<" ntrksub = "<<sumsubntrk<<" ntrkMax = "<<sumtranmaxntrk<<" ntrkMin = "<<sumtranminntrk<<endl;

        	leadjetntrkvsxaxis->Fill(xvariable,sumleadntrk);
        	hleadjetntrkvsxaxis->Fill(xvariable,sumleadntrk);
        	subjetntrkvsxaxis->Fill(xvariable,sumsubntrk);
        	hsubjetntrkvsxaxis->Fill(xvariable,sumsubntrk);
        	tranntrkvsxaxis->Fill(xvariable,(sumtranmaxntrk+sumtranminntrk)/2.);		// do not use tranntrk directly because it is int type and it is the sum divided by 2 
        	htranntrkvsxaxis->Fill(xvariable,(sumtranmaxntrk+sumtranminntrk)/2.);

		// for multiplicity, the largest one with ntrk will be tranmax
		MaxOrMin(sumtranmaxntrk,sumtranminntrk);		
        	tranmaxntrkvsxaxis->Fill(xvariable,sumtranmaxntrk);
        	htranmaxntrkvsxaxis->Fill(xvariable,sumtranmaxntrk);
        	tranminntrkvsxaxis->Fill(xvariable,sumtranminntrk);
        	htranminntrkvsxaxis->Fill(xvariable,sumtranminntrk);


		// for track <pT>, the largest region with sum pT will be tranmax
		sumtranpt=sumtranpt/2.;		// take the average 
        	leadjetptsumvsxaxis->Fill(xvariable,sumleadpt);
        	subjetptsumvsxaxis->Fill(xvariable,sumsubpt);
        	tranptsumvsxaxis->Fill(xvariable,sumtranpt);
        	hleadjetptsumvsxaxis->Fill(xvariable,sumleadpt);
        	hsubjetptsumvsxaxis->Fill(xvariable,sumsubpt);
        	htranptsumvsxaxis->Fill(xvariable,sumtranpt);

		// for sum pt, the largest one with sum pt will be tranmax
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
	if(what2fill.Contains("trackpt",TString::kIgnoreCase)||what2fill.Contains("maxtrack",TString::kIgnoreCase))  drawxmax=21;
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
	leg->DrawClone("same");
	if(savefig) c[1]->SaveAs(Form("figs/%sVs%s_%s.png","Ntrk",what2fill.Data(),filetag.Data()));

	c[2]->cd();
	htmp[2] = new TH2D("htmp2","",1000,0,drawxmax,1000,0,30);//20);
	htmp[2]->GetXaxis()->SetTitle(xvariablename);
	htmp[2]->GetYaxis()->SetTitle("Event #LT#sum p_{T}#GT");
	htmp[2]->Draw();
	lat->DrawLatex(0.5,0.94,"Event #LT#sum p_{T}#GT");
	leadjetptsumvsxaxis->Draw("same");
	subjetptsumvsxaxis->Draw("same");
	tranmaxptsumvsxaxis->Draw("same");
	tranminptsumvsxaxis->Draw("same");	
	tranptsumvsxaxis->Draw("same");	
	leg->DrawClone("same");
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
	leg->DrawClone("same");
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
		TString outfilepath = dir+what2fill+"hist4"+filetag+jettag+".root";
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

		if(what2fill.Contains("trackpt",TString::kIgnoreCase)||what2fill.Contains("maxtrack",TString::kIgnoreCase)) {
			hphidiff->Write();
			hphidiffVsJetpT->Write();
			hphidiffVsMaxTrackpT->Write();
		}
			
		fout->Close();
	}

}


