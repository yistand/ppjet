//===============================================================================================================
//
//	2016.01.20	Li Yi
//	read from Pico data, plot <pT> distribution as function of mulitplicity or selected variables
//
//===============================================================================================================


#include <iostream>

#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

#include "Jet.C"

using namespace std;


double getphi(double px, double py) {
	
	double phi = ((py==0)?0:atan(py/px)) ;

	if(px<0&&py<0) phi-=TMath::Pi();
	if(px<0&&py>0) phi+=TMath::Pi();
		
	return phi;
}


double geteta(double px, double py, double pz) {

	double p = sqrt(px*px+py*py+pz*pz);

	if(fabs(p-pz)<1e-12) return 1000;

	double e = 0.5*log((p+pz)/(p-pz));

	return e;
}

double getpt(double px, double py) {
	return sqrt(px*px+py*py);
}

float inveff_pion(float pt) {                                                         
                                                                                      
        if(pt<0.1) return 0;                                                          
        
        float eff = 0;
        TF1* feff=new TF1("feff","[0]*(exp(-pow([1]/x,[2])))", 0.1, 4.5);             
        feff->SetParameters(0.874739, 0.156624, 5.67316);                             
        eff = feff->Eval(pt);
                                                                                      
        if(pt>4.5) eff = feff->Eval(4.5);
                                                                                      
        delete feff;            // prevent memory leak                                
        
        if(eff>0) return 1./eff;                                                      
        else return 0;                                                                

}


void plotpT_fromJetTree() {

        TChain *chain = new TChain("JetTree");
	//chain->Add("/projects/rhig/ly247/stardata/pp200Y12picoMB_151113/sum*root");
	//chain->Add("/home/ly247/data/pp200Y12picoJP2_151030/*root");
	//chain->Add("/projects/rhig/ly247/stardata/pp200Y12picoMB_151207/*_*.root");
	chain->Add("/home/hep/caines/ly247/Scratch/pp12MBPico_151207/*.root");

	Jet *t = new Jet(chain);

	const int NoZdc = 3;
	const double ZDChigh[NoZdc] = {4000,9000,14000};
	const double ZDClow[NoZdc] = {2000,7000,12000};

	TH2D *hptVsmult[NoZdc]; 
	TH2D *hgptVsmult[NoZdc]; 
	TH2D *hrefmutVsNtracks[NoZdc]; 

	TH2D *hptVsrefmult[NoZdc]; 
	TH2D *hgptVsrefmult[NoZdc]; 

	TH2D *hptVsbemcmatchmult[NoZdc];
	TH2D *hptVsbemcmatchrefmult[NoZdc];
	TH2D *hbemcmatchrefmultVsrefmult[NoZdc];
	TH2D *hbemcmatchNtracksVsNtracks[NoZdc];

	TH2D *hptVsbtofmatchmult[NoZdc];
	TH2D *hptVsbtofmatchrefmult[NoZdc];
	TH2D *hbtofmatchrefmultVsrefmult[NoZdc];
	TH2D *hbtofmatchNtracksVsNtracks[NoZdc];

	TH2D *hptVsbothmatchmult[NoZdc];
	TH2D *hptVsbothmatchrefmult[NoZdc];
	TH2D *hbothmatchrefmultVsrefmult[NoZdc];
	TH2D *hbothmatchNtracksVsNtracks[NoZdc];

	TH2D *hptVsnotower[NoZdc];
	TH2D *hnotowerVsbemcmatchNtracks[NoZdc];
	TH2D *hnotowerVsbtofmatchNtracks[NoZdc];
	TH2D *hbemcmatchNtracksVsbtofmatchNtracks[NoZdc];
	TH2D *hbemcmatchNtracksVsbothmatchNtracks[NoZdc];
	TH2D *hbtofmatchNtracksVsbothmatchNtracks[NoZdc];


	for(int in = 0; in<NoZdc ; in++) {
		hptVsmult[in] = new TH2D(Form("hptVsmult%d",in),Form("pT vs Mulitplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hgptVsmult[in] = new TH2D(Form("hgptVsmult%d",in),Form("global pT vs Mulitplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hrefmutVsNtracks[in] = new TH2D(Form("hrefmultVsNtracks%d",in),Form("refmult vs Ntracks %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);

		hptVsrefmult[in] = new TH2D(Form("hptVsrefmult%d",in),Form("pT vs refMult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hgptVsrefmult[in] = new TH2D(Form("hgptVsrefmult%d",in),Form("global pT vs refMult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);

		hptVsbemcmatchmult[in] = new TH2D(Form("hptVsbemcmatchmult%d",in),Form("pT vs bemcmatch Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hptVsbemcmatchrefmult[in] = new TH2D(Form("hptVsbemcmatchrefmult%d",in),Form("pT vs bemcmatch refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hbemcmatchrefmultVsrefmult[in] = new TH2D(Form("hbemcmatchrefmultVsrefmult%d",in),Form("Bemcmatch refmult vs refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbemcmatchNtracksVsNtracks[in] = new TH2D(Form("hbemcmatchNtracksVsNtracks%d",in),Form("Bemcmatch Multiplicity vs Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);

		hptVsbtofmatchmult[in] = new TH2D(Form("hptVsbtofmatchmult%d",in),Form("pT vs btofmatch Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hptVsbtofmatchrefmult[in] = new TH2D(Form("hptVsbtofmatchrefmult%d",in),Form("pT vs btofmatch refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hbtofmatchrefmultVsrefmult[in] = new TH2D(Form("hbtofmatchrefmultVsrefmult%d",in),Form("Btofmatch refmult vs refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbtofmatchNtracksVsNtracks[in] = new TH2D(Form("hbtofmatchNtracksVsNtracks%d",in),Form("Btofmatch Multiplicity vs Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);

		hptVsbothmatchmult[in] = new TH2D(Form("hptVsbothmatchmult%d",in),Form("pT vs bothmatch Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hptVsbothmatchrefmult[in] = new TH2D(Form("hptVsbothmatchrefmult%d",in),Form("pT vs bothmatch refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hbothmatchrefmultVsrefmult[in] = new TH2D(Form("hbothmatchrefmultVsrefmult%d",in),Form("bothmatch refmult vs refmult %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbothmatchNtracksVsNtracks[in] = new TH2D(Form("hbothmatchNtracksVsNtracks%d",in),Form("bothmatch Multiplicity vs Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);

		hptVsnotower[in] = new TH2D(Form("hptVsnotower%d",in),Form("pT vs notower %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,1000,0,25);
		hnotowerVsbemcmatchNtracks[in] = new TH2D(Form("hnotowerVsbemcmatchNtracks%d",in),Form("Number of Towers vs Bemcmatch Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hnotowerVsbtofmatchNtracks[in] = new TH2D(Form("hnotowerVsbtofmatchNtracks%d",in),Form("Number of Towers vs Btofmatch Multiplicity %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbemcmatchNtracksVsbtofmatchNtracks[in] = new TH2D(Form("hbemcmatchNtracksVsbtofmatchNtracks%d",in),Form("Bemcmatch Vs Btofmatch %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbemcmatchNtracksVsbothmatchNtracks[in] = new TH2D(Form("hbemcmatchNtracksVsbothmatchNtracks%d",in),Form("Bemcmatch Vs Bothmatch %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
		hbtofmatchNtracksVsbothmatchNtracks[in] = new TH2D(Form("hbtofmatchNtracksVsbothmatchNtracks%d",in),Form("Btofmatch Vs Bothmatch %g<zdc<%g",ZDClow[in],ZDChigh[in]),50,0,50,50,0,50);
	}

	float zdc = 0;
	int refmult = 0;
	int bemcmatchrefmult = 0;
	int btofmatchrefmult = 0;
	int bothmatchrefmult = 0;
	int Ngoodtracks = 0; 
	int Nbemcmatchtracks = 0;
	int Nbtofmatchtracks = 0;
	int Nbothmatchtracks = 0;
	int runnumber = 0;
	double pt[1000] = {0};
	double gpt[1000] = {0};
	double meanpt;
	double peta[1000] = {0};
	double pphi[1000] = {0};
	double pdca[1000] = {0};

	TTree *tree = new TTree("MeanPtTree","recreate");
	tree->Branch("zdc",&zdc,"zdc/F");
	tree->Branch("btofmatch",&Nbtofmatchtracks,"btofmatch/I");
	tree->Branch("bemcmatch",&Nbemcmatchtracks,"bemcmatch/I");
	tree->Branch("bothmatch",&Nbothmatchtracks,"bothmatch/I");
	tree->Branch("refmult",&refmult,"refmult/I");
	tree->Branch("Ngoodtracks",&Ngoodtracks,"Ngoodtracks/I");
	tree->Branch("runnumber",&runnumber,"runnumber/I");
	tree->Branch("meanpt",&meanpt,"meanpt/D");
	tree->Branch("pt",&pt,"pt[Ngoodtracks]/D");
	tree->Branch("pphi",&pphi,"pphi[Ngoodtracks]/D");
	tree->Branch("peta",&peta,"peta[Ngoodtracks]/D");
	tree->Branch("pdca",&pdca,"pdca[Ngoodtracks]/D");
	

	cout<<"Events: "<<chain->GetEntries()<<endl;
	for(int ievt = 0 ; ievt<chain->GetEntries() ; ievt++) {
		t->GetEntry(ievt);
		if(ievt%100000==0) cout<<ievt<<endl;

		//if(ievt>1000000) break;	// test


		// event cuts
		//if(t->mZDCCointcidenceRate>5500) continue;

		if(fabs(t->fEventHeader_fPVz)>30) continue;
		if(fabs(t->fEventHeader_fPVz - t->fEventHeader_fvpdVz)>3) continue;
		if(t->fEventHeader_fRank<0) continue;
		int triggered = 0;
		for(int itrg = 0; itrg<t->fEventHeader_fNOfTriggerIds; itrg++) {
			if(t->fEventHeader_fTriggerIdArray[itrg]==370011) triggered = 1;		// VPDMB_nobsmd
			//if(t->mTrigId[itrg]==370341) triggered = 1;		// TOFMult4
			//if(t->mTrigId[itrg]==370361) triggered = 1;		// TOFMult3*VPD
		}
		if(triggered==0) continue;	

		runnumber = t->fEventHeader_fRunId;

		Ngoodtracks = 0;
		Nbemcmatchtracks = 0;
		Nbtofmatchtracks = 0;
		Nbothmatchtracks = 0;
		refmult = 0;
		bemcmatchrefmult = 0;
		btofmatchrefmult = 0;
		bothmatchrefmult = 0;
		meanpt = 0;
		for(int j = 0; j<t->fPrimaryTracks_; j++) {

			// track cuts
			double ieta = geteta(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j],t->fPrimaryTracks_fPz[j]);
			double ipt = getpt(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double iphi = getphi(t->fPrimaryTracks_fPx[j],t->fPrimaryTracks_fPy[j]);
			double idca = fabs(t->fPrimaryTracks_fDCA[j]);
			

			if((t->fPrimaryTracks_fFlag[j]<=0)) continue;	// test
			//cout<<j<<endl;

			if((t->fPrimaryTracks_fBemcMatchFlag[j]>0)) {
				Nbemcmatchtracks++;
				//cout<<"bemc :"<<Nbemcmatchtracks<<endl;

				if(fabs(ieta)<0.5) bemcmatchrefmult++;			
			}

			if((t->fPrimaryTracks_fTofMatchFlag[j]>0)) {
			//if((t->fPrimaryTracks_fTofBeta[j]>0)) {	// test
				Nbtofmatchtracks++;
				//cout<<"btof :"<<Nbtofmatchtracks<<endl;

				if(fabs(ieta)<0.5) btofmatchrefmult++;			
			}

			if((t->fPrimaryTracks_fTofMatchFlag[j]>0)&&(t->fPrimaryTracks_fBemcMatchFlag[j]>0)) {
			//if((t->fPrimaryTracks_fTofBeta[j]>0)&&(t->fPrimaryTracks_fBemcMatchFlag[j]>0)) {	// test
				Nbothmatchtracks++;
				//cout<<"both :"<<Nbothmatchtracks<<endl;

				if(fabs(ieta)<0.5) bothmatchrefmult++;
			}	
			//cout<<"btof = "<<Nbtofmatchtracks<<"\tbemc = "<<Nbemcmatchtracks<<"\tboth = "<<Nbothmatchtracks<<endl<<endl;

			if(fabs(ieta)>1) continue;			
			if(ipt<0.2 && ipt>20 ) continue;
			if(idca>1) continue;
			if(t->fPrimaryTracks_fNFittedHits[j]<25) continue;
			if(1.*t->fPrimaryTracks_fNFittedHits[j]/t->fPrimaryTracks_fNHitsPoss[j]<0.52) continue;
			//test if((t->fPrimaryTracks_fTofMatchFlag[j]<=0)&&(t->fPrimaryTracks_fBemcMatchFlag[j]<=0)) continue;
			if((t->fPrimaryTracks_fTofMatchFlag[j]<=0)) continue;	// test
			//if((t->fPrimaryTracks_fBemcMatchFlag[j]<=0)) continue;	// test

			if(fabs(ieta)<0.5) refmult++;

			//if(t->mPt[j]<0.5) continue;
			//pt[Ngoodtracks] = t->mPt[j];
			pt[Ngoodtracks] = ipt;
			gpt[Ngoodtracks] = ipt;
			pphi[Ngoodtracks] = iphi;
			peta[Ngoodtracks] = ieta;
			pdca[Ngoodtracks] = idca;
			meanpt+=ipt;
			Ngoodtracks++;
			//hptVsmult->Fill(t->mBTOFMatch,t->mPt[j],inveff_pion(t->mPt[j]));
			//hptVsmult->Fill(t->mrefMult,t->mPt[j],inveff_pion(t->mPt[j]));
			//hptVsmult->Fill(t->mrefMult,t->mPt[j]);
			//hptVsmult->Fill(t->mBTOFMatch,t->mPt[j]);
		}

		//if(fabs(Nbtofmatchtracks-Nbemcmatchtracks)>5) continue;	//test

		if(Ngoodtracks>0) meanpt/=Ngoodtracks;

		int notower = t->fEventHeader_fNOfTowers;		
		

		zdc = t->fEventHeader_fZdcCoincidenceRate;
		//if(zdc<5500) {	
		if(zdc<ZDChigh[0]&&zdc>ZDClow[0]) {	
			//cout<<"btof = "<<Nbtofmatchtracks<<"\tbemc = "<<Nbemcmatchtracks<<"\tboth = "<<Nbothmatchtracks<<endl;
			hrefmutVsNtracks[0]->Fill(Ngoodtracks,refmult);
			hbemcmatchrefmultVsrefmult[0]->Fill(refmult,bemcmatchrefmult);
			hbemcmatchNtracksVsNtracks[0]->Fill(Ngoodtracks,Nbemcmatchtracks);
			hnotowerVsbemcmatchNtracks[0]->Fill(Nbemcmatchtracks,notower);
			hbtofmatchrefmultVsrefmult[0]->Fill(refmult,btofmatchrefmult);
			hbtofmatchNtracksVsNtracks[0]->Fill(Ngoodtracks,Nbtofmatchtracks);
			hnotowerVsbtofmatchNtracks[0]->Fill(Nbtofmatchtracks,notower);
			hbemcmatchNtracksVsbtofmatchNtracks[0]->Fill(Nbtofmatchtracks,Nbemcmatchtracks);
			hbothmatchrefmultVsrefmult[0]->Fill(refmult,bothmatchrefmult);
			hbothmatchNtracksVsNtracks[0]->Fill(Ngoodtracks,Nbothmatchtracks);
			hbemcmatchNtracksVsbothmatchNtracks[0]->Fill(Nbothmatchtracks,Nbemcmatchtracks);
			hbtofmatchNtracksVsbothmatchNtracks[0]->Fill(Nbothmatchtracks,Nbtofmatchtracks);
		}
		//else if(zdc>10000) {
		else if(zdc>ZDClow[1]&&zdc<ZDChigh[1]) {
			hrefmutVsNtracks[1]->Fill(Ngoodtracks,refmult);
			hbemcmatchrefmultVsrefmult[1]->Fill(refmult,bemcmatchrefmult);
			hbemcmatchNtracksVsNtracks[1]->Fill(Ngoodtracks,Nbemcmatchtracks);
			hnotowerVsbemcmatchNtracks[1]->Fill(Nbemcmatchtracks,notower);
			hbtofmatchrefmultVsrefmult[1]->Fill(refmult,btofmatchrefmult);
			hbtofmatchNtracksVsNtracks[1]->Fill(Ngoodtracks,Nbtofmatchtracks);
			hnotowerVsbtofmatchNtracks[1]->Fill(Nbtofmatchtracks,notower);
			hbemcmatchNtracksVsbtofmatchNtracks[1]->Fill(Nbtofmatchtracks,Nbemcmatchtracks);
			hbothmatchrefmultVsrefmult[1]->Fill(refmult,bothmatchrefmult);
			hbothmatchNtracksVsNtracks[1]->Fill(Ngoodtracks,Nbothmatchtracks);
			hbemcmatchNtracksVsbothmatchNtracks[1]->Fill(Nbothmatchtracks,Nbemcmatchtracks);
			hbtofmatchNtracksVsbothmatchNtracks[1]->Fill(Nbothmatchtracks,Nbtofmatchtracks);
		}
		//else {
		else if(zdc>ZDClow[2]&&zdc<ZDChigh[2]) {
			hrefmutVsNtracks[2]->Fill(Ngoodtracks,refmult);
			hbemcmatchrefmultVsrefmult[2]->Fill(refmult,bemcmatchrefmult);
			hbemcmatchNtracksVsNtracks[2]->Fill(Ngoodtracks,Nbemcmatchtracks);
			hnotowerVsbemcmatchNtracks[2]->Fill(Nbemcmatchtracks,notower);
			hbtofmatchrefmultVsrefmult[2]->Fill(refmult,btofmatchrefmult);
			hbtofmatchNtracksVsNtracks[2]->Fill(Ngoodtracks,Nbtofmatchtracks);
			hnotowerVsbtofmatchNtracks[2]->Fill(Nbtofmatchtracks,notower);
			hbemcmatchNtracksVsbtofmatchNtracks[2]->Fill(Nbtofmatchtracks,Nbemcmatchtracks);
			hbothmatchrefmultVsrefmult[2]->Fill(refmult,bothmatchrefmult);
			hbothmatchNtracksVsNtracks[2]->Fill(Ngoodtracks,Nbothmatchtracks);
			hbemcmatchNtracksVsbothmatchNtracks[2]->Fill(Nbothmatchtracks,Nbemcmatchtracks);
			hbtofmatchNtracksVsbothmatchNtracks[2]->Fill(Nbothmatchtracks,Nbtofmatchtracks);
		}
		
		
		for(int j = 0; j<Ngoodtracks; j++) {		
			//hptVsmult->Fill(Ngoodtracks, pt[j]);
			//if(zdc<5500) {
			if(zdc<ZDChigh[0]&&zdc>ZDClow[0]) {	
				hptVsmult[0]->Fill(Ngoodtracks, pt[j]);
				hgptVsmult[0]->Fill(Ngoodtracks, gpt[j]);

				hptVsrefmult[0]->Fill(refmult, pt[j]);
				hgptVsrefmult[0]->Fill(refmult, gpt[j]);

				hptVsbemcmatchmult[0]->Fill(Nbemcmatchtracks, pt[j]);
				hptVsbemcmatchrefmult[0]->Fill(bemcmatchrefmult, pt[j]);

				hptVsbtofmatchmult[0]->Fill(Nbtofmatchtracks, pt[j]);
				hptVsbtofmatchrefmult[0]->Fill(btofmatchrefmult, pt[j]);

				hptVsbothmatchmult[0]->Fill(Nbothmatchtracks, pt[j]);
				hptVsbothmatchrefmult[0]->Fill(bothmatchrefmult, pt[j]);

				hptVsnotower[0]->Fill(notower,pt[j]);
			}
			//else if(zdc>10000)	{
			else if(zdc>ZDClow[1]&&zdc<ZDChigh[1]) {
				hptVsmult[1]->Fill(Ngoodtracks, pt[j]);
				hgptVsmult[1]->Fill(Ngoodtracks, gpt[j]);

				hptVsrefmult[1]->Fill(refmult, pt[j]);
				hgptVsrefmult[1]->Fill(refmult, gpt[j]);

				hptVsbemcmatchmult[1]->Fill(Nbemcmatchtracks, pt[j]);
				hptVsbemcmatchrefmult[1]->Fill(bemcmatchrefmult, pt[j]);

				hptVsbtofmatchmult[1]->Fill(Nbtofmatchtracks, pt[j]);
				hptVsbtofmatchrefmult[1]->Fill(btofmatchrefmult, pt[j]);

				hptVsbothmatchmult[1]->Fill(Nbothmatchtracks, pt[j]);
				hptVsbothmatchrefmult[1]->Fill(bothmatchrefmult, pt[j]);

				hptVsnotower[1]->Fill(notower,pt[j]);
			}
			//else {
			else if(zdc>ZDClow[2]&&zdc<ZDChigh[2]) {
				hptVsmult[2]->Fill(Ngoodtracks, pt[j]);
				hgptVsmult[2]->Fill(Ngoodtracks, gpt[j]);

				hptVsrefmult[2]->Fill(refmult, pt[j]);
				hgptVsrefmult[2]->Fill(refmult, gpt[j]);

				hptVsbemcmatchmult[2]->Fill(Nbemcmatchtracks, pt[j]);
				hptVsbemcmatchrefmult[2]->Fill(bemcmatchrefmult, pt[j]);

				hptVsbtofmatchmult[2]->Fill(Nbtofmatchtracks, pt[j]);
				hptVsbtofmatchrefmult[2]->Fill(btofmatchrefmult, pt[j]);

				hptVsbothmatchmult[2]->Fill(Nbothmatchtracks, pt[j]);
				hptVsbothmatchrefmult[2]->Fill(bothmatchrefmult, pt[j]);

				hptVsnotower[2]->Fill(notower,pt[j]);
			}
		}

		//cout<<"Nbtofmatchtracks = "<<Nbtofmatchtracks<<endl;	// test
		tree->Fill();

			
		//for(int k = 0; k<t->fEventHeader_fNOfTowers;k++) {
		//	hptVsmult->Fill(t->fEventHeader_fNOfMatchedTracks,t->fTowers_fEnergy[k]);
		//}
	}

	cout<<"loop done"<<endl;

	hptVsmult[0]->SetLineColor(3);
	hptVsmult[1]->SetLineColor(1);
	hptVsmult[2]->SetLineColor(2);
	hgptVsmult[0]->SetLineColor(3);
	hgptVsmult[1]->SetLineColor(1);
	hgptVsmult[2]->SetLineColor(2);
	hgptVsmult[0]->SetLineStyle(7);
	hgptVsmult[1]->SetLineStyle(7);
	hgptVsmult[2]->SetLineStyle(7);

	hptVsbemcmatchmult[0]->SetLineColor(3);
	hptVsbemcmatchmult[1]->SetLineColor(1);
	hptVsbemcmatchmult[2]->SetLineColor(2);
	hptVsbemcmatchmult[0]->SetMarkerStyle(8);
	hptVsbemcmatchmult[1]->SetMarkerStyle(8);
	hptVsbemcmatchmult[2]->SetMarkerStyle(8);

	hptVsbtofmatchmult[0]->SetLineColor(3);
	hptVsbtofmatchmult[1]->SetLineColor(1);
	hptVsbtofmatchmult[2]->SetLineColor(2);
	hptVsbtofmatchmult[0]->SetMarkerStyle(4);
	hptVsbtofmatchmult[1]->SetMarkerStyle(4);
	hptVsbtofmatchmult[2]->SetMarkerStyle(4);

	hptVsbothmatchmult[0]->SetLineColor(3);
	hptVsbothmatchmult[1]->SetLineColor(1);
	hptVsbothmatchmult[2]->SetLineColor(2);
	hptVsbothmatchmult[0]->SetMarkerStyle(28);
	hptVsbothmatchmult[1]->SetMarkerStyle(28);
	hptVsbothmatchmult[2]->SetMarkerStyle(28);

	hptVsnotower[0]->SetLineColor(3);
	hptVsnotower[1]->SetLineColor(1);
	hptVsnotower[2]->SetLineColor(2);
	hptVsnotower[0]->SetMarkerColor(3);
	hptVsnotower[1]->SetMarkerColor(1);
	hptVsnotower[2]->SetMarkerColor(2);
	hptVsnotower[0]->SetMarkerStyle(23);
	hptVsnotower[1]->SetMarkerStyle(23);
	hptVsnotower[2]->SetMarkerStyle(23);

	hptVsrefmult[0]->SetLineColor(3);
	hptVsrefmult[1]->SetLineColor(1);
	hptVsrefmult[2]->SetLineColor(2);
	hgptVsrefmult[0]->SetLineColor(3);
	hgptVsrefmult[1]->SetLineColor(1);
	hgptVsrefmult[2]->SetLineColor(2);
	hgptVsrefmult[0]->SetLineStyle(7);
	hgptVsrefmult[1]->SetLineStyle(7);
	hgptVsrefmult[2]->SetLineStyle(7);

	hptVsbemcmatchrefmult[0]->SetLineColor(3);
	hptVsbemcmatchrefmult[1]->SetLineColor(1);
	hptVsbemcmatchrefmult[2]->SetLineColor(2);
	hptVsbemcmatchrefmult[0]->SetMarkerStyle(8);
	hptVsbemcmatchrefmult[1]->SetMarkerStyle(8);
	hptVsbemcmatchrefmult[2]->SetMarkerStyle(8);

	hptVsbtofmatchrefmult[0]->SetLineColor(3);
	hptVsbtofmatchrefmult[1]->SetLineColor(1);
	hptVsbtofmatchrefmult[2]->SetLineColor(2);
	hptVsbtofmatchrefmult[0]->SetMarkerStyle(4);
	hptVsbtofmatchrefmult[1]->SetMarkerStyle(4);
	hptVsbtofmatchrefmult[2]->SetMarkerStyle(4);

	hptVsbothmatchrefmult[0]->SetLineColor(3);
	hptVsbothmatchrefmult[1]->SetLineColor(1);
	hptVsbothmatchrefmult[2]->SetLineColor(2);
	hptVsbothmatchrefmult[0]->SetMarkerStyle(28);
	hptVsbothmatchrefmult[1]->SetMarkerStyle(28);
	hptVsbothmatchrefmult[2]->SetMarkerStyle(28);

	TCanvas *c1 = new TCanvas();
	hptVsmult[0]->ProfileX("ptvsmult0")->Draw();
	hptVsmult[1]->ProfileX("ptvsmult1")->Draw("same");
	hptVsmult[2]->ProfileX("ptvsmult2")->Draw("same");
	//hgptVsmult[0]->ProfileX("gptvsmult0")->Draw("same");
	//hgptVsmult[1]->ProfileX("gptvsmult1")->Draw("same");
	//hgptVsmult[2]->ProfileX("gptvsmult2")->Draw("same");
	hptVsbemcmatchmult[0]->ProfileX("ptvsbemcmatchmult0")->Draw("same");
	hptVsbemcmatchmult[1]->ProfileX("ptvsbemcmatchmult1")->Draw("same");
	hptVsbemcmatchmult[2]->ProfileX("ptvsbemcmatchmult2")->Draw("same");
	hptVsnotower[0]->ProfileX("ptvsnotower0")->Draw("same");
	hptVsnotower[1]->ProfileX("ptvsnotower1")->Draw("same");
	hptVsnotower[2]->ProfileX("ptvsnotower2")->Draw("same");
	hptVsbtofmatchmult[0]->ProfileX("ptvsbtofmatchmult0")->Draw("same");
	hptVsbtofmatchmult[1]->ProfileX("ptvsbtofmatchmult1")->Draw("same");
	hptVsbtofmatchmult[2]->ProfileX("ptvsbtofmatchmult2")->Draw("same");
	hptVsbothmatchmult[0]->ProfileX("ptvsbothmatchmult0")->Draw("same");
	hptVsbothmatchmult[1]->ProfileX("ptvsbothmatchmult1")->Draw("same");
	hptVsbothmatchmult[2]->ProfileX("ptvsbothmatchmult2")->Draw("same");

	TCanvas *c2 = new TCanvas();
	hptVsrefmult[0]->ProfileX("ptvsrefmult0")->Draw();
	hptVsrefmult[1]->ProfileX("ptvsrefmult1")->Draw("same");
	hptVsrefmult[2]->ProfileX("ptvsrefmult2")->Draw("same");
	//hgptVsrefmult[0]->ProfileX("gptvsrefmult0")->Draw("same");
	//hgptVsrefmult[1]->ProfileX("gptvsrefmult1")->Draw("same");
	//hgptVsrefmult[2]->ProfileX("gptvsrefmult2")->Draw("same");
	hptVsbemcmatchrefmult[0]->ProfileX("ptvsbemcmatchrefmult0")->Draw("same");
	hptVsbemcmatchrefmult[1]->ProfileX("ptvsbemcmatchrefmult1")->Draw("same");
	hptVsbemcmatchrefmult[2]->ProfileX("ptvsbemcmatchrefmult2")->Draw("same");
	hptVsbtofmatchrefmult[0]->ProfileX("ptvsbtofmatchrefmult0")->Draw("same");
	hptVsbtofmatchrefmult[1]->ProfileX("ptvsbtofmatchrefmult1")->Draw("same");
	hptVsbtofmatchrefmult[2]->ProfileX("ptvsbtofmatchrefmult2")->Draw("same");
	hptVsbothmatchrefmult[0]->ProfileX("ptvsbothmatchrefmult0")->Draw("same");
	hptVsbothmatchrefmult[1]->ProfileX("ptvsbothmatchrefmult1")->Draw("same");
	hptVsbothmatchrefmult[2]->ProfileX("ptvsbothmatchrefmult2")->Draw("same");


	TFile *f = new TFile("ptVsmult_pp12MBfromJetTree.root","RECREATE");
	hptVsmult[0]->Write();
	hptVsmult[1]->Write();
	hptVsmult[2]->Write();
	hptVsbemcmatchmult[0]->Write();
	hptVsbemcmatchmult[1]->Write();
	hptVsbemcmatchmult[2]->Write();
	hptVsbtofmatchmult[0]->Write();
	hptVsbtofmatchmult[1]->Write();
	hptVsbtofmatchmult[2]->Write();
	hptVsbothmatchmult[0]->Write();
	hptVsbothmatchmult[1]->Write();
	hptVsbothmatchmult[2]->Write();
	hptVsnotower[0]->Write();
	hptVsnotower[1]->Write();
	hptVsnotower[2]->Write();
	hgptVsmult[0]->Write();
	hgptVsmult[1]->Write();
	hgptVsmult[2]->Write();
	hptVsrefmult[0]->Write();
	hptVsrefmult[1]->Write();
	hptVsrefmult[2]->Write();
	hptVsbemcmatchrefmult[0]->Write();
	hptVsbemcmatchrefmult[1]->Write();
	hptVsbemcmatchrefmult[2]->Write();
	hptVsbtofmatchrefmult[0]->Write();
	hptVsbtofmatchrefmult[1]->Write();
	hptVsbtofmatchrefmult[2]->Write();
	hptVsbothmatchrefmult[0]->Write();
	hptVsbothmatchrefmult[1]->Write();
	hptVsbothmatchrefmult[2]->Write();
	hgptVsrefmult[0]->Write();
	hgptVsrefmult[1]->Write();
	hgptVsrefmult[2]->Write();
	hrefmutVsNtracks[0]->Write();
	hrefmutVsNtracks[1]->Write();
	hrefmutVsNtracks[2]->Write();
	hbemcmatchrefmultVsrefmult[0]->Write();
	hbemcmatchrefmultVsrefmult[1]->Write();
	hbemcmatchrefmultVsrefmult[2]->Write();
	hbemcmatchNtracksVsNtracks[0]->Write();	
	hbemcmatchNtracksVsNtracks[1]->Write();	
	hbemcmatchNtracksVsNtracks[2]->Write();	
	hnotowerVsbemcmatchNtracks[0]->Write();
	hnotowerVsbemcmatchNtracks[1]->Write();
	hnotowerVsbemcmatchNtracks[2]->Write();
	hbtofmatchrefmultVsrefmult[0]->Write();
	hbtofmatchrefmultVsrefmult[1]->Write();
	hbtofmatchrefmultVsrefmult[2]->Write();
	hbtofmatchNtracksVsNtracks[0]->Write();	
	hbtofmatchNtracksVsNtracks[1]->Write();	
	hbtofmatchNtracksVsNtracks[2]->Write();	
	hnotowerVsbtofmatchNtracks[0]->Write();
	hnotowerVsbtofmatchNtracks[1]->Write();
	hnotowerVsbtofmatchNtracks[2]->Write();
	hbemcmatchNtracksVsbtofmatchNtracks[0]->Write();
	hbemcmatchNtracksVsbtofmatchNtracks[1]->Write();
	hbemcmatchNtracksVsbtofmatchNtracks[2]->Write();
	hbemcmatchNtracksVsbothmatchNtracks[0]->Write();
	hbemcmatchNtracksVsbothmatchNtracks[1]->Write();
	hbemcmatchNtracksVsbothmatchNtracks[2]->Write();
	hbtofmatchNtracksVsbothmatchNtracks[0]->Write();
	hbtofmatchNtracksVsbothmatchNtracks[1]->Write();
	hbtofmatchNtracksVsbothmatchNtracks[2]->Write();

	tree->Write();

	f->Close();

}




