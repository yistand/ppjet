//====================================================================================================
//
//		2016.11.29	Li Yi
//		Dihadron DeltaPhi-DeltaEta
//
//====================================================================================================



#include "AjParameters.hh"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TFile.h>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

#include "TClonesArray.h"

#include <utility>	// std::pair, std::make_pair
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <exception>
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>    // std::random_shuffle

#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"

#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"

#include "TStarJetPicoTriggerInfo.h"


using namespace std;

typedef std::pair<float, float> EtaPhiPair;

bool readinbadrunlist(std::set<int> & badrun, TString csvfile="./include/pp200Y12_badrun.list") {

	// open infile
	std::string line;
	std::ifstream inFile (csvfile );

	std::cout<<"Loading bad run id from "<< csvfile.Data()<<std::endl;;

	if ( !inFile.good() ) {
		std::cout<<"Can't open "<<csvfile.Data()<<std::endl;
		return false;
	}

	while (std::getline (inFile, line) ){
		if ( line.size()==0 ) continue; // skip empty lines
		if ( line[0] == '#' ) continue; // skip comments

		std::istringstream ss( line );
		while( ss ){
			std::string entry;
			std::getline( ss, entry, ',' );
			int ientry = atoi(entry.c_str());
			if (ientry) {
				badrun.insert( ientry );
				std::cout<<"Added bad runid "<<ientry<<std::endl;
			}
		}
	}

	return true;
}



// Helper to deal with repetitive stuff
TStarJetPicoReader SetupReader ( TChain* chain, TString TriggerString, const double RefMultCut ){
	//TStarJetPicoDefinitions::SetDebugLevel(10); // 10 for more output, 0 for less output

	TStarJetPicoReader reader;
	reader.SetInputChain (chain);

	// Event and track selection
	// -------------------------
	TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
	evCuts->SetTriggerSelection( TriggerString ); //All, MB, HT, pp, ppHT, ppJP
	// Additional cuts 
	evCuts->SetVertexZCut (AjParameters::VzCut);
	evCuts->SetRefMultCut ( RefMultCut );
	//pp evCuts->SetVertexZDiffCut( AjParameters::VzDiffCut );
	evCuts->SetVertexZDiffCut( 999999 );	// test pAu

	evCuts->SetMaxEventPtCut ( AjParameters::MaxEventPtCut );
	evCuts->SetMaxEventEtCut ( AjParameters::MaxEventEtCut );

	evCuts->SetPVRankingCut ( 0 );		// Vertex ranking > 0. Use SetPVRankingCutOff() to turn off vertex ranking cut.  default is OFF

	std::cout << "Exclude event with track > " << evCuts->GetMaxEventPtCut() << std::endl;
	std::cout << "Exclude event with tower > " << evCuts->GetMaxEventEtCut() << std::endl;

	// Tracks cuts
	TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
	trackCuts->SetDCACut( AjParameters::DcaCut );
	trackCuts->SetMinNFitPointsCut( AjParameters::NMinFit );
	trackCuts->SetFitOverMaxPointsCut( AjParameters::FitOverMaxPointsCut );
	trackCuts->SetMaxPtCut ( AjParameters::MaxTrackPt );

	std::cout << "Using these track cuts:" << std::endl;
	std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
	std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
	std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
	std::cout << " maxpt : " << trackCuts->GetMaxPtCut (  ) << std::endl;

	// Towers
	TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
	towerCuts->SetMaxEtCut(AjParameters::MaxEtCut);
	if(TriggerString.Contains("pAu")) {
		towerCuts->AddBadTowers("./include/pAu200Y15_hottower.list");		// #LY CHECK where is the bad tower list
	}
	else {		// Default is pp@200GeV Y12
		towerCuts->AddBadTowers("./include/pp200Y12_badtower.list");		// #LY CHECK where is the bad tower list
	}

	// Tower energy correction (subtract associated charged particle deposit energy). By default, it is MIP correction (comment out the following 3 lines)
	reader.SetApplyFractionHadronicCorrection(kTRUE);
	reader.SetFractionHadronicCorrection(0.9999);
	reader.SetRejectTowerElectrons( kFALSE );


	std::cout << "Using these tower cuts:" << std::endl;
	std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
	std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;

	// V0s: Turn off
	reader.SetProcessV0s(false);

	return reader;

}


int main ( int argc, const char** argv ) {

	// Set up some convenient default
	// ------------------------------
	const char *defaults[] = {"DPhiDEta","DiHadronpAutest.root","pAuVPDMB","~/Scratch/run15pAu/*.root", "11", "8", "500" };
	// {Code name, to be discard but needed since argv will use command name as the [0], output file name, triggername, intput file list, TPtNo, MultL, MultH}
	//


	if ( argc==1 ) {
		argv=defaults;
		argc=sizeof (defaults ) / sizeof (defaults[0] );
	}

	// Throw arguments in a vector
	// ---------------------------
	vector<string> arguments(argv + 1, argv + argc);

	TPtNo = atoi( arguments.at(3).data() );

	//Set Constants
	const float pi=TMath::Pi();
	const int NoAssocPt = 9+2+1;
	const float PtLow[NoAssocPt]={0.15,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,0.15,0.3,1};
	const float PtHigh[NoAssocPt]={0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,10.0,0.3,0.5,3};
	if(TPtNo>=NoAssocPt||TPtNo<0) { cout<<"TPtNo = "<<TPtNo<<" exceeds array"<<endl; exit(1);}
	const float TPtLow=PtLow[TPtNo];
	const float TPtHigh=PtHigh[TPtNo];


	// Initiate Histogram
	// --------------------
	TH2D *deta_dphi[NoAssocPt];                   // real events
	TH2D *mdeta_mdphi[NoAssocPt];                 // mix events
	TH1D *EtaDistribution[NoAssocPt];             // associate particle
	TH1D *pTDistribution[NoAssocPt];              // associate particle
	TH1D *TrigTot = new TH1D("TrigTot","Number of Trig.Particles",1,0,1);
	TH1D *TpTDistribution = new TH1D("TpTDistribution",Form("TpTDistribution pT %g-%g",TPtLow,TPtHigh),2000,0,20);
	TH1D *TEtaDistribution = new TH1D("TEtaDistribution",Form("TEtaDistribution pT %g-%g",TPtLow,TPtHigh),etabinning,-1,1);
	// for Vz Mix weight use
	TH1D *EventTriggerVz=new TH1D("EventTriggerVz","All Trigger Particle Vz distribution",10000,-50, 50);            // 2014.03.31 ly, will weighted mix event by the real/mix ratio later
	TH1D *MatchedTriggerVz=new TH1D("MatchedTriggerVz","Matched Trigger Particle Vz distribution",10000,-50,50);              // 2014.03.31 ly, will weighted mix event by the real/mix ratio later
	TH1D *VzDifference = new TH1D("VzDifference","Vz difference",1000,-10,10);  // 2014.04.01 ly for monitoring purpose

	// histogram format
	int phibinning = 96;//48;
	int etabinning = 80;//48;
	double detamin = -2, detamax = 2;
	double etamin = -1, etamax = 1;
	for(int i =0;i<NoAssocPt;i++) {
		deta_dphi[i]=new TH2D(Form("deta_dphi_pT%d",i),Form("Delta eta vs. delta phi: pT %.2f-%.2f ",PtLow[i],PtHigh[i]),etabinning,detamin,detamax,phibinning,-pi/2,3*pi/2);
		mdeta_mdphi[i]=new TH2D(Form("mdeta_mdphi_pT%d",i),Form("Mixed delta eta vs. delta phi: pT %.2f-%.2f ",PtLow[i],PtHigh[i]),etabinning,detamin,detamax,phibinning,-pi/2,3*pi/2);
		EtaDistribution[i]=new TH1D(Form("EtaDistribution_pT%d",i),Form("Eta Distribution for pT %.2f-%.2f ",PtLow[i],PtHigh[i]),etabinning,etamin,etamax);
		pTDistribution[i]=new TH1D(Form("pTDistribution_pT%d",i),Form("pT Distribution for pT %.2f-%.2f ",PtLow[i],PtHigh[i]),2000,0,20);
		deta_dphi[i]->Sumw2();
		mdeta_mdphi[i]->Sumw2();
		EtaDistribution[i]->Sumw2();
		pTDistribution[i]->Sumw2();
	}

	TrigTot->Sumw2();
	TpTDistribution->Sumw2();
	TEtaDistribution->Sumw2();

	// histogram for event info.

	int maxntrack = 200;
	int maxbbc = 50000, bbcbins = 200;

	TH2D *BBCEvsrefMult=new TH2D("BBCEvsrefMult","Total BBC East vs TPC refMult",maxntrack,0,maxntrack,bbcbins,0,maxbbc);
	BBCEvsrefMult->Sumw2();
	TH2D *BBCEvsNPrimaryTrk=new TH2D("BBCEvsNPrimaryTrk","Total BBC East vs TPC NPrimaryTrk",maxntrack,0,maxntrack,bbcbins,0,maxbbc);
	BBCEvsNPrimaryTrk->Sumw2();
	TH2D *BBCEvsNGlobalTrk=new TH2D("BBCEvsNGlobalTrk","Total BBC East vs TPC NGlobalTrk",maxntrack,0,5000,bbcbins,0,maxbbc);
	BBCEvsNGlobalTrk->Sumw2();

	TH2D *refMultvsNGlobalTrk=new TH2D("refMultvsNGlobalTrk","TPC refMult vs NGlobalTrk",maxntrack,0,5000,maxntrack,0,maxntrack);
	TH2D *refMultvsNtrk=new TH2D("refMultvsNtrk","TPC refMult vs Ntrk",maxntrack,0,maxntrack,maxntrack,0,maxntrack);
	refMultvsNPrimaryTrk->Sumw2();
	refMultvsNGlobalTrk->Sumw2();
	refMultvsNtrk->Sumw2();


	// Matching control
	//Max times matched for event
	Int_t MaxMatch = 10;          // 3;     
	Double_t VzMatch = 5;         // 1;     
	//Double_t VrMatch = 0.05;                      // 2014.09.29 ly
	int scratchy=0;                       // how many times overflow occurs and need to delte the old bin for current event recording

	//tirgger particle info in previous events for event mixing
	const int ArrayLength = 10000;         // max number events in history record
	const int MaxT = 100;                 // max number of trigger particle in one events
	Int_t HistTimes[ArrayLength] = {0};   // how many times been matched
	Int_t HistMult[ArrayLength] = {0};      // match Multiplicity for TPC, in ZDC difference less than 10
	Double_t HistVz[ArrayLength] = {0};      // will match Vz for mix events   2014.04.01 ly fix bugs. It was Int_t. 
	Double_t HistVx[ArrayLength] = {0};      // will match Vx for mix events   2014.09.30 ly for better mix event
	Double_t HistVy[ArrayLength] = {0};      // will match Vy for mix events   2014.09.30 ly for better mix event
	Int_t HistTotTrig[ArrayLength] = {0};
	Double_t HistTrigPhi[ArrayLength][MaxT] = {{0}};
	Double_t HistTrigEta[ArrayLength][MaxT] = {{0}};

	TString OutFileName = arguments.at(0);

	cout<<"TriggerName: "<<arguments.at(1)<<endl;
	TString TriggerName = arguments.at(1);

	// input tree
	cout<<"Chain data: "<<arguments.at(2).data()<<" for "<<ChainName<<endl;
	TChain* chain = new TChain( ChainName );
	chain->Add( arguments.at(2).data() );

	cout<<"SetupReader for pico"<<endl;
	double RefMultCut = 0;
	TStarJetPicoReader reader = SetupReader( chain, TriggerName, RefMultCut );			// #ly note: Events & Tracks & Towers cuts are set here
	//reader.SetTrackPileUpCut(kTRUE);		// #ly	tpc track matching to bemc or tof
	reader.SetTrackPileUpCut(1);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	if( OutFileName.Contains ("NoTofMatch") ) {
		reader.SetTrackPileUpCut(0);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcOrTofMatch") ) {
		reader.SetTrackPileUpCut(1);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcMatch") ) {
		reader.SetTrackPileUpCut(3);		// #ly	3: tpc track matching to bemc.		1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("TofMatch") ) {
		reader.SetTrackPileUpCut(2);		// #ly	3: tpc track matching to bemc.		1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	TStarJetPicoDefinitions::SetDebugLevel(0);

	// -------------------------
	// Cycle through events
	// --------------------
	TStarJetVectorContainer<TStarJetVector>* container;		// for underlying event loop
	TStarJetVector* sv; // TLorentzVector* would be sufficient. 


	Long64_t nEvents=-1; // -1 for all
	//nEvents=10000;	// test
	cout<<"init..."<<endl;
	reader.Init(nEvents);
	int count = 0;



	std::set<int>badrun;
	badrun.clear();

	if(OutFileName.Contains("pAu")) {
		readinbadrunlist(badrun,"./include/pAu200Y15_badrun.list");        
	}
	else {		// Default is pp@200GeV Y12
		readinbadrunlist(badrun);        
	}

	Int_t nTrigProcessed = 0;

	float MultL = atof(arguments.at(4).data());
	float MultH = atof(arguments.at(5).data());

	try{
		while ( reader.NextEvent() ) {
			reader.PrintStatus(10);
			if(count%10000==0) cout<<"event "<<count<<endl;
			count++;


			// event info
			// ----------
			//cout<<"load event header"<<endl;
			TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();

			// eventid = header->GetEventId();
			int runid   = header->GetRunId();
			if(badrun.count(runid)>0) continue;			// in bad run list
			//if(runid<13047000) continue;				// test beginning of MB run looks different

			//if(header->GetZdcCoincidenceRate()<6000 || header->GetZdcCoincidenceRate()>10000) continue;		// test


			Double_t Vz=header->GetPrimaryVertexZ();
			Double_t Vx=header->GetPrimaryVertexX();
			Double_t Vy=header->GetPrimaryVertexY();

			Double_t Phi=0;
			Double_t Pt=0;
			Double_t Eta=0;
			Int_t TotTrig=0;
			Double_t   TrigPt[MaxT]={0};
			Double_t   TrigPhi[MaxT]={0};
			Double_t   TrigEta[MaxT]={0};


			// "Centrality" cuts
			// ----------
			// Loop first time to get number of tracks matched to TOF or BEMC
			// ----------
			container = reader.GetOutputContainer();
			float tracks = 0;
			for (int ip = 0; ip<container->GetEntries() ; ++ip ){
				sv = container->Get(ip);  // Note that TStarJetVector contains more info, such as charge;

				if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

				if (sv->GetCharge()==0 ) continue;		// Charged particle only 
				tracks=tracks+1;
			}    

			if(tracks<MultL||tracks>MultH) continue;         // 2014.03.31 ly. so [MultL, MultH] now


			//============== record event info.
			int refMult = header->GetGReferenceMultiplicity();
			int NPrimaryTrk = header->GetNOfPrimaryTracks();
			int NGlobalTrk = header->GetNGlobalTracks();
			int BBCE = header->GetBbcAdcSumEast();

			BBCEvsrefMult->Fill(refMult,BBCE);
			BBCEvsNPrimaryTrk->Fill(NPrimaryTrk,BBCE);
			BBCEvsNGlobalTrk->Fill(NGlobalTrk,BBCE);

			refMultvsNtrk->Fill(Ntrk,refMult);
			refMultvsNPrimaryTrk->Fill(NPrimaryTrk,refMult);
			refMultvsNGlobalTrk->Fill(NGlobalTrk,refMult);

			// ----------
			// Loop second time to do the real analysis
			// ----------
			//=========================== real events
			for (int ip = 0; ip<container->GetEntries() ; ++ip ){
				sv = container->Get(ip);  // Note that TStarJetVector contains more info, such as charge;

				if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

				if (sv->GetCharge()==0 ) continue;		// Charged particle only 

				if (fabs(sv->eta())>1) continue;

				if (sv->perp()<TPtLow||sv->perp()>=TPtHigh)continue;

				TrigPt[TotTrig] = sv->perp();
				TrigPhi[TotTrig] = sv->phi_std();
				TrigEta[TotTrig] = sv->eta();

				//fill trigger particle histogram
				TEtaDistribution->Fill(TrigEta[TotTrig]);
				TpTDistribution->Fill(TrigPt[TotTrig]);


				////Begin Associate Loop for real events
				for (int ia = 0; ia<container->GetEntries() ; ++ia ){
					if(ip==ia) continue;		// self	

					sv = container->Get(ia);

					if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

					if (sv->GetCharge()==0 ) continue;		// Charged particle only 

					if (fabs(sv->eta())>1) continue;

					Pt = sv->perp();
					Phi = sv->phi_std();
					Eta = sv->eta();

					Float_t dPhi = Phi - TrigPhi[TotTrig];
					Float_t dEta = Eta - TrigEta[TotTrig];

					//Folding 
					if(dPhi<-pi/2) dPhi+=2*pi;            // ly 05.09
					if(dPhi>3*pi/2) dPhi-=2*pi;

					Int_t region=999;
					for(int iptr = 0;iptr<NoAssocPt;iptr++) {
						if(Pt<PtHigh[iptr]&&Pt>=PtLow[iptr]) {
							region = iptr;              // found pt bin

							if(TotTrig==0)        // fill associate particle info only once per event not per trigger particle
							{
								pTDistribution[region]->Fill(Pt);
								EtaDistribution[region]->Fill(Eta);

							}

							//fill real event \Delta\eta - \Delta\phi distribution
							deta_dphi[region]->Fill(dEta,dPhi);
						}
					}
				}//End of Associate Loop

				TotTrig++;

			}//End of trig loop

			// fill no. of trigger particle in event
			TrigTot->Fill(0.5,TotTrig);
			EventTriggerVz->Fill(Vz,TotTrig);               // 2014.03.31 ly  Trigger particle Vz dist. for all trigger particle read in
			nTrigProcessed+=TotTrig;

			//========================  mix events
			// Loop events in history
			for(int ihevt = 0; ihevt < ArrayLength; ihevt ++) {

				if(HistTimes[ihevt] > MaxMatch) continue;

				if(HistTotTrig[ihevt]<=0) continue;

				double diffVz = HistVz[ihevt] - Vz;
				if(fabs(diffVz) > VzMatch) continue;

				//cout<<"diffVz = "<<diffVz<<endl; // test
				VzDifference->Fill(diffVz);

				int flag = 0;           // 2014.03.31 ly whether any mixed pair

				// Begin History Trigger loop for mix events
				for(int iht = 0; iht<HistTotTrig[ihevt] ; iht++)
				{
					double weight4vz = 1;

					for (int ja = 0; ja<container->GetEntries() ; ++ja ){

						sv = container->Get(ja);

						if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

						if (sv->GetCharge()==0 ) continue;		// Charged particle only 

						if (fabs(sv->eta())>1) continue;

						Pt = sv->perp();
						Phi = sv->phi_std();
						Eta = sv->eta();

						Float_t mdPhi = Phi - HistTrigPhi[ihevt][iht];
						Float_t mdEta = Eta - HistTrigEta[ihevt][iht];

						//Folding 
						if(mdPhi<-pi/2) mdPhi+=2*pi;
						if(mdPhi>3*pi/2) mdPhi-=2*pi;


						Int_t region = 999;
						for(int iptr = 0 ;iptr<NoAssocPt; iptr++) {
							if(Pt<PtHigh[iptr]&&Pt>=PtLow[iptr]) {
								region = iptr;          // found pt bin

								// fill mix events


								mdeta_mdphi[region]->Fill(mdEta,mdPhi, weight4vz);             // 2014.04.01 ly fill by real/mix ratio
								flag++;
							}
						}
					}// End of Associate loop for mixing event
					if(flag>0) {
						MatchedTriggerVz->Fill(HistVz[ihevt],weight4vz );          // 2014.03.31 ly  Mixed Trigger Particle Vz dist. for all trigger particle in history found to be mix // 2014.04.01 ly 
					}
				}// End of trig loop for mix event
				if(flag>0) {
					HistTimes[ihevt]++;
				}
			}// End of Mix event array

			//record trigger particle for later mixing events
			if(TotTrig<=0) continue;
			//find a spot for current event
			int seat = -1;
			//cout<<"find seat?"<<endl;     // test
			for(int isp = 0 ;isp<ArrayLength; isp++) {
				//            cout<<isp<<"\t"<<HistTimes[isp]<<"\t"<<HistTotTrig[isp]<<endl;
				//            if( HistTimes[isp]>MaxMatch ) cout<<"history["<<isp<<"] matched enough"<<endl;
				if( (HistTimes[isp]==0&&HistTotTrig[isp]==0) || (HistTimes[isp]>MaxMatch)  ) {
					seat = isp;
					break;
				}
			}
			if(seat==-1) {            // overflow, no seat for current event trigger particles to record
				int DelFlag = 0;                        // test if successful delete one bin
				for(int imatch = MaxMatch; imatch >= 0  ; imatch --) {
					//                cout<<"match = "<<imatch<<endl;
					for(int dary = 0; dary < ArrayLength; dary++) {
						if(HistTimes[dary]==(imatch)) {          // delete the old bin which been matched 2 (MaxMatch-1) times, if not search for 1 times, 0 times
							//                        cout<<dary<<"\t"<<HistTimes[dary]<<endl;
							HistTimes[dary] = MaxMatch+1;            // set as MaxMatch times 
							seat = dary;                                // put in the active event here
							DelFlag=1;                                  // found
							break;
						}
					}
					if(DelFlag==1) break;
				}
				scratchy++;             // how many times delete the old bins to make room for the later one
				//            cout<<"delete "<<seat<<endl;
			}
			//cout<<"seat = "<<seat<<endl;  //test
			//store info.
			HistTimes[seat] = 1;
			HistMult[seat] = tracks;                // 2014.04.08 ly    match multiplicity now
			HistVz[seat]=header->GetPrimaryVertexZ();
			HistVx[seat]=header->GetPrimaryVertexX();
			HistVy[seat]=header->GetPrimaryVertexY();
			HistTotTrig[seat] = TotTrig;
			for(int i = 0; i<TotTrig; i++) {
				HistTrigPhi[seat][i] = TrigPhi[i];
				HistTrigEta[seat][i] = TrigEta[i];
			}

		} // while NextEvent
	} catch ( exception& e) {
		cerr << "Caught " << e.what() << endl;
		return -1;
	}
	cout << "##################################################################" << endl;

	//Long64_t nEventsUsed=reader.GetNOfEvents();  

	// Write output 
	// -------------
	TFile *fout = new TFile(OutFileName,"recreate");
	fout->cd();
	//Write Histograms
	for(int i =0;i<NoAssocPt;i++)
	{
		deta_dphi[i]->Write();
		mdeta_mdphi[i]->Write();
		EtaDistribution[i]->Write();
		pTDistribution[i]->Write();
	}


	TrigTot->Write();
	TpTDistribution->Write();
	TEtaDistribution->Write();


	EventTriggerVz->Write();
	MatchedTriggerVz->Write();

	BBCEvsrefMult->Write();
	BBCEvsNPrimaryTrk->Write();
	BBCEvsNGlobalTrk->Write();

	refMultvsNtrk->Write();
	refMultvsNPrimaryTrk->Write();
	refMultvsNGlobalTrk->Write();

	fout->Close();

	cout<<"Times delete events without mix "<<MaxMatch<<" x: "<<scratchy<<endl;


	vector<int> timesummary (MaxMatch+3,0);
	for(int i = 0; i<ArrayLength; i++) {
		//        cout<<i<<"\t"<<HistTimes[i]<<"\t"<<HistTotTrig[i]<<endl;
		if(HistTotTrig[i]!=0) {
			if(HistTimes[i]<=MaxMatch+1 && HistTimes[i]>=0) {
				timesummary[HistTimes[i]]++;
			}
			else {
				timesummary[MaxMatch+2]++;              // overflow or underflow
			}
		}
	}
	cout<<endl;
	cout<<"Total Trig Particle = "<<nTrigProcessed<<endl;
	for(int ib = 0; ib<MaxMatch+2; ib++) {
		cout<<timesummary[ib]<<" events mixed "<<ib<<" times in buffer"<<endl;
	}
	if(timesummary[MaxMatch+2]) cout<<" and addtional "<<timesummary[MaxMatch+2]<<" are overflow or underflow mixed"<<endl;
	cout<<endl;


	cout << "Bye." << endl;

	return 0;
}

