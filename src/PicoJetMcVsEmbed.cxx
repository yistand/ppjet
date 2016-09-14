//==============================================================================================================
//
//	2016.06.16	Li YI
//	Prepare for unfolding jet, we need to get MC particle level and Rc detector level jets from embedding
//
//
//	2016.07.20	Li YI
//	Need to do it for each simulated parton pt bin. Then weigthed by cross section/generated events
//
//==============================================================================================================



#include "jetpt_McVsEmbed.hh"
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
using namespace fastjet;

typedef std::pair<float, float> EtaPhiPair;



float CorrectBemcVzEta(float geoEta, float PrimVertexZ, float radius = 224){	// bemc r_int  = 224 cm
// eta in trigger info for bemc is just geometry eta. Need to shift according to the vertex position of each collisions
	float corrEta = 0;
	float T1;
        T1=2*atan(exp(-geoEta));
        Double_t zNew;
        if(geoEta!=0){zNew=radius/tan(T1);} //radius in cm 
        if(geoEta==0){zNew=0;}
        double zNom=zNew-PrimVertexZ;
        double THETA=atan2(radius,zNom);
        corrEta=-log(tan(THETA/2));

	return corrEta;

}

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
// Reader for Rc (reconstructed)
TStarJetPicoReader SetupReader ( TChain* chain, TString TriggerString, const double RefMultCut ){
	//TStarJetPicoDefinitions::SetDebugLevel(10); // 10 for more output, 0 for less output

	TStarJetPicoReader reader;
	reader.SetInputChain (chain);

	// Event and track selection
	// -------------------------
	TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
	evCuts->SetTriggerSelection( TriggerString ); //All, MB, HT, pp, ppHT, ppJP
	//evCuts->SetTriggerSelection( "All" ); //All, MB, HT, pp, ppHT, ppJP	// test
	// Additional cuts 
	evCuts->SetVertexZCut (AjParameters::VzCut);
	evCuts->SetRefMultCut ( RefMultCut );
	evCuts->SetVertexZDiffCut( 9999 ); //AjParameters::VzDiffCut );		// vpd not available for embedding we have at 2016.07...

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
	towerCuts->AddBadTowers("./include/pp200Y12_badtower.list");		// #LY CHECK where is the bad tower list

	// Tower energy correction (subtract associated charged particle deposit energy). By default, it is MIP correction (comment out the following 3 lines)
	reader.SetApplyFractionHadronicCorrection(kTRUE);
	reader.SetFractionHadronicCorrection(0.9999);
	reader.SetRejectTowerElectrons( kFALSE );


	std::cout << "Using these tower cuts:" << std::endl;
	std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;

	// V0s: Turn off
	reader.SetProcessV0s(false);

	return reader;

}

// Reader for Mc
TStarJetPicoReader SetupMcReader ( TChain* chain){

	TStarJetPicoReader reader;
	reader.SetInputChain (chain);

	// Event and track selection
	// -------------------------
	TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
	evCuts->SetTriggerSelection( "All" ); //All, MB, HT, pp, ppHT, ppJP
	// No Additional cuts 
	evCuts->SetVertexZCut (AjParameters::VzCut);	//test
	//test evCuts->SetVertexZCut (99999);
	evCuts->SetRefMultCut (0);
	evCuts->SetVertexZDiffCut(999999);

	evCuts->SetMaxEventPtCut (99999);
	evCuts->SetMaxEventEtCut (99999);

	evCuts->SetPVRankingCutOff();		//  Use SetPVRankingCutOff() to turn off vertex ranking cut.  default is OFF

	// Tracks cuts
	TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
	trackCuts->SetDCACut(99999);
	trackCuts->SetMinNFitPointsCut(-1);
	trackCuts->SetFitOverMaxPointsCut(-1);
	trackCuts->SetMaxPtCut (99999);

	// Towers: should be no tower in MC. All (charged or neutral) are handled in track
	TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
	towerCuts->SetMaxEtCut(99999);
	// Tower energy correction (subtract associated charged particle deposit energy). By default, it is MIP correction (comment out the following 3 lines)
	//reader.SetApplyFractionHadronicCorrection(kTRUE);
	//reader.SetFractionHadronicCorrection(0.9999);
	//reader.SetRejectTowerElectrons( kFALSE );

	// V0s: Turn off
	reader.SetProcessV0s(false);

	return reader;

}


int main ( int argc, const char** argv ) {

	const char *defaults[] = {"PicoJetMcVsEmbed","/home/fas/caines/ly247/Scratch/embedPythia/pt2_3_JetMcVsEmbedMatchTrig.root","ppJP2","/home/fas/caines/ly247/Scratch/embedPythia/160726/pt2_3*.root"};
	//const char *defaults[] = {"PicoJetMcVsEmbed","test.root","ppJP2","~/Scratch/embedPythia/160726/pp12Pico_pt35_-1*.root"};


	if ( argc==1 ) {
		argv=defaults;
		argc=sizeof (defaults ) / sizeof (defaults[0] );
	}

	// Throw arguments in a vector
	// ---------------------------
	vector<string> arguments(argv + 1, argv + argc);

	// Load and set up tree
	// --------------------
	TString McChainName  = "JetTreeMc";
	TString RcChainName  = "JetTree";
	TString OutFileName = arguments.at(0);


	// jet resolution parameter
	// ------------------------
	float R = 0.6;

	if ( OutFileName.Contains ("R0.2") ){
		R=0.2;
	}
	if ( OutFileName.Contains ("R0.4") ){
		R=0.4;
	}
	if ( OutFileName.Contains ("R0.6") ){
		R=0.6;
	}




	cout << " ################################################### " << endl;
	cout << "Triggering with R=" << R << endl;
	cout << " ################################################### " << endl;

	cout<<"TriggerName: "<<arguments.at(1)<<endl;
	TString TriggerName = arguments.at(1);

	int TrigFlagId = 0;
	if(TriggerName.EqualTo("ppJP2")) TrigFlagId = 1236;		//// JP2               HERE NEED TO IMPROVE, NOW IT IS PUT IN BY HAND


	cout<<"Chain data: "<<arguments.at(2).data()<<" for "<<RcChainName<<" and "<<McChainName<<endl;
	TChain* chain = new TChain( RcChainName );
	chain->Add( arguments.at(2).data() );
	TChain* Mcchain = new TChain( McChainName );
	Mcchain->Add( arguments.at(2).data() );





	cout<<"SetupReader for RcPico"<<endl;
	double RefMultCut = 0;
	TStarJetPicoReader reader = SetupReader( chain, TriggerName,RefMultCut );			// #ly note: Events & Tracks & Towers cuts are set here
	reader.SetTrackPileUpCut(0);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	if( OutFileName.Contains ("NoTofMatch") ) {
	  reader.SetTrackPileUpCut(0);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcOrTofMatch") ) {		// NOT IMPLEMMENTED YET!!!!!!!!!!!!!!!!!!
	  reader.SetTrackPileUpCut(1);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcMatch") ) {
	  reader.SetTrackPileUpCut(3);		// #ly	3: tpc track matching to bemc.		1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}

	TStarJetPicoDefinitions::SetDebugLevel(0);

	TStarJetPicoReader Mcreader = SetupMcReader( Mcchain); 

	// Initialize analysis class
	// -------------------------

	string jetalgorithm = "antikt";		
	if(OutFileName.Contains ("kT") && (!(OutFileName.Contains ("AntikT")))) {
		jetalgorithm = "kt";
	}

	std::cout<<"OutFileName="<<OutFileName<<std::endl;

	jetpt_McVsEmbed *jme = new jetpt_McVsEmbed(	R,
							AjParameters::max_track_rap,	
							0.2,                    // pt min for const.
							jetalgorithm,
							OutFileName
						  );

	// initial ttree & histograms in jme
	jme->Init();

	if(OutFileName.Contains ("pt2_3")) jme->SetOutlierMcpTCut(20);	// No event with Mc pT > 20 for 2<pt<3 bin

	
	if(OutFileName.Contains ("MatchTrig")&&TriggerName.EqualTo("ppJP2") ) jme->SetToMatchJetTrigger(true);			// whether match jet found with fastjet with the location which fired the trigger, NEED TO CHECK TrigFlagId
	else jme->SetToMatchJetTrigger(false);


	// Cycle through events
	// --------------------
	vector<PseudoJet> Mcparticles;		
	vector<PseudoJet> particles;		
	TStarJetVectorContainer<TStarJetVector>* Mccontainer;		
	TStarJetVectorContainer<TStarJetVector>* container;		
	TStarJetVector* Mcsv; // TLorentzVector* would be sufficient. 
	TStarJetVector* sv; // TLorentzVector* would be sufficient. 
	PseudoJet Mcpj;
	PseudoJet pj;

	std::vector<EtaPhiPair> TrigLoc2Match;		// trigger of High Tower or Jet Patch

	//std::set<int>badrun;
	//badrun.clear();

	//readinbadrunlist(badrun);        

	jme->SetVerbose(0);

	Long64_t nEvents=-1; // -1 for all
	//nEvents=1000;	// test
	cout<<"init..."<<nEvents<<endl;
	Mcreader.Init(nEvents);
	reader.Init(nEvents);
	int count = 0;

	try{
		cout<<"START:"<<endl;
		while ( Mcreader.NextEvent() ) {
			Mcreader.PrintStatus(10);
			if(count%10000==0) cout<<"event "<<count<<endl;
			count++;
				
			
			// clean-up
			Mcparticles.clear();
			particles.clear();
			TrigLoc2Match.clear();

			// First do Mc event to see whether this event passed. 

			// event info
			// ----------
			//cout<<"load event header"<<endl;
			TStarJetPicoEventHeader* Mcheader = Mcreader.GetEvent()->GetHeader();
			int Mcrunid   = Mcheader->GetRunId();
			int Mceventid   = Mcheader->GetEventId();
			//cout<<"Mc runid="<<Mcrunid<<" Mceventid="<<Mceventid<<endl;		


			// Load event particles
			// ----------
			Mccontainer = Mcreader.GetOutputContainer();

			// Make particle vector
			// --------------------

			for (int mcip = 0; mcip<Mccontainer->GetEntries() ; ++mcip ){
				Mcsv = Mccontainer->Get(mcip);  // Note that TStarJetVector contains more info

				Mcpj=MakePseudoJet( Mcsv );
				Mcpj.set_user_info ( new JetAnalysisUserInfo( 3*Mcsv->GetCharge(), Mcsv->GetFeatureD(TStarJetVector::_DEDX), 0)  );	// for Mc track, DEDX is actually pdg id. TOF one is set to 0.
				Mcpj.set_user_index(mcip);		
				//cout<<"input "<<Mcsv->GetCharge() <<" -> "<<Mcpj.user_info<JetAnalysisUserInfo>().GetQuarkCharge()<<endl;	 

				Mcparticles.push_back ( Mcpj );
			}    



			// Then loop over Rc event
			// =========================
			if(reader.ReadEvent(Mcreader.GetNOfCurrentEvent())) { 	// try to load the same event as Mc one
				// If Rc event passed the cuts:
				// Need to check whether the same event. It SHOULD!
				TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();
				int runid   = header->GetRunId();
				int eventid   = header->GetEventId();


				//test #ly somehow McEvent doesnot have correct runid. they are all 0... if(runid!=Mcrunid || eventid!=Mceventid) {
				if(eventid!=Mceventid) {	// test #ly only check eventid for now
					cout<<"ERROR!!!!!!!!!!!!!!!!!! MC and RC not consistent!!!!!!!"<<endl;
					cout<<"runid Rc = "<<runid<<" VS  Mc = "<<Mcrunid<<endl;
					cout<<"eventid Rc = "<<eventid<<" VS  Mc = "<<Mceventid<<endl;
					return -1;
				}

				//cout<<"eventid Rc = "<<eventid<<endl<<endl<<endl;	

				// event info
				// ----------
				//cout<<"load event header"<<endl;

				// eventid = header->GetEventId();
				//if(badrun.count(runid)>0) continue;			// in bad run list
				//if(header->GetZdcCoincidenceRate()>6000) continue;		// test
	


				// Load event ht/jetpatch trigger objs
				// ----------
				//std::cout<<"load trigger objs"<<endl;	
				//test if(jme->GetToMatchJetTrigger()) {
					TClonesArray *trigobj = reader.GetEvent()->GetTrigObjs();
					for(int itrg = 0; itrg<trigobj->GetEntries(); itrg++) {
						if( ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetTriggerFlag()==TrigFlagId )	 { 
							EtaPhiPair itrigloc =std::make_pair(CorrectBemcVzEta(((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetEta(),header->GetPrimaryVertexZ()), ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetPhi()) ;
							TrigLoc2Match.push_back(itrigloc);
						}
					}
				//test }


				// Load event particles
				// ----------
				container = reader.GetOutputContainer();

				// Make particle vector
				// --------------------

				for (int ip = 0; ip<container->GetEntries() ; ++ip ){
					sv = container->Get(ip);  // Note that TStarJetVector contains more info, such as charge;


					//if (sv->GetCharge()==0 ) (*sv) *= fTowScale; // for systematics
					pj=MakePseudoJet( sv );
					pj.set_user_info ( new JetAnalysisUserInfo( 3*sv->GetCharge(), sv->GetFeatureD(TStarJetVector::_DEDX), sv->GetFeatureD(TStarJetVector::_TOFBETA) ) );
					pj.set_user_index(ip);		// #ly	link fastjet::PseudoJet to TStarJetVector class	--> NEED TO FIX THIS, NOT SURE WHY USER_INFO IS NOT PASSED TO JAResult.at(0).constituents() in UnderlyingAna.cxx
					//cout<<"input "<<sv->GetCharge() <<" -> "<<pj.user_info<JetAnalysisUserInfo>().GetQuarkCharge()<<endl;	 

					particles.push_back ( pj );
					//}	      
				}    
			}


			// Run jet finding analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	

			jme->Make( 	Mcparticles, 
					particles, 
					//*container, 
					reader.GetEvent()->GetHeader()->GetEventId(),
					reader.GetEvent()->GetHeader()->GetRunId(),
					Mcreader.GetEvent()->GetHeader()->GetGReferenceMultiplicity(),
					Mcreader.GetEvent()->GetHeader()->GetPrimaryVertexZ(),
					reader.GetEvent()->GetHeader()->GetGReferenceMultiplicity(),
					reader.GetEvent()->GetHeader()->GetPrimaryVertexZ(),
					TrigLoc2Match
			);



		} // while NextEvent
	} catch ( exception& e) {
		cerr << "Caught " << e.what() << endl;
		return -1;
	}
	cout << "##################################################################" << endl;

	//Long64_t nEventsUsed=reader.GetNOfEvents();  


	if(jme->IsOutlierMcpTCutApplied() && jme->GetNEventOutlierMcpTCut()> 3) { cout << " ERROR: there is "<<jme->GetNEventOutlierMcpTCut()<<" events were rejected due to MC pT>"<<jme->GetOutlierMcpTCut()<<" cut. NEED TO CHECK!!!!!! "<<endl; }

	// Close up shop
	// -------------
	jme->Finish();


	//delete jme;
	cout << "Bye." << endl;

	return 0;
}


