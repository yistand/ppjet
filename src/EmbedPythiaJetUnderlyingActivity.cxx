//====================================================================================================
//
//	2015.09.22	Li Yi
//	modified from Kolja Kauder's Aj analysis to study underlying event activity dependence on jet
//	energy in pp 200 GeV
//
//====================================================================================================


/** 
  @author Kolja Kauder
  @version Revision 0.2
  @brief Aj analysis and embedding prep in p+p.
  @details Perform Aj analysis in a given TStarPicoJetTree chain. Can also save events with 10 GeV jets for embedding.
  @date Mar 04, 2015
  */

#include "UnderlyingAna.hh"
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


/** 
  - Set up input tree
  - Set up output histos and tree
  - Initialize UnderlyingAna object
  - Loop through events
  \arg argv: 
  - [1] : output file
  - [2] : trigger name
  - [3] : Input file pattern. Let TChain handle the globbing, i.e. use for example
  <BR><tt>% PicoJetUnderlyingActivity '~putschke/Data/ppHT/&lowast;.root'</tt>
  <BR>For cases like 
  <BR><tt>% PicoJetUnderlyingActivity ~putschke/Data/&lowast;/&lowast;.root</tt>
  <BR>change this macro
  - [4] : tower uncertainty switch ( -1/0/1 )
  - [5] : efficiency correction uncertainty switch ( -1/0/1 )
  <BR> Note that we use a vector for argv and shift the index down by one.
  */

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
TStarJetPicoReader SetupReader ( TChain* chain, TString TriggerString, const double RefMultCut ){
	TStarJetPicoDefinitions::SetDebugLevel(10); // 10 for more output, 0 for less output

	TStarJetPicoReader reader;
	reader.SetInputChain (chain);

	// Event and track selection
	// -------------------------
	TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
	evCuts->SetTriggerSelection( TriggerString ); //All, MB, HT, pp, ppHT, ppJP
	// Additional cuts 
	evCuts->SetVertexZCut (AjParameters::VzCut);
	evCuts->SetRefMultCut ( RefMultCut );
	evCuts->SetVertexZDiffCut( AjParameters::VzDiffCut );

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
	std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;

	// V0s: Turn off
	reader.SetProcessV0s(false);

	return reader;

}


int main ( int argc, const char** argv ) {

	// Set up some convenient default
	// ------------------------------
	const char *defaults[] = {"EmbedPythiaJetUnderlyingActivity","/home/fas/caines/ly247/Scratch/embedPythia/FullJet_TransCharged_embedPythia.root","All","/home/fas/caines/ly247/Scratch/embedPythia/embedpp200Y12Pico_160509.root", "0", "0" };
	// {Code name, to be discard but needed since argv will use command name as the [0], output file name, triggername, intput file list, for variable IntTowScale to scale tower as systematics study, which effiencey file to use }
	//
	// output file name can include(optional): 
	// 	"R0.6" (default) OR "R0.4" OR "R0.2";
	// 	"FullJet"(default) OR "ChargeJet" OR "NeutralJet";
	// 	"TransCharged" (default) particle only OR "TransNeutral" particle only;
	// 	"AntikT" (default) OR "kT"
	// 	"MatchTrig" together with "ppJP2" will match leading jet with the trigger jet phi, eta
	// 	"Monojet" (default) OR "Dijet"
	//	"TranPhi60" (default) OR "TranPhi30"
	

	if ( argc==1 ) {
		argv=defaults;
		argc=sizeof (defaults ) / sizeof (defaults[0] );
	}

	// Throw arguments in a vector
	// ---------------------------
	vector<string> arguments(argv + 1, argv + argc);

	// Load and set up tree
	// --------------------
	TString ChainName  = "JetTree";
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


	cout<<"Chain data: "<<arguments.at(2).data()<<" for "<<ChainName<<endl;
	TChain* chain = new TChain( ChainName );
	chain->Add( arguments.at(2).data() );


	cout<<"SetupReader for pico for embedded pythia"<<endl;
	double RefMultCut = 0;
	TStarJetPicoReader reader = SetupReader( chain, TriggerName, RefMultCut );			// #ly note: Events & Tracks & Towers cuts are set here
	reader.SetTrackPileUpCut(0);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	TStarJetPicoDefinitions::SetDebugLevel(0);

	// Initialize analysis class
	// -------------------------
	//cout<<"initialize analysis"<<endl;		

	//TString OutFileName = "test.root"; //test 
	//float R = 0.6;	// test
	string jetalgorithm = "antikt";		
	if(OutFileName.Contains ("kT") && (!(OutFileName.Contains ("AntikT")))) {
		jetalgorithm = "kt";
	}
	UnderlyingAna *ula = new UnderlyingAna( R,
			AjParameters::max_track_rap,
			0.2,			// pt min for const.
			AjParameters::dPhiCut,
			OutFileName,
			jetalgorithm
	);  


	// Charge selection for jet and underlying particle
	int jetchargecode = 2;
	if ( OutFileName.Contains ("ChargeJet") ){
		jetchargecode = 1;
	}
	else if ( OutFileName.Contains ("NeutralJet") ){
		jetchargecode = 0;
	}
	else if ( OutFileName.Contains ("FullJet") ){
		jetchargecode = 2;
	}

	int underlyingchargecode = 2;
	if ( OutFileName.Contains ("TransCharged") ){
		underlyingchargecode = 1;
	}
	else if ( OutFileName.Contains ("TransNeutral") ){
		underlyingchargecode = 0;
	}


	cout << " ################################################### " << endl;
	cout << " jetchargecode = " << jetchargecode <<endl; 
	cout << " underlyingchargecode = " << underlyingchargecode <<endl; 
	cout << " ################################################### " << endl;

	// initial ttree & histograms in ula
	ula->Init();

	//if(OutFileName.Contains ("MatchTrig")&&TriggerName.EqualTo("ppJP2") ) ula->SetToMatchJetTrigger(true);			// whether match jet found with fastjet with the location which fired the trigger, NEED TO CHECK TrigFlagId
	//else 
	ula->SetToMatchJetTrigger(false);	// No jet matching for MC

	if(jetchargecode==2) ula->SetNetraulJetFracCut(true);			// whether apply neutral energy fraction in jet cut
	else ula->SetNetraulJetFracCut(false);

	ula->SetJetCharge(jetchargecode);			// Jet charge

	ula->SetUnderlyingParticleCharge(underlyingchargecode);			// underlying event charge: 0 for netural, 1 for charged, 2 for all

	if(OutFileName.Contains ("Dijet")) {
		ula->SetDiJetAngle(1);					// Use Dijet angle (1) or Monojet angle (0) 
	}
	else {
		ula->SetDiJetAngle(0);					// Use Dijet angle (1) or Monojet angle (0) 
	}
	
	if ( OutFileName.Contains ("TranPhi30") ) {
		ula->SetTransversePhiSize(30);	
	}

	// Cycle through events
	// --------------------
	vector<PseudoJet> particles;		// for jet finding
	TStarJetVectorContainer<TStarJetVector>* container;		// for underlying event loop
	TStarJetVector* sv; // TLorentzVector* would be sufficient. 
	PseudoJet pj;

	std::vector<EtaPhiPair> TrigLoc2Match;		// trigger of High Tower or Jet Patch


	Long64_t nEvents=-1; // -1 for all
	//nEvents=10000;	// test
	cout<<"init..."<<endl;
	reader.Init(nEvents);
	int count = 0;


	try{
		while ( reader.NextEvent() ) {
			reader.PrintStatus(10);
			if(count%10000==0) cout<<"event "<<count<<endl;
			count++;

			// event info
			// ----------
			//cout<<"load event header"<<endl;
			TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();

			//if(header->GetZdcCoincidenceRate()>6000) continue;		// test


			// Load event particles
			// ----------
			container = reader.GetOutputContainer();

			// Make particle vector
			// --------------------
			particles.clear();
			TrigLoc2Match.clear();

			for (int ip = 0; ip<container->GetEntries() ; ++ip ){
				sv = container->Get(ip);  // Note that TStarJetVector contains more info, such as charge;

				if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

				Float_t fTowScale=1;
				if (sv->GetCharge()==0 ) (*sv) *= fTowScale; // for systematics
				pj=MakePseudoJet( sv );
				pj.set_user_info ( new JetAnalysisUserInfo( 3*sv->GetCharge(), sv->GetFeatureD(TStarJetVector::_DEDX), sv->GetFeatureD(TStarJetVector::_TOFBETA) ) );
				pj.set_user_index(ip);		// #ly	link fastjet::PseudoJet to TStarJetVector class	--> NEED TO FIX THIS, NOT SURE WHY USER_INFO IS NOT PASSED TO JAResult.at(0).constituents() in UnderlyingAna.cxx
				//cout<<"input "<<sv->GetCharge() <<" -> "<<pj.user_info<JetAnalysisUserInfo>().GetQuarkCharge()<<endl;	// test 

				particles.push_back ( pj );
				//}	      
			}    
			// Run analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	// test

			ula->AnalyzeAndFill( particles, 
					//*container, 
					reader.GetEvent()->GetHeader()->GetEventId(),
					reader.GetEvent()->GetHeader()->GetRunId(),
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

	// Close up shop
	// -------------
	ula->Finish();


	//delete ula;
	cout << "Bye." << endl;

	return 0;
}

