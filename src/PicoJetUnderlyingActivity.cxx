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
	const char *defaults[] = {"PicoJetUnderlyingActivity","/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/FullJet_TransCharged_MatchTrig_ppJP2.root","ppJP2","/home/hep/caines/ly247/Scratch/pp12JP2Pico_151018/*.root", "0", "0" };
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

	int TrigFlagId = 0;
	if(TriggerName.EqualTo("ppJP2")) TrigFlagId = 1236;		//// JP2               HERE NEED TO IMPROVE, NOW IT IS PUT IN BY HAND
	if(TriggerName.EqualTo("ppJP1")) TrigFlagId = 1228;		//// JP1               HERE NEED TO IMPROVE, NOW IT IS PUT IN BY HAND
	if(TriggerName.EqualTo("ppJP0")) TrigFlagId = 1220;		//// JP0               HERE NEED TO IMPROVE, NOW IT IS PUT IN BY HAND


	cout<<"Chain data: "<<arguments.at(2).data()<<" for "<<ChainName<<endl;
	TChain* chain = new TChain( ChainName );
	chain->Add( arguments.at(2).data() );

	// Au+Au?
	// ------
	bool isAuAu=false;
	if (arguments.at(2).find("AuAu") != std::string::npos ) isAuAu=true; // Quick and dirty...
	if (arguments.at(2).find("auau") != std::string::npos ) isAuAu=true; // Quick and dirty...

	// for systematics
	// ---------------
	Int_t IntTowScale=atoi( arguments.at(3).data() ); // +/- 2%
	Float_t fTowScale = 1.0 + IntTowScale*0.02;
	Int_t mEffUn=atoi( arguments.at(4).data() ) ;

	switch ( mEffUn ){
		case 0 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Eff0_")+ gSystem->BaseName(OutFileName));
			break;
		case 1 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Eff1_")+ gSystem->BaseName(OutFileName));
			break;
		case -1 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Eff-1_")+ gSystem->BaseName(OutFileName));
			break;
		default :
			cerr << "mEffUn = " << mEffUn << " not supported." <<endl;
			return -1;
	}

	switch ( IntTowScale ){
		case 0 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Tow0_")+ gSystem->BaseName(OutFileName));
			break;
		case 1 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Tow1_")+ gSystem->BaseName(OutFileName));
			break;
		case -1 :
			if ( !isAuAu ) OutFileName.ReplaceAll ( gSystem->BaseName(OutFileName), TString ("Tow-1_")+ gSystem->BaseName(OutFileName));
			break;
		default :
			cerr << "IntTowScale = " << IntTowScale << " not supported." <<endl;
			return -1;
	}
	if ( isAuAu && IntTowScale ){
		cerr << "IntTowScale = " << IntTowScale << " not supported in AuAu." <<endl;
		return -1;    
	}

	cout<<"SetupReader for pico"<<endl;
	double RefMultCut = 0;
	TStarJetPicoReader reader = SetupReader( chain, TriggerName, RefMultCut );			// #ly note: Events & Tracks & Towers cuts are set here
	//reader.SetTrackPileUpCut(kTRUE);		// #ly	tpc track matching to bemc or tof
	reader.SetTrackPileUpCut(2);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	if( OutFileName.Contains ("NoTofMatch") ) {
	  reader.SetTrackPileUpCut(0);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcOrTofMatch") ) {
	  reader.SetTrackPileUpCut(1);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	if( OutFileName.Contains ("BemcMatch") ) {
	  reader.SetTrackPileUpCut(3);		// #ly	3: tpc track matching to bemc.		1: tpc track matching to bemc or tof. 	2: tof match only.    0: no requirement for fast detector matching
	}
	TStarJetPicoDefinitions::SetDebugLevel(0);

	// Initialize analysis class
	// -------------------------
	//cout<<"initialize analysis"<<endl;		

	//TString OutFileName = "test.root"; 
	//float R = 0.6;	
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
	int jetchargecode = 2;			// default one is to do full jet (charged + neutral) for jet finding
	if ( OutFileName.Contains ("ChargeJet") ){
		jetchargecode = 1;
	}
	else if ( OutFileName.Contains ("NeutralJet") ){
		jetchargecode = 0;
	}
	else if ( OutFileName.Contains ("FullJet") ){
		jetchargecode = 2;
	}

	int underlyingchargecode = 2;		// default one is takeing both charged + neutral particles for underlying event
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

	if(OutFileName.Contains ("MatchTrig")&&TriggerName.Contains("ppJP") ) ula->SetToMatchJetTrigger(true);			// whether match jet found with fastjet with the location which fired the trigger, NEED TO CHECK TrigFlagId
	else ula->SetToMatchJetTrigger(false);

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


        // problematic runs, need future investigation			<--------- Moved to include/pp200Y12_badrun.list
        //const int NoBadRun = 185;
        //int badrun[NoBadRun] = {13044118, 13044123, 13044124, 13044125, 13045001, 13045003, 13045005, 13045006, 13045007, 13045012, 13045029, 13046002, 13046008, 13046010, 13046029, 13046118, 13046119, 13046120, 13047004, 13047014, 13047018, 13047036, 13047037, 13047039, 13047040, 13047041, 13047042, 13047043, 13047044, 13047045, 13047046, 13047047, 13047048, 13047049, 13047050, 13047051, 13047052, 13047053, 13047054, 13047055, 13048007, 13048022, 13048046, 13049004, 13049005, 13049050, 13049052, 13049075, 13049086, 13049087, 13049088, 13049089, 13050007, 13050025, 13050026, 13050027, 13050033, 13050039, 13050043, 13050044, 13050046, 13050047, 13050049, 13050050, 13051068, 13051080, 13051088, 13051095, 13051102, 13052021, 13052022, 13052054, 13052063, 13052068, 13053010, 13053021, 13054004, 13054005, 13054006, 13054007, 13054008, 13054009, 13054011, 13054012, 13054013, 13054014, 13054015, 13054016, 13054017, 13054018, 13054019, 13054020, 13054022, 13054042, 13054045, 13054046, 13054057, 13055015, 13055072, 13055081, 13055082, 13055086, 13055087, 13055088, 13055089, 13055090, 13056011, 13056012, 13056034, 13056035, 13056037, 13056038, 13056039, 13057038, 13057039, 13058019, 13058030, 13058047, 13058048, 13059003, 13059004, 13059005, 13059006, 13059007, 13059008, 13059009, 13059010, 13059019, 13059035, 13059082, 13059083, 13059084, 13059085, 13059086, 13059087, 13060001, 13060002, 13060003, 13060009, 13060012, 13061026, 13063033, 13064030, 13064057, 13064059, 13064074, 13066035, 13066036, 13066101, 13066102, 13066104, 13066109, 13066110, 13067001, 13067002, 13067003, 13067004, 13067005, 13067006, 13067007, 13067008, 13067009, 13067010, 13067011, 13067012, 13067013, 13067014, 13067015, 13067017, 13068017, 13068022, 13068027, 13068029, 13068034, 13068036, 13068037, 13069006, 13069009, 13069029, 13070030, 13070056, 13071034, 13071037, 13071038, 13071040};
	
	std::set<int>badrun;
	badrun.clear();

	readinbadrunlist(badrun);        


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
			//if(runid>=13058000&& runid<13061000) continue;          // a dip in TPC primary tracks. problematic runs, need future investigation	<--------- Moved to include/pp200Y12_badrun.list
			//for(int i = 0; i<NoBadRun; i++) {
                        //	if(runid==badrun[i]) continue;
                	//}

			//if(header->GetZdcCoincidenceRate()>6000) continue;		// test
	


			// Load event ht/jetpatch trigger objs
			// ----------
			//std::cout<<"load trigger objs"<<endl;	
			TrigLoc2Match.clear();
			if(ula->GetToMatchJetTrigger()) {
				TClonesArray *trigobj = reader.GetEvent()->GetTrigObjs();
				for(int itrg = 0; itrg<trigobj->GetEntries(); itrg++) {
					if( ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetTriggerFlag()==TrigFlagId )	 { 
						EtaPhiPair itrigloc =std::make_pair(CorrectBemcVzEta(((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetEta(),header->GetPrimaryVertexZ()), ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetPhi()) ;
						TrigLoc2Match.push_back(itrigloc);
					}
				}
			}


			// Load event particles
			// ----------
			container = reader.GetOutputContainer();

			// Make particle vector
			// --------------------
			particles.clear();

			for (int ip = 0; ip<container->GetEntries() ; ++ip ){
				sv = container->Get(ip);  // Note that TStarJetVector contains more info, such as charge;

				if(fabs(sv->perp())<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

				if (sv->GetCharge()==0 ) (*sv) *= fTowScale; // for systematics
				pj=MakePseudoJet( sv );
				pj.set_user_info ( new JetAnalysisUserInfo( 3*sv->GetCharge(), sv->GetFeatureD(TStarJetVector::_DEDX), sv->GetFeatureD(TStarJetVector::_TOFBETA) ) );
				pj.set_user_index(ip);		// #ly	link fastjet::PseudoJet to TStarJetVector class	--> NEED TO FIX THIS, NOT SURE WHY USER_INFO IS NOT PASSED TO JAResult.at(0).constituents() in UnderlyingAna.cxx
				//cout<<"input "<<sv->GetCharge() <<" -> "<<pj.user_info<JetAnalysisUserInfo>().GetQuarkCharge()<<endl;	 

				particles.push_back ( pj );

				//}	      
			}    
			// Run analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	

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

