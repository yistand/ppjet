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
#include <cmath>
#include <exception>
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>    // std::random_shuffle

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

int main ( int argc, const char** argv ) {

	// Set up some convenient default
	// ------------------------------
	//const char *defaults[] = {"PicoJetUnderlyingActivity","/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/test.root","ppHT","/home/hep/caines/ly247/Scratch/pp12Pico_150407/*root", "0", "0" };
	//const char *defaults[] = {"PicoJetUnderlyingActivity","/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/MatchTrig_ppJP2.root","ppJP2","/home/hep/caines/ly247/Scratch/pp12Pico_150407/*root", "0", "0" };
	//const char *defaults[] = {"PicoJetUnderlyingActivity","/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/TransCharge0_MatchTrig_ppJP2.root","ppJP2","/home/hep/caines/ly247/Scratch/pp12JP2Pico_151018/*.root", "0", "0" };
	const char *defaults[] = {"PicoJetUnderlyingActivity","/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/TransCharged_MatchTrig_ppJP2.root","ppJP2","/home/hep/caines/ly247/Scratch/pp12JP2Pico_151018/*.root", "0", "0","ChargeJet","Charge" };
	// {Code name, to be discard but needed since argv will use command name as the [0], output file name, triggername, intput file list, for variable IntTowScale to scale tower as systematics study, which effiencey file to use }
	// output file name can include "R0.6" OR "R0.4" OR "R0.2", "ChargeJet" OR "FullJet" OR "NeutralJet", "TransCharged" particle only OR "TransNeutral" particle only 
	

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

	//#ly // Also pull a random offset from the file name
	//#ly // --------------------------------------------
	//#ly // To seed the backgrounder in different ways
	//#ly // allows for patterns like output/rndm0/test.root
	//#ly // ONLY pulls one digit
	//#ly int randomoff = TString( OutFileName( OutFileName.Index("rndm") + 4,  1 )).Atoi(); // defaults to zero
	//#ly // eventid: typically < 1M --> shift by 10M
	//#ly // runid: typically 8XXYYYY --> shift by 10M
	//#ly randomoff *= 10000000;
	//#ly cout << " ################################################### " << endl;
	//#ly cout << "   FastJet random seeds offset by " << randomoff << endl;
	//#ly cout << " ################################################### " << endl;


	cout << " ################################################### " << endl;
	cout << "Triggering with R=" << R << endl;
	cout << " ################################################### " << endl;

	cout<<"TriggerName: "<<arguments.at(1)<<endl;
	TString TriggerName = arguments.at(1);

	int TrigFlagId = 0;
	if(TriggerName.EqualTo("JP2")) TrigFlagId = 1236;		//// JP2               HERE NEED TO IMPROVE, NOW IT IS PUT IN BY HAND


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
	reader.SetTrackPileUpCut(2);		// #ly	1: tpc track matching to bemc or tof. 	2: tof match ly.    0: no requirement for fast detector matching
	TStarJetPicoDefinitions::SetDebugLevel(2);

	// // Files and histograms
	// // --------------------
	// TFile* fout = new TFile( OutFileName, "RECREATE");
	// assert ( fout->IsOpen() );
	// cout << " ################################################### " << endl;
	// cout << "Writing to: " << fout->GetName() << endl;
	// cout << " ################################################### " << endl;

	// TH1::SetDefaultSumw2(true);
	// TH2::SetDefaultSumw2(true);
	// TH3::SetDefaultSumw2(true);

	// TH1D* LeadJetPt = new TH1D("LeadingJetPt","Leading Jet Pt",100,0,100);
	// TH1D* SubJetPt = new TH1D("SubLeadingJetPt","SubLeading Jet Pt",100,0,100);

	// TProfile* LeadJetNtrkvsLeadJetPt = new TProfile("LeadingJetNtrkvsLeadJetPt","Leading Jet Ntrk vs Leading Jet Pt",100,0,100);
	// TProfile* SubJetNtrkvsLeadJetPt = new TProfile("SubJetNtrkvsLeadJetPt","SubLeading Jet Ntrk vs Leading Jet Pt",100,0,100);
	// TProfile* TranMaxNtrkvsLeadJetPt = new TProfile("TranMaxNtrkvsLeadJetPt","Transverse Max Ntrk vs Leading Jet Pt",100,0,100);
	// TProfile* TranMinNtrkvsLeadJetPt = new TProfile("TranMinNtrkvsLeadJetPt","Transverse Min Ntrk vs Leading Jet Pt",100,0,100);
	// TProfile* TranNtrkvsLeadJetPt = new TProfile("TranNtrkvsLeadJetPt","Transverse Ntrk vs Leading Jet Pt",100,0,100);


	// TProfile* LeadJetPtvsLeadJetPt = new TProfile("LeadingJetPtvsLeadJetPt","Leading Jet <Pt> vs Leading Jet Pt",100,0,100);
	// TProfile* SubJetPtvsLeadJetPt = new TProfile("SubJetPtvsLeadJetPt","SubLeading Jet <Pt> vs Leading Jet Pt",100,0,100);
	// TProfile* TranMaxPtvsLeadJetPt = new TProfile("TranMaxPtvsLeadJetPt","Transverse Max <Pt> vs Leading Jet Pt",100,0,100);
	// TProfile* TranMinPtvsLeadJetPt = new TProfile("TranMinPtvsLeadJetPt","Transverse Min <Pt> vs Leading Jet Pt",100,0,100);
	// TProfile* TranPtvsLeadJetPt = new TProfile("TranPtvsLeadJetPt","Transverse <Pt> vs Leading Jet Pt",100,0,100);

	// TH2D* Spectrum_LeadJetPtvsLeadJetPt = new TH2D("Spectrum_LeadingJetPtvsLeadJetPt","Leading Jet Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,10);
	// TH2D* Spectrum_SubJetPtvsLeadJetPt = new TH2D("Spectrum_SubJetPtvsLeadJetPt","SubLeading Jet Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,10);
	// TH2D* Spectrum_TranMaxPtvsLeadJetPt = new TH2D("Spectrum_TranMaxPtvsLeadJetPt","Transverse Max Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,10);
	// TH2D* Spectrum_TranMinPtvsLeadJetPt = new TH2D("Spectrum_TranMinPtvsLeadJetPt","Transverse Min Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,10);
	// TH2D* Spectrum_TranPtvsLeadJetPt = new TH2D("Spectrum_TranPtvsLeadJetPt","Transverse Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,10);


	// TProfile* TranPionNtrkvsLeadJetPt = new TProfile("TranPionNtrkvsLeadJetPt","Transverse Pion Ntrk vs Leading Jet Pt",100,0,100);

	// TProfile* TranPionPtvsLeadJetPt = new TProfile("TranPionPtvsLeadJetPt","Transverse Pion <Pt> vs Leading Jet Pt",100,0,100);

	// TProfile* TranProtonNtrkvsLeadJetPt = new TProfile("TranProtonNtrkvsLeadJetPt","Transverse Proton Ntrk vs Leading Jet Pt",100,0,100);
	// TProfile* TranProtonPtvsLeadJetPt = new TProfile("TranProtonPtvsLeadJetPt","Transverse Proton <Pt> vs Leading Jet Pt",100,0,100);



	//#ly // DEBUG histos
	//#ly // ------------
	//#ly TH3D* ptphieta = new TH3D("ptphieta","",500, 0.2, 50.2, 100, 0, TMath::TwoPi(), 100, -1, 1);
	//#ly TH1D* csize = new TH1D("csize","",5000, -0.5, 4999.5 );
	//#ly // Find some info on the background
	//#ly TH1D* hrho = new TH1D( "hrho","#rho", 240, 0, 120 );
	//#ly TH1D* hrhoerr = new TH1D( "hrhoerr","#sigma_{#rho}", 100, 0, 50 );

	//#ly // can't compute area_error(), take rms of this instead
	//#ly TH1D* hj1area    = new TH1D( "hj1area","j1 area", 320, 0.0, 0.8 );
	//#ly TH1D* hj2area    = new TH1D( "hj2area","j2 area", 320, 0.0, 0.8 );
	//#ly TH1D* hrestarea    = new TH1D( "hrestarea","restarea", 160, 0.0, 0.8 );

	//// Save results
	//// ------------
	//TTree* ResultTree=new TTree("ResultTree","Result Jets");
	//TLorentzVector j1, j2;
	//double leadpt=0, subpt=0;	// <pT>
	//int leadntrk=0, subntrk=0;	// Ntrk in the leading & subleading jet angle regions
	//double tranpt=0, tranmaxpt, tranminpt=0;
	//int tranntrk=0, tranmaxntrk, tranminntrk=0;
	//double tranpionpt=0, tranprotonpt=0;
	//int tranpionntrk=0, tranprotonntrk=0;
	//ResultTree->Branch("LeadJet",&j1);
	//ResultTree->Branch("SubLeadJet",&j2);
	//ResultTree->Branch("LeadAreaPt",&leadpt);
	//ResultTree->Branch("SubLeadAreaPt",&subpt);
	//ResultTree->Branch("TranPt",&tranpt);
	//ResultTree->Branch("TranMaxPt",&tranmaxpt);
	//ResultTree->Branch("TranMinPt",&tranminpt);
	//ResultTree->Branch("LeadAreaNtrk",&leadntrk);
	//ResultTree->Branch("SubLeadAreaNtrk",&subntrk);
	//ResultTree->Branch("TranNtrk",&tranntrk);
	//ResultTree->Branch("TranMaxNtrk",&tranmaxntrk);
	//ResultTree->Branch("TranMinNtrk",&tranminntrk);
	//ResultTree->Branch("TranPionNtrk",&tranpionntrk);
	//ResultTree->Branch("TranProtonNtrk",&tranprotonntrk);
	//ResultTree->Branch("TranPionPt",&tranpionpt);
	//ResultTree->Branch("TranProtonPt",&tranprotonpt);

	//int eventid=0;
	//int runid=0;
	//double refmult=0; // Really an int, but may change when using refmultcorr
	//float rho=0;
	//float rhoerr=0;
	//float j1area=0;
	//float j2area=0;

	//ResultTree->Branch("eventid",&eventid, "eventid/i");
	//ResultTree->Branch("runid",&runid, "runid/i");
	//ResultTree->Branch("refmult",&refmult, "refmult/d");
	//ResultTree->Branch("rho",&rho, "rho/f");
	//ResultTree->Branch("rhoerr",&rhoerr, "rhoerr/f");
	//ResultTree->Branch("j1area",&j1area, "j1area/f");
	//ResultTree->Branch("j2area",&j2area, "j2area/f");

	// #ly // area and pT of all remaining jets (those used for rho)
	// #ly static const Int_t kmaxJ=500; // max # of jets
	// #ly int nRestJ=0;
	// #ly float RestArea[kmaxJ];
	// #ly float RestPt[kmaxJ];
	// #ly ResultTree->Branch("nRestJ",&nRestJ, "nRestJ/i");
	// #ly ResultTree->Branch("RestArea",RestArea, "RestArea[nRestJ]/f");
	// #ly ResultTree->Branch("RestPt",  RestPt,   "RestPt[nRestJ]/f");

	// Initialize tracking efficiency
	// ------------------------------
	//#ly if ( !isAuAu ) {
	//#ly 	tEff = new ktTrackEff();
	//#ly 	tEff->SetSysUncertainty(mEffUn);
	//#ly 	cout<<endl;
	//#ly 	tEff->PrintInfo();
	//#ly 	cout<<endl;
	//#ly }

	// Initialize analysis class
	// -------------------------
	//cout<<"initialize analysis"<<endl;		// test

	//TString OutFileName = "test.root"; //test 
	//float R = 0.6;	// test
	UnderlyingAna *ula = new UnderlyingAna( R,
			AjParameters::max_track_rap,
			AjParameters::dPhiCut,
			OutFileName
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

	if(jetchargecode!=1) ula->SetToMatchJetTrigger(true);			// whether match jet found with fastjet with the location which fired the trigger, NEED TO CHECK TrigFlagId
	if(jetchargecode==2) ula->SetNetraulJetFracCut(true);			// whether apply neutral energy fraction in jet cut
	ula->SetUnderlyingParticleCharge(underlyingchargecode);			// underlying event charge: 0 for netural, 1 for charged, 2 for all
	ula->SetDiJetAngle(0);					// Use Dijet angle (1) or Monojet angle (0) 

	std::vector<EtaPhiPair> TrigLoc2Match;		// trigger of High Tower or Jet Patch

	// Cycle through events
	// --------------------
	vector<PseudoJet> particles;
	TStarJetVectorContainer<TStarJetVector>* container;
	TStarJetVector* sv; // TLorentzVector* would be sufficient. 
	PseudoJet pj;


	//int nHardDijets = 0;
	//int nCorrespondingLowDijets = 0;
	//int nMatchedDijets=0;

	Long64_t nEvents=-1; // -1 for all
	//nEvents=10000;	// test
	cout<<"init..."<<endl;
	reader.Init(nEvents);
	int count = 0;

	//#ly fastjet::GhostedAreaSpec TmpArea; // for access to static random seed
	//#ly vector<int> SeedStatus;

	//#ly fastjet::Selector GrabCone = fastjet::SelectorCircle( R );    

	try{
		while ( reader.NextEvent() ) {
			reader.PrintStatus(10);
			if(count%1000==0) cout<<"event "<<count<<endl;
			count++;

			// event info
			// ----------
			//cout<<"load event header"<<endl;
			TStarJetPicoEventHeader* header = reader.GetEvent()->GetHeader();

			// eventid = header->GetEventId();
			// runid   = header->GetRunId();

			//#ly // Let's use the eventid as random seed.
			//#ly // that way things stay reproducible between different trees
			//#ly // but at the same time there's enough randomness
			//#ly gRandom->SetSeed( eventid ); 

			//#ly // NEW 05/07/15: For repeatability across different picoDSTs, set random seed
			//#ly // Static member, so we can set it here
			//#ly // Annoyingly, the getter and setter isn't static, so we need to instantiate
			//#ly // Apparently, the seed is always an int[2], so it's natural to seed it with runid and eventid      
			//#ly TmpArea.get_random_status(SeedStatus);
			//#ly if ( SeedStatus.size() !=2 ) {
			//#ly 	throw std::string("SeedStatus.size() !=2");
			//#ly 	return -1;
			//#ly } 
			//#ly SeedStatus.at(0) = runid   + randomoff;
			//#ly SeedStatus.at(1) = eventid + randomoff;
			//#ly TmpArea.set_random_status(SeedStatus);

			//#ly refmult=0;
			//#ly if ( isAuAu ) refmult=header->GetGReferenceMultiplicity();

			// Load event ht/jetpatch trigger objs
			// ----------
			//std::cout<<"load trigger objs"<<endl;
			TClonesArray *trigobj = reader.GetEvent()->GetTrigObjs();
			for(int itrg = 0; itrg<trigobj->GetEntries(); itrg++) {
				if( ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetTriggerFlag()==TrigFlagId )	 { 
					EtaPhiPair itrigloc =std::make_pair(CorrectBemcVzEta(((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetEta(),header->GetPrimaryVertexZ()), ((TStarJetPicoTriggerInfo *)((*trigobj)[itrg]))->GetPhi()) ;
					TrigLoc2Match.push_back(itrigloc);
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
				if(jetchargecode==1&&sv->GetCharge()==0) continue;		// ChargeJet
				if(jetchargecode==0&&sv->GetCharge()==1) continue;		// NeutralJet
				if (sv->GetCharge()==0 ) (*sv) *= fTowScale; // for systematics
				pj=MakePseudoJet( sv );
				pj.set_user_info ( new JetAnalysisUserInfo( 3*sv->GetCharge() ) );
				pj.set_user_index(ip);		// #ly	link fastjet::PseudoJet to TStarJetVector class

				//if ( sv->GetCharge()!=0 && tEff ) {		//#ly not apply Au+Au eff to pp, we are now calculate pp only
				//  Double_t reff=tEff->EffRatio_20(sv->Eta(),sv->Pt());
				//  Double_t mran=gRandom->Uniform(0,1);
				//  // cout << reff << "  " << mran << endl;
				//  if (mran<reff)  {
				//    particles.push_back ( pj );
				//  }
				//} else { // no charge or no efficiency class
				particles.push_back ( pj );
				//}	      
			}    

			// Run analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	// test
			//int ret;      
			//ret =

			ula->AnalyzeAndFill( particles, 
					*container, 
					reader,
					TrigLoc2Match
					//refmult
					//LeadJetPt, SubJetPt
					//LeadJetNtrkvsLeadJetPt, SubJetNtrkvsLeadJetPt, TranMaxNtrkvsLeadJetPt, TranMinNtrkvsLeadJetPt, TranNtrkvsLeadJetPt, 
					//LeadJetPtvsLeadJetPt, SubJetPtvsLeadJetPt, TranMaxPtvsLeadJetPt, TranMinPtvsLeadJetPt, TranPtvsLeadJetPt, 
					//TranPionNtrkvsLeadJetPt, TranPionPtvsLeadJetPt, TranProtonNtrkvsLeadJetPt, TranProtonPtvsLeadJetPt,
					//Spectrum_LeadJetPtvsLeadJetPt, Spectrum_SubJetPtvsLeadJetPt, Spectrum_TranMaxPtvsLeadJetPt, Spectrum_TranMinPtvsLeadJetPt, Spectrum_TranPtvsLeadJetPt,
			);


			//cout<<" -- "<<ret<<endl;		// test
			
			//#ly no need anymore. the rest are moved to UnderlyingAna.cxx
			//#ly switch ( ret ){
			//#ly 	case 1 : nHardDijets++;
			//#ly 		 // FALLTHROUGH
			//#ly 	case 0 : // Nothing found 
			//#ly 		 break;
			//#ly 	default :
			//#ly 		 cerr << "Unrecognized return value!" << endl;
			//#ly 		 throw(-1);
			//#ly 		 return -1;
			//#ly 		 break;      
			//#ly }


			// save event info
			//#ly if ( ret>0 ) {
			//#ly 	csize->Fill(ula.GetConstituents().size());
			//#ly 	for (int i =0; i< ula.GetConstituents().size(); ++i ){
			//#ly 		ptphieta->Fill( ula.GetConstituents().at(i).pt(), ula.GetConstituents().at(i).phi(), ula.GetConstituents().at(i).eta() );
			//#ly 	}
			//#ly }

			// #ly now this part is moved into UnderlyingAna.cxx
			// Save results
			//cout<<"Save result"<<endl; 		// test
			//vector<PseudoJet> DiJets = ula.GetDiJets();
			//vector<PseudoJet> JetsResult = ula.GetJAResult();

			//// This Part needs to be rewrite. 
			//// Currently if there is not dijet, rho & area info are not recorded properly
			//if ( JetsResult.size()>0 ) {
			//	if ( DiJets.size()==2 ){
			//		//cout<<"Dijet found"<<endl;		//test
			//		j1 = MakeTLorentzVector( DiJets.at(0) );
			//		j2 = MakeTLorentzVector( DiJets.at(1) );

			//		// DEBUG
			//		rho=rhoerr=j1area=j2area=0;
			//		fastjet::Selector selector_bkgd = fastjet::SelectorAbsRapMax( 0.6 ) * (!fastjet::SelectorNHardest(2));
			//		if ( refmult >= 0 ){
			//			//cout<<"cal rho"<<endl;		//test
			//			rho=ula.GetJA()->GetBackgroundEstimator()->rho() ;
			//			rhoerr=ula.GetJA()->GetBackgroundEstimator()->sigma() ;

			//			hrho->Fill(ula.GetJA()->GetBackgroundEstimator()->rho()) ;
			//			hrhoerr->Fill(ula.GetJA()->GetBackgroundEstimator()->sigma()) ;

			//			j1area=DiJets.at(0).area();
			//			j2area=DiJets.at(1).area();
			//			hj1area->Fill(DiJets.at(0).area());
			//			hj2area->Fill(DiJets.at(1).area());

			//			vector<PseudoJet> lojets = fastjet::sorted_by_pt( selector_bkgd ( ula.GetJA()->inclusive_jets() ) );
			//			for ( int i=2; i <lojets.size() ; ++ i ){
			//				hrestarea->Fill(lojets.at(i).area());
			//			} 

			//			//cout << DiJets.at(0).pt()<< "  "  << j1.Pt() << "  " << lojets.at(0).pt() << endl;	// test
			//			nRestJ=0;
			//			for ( int i=0; i<lojets.size()-2 ; ++i ){
			//				if ( lojets.at(i+2).pt() < 1e-4 ) continue;
			//				RestArea[i] = lojets.at(i+2).area();
			//				RestPt[i] = lojets.at(i+2).pt();
			//				nRestJ++;
			//			}

			//		}
			//	}
			//	else {
			//		j1 = MakeTLorentzVector( JetsResult.at(0) );
			//	}

			//	//cout<<"get underlying info"<<endl;		//test
			//	leadpt = ula.GetLeadJetPt();
			//	subpt = ula.GetSubJetPt();
			//	tranpt = ula.GetTranPt();
			//	tranmaxpt = ula.GetTranMaxPt();
			//	tranminpt = ula.GetTranMinPt();

			//	leadntrk = ula.GetLeadJetNtrk();
			//	subntrk = ula.GetLeadJetNtrk();
			//	tranntrk = ula.GetTranNtrk();
			//	tranmaxntrk = ula.GetTranMaxNtrk();
			//	tranminntrk = ula.GetTranMinNtrk();

			//	//cout<<"fill tree"<<endl;			//test
			//	ResultTree->Fill();
			//}


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


	//fout->Write();

	//cout << "In " << nEventsUsed << " events, found " << endl
	//	<< nHardDijets << " dijets"  << endl ; 


	//cout << "Wrote to " << fout->GetName() << endl;

	//delete ula;
	cout << "Bye." << endl;

	// cout << "Aj->GetEntries() = " << AJ_hi->GetEntries() << endl;
	// cout << "Unmatched GetEntries() = " << UnmatchedAJ_hi->GetEntries() << endl;

	//cout << "TranNtrkvsLeadJetPt->GetEntries() = " << TranNtrkvsLeadJetPt->GetEntries() << endl;
	return 0;

}

