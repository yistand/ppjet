//============================================================================================================================================
//
//		2016.08.31		Li YI
//		Modified code from leading jet pt distribution only to underlying activity vs leading jet pt for 2D unfolding later on
//		Events will be recorded when either Mc *OR* Embed has a jet.
//
//		The following flag should be stored in the tree for later unfolding purpose:
//		event flag:
//		1. Rc is JP2 event? 		// input from main macro
//		2. Leading Rc jet is matched to jet patch eta, phi?
//		jet flag
//		1. One-to-one good match for leading to leading
//		2. Embed (Rc) leading jet matched to a Mc subleading jet (no eta requirement) and its Mc leading jet is in good eta.
//		3. Rc jet matched to any Mc jet (no eta requirement)
//
//
//
//============================================================================================================================================
//
//		2016.07.14		Li YI
//		jetpt_McVsEmbed.cxx 
//		read pythia Embedded zerobias event
//		obtain jet_particle (MC) and jet_detector (embedding) covariation matrix for furture unfolding jet pT
//
//		first do no fast detector matching, but try to add the option possible
//
//
//		2016.07.25	Li Yi
//		pt2_3 has one entry with extremely high Mc pt. This may cause trouble when merging
//		because low pT bin will have huge weight, therefore we delete this bin
//
//
//============================================================================================================================================

#include "underlying_McVsEmbed.hh"
#include <string>

// Standard ctor
underlying_McVsEmbed::underlying_McVsEmbed (	double R, 
					double max_const_rap,
					double min_const_pt,
					std::string jetalgo,
					TString name

					) 
					: R(R),
					max_const_rap (max_const_rap),
					min_const_pt (min_const_pt),
					dPhiCut(0.4),		// Dijet opening angle requirement. Accept only  |&phi;1 - &phi;2 - &pi;| < dPhiCut
					OutFileName (name),
					mVerbose(0)
{

	// derived rapidity cuts
	// ---------------------
	max_rap      = max_const_rap-R;
	ghost_maxrap = max_rap + 2.0 * R;

	// Constituent selectors
	// ---------------------
	select_const_rap = fastjet::SelectorAbsRapMax(max_const_rap);
	select_const_ptmin = fastjet::SelectorPtMin(min_const_pt);

	// Constituent Selector for jet finding
	// ---------------------
	//2017.01.10  Mcsconst     = fastjet::SelectorIdentity(); //select_const_ptmin;			
	Mcsconst     = select_const_ptmin;			
	Rcsconst     = select_const_rap && select_const_ptmin;				// detector eta acceptance

	// Constituent Selector for underlying events
	// ---------------------
	//2017.01.10  McsUconst     = select_const_rap; // fastjet::SelectorIdentity(); //select_const_ptmin;		For underlying event quantities, still apply detector eta cut
	McsUconst     = select_const_rap && select_const_ptmin;		 
	RcsUconst     = select_const_rap && select_const_ptmin;				// detector eta acceptance


	// Jet candidate selectors
	// ---------------------
	select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);			// select objects with |rap| <= absrapmax 
	sjet = select_jet_rap;
	NoGhosts = !fastjet::SelectorIsPureGhost();

	// Choose a jet and area definition
	// --------------------------------
	std::cout<<"INFO: jetalgo = "<<jetalgo<<std::endl;
	jet_def = fastjet::JetDefinition(AlgoFromString(jetalgo), R);	// in JetAnalyzer.cxx

	// create an area definition for the clustering
	//----------------------------------------------------------
	// ghosts should go up to the acceptance of the detector or
	// (with infinite acceptance) at least 2R beyond the region
	// where you plan to investigate jets.
	area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, AjParameters::ghost_repeat, AjParameters::ghost_area );
	area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

	// set up output root file
	//----------------------------------------------------------
	fout = new TFile(OutFileName,"RECREATE");	// just the name from input 
	assert ( fout->IsOpen() );
	std::cout << " ################################################### " << std::endl;
	std::cout << "Writing to: " << fout->GetName() << std::endl;
	std::cout << " ################################################### " << std::endl;


	// event counter
	nEventPassed = 0;  
	nEventOutlierMcpT = 0;  		// how many events were cut off because too large Mc pT
	nEventOutlierRcpT = 0;  		// how many events were cut off because too large Rc pT

	// detector level jet parameters
	mNeedToMatchTrig = true;
	mNeutralJetFracCut = true;

	// OutlierMcpTCut flag
	DoOutlierMcpTCut = false;
	mOutlierMcpTCut = 99999;
	DoOutlierRcpTCut = false;
	mOutlierRcpTCut = 99999;

	// Jet Charge
	mJetCharge = 2;			// Full jet
	// Underlying event parameters
	mUnderlyingParticleCharge = 1;	// Charged only
	mTranPhiSize = 60;


	// output
	//----------------------------------------------------------
	// Tree
	ResultTree = NULL;				// event will be filled if there is >=1 Mc jet found

	// Histogram
	McMatchedLeadJetPt = NULL;			// Mc Jet, when Detector Leading and Particle level Leading Jet geometrical matched, take the leading Mc/Particle-level jet
	McMatchedLeadOrSubJetPt = NULL;			// Mc Jet, when Detector Leading and Particle level Leading or subleading Jet geometrical matched, take the leading Mc/Particle-level jet (Leading particle-level jet is required to be inside good acceptance, subleading is not as long as it matched to detector-level)
	CovMcMatchedLeadJetVsRcJet = NULL;
	CovMcMatchedLeadOrSubJetVsRcJet = NULL;
	McLeadJetPt = NULL;				// Mc Leading Jet, no matter whether there is a detector level jet matched to it or not
	McSubLeadJetPt = NULL;				// Mc SubLeading Jet, no matter whether there is a detector level jet matched to it or not
	RcLeadJetPt = NULL;				// Rc Leading Jet, no matter whether there is a particle level jet matched to it or not
	RcSubLeadJetPt = NULL;				// Rc SubLeading Jet, no matter whether there is a particle level jet matched to it or not
	
}


underlying_McVsEmbed::~underlying_McVsEmbed() 
{
	delete ResultTree;
	delete McMatchedLeadJetPt;
	delete McMatchedLeadOrSubJetPt;
	delete CovMcMatchedLeadJetVsRcJet;
	delete CovMcMatchedLeadOrSubJetVsRcJet;
	delete McLeadJetPt;
	delete McSubLeadJetPt;
	delete RcLeadJetPt;
	delete RcSubLeadJetPt;
}


int underlying_McVsEmbed::Init() 
{

	// set up TTree
	ResultTree=new TTree("ResultTree","Result Jets");

	ResultTree->Branch("eventid",&eventid, "eventid/I");
	ResultTree->Branch("runid",&runid, "runid/I");

	//ResultTree->Branch("weight",&weight, "weight/D");


	// Match flag
	// event fired
	ResultTree->Branch("flagIsTrigger",&flagIsTrigger, "flagIsTrigger/O");	// caps letter o for Bool_t
	// Match to trig
	ResultTree->Branch("trigmatch",&trigmatch, "trigmatch/O");   
	// matched to leading Mc
	ResultTree->Branch("flagMatch2Lead",&flagMatch2Lead, "flagMatch2Lead/O");	// caps letter o for Bool_t
	// matched to subleading Mc and leading Mc inside good eta
	ResultTree->Branch("flagMatch2Sub",&flagMatch2Sub, "flagMatch2Sub/O");	// o for Bool_t
	// matched to any Mc jet
	ResultTree->Branch("flagMatch2McJet",&flagMatch2McJet, "flagMatch2McJet/O");	// o for Bool_t

	// Particle level
	ResultTree->Branch("Mcrefmult",&Mcrefmult, "Mcrefmult/D");
	ResultTree->Branch("Mcvz",&Mcvz, "Mcvz/D");

	ResultTree->Branch("Mcrho",&Mcrho, "Mcrho/F");
	ResultTree->Branch("Mcrhoerr",&Mcrhoerr, "Mcrhoerr/F");

	ResultTree->Branch("Mcj1pt",&Mcj1pt, "Mcj1pt/F");
	ResultTree->Branch("Mcjaspt",&Mcjaspt, "Mcjaspt/F");
	ResultTree->Branch("Mcj2pt",&Mcj2pt, "Mcj2pt/F");
	ResultTree->Branch("Mcj1phi",&Mcj1phi, "Mcj1phi/F");
	ResultTree->Branch("Mcjasphi",&Mcjasphi, "Mcjasphi/F");
	ResultTree->Branch("Mcj2phi",&Mcj2phi, "Mcj2phi/F");
	ResultTree->Branch("Mcj1eta",&Mcj1eta, "Mcj1eta/F");
	ResultTree->Branch("Mcjaseta",&Mcjaseta, "Mcjaseta/F");
	ResultTree->Branch("Mcj2eta",&Mcj2eta, "Mcj2eta/F");
	ResultTree->Branch("Mcj1area",&Mcj1area, "Mcj1area/F");
	ResultTree->Branch("Mcjasarea",&Mcjasarea, "Mcjasarea/F");
	ResultTree->Branch("Mcj2area",&Mcj2area, "Mcj2area/F");
	ResultTree->Branch("Mcj1area_err",&Mcj1area_err, "Mcj1area_err/F");
	ResultTree->Branch("Mcjasarea_err",&Mcjasarea_err, "Mcjasarea_err/F");
	ResultTree->Branch("Mcj2area_err",&Mcj2area_err, "Mcj2area_err/F");

	ResultTree->Branch("Mcj1neutralfrac",&Mcj1neutralfrac, "Mcj1neutralfrac/F");


	ResultTree->Branch("Mcj3pt",&Mcj3pt, "Mcj3pt/F");
	ResultTree->Branch("Mcj4pt",&Mcj4pt, "Mcj4pt/F");
	ResultTree->Branch("Mcj3phi",&Mcj3phi, "Mcj3phi/F");
	ResultTree->Branch("Mcj4phi",&Mcj4phi, "Mcj4phi/F");
	ResultTree->Branch("Mcj3eta",&Mcj3eta, "Mcj3eta/F");
	ResultTree->Branch("Mcj4eta",&Mcj4eta, "Mcj4eta/F");


	ResultTree->Branch("Mcj1constntrk",&Mcj1constntrk, "Mcj1constntrk/I");
	ResultTree->Branch("Mcj1constpt",Mcj1constpt, "Mcj1constpt[Mcj1constntrk]/F");
	ResultTree->Branch("Mcj1constphi",Mcj1constphi, "Mcj1constphi[Mcj1constntrk]/F");
	ResultTree->Branch("Mcj1consteta",Mcj1consteta, "Mcj1consteta[Mcj1constntrk]/F");

	// Mc underlying 
	ResultTree->Branch("McLeadAreaPtSum",&McLeadAreaPt,"McLeadAreaPtSum/F");
	ResultTree->Branch("McSubLeadAreaPtSum",&McSubAreaPt,"McSubLeadAreaPtSum/F");
	ResultTree->Branch("McTranPtSum",&McTranPt,"McTranPtSum/F");
	ResultTree->Branch("McTranMaxPtSum",&McTranMaxPt,"McTranMaxPtSum/F");
	ResultTree->Branch("McTranMinPtSum",&McTranMinPt,"McTranMinPtSum/F");
	ResultTree->Branch("McLeadAreaNtrk",&McLeadAreaNtrk,"McLeadAreaNtrk/I");
	ResultTree->Branch("McSubAreaNtrk",&McSubAreaNtrk,"McSubAreaNtrk/I");
	ResultTree->Branch("McTranMaxNtrk",&McTranMaxNtrk,"McTranMaxNtrk/I");
	ResultTree->Branch("McTranMinNtrk",&McTranMinNtrk,"McTranMinNtrk/I");
	ResultTree->Branch("McTrkTranMaxPt",McTrkTranMaxPt,"McTrkTranMaxPt[McTranMaxNtrk]/F");
	ResultTree->Branch("McTrkTranMaxPhi",McTrkTranMaxPhi,"McTrkTranMaxPhi[McTranMaxNtrk]/F");
	ResultTree->Branch("McTrkTranMaxEta",McTrkTranMaxEta,"McTrkTranMaxEta[McTranMaxNtrk]/F");
	ResultTree->Branch("McTrkTranMaxId",McTrkTranMaxId,"McTrkTranMaxId[McTranMaxNtrk]/I");
	ResultTree->Branch("McTrkTranMinPt",McTrkTranMinPt,"McTrkTranMinPt[McTranMinNtrk]/F");
	ResultTree->Branch("McTrkTranMinPhi",McTrkTranMinPhi,"McTrkTranMinPhi[McTranMinNtrk]/F");
	ResultTree->Branch("McTrkTranMinEta",McTrkTranMinEta,"McTrkTranMinEta[McTranMinNtrk]/F");
	ResultTree->Branch("McTrkTranMinId",McTrkTranMinId,"McTrkTranMinId[McTranMinNtrk]/I");
	ResultTree->Branch("McTrkLeadAreaPt",McTrkLeadAreaPt,"McTrkLeadAreaPt[McLeadAreaNtrk]/F");
	ResultTree->Branch("McTrkLeadAreaPhi",McTrkLeadAreaPhi,"McTrkLeadAreaPhi[McLeadAreaNtrk]/F");
	ResultTree->Branch("McTrkLeadAreaEta",McTrkLeadAreaEta,"McTrkLeadAreaEta[McLeadAreaNtrk]/F");
	ResultTree->Branch("McTrkLeadAreaId",McTrkLeadAreaId,"McTrkLeadAreaId[McLeadAreaNtrk]/I");
	ResultTree->Branch("McTrkSubAreaPt",McTrkSubAreaPt,"McTrkSubAreaPt[McSubAreaNtrk]/F");
	ResultTree->Branch("McTrkSubAreaPhi",McTrkSubAreaPhi,"McTrkSubAreaPhi[McSubAreaNtrk]/F");
	ResultTree->Branch("McTrkSubAreaEta",McTrkSubAreaEta,"McTrkSubAreaEta[McSubAreaNtrk]/F");
	ResultTree->Branch("McTrkSubAreaId",McTrkSubAreaId,"McTrkSubAreaId[McSubAreaNtrk]/I");
 

	// Detector level
	ResultTree->Branch("Rcrefmult",&Rcrefmult, "Rcrefmult/D");
	ResultTree->Branch("Rcvz",&Rcvz, "Rcvz/D");

	ResultTree->Branch("Rcrho",&Rcrho, "Rcrho/F");
	ResultTree->Branch("Rcrhoerr",&Rcrhoerr, "Rcrhoerr/F");

	ResultTree->Branch("Rcj1pt",&Rcj1pt, "Rcj1pt/F");
	ResultTree->Branch("Rcjaspt",&Rcjaspt, "Rcjaspt/F");
	ResultTree->Branch("Rcj2pt",&Rcj2pt, "Rcj2pt/F");
	ResultTree->Branch("Rcj1phi",&Rcj1phi, "Rcj1phi/F");
	ResultTree->Branch("Rcjasphi",&Rcjasphi, "Rcjasphi/F");
	ResultTree->Branch("Rcj2phi",&Rcj2phi, "Rcj2phi/F");
	ResultTree->Branch("Rcj1eta",&Rcj1eta, "Rcj1eta/F");
	ResultTree->Branch("Rcjaseta",&Rcjaseta, "Rcjaseta/F");
	ResultTree->Branch("Rcj2eta",&Rcj2eta, "Rcj2eta/F");
	ResultTree->Branch("Rcj1area",&Rcj1area, "Rcj1area/F");
	ResultTree->Branch("Rcjasarea",&Rcjasarea, "Rcjasarea/F");
	ResultTree->Branch("Rcj2area",&Rcj2area, "Rcj2area/F");
	ResultTree->Branch("Rcj1area_err",&Rcj1area_err, "Rcj1area_err/F");
	ResultTree->Branch("Rcjasarea_err",&Rcjasarea_err, "Rcjasarea_err/F");
	ResultTree->Branch("Rcj2area_err",&Rcj2area_err, "Rcj2area_err/F");

	ResultTree->Branch("Rcj1neutralfrac",&Rcj1neutralfrac, "Rcj1neutralfrac/F");


	ResultTree->Branch("Rcj3pt",&Rcj3pt, "Rcj3pt/F");
	ResultTree->Branch("Rcj4pt",&Rcj4pt, "Rcj4pt/F");
	ResultTree->Branch("Rcj3phi",&Rcj3phi, "Rcj3phi/F");
	ResultTree->Branch("Rcj4phi",&Rcj4phi, "Rcj4phi/F");
	ResultTree->Branch("Rcj3eta",&Rcj3eta, "Rcj3eta/F");
	ResultTree->Branch("Rcj4eta",&Rcj4eta, "Rcj4eta/F");


	ResultTree->Branch("Rcj1constntrk",&Rcj1constntrk, "Rcj1constntrk/I");
	ResultTree->Branch("Rcj1constpt",Rcj1constpt, "Rcj1constpt[Rcj1constntrk]/F");
	ResultTree->Branch("Rcj1constphi",Rcj1constphi, "Rcj1constphi[Rcj1constntrk]/F");
	ResultTree->Branch("Rcj1consteta",Rcj1consteta, "Rcj1consteta[Rcj1constntrk]/F");


	// Rc underlying events
	ResultTree->Branch("RcLeadAreaPtSum",&RcLeadAreaPt,"RcLeadAreaPtSum/F");
	ResultTree->Branch("RcSubLeadAreaPtSum",&RcSubAreaPt,"RcSubLeadAreaPtSum/F");
	ResultTree->Branch("RcTranPtSum",&RcTranPt,"RcTranPtSum/F");
	ResultTree->Branch("RcTranMaxPtSum",&RcTranMaxPt,"RcTranMaxPtSum/F");
	ResultTree->Branch("RcTranMinPtSum",&RcTranMinPt,"RcTranMinPtSum/F");
	ResultTree->Branch("RcLeadAreaNtrk",&RcLeadAreaNtrk,"RcLeadAreaNtrk/I");
	ResultTree->Branch("RcSubAreaNtrk",&RcSubAreaNtrk,"RcSubAreaNtrk/I");
	ResultTree->Branch("RcTranMaxNtrk",&RcTranMaxNtrk,"RcTranMaxNtrk/I");
	ResultTree->Branch("RcTranMinNtrk",&RcTranMinNtrk,"RcTranMinNtrk/I");
	ResultTree->Branch("RcTrkTranMaxPt",RcTrkTranMaxPt,"RcTrkTranMaxPt[RcTranMaxNtrk]/F");
	ResultTree->Branch("RcTrkTranMaxPhi",RcTrkTranMaxPhi,"RcTrkTranMaxPhi[RcTranMaxNtrk]/F");
	ResultTree->Branch("RcTrkTranMaxEta",RcTrkTranMaxEta,"RcTrkTranMaxEta[RcTranMaxNtrk]/F");
	ResultTree->Branch("RcTrkTranMaxMcId",RcTrkTranMaxMcId,"RcTrkTranMaxMcId[RcTranMaxNtrk]/I");
	ResultTree->Branch("RcTrkTranMinPt",RcTrkTranMinPt,"RcTrkTranMinPt[RcTranMinNtrk]/F");
	ResultTree->Branch("RcTrkTranMinPhi",RcTrkTranMinPhi,"RcTrkTranMinPhi[RcTranMinNtrk]/F");
	ResultTree->Branch("RcTrkTranMinEta",RcTrkTranMinEta,"RcTrkTranMinEta[RcTranMinNtrk]/F");
	ResultTree->Branch("RcTrkTranMinMcId",RcTrkTranMinMcId,"RcTrkTranMinMcId[RcTranMinNtrk]/I");
	ResultTree->Branch("RcTrkLeadAreaPt",RcTrkLeadAreaPt,"RcTrkLeadAreaPt[RcLeadAreaNtrk]/F");
	ResultTree->Branch("RcTrkLeadAreaPhi",RcTrkLeadAreaPhi,"RcTrkLeadAreaPhi[RcLeadAreaNtrk]/F");
	ResultTree->Branch("RcTrkLeadAreaEta",RcTrkLeadAreaEta,"RcTrkLeadAreaEta[RcLeadAreaNtrk]/F");
	ResultTree->Branch("RcTrkLeadAreaMcId",RcTrkLeadAreaMcId,"RcTrkLeadAreaMcId[RcLeadAreaNtrk]/I");
	ResultTree->Branch("RcTrkSubAreaPt",RcTrkSubAreaPt,"RcTrkSubAreaPt[RcSubAreaNtrk]/F");
	ResultTree->Branch("RcTrkSubAreaPhi",RcTrkSubAreaPhi,"RcTrkSubAreaPhi[RcSubAreaNtrk]/F");
	ResultTree->Branch("RcTrkSubAreaEta",RcTrkSubAreaEta,"RcTrkSubAreaEta[RcSubAreaNtrk]/F");
	ResultTree->Branch("RcTrkSubAreaMcId",RcTrkSubAreaMcId,"RcTrkSubAreaMcId[RcSubAreaNtrk]/I");
 

	// Detector-level jet matched to any Particle-level jet
	ResultTree->Branch("MatchedNthMcj",&MatchedNthMcj, "MatchedNthMcj/I");
	ResultTree->Branch("MatchedMcjpt",&MatchedMcjpt, "MatchedMcjpt/F");
	ResultTree->Branch("MatchedMcjphi",&MatchedMcjphi, "MatchedMcjphi/F");
	ResultTree->Branch("MatchedMcjeta",&MatchedMcjeta, "MatchedMcjeta/F");
	ResultTree->Branch("LeadInMatchedMcjpt",&LeadInMatchedMcjpt, "LeadInMatchedMcjpt/F");
	ResultTree->Branch("LeadInMatchedMcjphi",&LeadInMatchedMcjphi, "LeadInMatchedMcjphi/F");
	ResultTree->Branch("LeadInMatchedMcjeta",&LeadInMatchedMcjeta, "LeadInMatchedMcjeta/F");
	// Particle-level jet matched to any Detector-level jet
	ResultTree->Branch("MatchedNthRcj",&MatchedNthRcj, "MatchedNthRcj/I");
	ResultTree->Branch("MatchedRcjpt",&MatchedRcjpt, "MatchedRcjpt/F");
	ResultTree->Branch("MatchedRcjphi",&MatchedRcjphi, "MatchedRcjphi/F");
	ResultTree->Branch("MatchedRcjeta",&MatchedRcjeta, "MatchedRcjeta/F");


	// set up histograms
	TH1::SetDefaultSumw2(true);
	TH2::SetDefaultSumw2(true);


	McMatchedLeadJetPt = new TH1D("McMatchedLeadJetPt","MC Matched Event Leading Jet p_{T}",100,0,100);
	McMatchedLeadOrSubJetPt = new TH1D("McMatchedLeadOrSubJetPt","MC Matched Event LeadOrSubing Jet p_{T}",100,0,100);		// if matched to subjet, still fill the leading jet pt
	CovMcMatchedLeadJetVsRcJet = new TH2D("CovMcMatchedLeadJetVsRcJet","MC Leading Jet p_{T} vs Smearing Matrix Matched Event Reconstructed Leading Jet p_{T} ",100,0,100,100,0,100 );
	CovMcMatchedLeadOrSubJetVsRcJet = new TH2D("CovMcMatchedLeadOrSubJetVsRcJet","MC Leading Or Subleading Jet p_{T} vs Smearing Matrix Matched Event Reconstructed Leading Jet p_{T} ",100,0,100,100,0,100 );
	McLeadJetPt = new TH1D("McLeadJetPt","MC Leading Jet p_{T} (Matched or Not)",100,0,100);
	McSubLeadJetPt = new TH1D("McSubLeadJetPt","MC Subleading Jet p_{T} (Matched or Not)",100,0,100);
	RcLeadJetPt = new TH1D("RcLeadJetPt","Rc Leading Jet p_{T} (Matched or Not)",100,0,100);
	RcSubLeadJetPt = new TH1D("RcSubLeadJetPt","Rc Subleading Jet p_{T} (Matched or Not)",100,0,100);

	return true;
}


// Main analysis method
int underlying_McVsEmbed::Make (	const std::vector<fastjet::PseudoJet>& Mcparticles,
				const std::vector<fastjet::PseudoJet>& Rcparticles,
				int ineventid, int inrunid, 
				double inMcrefmult, double inMcvz,
				double inRcrefmult, double inRcvz,
				const std::vector<std::pair<float,float> > &ToMatch,			// trigger 
				Bool_t is_trigger
				//double weightbyXsec
				) 
{

	//std::cout<<"Make called"<<std::endl;

	McJconstituents.clear();
	RcJconstituents.clear();

	McUconstituents.clear();
	RcUconstituents.clear();


	eventid = ineventid;
	runid = inrunid;

	Mcrefmult = inMcrefmult;	
	Mcvz = inMcvz;

	Rcrefmult = inRcrefmult;	
	Rcvz = inRcvz;


	Mcj1pt=0, Mcjaspt=0, Mcj2pt=0;
	Mcj1phi=-999, Mcjasphi=-999, Mcj2phi=-999;
	Mcj1eta=-999, Mcjaseta=-999, Mcj2eta=-999;
	Mcj1area=0, Mcjasarea=0, Mcj2area=0;
	Mcj1area_err=0, Mcjasarea_err=0, Mcj2area_err=0;
	Mcj1neutralfrac=999;

	Mcrho=0, Mcrhoerr=0;

	Mcj3pt=0, Mcj4pt=0;
	Mcj3phi=-999, Mcj4phi=-999;
	Mcj3eta=-999, Mcj4eta=-999;

	Mcj1constntrk=0;


	Rcj1pt=0, Rcjaspt=0, Rcj2pt=0;
	Rcj1phi=-999, Rcjasphi=-999, Rcj2phi=-999;
	Rcj1eta=-999, Rcjaseta=-999, Rcj2eta=-999;
	Rcj1area=0, Rcjasarea=0, Rcj2area=0;
	Rcj1area_err=0, Rcjasarea_err=0, Rcj2area_err=0;
	Rcj1neutralfrac=999;

	Rcrho=0, Rcrhoerr=0;

	Rcj3pt=0, Rcj4pt=0;
	Rcj3phi=-999, Rcj4phi=-999;
	Rcj3eta=-999, Rcj4eta=-999;

	Rcj1constntrk=0;

	MatchedNthMcj=-1, MatchedMcjpt=0, MatchedMcjphi=-999, MatchedMcjeta=-999;
	LeadInMatchedMcjpt=0, LeadInMatchedMcjphi=-999, LeadInMatchedMcjeta=-999;

	MatchedNthRcj=-1, MatchedRcjpt=0, MatchedRcjphi=-999, MatchedRcjeta=-999;

	flagIsTrigger = is_trigger;
	//weight = weightbyXsec;

	// underlying event info zero
	McLeadAreaPt = 0;
	McSubAreaPt = 0;
	McTranMaxPt = 0;
	McTranMinPt = 0;
	McTranPt = 0;
	
	McLeadAreaNtrk = 0;
	McSubAreaNtrk = 0;
	McTranMaxNtrk = 0;
	McTranMinNtrk = 0;
		
	RcLeadAreaPt = 0;
	RcSubAreaPt = 0;
	RcTranMaxPt = 0;
	RcTranMinPt = 0;
	RcTranPt = 0;
	
	RcLeadAreaNtrk = 0;
	RcSubAreaNtrk = 0;
	RcTranMaxNtrk = 0;
	RcTranMinNtrk = 0;
		
	// Select particles to perform analysis on
	// ---------------------------------------
	// Constituent Selector for jet finding
	McJconstituents = Mcsconst( Mcparticles );
	RcJconstituents = Rcsconst( Rcparticles );
	// Constituent Selector for underlying events 
	McUconstituents = McsUconst( Mcparticles );
	RcUconstituents = RcsUconst( Rcparticles );


	// Background selector
	// -------------------
	fastjet::Selector selector_bkgd = fastjet::SelectorAbsRapMax( max_rap ) * (!fastjet::SelectorNHardest(2));	//selector for fastjet::JetMedianBackgroundEstimator: the selection criterion is typically a geometrical one (e.g. all jets with |y|<2) sometimes supplemented with some kinematical restriction (e.g. exclusion of the two hardest jets). It is passed to the class through a Selector.


	// find jets
	// -----------------------------
	// NO background subtraction && NO jet eta constrain 
	McnosjetJA = new JetAnalyzer( McJconstituents, jet_def, area_def);
	McnosjetJAResult = fastjet::sorted_by_pt(NoGhosts(McnosjetJA->inclusive_jets()));


	flagGoodEtaMcJet = false; 
	if(McnosjetJAResult.size()>=1 && fabs(McnosjetJAResult.at(0).eta())<=max_rap) { // Mc leading jet inside good eta acceptance
		flagGoodEtaMcJet = true; 	
		McJAResult = fastjet::sorted_by_pt(sjet(McnosjetJAResult));
	}
	else {
		McJAResult.clear();
	}

	// Outlier Mc Pt Cut before proceeds to next step 
	if(DoOutlierMcpTCut && McJAResult.size()>=1 && McJAResult.at(0).pt()>mOutlierMcpTCut) {		// if we want to cut off the outlier MC jt (it can cause big effect in low pT bin as they get huge weigth)
// this should not be performed for more than 1 or 2 events, so we count and should be checked when finished
		std::cout<<"INFO: Mc leading pt = "<<McJAResult.at(0).pt()<<" does not pass mOutlierMcpTCut "<<mOutlierMcpTCut<<": Event Reject!!!!"<<std::endl;
		nEventOutlierMcpT++;
		return 1;	
	}





	RcnosjetJA = new JetAnalyzer( RcJconstituents, jet_def, area_def);			// no jet eta selection. but constituents still have eta cut (TPC acceptance) 
	RcnosjetJAResult = fastjet::sorted_by_pt(NoGhosts(RcnosjetJA->inclusive_jets()));
	flagGoodEtaRcJet = false;	// whether Rc jet the hardest one is inside jet eta cut
	if(RcnosjetJAResult.size()>=1 && fabs(RcnosjetJAResult.at(0).eta())<=max_rap) { 
		flagGoodEtaRcJet = true;
		RcJAResult = fastjet::sorted_by_pt(sjet(RcnosjetJAResult));
	}
	else {
		RcJAResult.clear();
	}

	// Outlier Rc Pt Cut before proceeds to next step 
	if(DoOutlierRcpTCut && RcJAResult.size()>=1 && RcJAResult.at(0).pt()>mOutlierRcpTCut) {		// if we want to cut off the outlier Rc jt (it can cause big effect in low pT bin as they get huge weigth)
// this should not be performed for more than 1 or 2 events, so we count and should be checked when finished
		std::cout<<"INFO: Rc leading pt = "<<RcJAResult.at(0).pt()<<" does not pass mOutlierRcpTCut "<<mOutlierRcpTCut<<": Event Reject!!!!"<<std::endl;
		nEventOutlierRcpT++;
		return 1;	
	}


	// WITH subtract background for rho estimation
	McJA_bkgsub = new JetAnalyzer( McJconstituents, jet_def , area_def, selector_bkgd);
	fastjet::Subtractor* McBackgroundSubtractor =  McJA_bkgsub->GetBackgroundSubtractor();
	McJAResult_bkgsub = fastjet::sorted_by_pt( sjet ((*McBackgroundSubtractor) ( McJA_bkgsub->inclusive_jets()) ) ); // with background subtraction

	RcJA_bkgsub = new JetAnalyzer( RcJconstituents, jet_def , area_def, selector_bkgd);
	fastjet::Subtractor* RcBackgroundSubtractor =  RcJA_bkgsub->GetBackgroundSubtractor();
	RcJAResult_bkgsub = fastjet::sorted_by_pt( sjet ((*RcBackgroundSubtractor) ( RcJA_bkgsub->inclusive_jets()) ) ); // with background subtraction


	// Any Match? 
	// ---------------------
	flagMatch2Lead = false;
	flagMatch2Sub = false;
	flagMatch2McJet = false;
	MatchedNthMcj = -1;
	bool flagMatch2RcJet = false;
	MatchedNthRcj = -1;

	if(flagGoodEtaMcJet && flagGoodEtaRcJet) {			// only matching Rc and Mc both in good eta acceptance
		if(RcJAResult.at(0).delta_R(McJAResult.at(0))<R) {		// Matched to leading one
			flagMatch2Lead = true;
			flagMatch2McJet = true;
			MatchedNthMcj = 0;	// the first hardest matched
		}
		//else if(McnosjetJAResult.size()>=2 && fabs(McnosjetJAResult.at(1).eta())<=max_rap && RcJAResult.at(0).delta_R(McnosjetJAResult.at(1))<R) {		// Matched to the subleading one who also inside good eta acceptance
		else if(McnosjetJAResult.size()>=2 && RcJAResult.at(0).delta_R(McnosjetJAResult.at(1))<R) {		// Matched to the subleading one No need to be inside good eta acceptance
			flagMatch2Sub = true;
			flagMatch2McJet = true;
			MatchedNthMcj = 1;	// the second hardest matched
		}

	}
	else if(flagGoodEtaRcJet){		// check whether Rc jet has any matched Mc jet, no need to be in good acceptance. This is for QA purpose
		for(int ijet = 0; ijet<McnosjetJAResult.size(); ijet++) {
			if(RcJAResult.at(0).delta_R(McnosjetJAResult.at(ijet))<R ) {		// match to any jet? 
				flagMatch2McJet = true;
				MatchedNthMcj = ijet;
			}
		}
	}
	// check whether good Mc jet (good for inside eta acceptance) has any matched Rc jet, no need to be in good acceptance. This is for QA purpose
	if(flagGoodEtaMcJet) {
		for(int ijet = 0; ijet<RcnosjetJAResult.size(); ijet++) {
			if(McJAResult.at(0).delta_R(RcnosjetJAResult.at(ijet))<R ) {		// match to any jet? 
				flagMatch2RcJet = true;
				MatchedNthRcj = ijet;
			}
		}
	}

	if(flagMatch2McJet) {
		MatchedMcjpt = McnosjetJAResult.at(MatchedNthMcj).pt();
		MatchedMcjeta = McnosjetJAResult.at(MatchedNthMcj).eta();
		MatchedMcjphi = McnosjetJAResult.at(MatchedNthMcj).phi();
		LeadInMatchedMcjpt = McnosjetJAResult.at(0).pt();
		LeadInMatchedMcjeta = McnosjetJAResult.at(0).eta();
		LeadInMatchedMcjphi = McnosjetJAResult.at(0).phi();
	}

	if(flagMatch2RcJet) {
		MatchedRcjpt = RcnosjetJAResult.at(MatchedNthRcj).pt();
		MatchedRcjeta = RcnosjetJAResult.at(MatchedNthRcj).eta();
		MatchedRcjphi = RcnosjetJAResult.at(MatchedNthRcj).phi();
	}

	// Do any Reconstructed (Rc) jets match to the one fired the trigger?
	// ---------------------------------------------------------
	bool flagtrigmatch = true;		// assuming no need to match. set to 1 for not care. 
	trigmatch = false;			// this one is the real value for matching to record in ResultTree 

	//std::cout<<std::endl<<"eventid = "<<eventid<<" Vz = "<<Rcvz<<std::endl;	
	//if(RcnosjetJAResult.size()>0) std::cout<<"nosjet jet pt "<<RcnosjetJAResult.at(0).pt()<<"at (eta,phi)=("<<RcnosjetJAResult.at(0).eta()<<","<<RcnosjetJAResult.at(0).phi()<<")"<<std::endl;
	if(  ToMatch.size()>0 && flagGoodEtaRcJet ) {
		for(unsigned int ito = 0; ito<ToMatch.size() ; ito++) {		// note: if using iteractor for vector, need 'const_iteractor' for const vector
			if(sqrt(pow(RcJAResult.at(0).eta()-ToMatch.at(ito).first,2)+pow(JetAnalyzer::phimod2pi(RcJAResult.at(0).phi()-ToMatch.at(ito).second),2))<R)    {
				trigmatch=true;
				if(mVerbose) std::cout<<"Jet pt "<<RcJAResult.at(0).pt()<<" at (eta,phi)=("<<RcJAResult.at(0).eta()<<","<<RcJAResult.at(0).phi()<<") in R="<<R<<" matched to trigger ("<<ToMatch.at(ito).first<<","<<ToMatch.at(ito).second<<")"<<std::endl;	
				break;
			}
			else if(mVerbose) std::cout<<"not for "<<ToMatch.at(ito).first<<","<<ToMatch.at(ito).second<<std::endl;
		}
	}


	if ( mNeedToMatchTrig && flagGoodEtaRcJet && !trigmatch){
		flagtrigmatch = false;		
		if(mVerbose) std::cout<<"No Match for Jet at (eta,phi)=("<<RcJAResult.at(0).eta()<<","<<RcJAResult.at(0).phi()<<") "<<std::endl;	
	}

	// Neutral fraction cut?
	// ---------------------------------------------------------
	bool flagneutralfrac = true;		// set 1 for not care
	if(flagGoodEtaRcJet) {
		//std::cout<<"Neutral/Total Pt of Jet < 90%"<<std::endl;	
		fastjet::PseudoJet NeutralPart  = fastjet::PseudoJet();
		fastjet::PseudoJet TotalPart  = fastjet::PseudoJet();	
		//fastjet::Selector NoGhosts = !fastjet::SelectorIsPureGhost();
		std::vector<fastjet::PseudoJet> constituents = NoGhosts(RcJAResult.at(0).constituents());
		int charge=-99;
		for(unsigned int jco = 0; jco<constituents.size(); jco++) {
			if ( constituents[jco].is_pure_ghost() ){
			} else {
				charge = (constituents[jco]).user_info<JetAnalysisUserInfo>().GetQuarkCharge();
				//std::cout<<"jet const charge = "<<charge<<std::endl;
				//if(constituents[jco].pt()<0.2) std::cout<<"pt = "<<constituents[jco].pt()<<std::endl;
				//if(fabs(constituents[jco].eta())>1) std::cout<<"eta = "<<constituents[jco].eta()<<std::endl;
			}
		    	TotalPart+=constituents[jco];			
		    	if( charge == 0 ) NeutralPart+=constituents[jco];	 
		}
		if(TotalPart.perp()!=0) {
			Rcj1neutralfrac = fabs(NeutralPart.perp()/TotalPart.perp());		// pt
		}
		else {
			Rcj1neutralfrac = 0;
		}


		if( mNeutralJetFracCut ) {	
			flagneutralfrac = false;		// if we want to cut on neutral energy fraction
			if( Rcj1neutralfrac <= AjParameters::JetNeutralPertMax )  {
				flagneutralfrac = true;		
			//std::cout<<":D Passed neutral cut~~~"<<std::endl;
			}
		}
	}




	if(flagGoodEtaMcJet) {
		fastjet::PseudoJet McNeutralPart  = fastjet::PseudoJet();
		fastjet::PseudoJet McTotalPart  = fastjet::PseudoJet();	
		fastjet::Selector McNoGhosts = !fastjet::SelectorIsPureGhost();
		std::vector<fastjet::PseudoJet> Mcconstituents = McNoGhosts(McJAResult.at(0).constituents());
		int charge=-99;
		for(unsigned int jco = 0; jco<Mcconstituents.size(); jco++) {
			if ( Mcconstituents[jco].is_pure_ghost() ){
			} else {
				charge = (Mcconstituents[jco]).user_info<JetAnalysisUserInfo>().GetQuarkCharge();
				//std::cout<<"Mc jet const charge = "<<charge<<std::endl;
				//if(Mcconstituents[jco].pt()<0.2) std::cout<<"pt = "<<Mcconstituents[jco].pt()<<std::endl;
				//if(fabs(Mcconstituents[jco].eta())>1) std::cout<<"eta = "<<Mcconstituents[jco].eta()<<std::endl;
			}
		    	McTotalPart+=Mcconstituents[jco];			
		    	if( charge == 0 ) McNeutralPart+=Mcconstituents[jco];	 
		}
		if(McTotalPart.perp()!=0) {
			Mcj1neutralfrac = fabs(McNeutralPart.perp()/McTotalPart.perp());		// pt
		}
		else {
			Mcj1neutralfrac = 0;
		}

	}



	// GoodMatch?
	// --------------------
	flagMatch2LeadGood = false;		// Rc jet matched 	(Good for optional: & matched to trigger  & pass netural fraction)
	flagMatch2SubGood = false;		// Rc jet matched 	(Good for optional: & matched to trigger  & pass netural fraction)
	if(flagMatch2Lead && flagtrigmatch==1 && flagneutralfrac==1) flagMatch2LeadGood = true;
	else if(flagMatch2Sub && flagtrigmatch==1 && flagneutralfrac==1) flagMatch2SubGood = true;


	// Dijet
	std::vector<fastjet::PseudoJet> McDiJets;  
	std::vector<fastjet::PseudoJet> RcDiJets;  

	McDiJets = SelectorDijets( dPhiCut ) ( McJAResult );		// McJAResult is non zero only when the hardest jet insdie good eta acceptance
	if(McDiJets.size()==2 && (McnosjetJAResult.size()>1&&fabs(McnosjetJAResult.at(1).eta())<=max_rap))  {		// has dijet and the recoil jet is the 2nd hardest one in whole event
		Mcjaspt = McDiJets.at(1).pt();
		Mcjasphi = McDiJets.at(1).phi();
		Mcjaseta = McDiJets.at(1).eta();
	}


	RcDiJets = SelectorDijets( dPhiCut ) ( RcJAResult );
	if(RcDiJets.size()==2 && (RcnosjetJAResult.size()>1&&fabs(RcnosjetJAResult.at(1).eta())<=max_rap) )  {		// has dijet and the recoil jet is the 2nd hardest one in whole event
		Rcjaspt = RcDiJets.at(1).pt();
		Rcjasphi = RcDiJets.at(1).phi();
		Rcjaseta = RcDiJets.at(1).eta();
	}



	// variables for tree
	// Tree only record Mc jet which are in good eta acceptance 
	if(flagGoodEtaMcJet) {
		Mcj1pt=McJAResult.at(0).pt();
		Mcj1phi=McJAResult.at(0).phi();
		Mcj1eta=McJAResult.at(0).eta();
		
		fastjet::Selector McNoGhosts = !fastjet::SelectorIsPureGhost();
		std::vector<fastjet::PseudoJet> Mcconstituents = McNoGhosts(McJAResult.at(0).constituents());
		Mcj1constntrk = Mcconstituents.size();
		for(unsigned int jco = 0; jco<Mcj1constntrk; jco++) {
			Mcj1constpt[jco] = Mcconstituents[jco].perp();
			Mcj1consteta[jco] = Mcconstituents[jco].eta();
			Mcj1constphi[jco] = Mcconstituents[jco].phi_std();
		}
			
		
		if(McnosjetJAResult.size()>=2 && fabs(McnosjetJAResult.at(1).eta())<=max_rap) {				// Mc subleading jet insdie good eta acceptance
			Mcj2pt=McJAResult.at(1).pt();
			Mcj2phi=McJAResult.at(1).phi();
			Mcj2eta=McJAResult.at(1).eta();
			if(McJAResult.size()>=3) {
				Mcj3pt=McJAResult.at(2).pt();
				Mcj3phi=McJAResult.at(2).phi();
				Mcj3eta=McJAResult.at(2).eta();
				if(McJAResult.size()>=4) {
					Mcj4pt=McJAResult.at(3).pt();
					Mcj4phi=McJAResult.at(3).phi();
					Mcj4eta=McJAResult.at(3).eta();
				}
			}
		}
	}



	// Tree only record Rc jet which are in good eta acceptance + (optionally match to trig and neutral fraction cut)
	if(flagGoodEtaRcJet && flagtrigmatch && flagneutralfrac) {
		Rcj1pt=RcJAResult.at(0).pt();
		Rcj1phi=RcJAResult.at(0).phi();
		Rcj1eta=RcJAResult.at(0).eta();

		fastjet::Selector RcNoGhosts = !fastjet::SelectorIsPureGhost();
		std::vector<fastjet::PseudoJet> Rcconstituents = RcNoGhosts(RcJAResult.at(0).constituents());
		Rcj1constntrk = Rcconstituents.size();
		for(unsigned int jco = 0; jco<Rcj1constntrk; jco++) {
			Rcj1constpt[jco] = Rcconstituents[jco].perp();
			Rcj1consteta[jco] = Rcconstituents[jco].eta();
			Rcj1constphi[jco] = Rcconstituents[jco].phi_std();
		}
			

		if(RcnosjetJAResult.size()>=2 && fabs(RcnosjetJAResult.at(1).eta())<=max_rap) {				// Rc subleading jet insdie good eta acceptance
			Rcj2pt=RcJAResult.at(1).pt();
			Rcj2phi=RcJAResult.at(1).phi();
			Rcj2eta=RcJAResult.at(1).eta();
			if(RcJAResult.size()>=3) {
				Rcj3pt=RcJAResult.at(2).pt();
				Rcj3phi=RcJAResult.at(2).phi();
				Rcj3eta=RcJAResult.at(2).eta();
				if(RcJAResult.size()>=4) {
					Rcj4pt=RcJAResult.at(3).pt();
					Rcj4phi=RcJAResult.at(3).phi();
					Rcj4eta=RcJAResult.at(3).eta();
				}
			}
		}
	}


	Mcrho = McJA_bkgsub->GetBackgroundEstimator()->rho() ;
	Mcrhoerr = McJA_bkgsub->GetBackgroundEstimator()->sigma() ;

	Rcrho = RcJA_bkgsub->GetBackgroundEstimator()->rho() ;
	Rcrhoerr = RcJA_bkgsub->GetBackgroundEstimator()->sigma() ;

	if(flagGoodEtaMcJet) {
		Mcj1area = McJAResult.at(0).area();
		Mcj1area_err = McJAResult.at(0).area_error();
		if(McJAResult.size()>=2) {
			Mcj2area = McJAResult.at(1).area();
			Mcj2area_err = McJAResult.at(1).area_error();
		}
	}

	if(flagGoodEtaRcJet && flagtrigmatch && flagneutralfrac) {
		Rcj1area = RcJAResult.at(0).area();
		Rcj1area_err = RcJAResult.at(0).area_error();
		if(RcJAResult.size()>=2) {
			Rcj2area = RcJAResult.at(1).area();
			Rcj2area_err = RcJAResult.at(1).area_error();
		}
	}



	// Fill histogram for jets
	// -----------------------------
	if(flagGoodEtaMcJet) {
		McLeadJetPt->Fill(Mcj1pt);
		McSubLeadJetPt->Fill(Mcj2pt);
	}

	if(flagGoodEtaRcJet && flagtrigmatch && flagneutralfrac) {
		RcLeadJetPt->Fill(Rcj1pt);
		RcSubLeadJetPt->Fill(Rcj2pt);
	}

	if(flagMatch2LeadGood) {
		McMatchedLeadJetPt->Fill(Mcj1pt);
		McMatchedLeadOrSubJetPt->Fill(Mcj1pt);

		CovMcMatchedLeadJetVsRcJet->Fill(Rcj1pt,Mcj1pt);
		CovMcMatchedLeadOrSubJetVsRcJet->Fill(Rcj1pt,Mcj1pt);
	}
	if(flagMatch2SubGood) {
		McMatchedLeadOrSubJetPt->Fill(Mcj1pt);
		CovMcMatchedLeadOrSubJetVsRcJet->Fill(Rcj1pt,Mcj1pt);
	}


	if(flagGoodEtaMcJet) {
		LoopUnderlying(Mcj1phi, McUconstituents, McLeadAreaNtrk, McSubAreaNtrk, McTranMaxNtrk, McTranMinNtrk, McLeadAreaPt, McSubAreaPt, McTranMaxPt, McTranMinPt, McTrkLeadAreaPt, McTrkLeadAreaPhi, McTrkLeadAreaEta, McTrkSubAreaPt, McTrkSubAreaPhi, McTrkSubAreaEta, McTrkTranMaxPt, McTrkTranMaxPhi, McTrkTranMaxEta, McTrkTranMinPt, McTrkTranMinPhi, McTrkTranMinEta, McTrkLeadAreaId, McTrkSubAreaId, McTrkTranMaxId, McTrkTranMinId);	
		McTranPt = (McTranMaxPt+McTranMinPt)/2.;
	}


	if(flagGoodEtaRcJet && flagtrigmatch && flagneutralfrac) {
		LoopUnderlying(Rcj1phi, RcUconstituents, RcLeadAreaNtrk, RcSubAreaNtrk, RcTranMaxNtrk, RcTranMinNtrk, RcLeadAreaPt, RcSubAreaPt, RcTranMaxPt, RcTranMinPt, RcTrkLeadAreaPt, RcTrkLeadAreaPhi, RcTrkLeadAreaEta, RcTrkSubAreaPt, RcTrkSubAreaPhi, RcTrkSubAreaEta, RcTrkTranMaxPt, RcTrkTranMaxPhi, RcTrkTranMaxEta, RcTrkTranMinPt, RcTrkTranMinPhi, RcTrkTranMinEta, RcTrkLeadAreaMcId, RcTrkSubAreaMcId, RcTrkTranMaxMcId, RcTrkTranMinMcId);	
		RcTranPt = (RcTranMaxPt+RcTranMinPt)/2.;
	}



	ResultTree->Fill();

	nEventPassed++;  


	// We want to hold onto the jetanalyzer objects, so they're created dynamically
	// Need to delete them by hand
	if (McnosjetJA){     delete McnosjetJA;    McnosjetJA=0;  }
	if (McJA_bkgsub){     delete McJA_bkgsub;    McJA_bkgsub=0;  }

	if (RcnosjetJA){     delete RcnosjetJA;    RcnosjetJA=0;  }
	if (RcJA_bkgsub){     delete RcJA_bkgsub;    RcJA_bkgsub=0;  }


	return 1;
}




int underlying_McVsEmbed::LoopUnderlying (float RefPhi, std::vector<fastjet::PseudoJet> Uconstituents, int &LeadAreaNtrk, int &SubAreaNtrk, int &TranMaxNtrk, int &TranMinNtrk, float &LeadAreaPt, float &SubAreaPt, float &TranMaxPt, float &TranMinPt, float *TrkLeadAreaPt, float *TrkLeadAreaPhi, float *TrkLeadAreaEta, float *TrkSubAreaPt, float *TrkSubAreaPhi, float *TrkSubAreaEta, float *TrkTranMaxPt, float *TrkTranMaxPhi, float *TrkTranMaxEta, float *TrkTranMinPt, float *TrkTranMinPhi, float *TrkTranMinEta, int *TrkLeadAreaId, int *TrkSubAreaId, int *TrkTranMaxId, int *TrkTranMinId)  {

	// Calculate underlying event and fill histos 
	std::vector<float> tmp1Pt;
	std::vector<float> tmp1Phi;
	std::vector<float> tmp1Eta;
	std::vector<int> tmp1Id;
	std::vector<float> tmp2Pt;
	std::vector<float> tmp2Phi;
	std::vector<float> tmp2Eta;
	std::vector<int> tmp2Id;

	//==================== Loop over TStarJetVectorContainer for underlying info ====================================
	int ntrklead = 0, ntrksublead = 0, ntrktran = 0, ntrktranmax = 0, ntrktranmin = 0;
	float ptlead = 0, ptsublead = 0, pttran = 0, pttranmax = 0, pttranmin = 0;

	//std::cout<<"# of particles = "<<particles.size()<<std::endl;

	for(std::vector<fastjet::PseudoJet>::const_iterator pj = Uconstituents.begin(); pj!=Uconstituents.end(); pj++) {
		double iphi = pj->phi_std();	// -pi - pi
		double ieta = pj->eta();
		double ipt = pj->perp();	
		int id = pj->user_info<JetAnalysisUserInfo>().GetId();			// MC is MC id; RC is matched MC id

		if(fabs(JetAnalyzer::phimod2pi(iphi-RefPhi))<((180.-mTranPhiSize)/2.)/180.*TMath::Pi()) {		// leading		phimod2pi(phi) gives -pi ->pi
			TrkLeadAreaPt[ntrklead] = ipt;
			TrkLeadAreaPhi[ntrklead] = iphi;
			TrkLeadAreaEta[ntrklead] = ieta;
			TrkLeadAreaId[ntrklead] = id;
			ntrklead++;
			ptlead+=ipt;		// scalar sum
			//std::cout<<"leading"<<std::endl;		
		}
		if(fabs(JetAnalyzer::phimod2pi(iphi-RefPhi))>((180.+mTranPhiSize)/2.)/180.*TMath::Pi()) {	//sub-leading
			TrkSubAreaPt[ntrksublead] = ipt;
			TrkSubAreaPhi[ntrksublead] = iphi;
			TrkSubAreaEta[ntrksublead] = ieta;
			TrkSubAreaId[ntrksublead] = id;
			ntrksublead++;
			ptsublead+=ipt;		// scalar sum
			//std::cout<<"subleading"<<std::endl;		
		}
		if(JetAnalyzer::phimod2pi(iphi-RefPhi)<=((180.+mTranPhiSize)/2.)/180.*TMath::Pi() && JetAnalyzer::phimod2pi(iphi-RefPhi)>=((180.-mTranPhiSize)/2.)/180.*TMath::Pi()) {
			tmp1Pt.push_back(ipt);
			tmp1Phi.push_back(iphi);
			tmp1Eta.push_back(ieta);
			tmp1Id.push_back(id);
			ntrktranmax++;		// will decide which one is max/min later and switch if needed
			pttranmax+=ipt;		// scalar sum
			//std::cout<<"transverse"<<std::endl;		
		}
		if(JetAnalyzer::phimod2pi(iphi-RefPhi)>=-((180.+mTranPhiSize)/2.)/180.*TMath::Pi() && JetAnalyzer::phimod2pi(iphi-RefPhi)<=-((180.-mTranPhiSize)/2.)/180.*TMath::Pi()) {
			tmp2Pt.push_back(ipt);
			tmp2Phi.push_back(iphi);
			tmp2Eta.push_back(ieta);
			tmp2Id.push_back(id);
			ntrktranmin++;
			pttranmin+=ipt;		// scalar sum
			//std::cout<<"transverse"<<std::endl;	
		}
	}

	ntrktran = ntrktranmax+ntrktranmin;
	pttran = pttranmax+pttranmin;

	// Check which is TranMax / TranMin
	if(pttranmax<pttranmin) {
		double tmpp = pttranmin;
		pttranmin = pttranmax;
		pttranmax = tmpp;

		int tmp = ntrktranmin;
		ntrktranmin = ntrktranmax;
		ntrktranmax = tmp;

		std::copy(tmp2Pt.begin(),tmp2Pt.end(),TrkTranMaxPt);
		std::copy(tmp2Phi.begin(),tmp2Phi.end(),TrkTranMaxPhi);
		std::copy(tmp2Eta.begin(),tmp2Eta.end(),TrkTranMaxEta);
		std::copy(tmp2Id.begin(),tmp2Id.end(),TrkTranMaxId);

		std::copy(tmp1Pt.begin(),tmp1Pt.end(),TrkTranMinPt);
		std::copy(tmp1Phi.begin(),tmp1Phi.end(),TrkTranMinPhi);
		std::copy(tmp1Eta.begin(),tmp1Eta.end(),TrkTranMinEta);
		std::copy(tmp1Id.begin(),tmp1Id.end(),TrkTranMinId);

	}
	else {
		std::copy(tmp1Pt.begin(),tmp1Pt.end(),TrkTranMaxPt);
		std::copy(tmp1Phi.begin(),tmp1Phi.end(),TrkTranMaxPhi);
		std::copy(tmp1Eta.begin(),tmp1Eta.end(),TrkTranMaxEta);
		std::copy(tmp1Id.begin(),tmp1Id.end(),TrkTranMaxId);

		std::copy(tmp2Pt.begin(),tmp2Pt.end(),TrkTranMinPt);
		std::copy(tmp2Phi.begin(),tmp2Phi.end(),TrkTranMinPhi);
		std::copy(tmp2Eta.begin(),tmp2Eta.end(),TrkTranMinEta);
		std::copy(tmp2Id.begin(),tmp2Id.end(),TrkTranMinId);

	}


	LeadAreaPt=ptlead;
	SubAreaPt=ptsublead;
	TranMaxPt=pttranmax;
	TranMinPt=pttranmin;

	LeadAreaNtrk=ntrklead;
	SubAreaNtrk=ntrksublead;
	TranMaxNtrk=ntrktranmax;
	TranMinNtrk=ntrktranmin;

}


int underlying_McVsEmbed::Finish () 
{
	
	fout->cd();
	ResultTree->Write();
	fout->Write();	// write histograms also

	std::cout << nEventPassed << " events passed " << std::endl;

	if(DoOutlierMcpTCut) {
		std::cout << nEventOutlierMcpT << " events does not pass due to large Mc pT > "<<mOutlierMcpTCut<<". " << std::endl;
	}
	if(DoOutlierRcpTCut) {
		std::cout << nEventOutlierRcpT << " events does not pass due to large Rc pT > "<<mOutlierRcpTCut<<". " << std::endl;
	}

	std::cout << "Wrote to " << OutFileName << std::endl;

	fout->Close();

	return 1;


}

	

void underlying_McVsEmbed::SetUnderlyingParticleCharge(int val) {

	mUnderlyingParticleCharge = val;

	fastjet::Selector select_Uconst_charge= fastjet::SelectorIdentity();
	//std::cout<<"UnderlyingParticleChargeCode = "<<mUnderlyingParticleCharge<<std::endl;	
	if (mUnderlyingParticleCharge==1) {
		select_Uconst_charge = !SelectorIsNeutralCharge();		// Charged UnderlyingParticle only
	}
	if (mUnderlyingParticleCharge==0) {
		select_Uconst_charge = SelectorIsNeutralCharge();		// Neutral UnderlyingParticle only
	}
	McsUconst     = McsUconst && select_Uconst_charge;			
	RcsUconst     = RcsUconst && select_Uconst_charge;			
}


void underlying_McVsEmbed::SetJetCharge(int val) {

	mJetCharge = val;

	fastjet::Selector select_const_charge= fastjet::SelectorIdentity();
	//std::cout<<"JetChargeCode = "<<mJetCharge<<std::endl;	
	if (mJetCharge==1) {
		select_const_charge = !SelectorIsNeutralCharge();		// Charged Jet only
	}
	if (mJetCharge==0) {
		select_const_charge = SelectorIsNeutralCharge();		// Neutral Jet only
	}
	Mcsconst     = Mcsconst && select_const_charge;			
	Rcsconst     = Rcsconst && select_const_charge;			
}


