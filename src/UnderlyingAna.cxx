/** @file UnderlyingAna.cxx
  @author Kolja Kauder
  @version Revision 0.1
  @brief Class for A<SUB>J</SUB> analysis
  @details Uses JetAnalyzer objects to perform A<SUB>J</SUB> analysis.
  @date Mar 02, 2015
  */

#include "UnderlyingAna.hh"
#include <algorithm>    // std::copy
#include <vector>

using std::endl;


// Standard ctor
UnderlyingAna::UnderlyingAna ( double R,
		//double jet_ptmin, double jet_ptmax,
		//double LeadPtMin, double SubLeadPtMin, 
		double max_track_rap, //double PtConsLo, double PtConsHi,
		double dPhiCut,
		TString name
		)
		: R(R),
	//jet_ptmin(jet_ptmin), jet_ptmax(jet_ptmax),
	//LeadPtMin(LeadPtMin), SubLeadPtMin(SubLeadPtMin),
	max_track_rap (max_track_rap), //PtConsLo (PtConsLo), PtConsHi (PtConsHi),
	dPhiCut (dPhiCut),
	OutFileName (name),
	pJA (0)
	 //pJAhi (0), pJAlo(0), pOtherJAlo(0)
{
	// derived rapidity cuts
	// ---------------------
	max_rap      = max_track_rap-R;
	ghost_maxrap = max_rap + 2.0 * R;

	// Constituent selectors
	// ---------------------
	select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);

	// Provide but turn off charge selector
	// fastjet::Selector select_track_charge= SelectorChargeRange( -3, 3);
	fastjet::Selector select_track_charge= fastjet::SelectorIdentity();
	sconst     = select_track_rap * select_track_charge;

	// Jet candidate selectors
	// -----------------------
	select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
	sjet = select_jet_rap ;

	// Choose a jet and area definition
	// --------------------------------
	jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, R);

	// create an area definition for the clustering
	//----------------------------------------------------------
	// ghosts should go up to the acceptance of the detector or
	// (with infinite acceptance) at least 2R beyond the region
	// where you plan to investigate jets.
	area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, AjParameters::ghost_repeat, AjParameters::ghost_area );
	// // DEBUG
	// // area_spec.set_grid_scatter(0);
	// // area_spec.set_pt_scatter(0);
	area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

	// DEBUG
	// area_spec = fastjet::VoronoiAreaSpec( 0.9 );
	// area_def = fastjet::AreaDefinition( fastjet::VoronoiAreaSpec( 0.9 ) );
	
	// eff
	//tEff = new ktTrackEff();

	std::cout << " ################################################### " << std::endl;
	//std::cout << "Leading jet above " << LeadPtMin << std::endl;
	//std::cout << "Sub-Leading jet above " << SubLeadPtMin << std::endl;
	std::cout << "Clustered with " << jet_def.description() << std::endl;
	std::cout << "Area Spec " << area_spec.description() << std::endl;
	std::cout << "Area Def  " << area_def.description() << std::endl;
	std::cout << " ################################################### " << std::endl;

	// std::cout << slo.applies_jet_by_jet() << std::endl;
	// std::cout << shi.applies_jet_by_jet() << std::endl;
	// std::cout << sjet.applies_jet_by_jet() << std::endl;  
	// std::cout << " ################################################### " << std::endl;

	// set up output root file
	fout = new TFile(OutFileName,"RECREATE");	// just the name from input 
	assert ( fout->IsOpen() );
	std::cout << " ################################################### " << std::endl;
	std::cout << "Writing to: " << fout->GetName() << std::endl;
	std::cout << " ################################################### " << std::endl;


  	nEventProcessed = 0;		// reset counter

	mNeedToMatchTrig = true;
	mNeutralJetFracCut = true;
	mUseDijetAngle = 0;
	mUnderlyingParticleCharge = 1;
	
	ResultTree = NULL;

	LeadJetPt = NULL;
	SubJetPt = NULL;

	LeadJetNtrkvsLeadJetPt = NULL;
	SubJetNtrkvsLeadJetPt = NULL;
	TranMaxNtrkvsLeadJetPt = NULL;
	TranMinNtrkvsLeadJetPt = NULL;
	TranNtrkvsLeadJetPt = NULL;
	
	LeadJetPtSumvsLeadJetPt = NULL;
	SubJetPtSumvsLeadJetPt = NULL;
	TranMaxPtSumvsLeadJetPt = NULL;
	TranMinPtSumvsLeadJetPt = NULL;
	TranPtSumvsLeadJetPt = NULL;
	  
	Spectrum_LeadJetPtvsLeadJetPt = NULL;
	Spectrum_SubJetPtvsLeadJetPt = NULL;
	Spectrum_TranMaxPtvsLeadJetPt = NULL;
	Spectrum_TranMinPtvsLeadJetPt = NULL;
	Spectrum_TranPtvsLeadJetPt = NULL;

}

UnderlyingAna::~UnderlyingAna() {
	
	delete ResultTree;

	delete LeadJetPt;
	delete SubJetPt;

	delete LeadJetNtrkvsLeadJetPt;
	delete SubJetNtrkvsLeadJetPt;
	delete TranMaxNtrkvsLeadJetPt;
	delete TranMinNtrkvsLeadJetPt;
	delete TranNtrkvsLeadJetPt;
	 
	delete LeadJetPtSumvsLeadJetPt;
	delete SubJetPtSumvsLeadJetPt;
	delete TranMaxPtSumvsLeadJetPt;
	delete TranMinPtSumvsLeadJetPt;
	delete TranPtSumvsLeadJetPt;
	   
	delete Spectrum_LeadJetPtvsLeadJetPt;
	delete Spectrum_SubJetPtvsLeadJetPt;
	delete Spectrum_TranMaxPtvsLeadJetPt;
	delete Spectrum_TranMinPtvsLeadJetPt;
	delete Spectrum_TranPtvsLeadJetPt;


}

int UnderlyingAna::Init() {

	// set up TTree
	ResultTree=new TTree("ResultTree","Result Jets");

	ResultTree->Branch("eventid",&eventid, "eventid/I");
	ResultTree->Branch("runid",&runid, "runid/I");
	ResultTree->Branch("refmult",&refmult, "refmult/D");
	ResultTree->Branch("rho",&rho, "rho/F");
	ResultTree->Branch("rhoerr",&rhoerr, "rhoerr/F");

	ResultTree->Branch("j1pt",&j1pt, "j1pt/F");
	ResultTree->Branch("j2pt",&j2pt, "j2pt/F");
	ResultTree->Branch("j1phi",&j1phi, "j1phi/F");
	ResultTree->Branch("j2phi",&j2phi, "j2phi/F");
	ResultTree->Branch("j1eta",&j1eta, "j1eta/F");
	ResultTree->Branch("j2eta",&j2eta, "j2eta/F");

	ResultTree->Branch("j1neutralfrac",&j1neutralfrac, "j1neutralfrac/F");

	ResultTree->Branch("LeadAreaPtSum",&mLeadAreaPt,"LeadAreaPtSum/F");
	ResultTree->Branch("SubLeadAreaPtSum",&mSubAreaPt,"SubLeadAreaPtSum/F");
	ResultTree->Branch("TranPtSum",&mTranPt,"TranPtSum/F");
	ResultTree->Branch("TranMaxPtSum",&mTranMaxPt,"TranMaxPtSum/F");
	ResultTree->Branch("TranMinPtSum",&mTranMinPt,"TranMinPtSum/F");
	ResultTree->Branch("LeadAreaNtrk",&mLeadAreaNtrk,"LeadAreaNtrk/I");
	ResultTree->Branch("SubAreaNtrk",&mSubAreaNtrk,"SubAreaNtrk/I");
	ResultTree->Branch("TranNtrk",&mTranNtrk,"TranNtrk/I");
	ResultTree->Branch("TranMaxNtrk",&mTranMaxNtrk,"TranMaxNtrk/I");
	ResultTree->Branch("TranMinNtrk",&mTranMinNtrk,"TranMinNtrk/I");

	ResultTree->Branch("TrkTranMaxdEdx",TrkTranMaxdEdx,"TrkTranMaxdEdx[TranMaxNtrk]/F");
	ResultTree->Branch("TrkTranMaxTofbeta",TrkTranMaxTofbeta,"TrkTranMaxTofbeta[TranMaxNtrk]/F");
	ResultTree->Branch("TrkTranMaxPt",TrkTranMaxPt,"TrkTranMaxPt[TranMaxNtrk]/F");
	ResultTree->Branch("TrkTranMindEdx",TrkTranMindEdx,"TrkTranMindEdx[TranMinNtrk]/F");
	ResultTree->Branch("TrkTranMinTofbeta",TrkTranMinTofbeta,"TrkTranMinTofbeta[TranMinNtrk]/F");
	ResultTree->Branch("TrkTranMinPt",TrkTranMinPt,"TrkTranMinPt[TranMinNtrk]/F");
	ResultTree->Branch("TrkLeadAreadEdx",TrkLeadAreadEdx,"TrkLeadAreadEdx[LeadAreaNtrk]/F");
	ResultTree->Branch("TrkLeadAreaTofbeta",TrkLeadAreaTofbeta,"TrkLeadAreaTofbeta[LeadAreaNtrk]/F");
	ResultTree->Branch("TrkLeadAreaPt",TrkLeadAreaPt,"TrkLeadAreaPt[LeadAreaNtrk]/F");
	ResultTree->Branch("TrkSubAreadEdx",TrkSubAreadEdx,"TrkSubAreadEdx[SubAreaNtrk]/F");
	ResultTree->Branch("TrkSubAreaTofbeta",TrkSubAreaTofbeta,"TrkSubAreaTofbeta[SubAreaNtrk]/F");
	ResultTree->Branch("TrkSubAreaPt",TrkSubAreaPt,"TrkSubAreaPt[SubAreaNtrk]/F");
 
	//ResultTree->Branch("nRestJ",&nRestJ, "nRestJ/i");
	//ResultTree->Branch("RestArea",RestArea, "RestArea[nRestJ]/f");
	//ResultTree->Branch("RestPt",  RestPt,   "RestPt[nRestJ]/f");

	// set up histograms
	TH1::SetDefaultSumw2(true);
	TH2::SetDefaultSumw2(true);
	TH3::SetDefaultSumw2(true);

	LeadJetPt = new TH1D("LeadingJetPt","Leading Jet Pt",100,0,100);
	SubJetPt = new TH1D("SubLeadingJetPt","SubLeading Jet Pt",100,0,100);

 
	LeadJetNtrkvsLeadJetPt = new TProfile("LeadJetAreaNtrkvsLeadJetPt","Leading Jet Area Ntrk vs Leading Jet Pt",100,0,100);
	SubJetNtrkvsLeadJetPt = new TProfile("SubJetAreaNtrkvsLeadJetPt","SubLeading Jet Area Ntrk vs Leading Jet Pt",100,0,100);
	TranMaxNtrkvsLeadJetPt = new TProfile("TranMaxNtrkvsLeadJetPt","Transverse Max Ntrk vs Leading Jet Pt",100,0,100);
	TranMinNtrkvsLeadJetPt = new TProfile("TranMinNtrkvsLeadJetPt","Transverse Min Ntrk vs Leading Jet Pt",100,0,100);
	TranNtrkvsLeadJetPt = new TProfile("TranNtrkvsLeadJetPt","Transverse Ntrk vs Leading Jet Pt",100,0,100);

	LeadJetPtSumvsLeadJetPt = new TProfile("LeadJetAreaPtSumvsLeadJetPt","Leading Jet Area Sum Pt vs Leading Jet Pt",100,0,100);
	SubJetPtSumvsLeadJetPt = new TProfile("SubJetAreaPtSumvsLeadJetPt","SubLeading Jet Area Sum Pt vs Leading Jet Pt",100,0,100);
	TranMaxPtSumvsLeadJetPt = new TProfile("TranMaxPtSumvsLeadJetPt","Transverse Max Sum Pt vs Leading Jet Pt",100,0,100);
	TranMinPtSumvsLeadJetPt = new TProfile("TranMinPtSumvsLeadJetPt","Transverse Min Sum Pt vs Leading Jet Pt",100,0,100);
	TranPtSumvsLeadJetPt = new TProfile("TranPtSumvsLeadJetPt","Transverse Sum Pt vs Leading Jet Pt",100,0,100);


	Spectrum_LeadJetPtvsLeadJetPt = new TH2D("Spectrum_LeadJetAreaPtvsLeadJetPt","Leading Jet Area Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,30);
	Spectrum_SubJetPtvsLeadJetPt = new TH2D("Spectrum_SubJetAreaPtvsLeadJetPt","SubLeading Jet Area Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,30);
	Spectrum_TranMaxPtvsLeadJetPt = new TH2D("Spectrum_TranMaxPtvsLeadJetPt","Transverse Max Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,30);
	Spectrum_TranMinPtvsLeadJetPt = new TH2D("Spectrum_TranMinPtvsLeadJetPt","Transverse Min Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,30);
	Spectrum_TranPtvsLeadJetPt = new TH2D("Spectrum_TranPtvsLeadJetPt","Transverse Pt Spectrum vs Leading Jet Pt",100,0,100,100,0,30);


	return 1;
}




// Main analysis method
// ====================
//int UnderlyingAna::AnalyzeAndFill ( std::vector<fastjet::PseudoJet>& particles, std::vector<fastjet::PseudoJet>& ToMatch,
int UnderlyingAna::AnalyzeAndFill ( const std::vector<fastjet::PseudoJet>& particles, 
		TStarJetVectorContainer<TStarJetVector>& container,
		TStarJetPicoReader &reader,
	 	//Int_t mEffUn,	
		const std::vector<std::pair<float,float> > &ToMatch
		//Double_t EventClassifier

		//TH1D* LeadJetPt, TH1D* SubJetPt
		//TProfile* LeadJetNtrkvsLeadJetPt, TProfile* SubJetNtrkvsLeadJetPt, TProfile* TranMaxNtrkvsLeadJetPt, TProfile* TranMinNtrkvsLeadJetPt, TProfile* TranNtrkvsLeadJetPt,
		//TProfile* LeadJetPtvsLeadJetPt, TProfile* SubJetPtvsLeadJetPt, TProfile* TranMaxPtvsLeadJetPt, TProfile* TranMinPtvsLeadJetPt, TProfile* TranPtvsLeadJetPt ,
		//TProfile* TranPionNtrkvsLeadJetPt, TProfile* TranPionPtvsLeadJetPt, TProfile* TranProtonNtrkvsLeadJetPt, TProfile* TranProtonPtvsLeadJetPt,
		//TH2D* Spectrum_LeadJetPtvsLeadJetPt, TH2D* Spectrum_SubJetPtvsLeadJetPt, TH2D* Spectrum_TranMaxPtvsLeadJetPt, TH2D* Spectrum_TranMinPtvsLeadJetPt, TH2D* Spectrum_TranPtvsLeadJetPt 
		){

	// We want to hold onto the jetanalyzer objects, so they're created dynamically
	// Need to delete them by hand
	if (pJA){    delete pJA;    pJA=0;  }

	DiJets.clear();

	Jconstituents.clear();

	Has10Gev=false;
	HasDijet=false;

	eventid = reader.GetEvent()->GetHeader()->GetEventId();	
	runid = reader.GetEvent()->GetHeader()->GetRunId();
	refmult = reader.GetEvent()->GetHeader()->GetGReferenceMultiplicity();

	j1pt=0, j2pt=0;
	j1phi=-999, j2phi=-999;
	j1eta=-999, j2eta=-999;
	j1area=0,j2area=0;
	rho=0, rhoerr=0;

	j1neutralfrac=0;

	mLeadAreaPt=0;
	mSubAreaPt=0;
	mTranMaxPt=0;
	mTranMinPt=0;
	mTranPt=0;

	mLeadAreaNtrk=0;
	mSubAreaNtrk=0;
	mTranMaxNtrk=0;
	mTranMinNtrk=0;
	mTranNtrk=0;

	for(int i = 0; i<MAXARRAYLENGTH; i++) {
		TrkTranMaxdEdx[i]=0;
		TrkTranMaxTofbeta[i]=0;
		TrkTranMaxPt[i]=0;
		TrkTranMindEdx[i]=0;
		TrkTranMinTofbeta[i]=0;
		TrkTranMinPt[i]=0;
		TrkLeadAreadEdx[i]=0;
		TrkLeadAreaTofbeta[i]=0;
		TrkLeadAreaPt[i]=0;
		TrkSubAreadEdx[i]=0;
		TrkSubAreaTofbeta[i]=0;
		TrkSubAreaPt[i]=0;
	}
	
	// Set eff
	//tEff->SetSysUncertainty(mEffUn);


	// Select particles to perform analysis on
	// ---------------------------------------
	//std::cout<<"selector for consti"<<std::endl; 		// test
	Jconstituents = sconst( particles );

	// Background selector
	// -------------------
	// It is unclear to me why, but it leads to segfaults if only a once-initialized member :-/
	//std::cout<<"selector for background"<<std::endl;		// test
	fastjet::Selector selector_bkgd = fastjet::SelectorAbsRapMax( max_rap ) * (!fastjet::SelectorNHardest(2));	//selector for fastjet::JetMedianBackgroundEstimator: the selection criterion is typically a geometrical one (e.g. all jets with |y|<2) sometimes supplemented with some kinematical restriction (e.g. exclusion of the two hardest jets). It is passed to the class through a Selector.
	// selector_bkgd=fastjet::SelectorAbsRapMax( max_rap );

	// find jets
	// -----------------------------
	pJA = new JetAnalyzer( Jconstituents, jet_def , area_def, selector_bkgd); // with background subtraction
	JetAnalyzer& JA = *pJA;
	//std::cout<<"selector on jet"<<std::endl; 			// test
	fastjet::Subtractor* BackgroundSubtractor =  JA.GetBackgroundSubtractor();
	JAResult = fastjet::sorted_by_pt( (*BackgroundSubtractor) ( JA.inclusive_jets() ) );

	if ( JAResult.size() < 1 )                 {     return 0; }
	if ( JAResult.at(0).pt() > 10 )            { Has10Gev=true; }

	//if ( JAResult.size() < 2 )                 {     return 0; }

	// back to back? Answer this question with a selector
	// ---------------------------------------------------
	DiJets = SelectorDijets( dPhiCut ) ( JAResult );
	if(DiJets.size()==2) HasDijet=true;
	//if ( DiJets.size() == 0 ) {
	// std::cout << " NO dijet found" << std::endl;
	// return 0;
	//}
	//assert ( DiJets.size() == 2 && "SelectorDijets returned impossible number of Dijets." );  

	// ---------------------------------------------------------
	// Do any jets match to the one fired the trigger?
	if ( mNeedToMatchTrig &&  ToMatch.size()>0 ){
		int flagtrigmatch = 0;
		for(unsigned int ito = 0; ito<ToMatch.size() ; ito++) {		// note: if using iteractor for vector, need 'const_iteractor' for const vector
			//if(IsMatched(JAResult.at(0), ToMatch.at(i), R)) {
			if(sqrt(pow(JAResult.at(0).eta()-ToMatch.at(ito).first,2)+pow(JAResult.at(0).phi()-ToMatch.at(ito).second,2))<R)    {
				flagtrigmatch=1;
				//std::cout<<"Jet at (eta,phi)=("<<JAResult.at(0).eta()<<","<<JAResult.at(0).phi()<<") in R="<<R<<" with ("<<ToMatch.at(ito).first<<","<<ToMatch.at(ito).second<<")"<<std::endl;	
				break;
			}
		}
		if(flagtrigmatch==0) return 0;
	}
	//if ( !IsMatched( DiJets, *ToMatch, R ) ) return 0;

	// ---------------------------------------------------------
	// Neutral/Total Pt of Jet < 90% cut
	//std::cout<<"mNeutralJetFracCut = "<<mNeutralJetFracCut<<std::endl;	//test
	//std::cout<<"mNeedToMatchTrig = "<<mNeedToMatchTrig<<std::endl;	//test
	//std::cout<<"mUseDijetAngle = "<<mUseDijetAngle<<std::endl;	//test
	//std::cout<<"mUnderlyingParticleCharge = "<<mUnderlyingParticleCharge<<std::endl;	//test

	if( mNeutralJetFracCut ) {
		//std::cout<<"Neutral/Total Pt of Jet < 90%"<<std::endl;	//test
		fastjet::PseudoJet NeutralPart  = fastjet::PseudoJet();
		fastjet::PseudoJet TotalPart  = fastjet::PseudoJet();		//test
		fastjet::Selector NoGhosts = !fastjet::SelectorIsPureGhost();
		std::vector<fastjet::PseudoJet> constituents = NoGhosts(JAResult.at(0).constituents());
		int charge=-99;
		for(unsigned int jco = 0; jco<constituents.size(); jco++) {
			if ( constituents[jco].is_pure_ghost() ){
				std::cout << "is a ghost" << endl;
			} else {
				charge = (constituents[jco]).user_info<JetAnalysisUserInfo>().GetQuarkCharge();
			}
		    	TotalPart+=constituents[jco];			// test
		    	if( charge == 0 ) NeutralPart+=constituents[jco];	//test 
		}
		//#ly NOTE: not sure why JAResult.at(0) is not equal to sum of its constituents: because particle sum by weight for jet
		//double frac = NeutralPart.perp2()/JAResult.at(0).perp2();
		j1neutralfrac = fabs(NeutralPart.perp()/TotalPart.perp());
		if( j1neutralfrac > AjParameters::JetNeutralPertMax )  {
		//if( (NeutralPart.perp2()/TotalPart.perp2()) > AjParameters::JetNeutralPertMax )  {
			//std::cout<<"Neutral Jet .. Pass: "<<frac<<" > "<<AjParameters::JetNeutralPertMax <<std::endl;
			//test #ly return 0;
		}
		//else std::cout<<":D Passed neutral cut~~~"<<std::endl;	// test
	}


	// ---------------------------------------------------------
	// Calculate underlying event and fill histos 
	double Ptleadingjet = JAResult.at(0).pt();

	LeadJetPt->Fill(JAResult.at(0).pt());
	j1 = MakeTLorentzVector(JAResult.at(0));
	j1pt = j1.Pt();	
	j1phi = j1.Phi();
	j1eta = j1.Eta();
	j1area = JAResult.at(0).area();
	if(HasDijet) { 
		SubJetPt->Fill(JAResult.at(1).pt());
		j2 = MakeTLorentzVector(JAResult.at(1));
		j2pt = j2.Pt();	
		j2phi = j2.Phi();
		j2eta = j2.Eta();
		j2area = JAResult.at(1).area();
	}
	else {
		SubJetPt->Fill(0);
		j2 = TLorentzVector(0,0,0,0);
	}
	rho = pJA->GetBackgroundEstimator()->rho() ;
	rhoerr = pJA->GetBackgroundEstimator()->sigma() ;

	double DiJetPhi;
	double DiJetEta; 

	if(HasDijet && mUseDijetAngle) {
		DiJetPhi = JetAnalyzer::phimod2pi( (JAResult.at(0).phi() + (TMath::Pi() - JAResult.at(1).phi()))/2. ); // -pi -> pi
		DiJetEta = (JAResult.at(0).eta() + JAResult.at(1).eta())/2.; 
	}
	else {
		DiJetPhi = JetAnalyzer::phimod2pi(JAResult.at(0).phi());	// -pi -> pi
		DiJetEta = JAResult.at(0).eta();
	}
	//std::cout<<"DiJetPhi = "<<JAResult.at(0).phi()<<"-"<<JAResult.at(1).phi()<<"/2="<<DiJetPhi<<std::endl;	



	TH1D *htmp1 = (TH1D*)Spectrum_TranMaxPtvsLeadJetPt->Clone("htmp1");
	htmp1->Reset();
	TH1D *htmp2 = (TH1D*)Spectrum_TranMinPtvsLeadJetPt->Clone("htmp2");
	htmp2->Reset();
	
	std::vector<float> tmp1dEdx;		
	std::vector<float> tmp1Tofbeta;
	std::vector<float> tmp1Pt;
	std::vector<float> tmp2dEdx;
	std::vector<float> tmp2Tofbeta;
	std::vector<float> tmp2Pt;

	//==================== Loop over TStarJetVectorContainer for underlying info ====================================
	TStarJetVector* sv; 		// derived from TLorentzVector
	int ntrklead = 0, ntrksublead = 0, ntrktran = 0, ntrktranmax = 0, ntrktranmin = 0;
	float ptlead = 0, ptsublead = 0, pttran = 0, pttranmax = 0, pttranmin = 0;

	for(int i = 0; i<container.GetEntries(); i++) {
		sv = container.Get(i);
		double iphi = sv->phi_std();	// -pi - pi
		double ieta = sv->eta();
		double ipt = sv->perp();	
		if(mUnderlyingParticleCharge==0 && sv->GetCharge()!=0) continue;		// neutral particle only
		if(mUnderlyingParticleCharge==1 && sv->GetCharge()==0) continue;		// charged particle only
		if(fabs(ieta)>max_track_rap) continue;	// rapidity cut

		float idedx = sv->GetFeatureD(TStarJetVector::_DEDX);
		float itofbeta = sv->GetFeatureD(TStarJetVector::_TOFBETA);
		
		//std::cout<<DiJetPhi<<"-"<<iphi<<" = "<<JetAnalyzer::phimod2pi(iphi-DiJetPhi)<<std::endl;		// test
		if(fabs(JetAnalyzer::phimod2pi(iphi-DiJetPhi))<60./180.*TMath::Pi()) {		// leading		phimod2pi(phi) gives -pi ->pi
			TrkLeadAreadEdx[ntrklead] = idedx;
			TrkLeadAreaTofbeta[ntrklead] = itofbeta;
			TrkLeadAreaPt[ntrklead] = ipt;
			ntrklead++;
			ptlead+=ipt;		// scalar sum
			Spectrum_LeadJetPtvsLeadJetPt->Fill(Ptleadingjet,ipt);	
			//std::cout<<"leading"<<std::endl;	// test
		}
		if(fabs(JetAnalyzer::phimod2pi(iphi-DiJetPhi))>120/180.*TMath::Pi()) {	//sub-leading
			TrkSubAreadEdx[ntrksublead] = idedx;
			TrkSubAreaTofbeta[ntrksublead] = itofbeta;
			TrkSubAreaPt[ntrksublead] = ipt;
			ntrksublead++;
			ptsublead+=ipt;		// scalar sum
			Spectrum_SubJetPtvsLeadJetPt->Fill(Ptleadingjet,ipt);	
			//std::cout<<"subleading"<<std::endl;	//test
		}
		if(JetAnalyzer::phimod2pi(iphi-DiJetPhi)<=120/180.*TMath::Pi() && JetAnalyzer::phimod2pi(iphi-DiJetPhi)>=60/180.*TMath::Pi()) {
			tmp1dEdx.push_back(idedx);
			tmp1Tofbeta.push_back(itofbeta);
			tmp1Pt.push_back(ipt);
			ntrktranmax++;		// will decide which one is max/min later and switch if needed
			pttranmax+=ipt;		// scalar sum
			Spectrum_TranPtvsLeadJetPt->Fill(Ptleadingjet,ipt);	
			htmp1->Fill(Ptleadingjet,ipt);
			//std::cout<<"transverse"<<std::endl;	//test
		}
		if(JetAnalyzer::phimod2pi(iphi-DiJetPhi)>=-120/180.*TMath::Pi() && JetAnalyzer::phimod2pi(iphi-DiJetPhi)<=-60/180.*TMath::Pi()) {
			tmp2dEdx.push_back(idedx);
			tmp2Tofbeta.push_back(itofbeta);
			tmp2Pt.push_back(ipt);
			ntrktranmin++;
			pttranmin+=ipt;		// scalar sum
			Spectrum_TranPtvsLeadJetPt->Fill(Ptleadingjet,ipt);	
			htmp2->Fill(Ptleadingjet,ipt);
			//std::cout<<"transverse"<<std::endl;	//test
		}
	}

	ntrktran = ntrktranmax+ntrktranmin;
	pttran = pttranmax+pttranmin;
	//if(ntrktran>0)  pttran/=ntrktran;

	// Check which is TranMax / TranMin
	if(pttranmax<pttranmin) {
		double tmpp = pttranmin;
		pttranmin = pttranmax;
		pttranmax = tmpp;

		int tmp = ntrktranmin;
		ntrktranmin = ntrktranmax;
		ntrktranmax = tmp;

		std::copy(tmp2dEdx.begin(),tmp2dEdx.end(),TrkTranMaxdEdx);	// #ly template <class InputIterator, class OutputIterator>  OutputIterator copy (InputIterator first, InputIterator last, OutputIterator result);
		std::copy(tmp2Tofbeta.begin(),tmp2Tofbeta.end(),TrkTranMaxTofbeta);
		std::copy(tmp2Pt.begin(),tmp2Pt.end(),TrkTranMaxPt);

		std::copy(tmp1dEdx.begin(),tmp1dEdx.end(),TrkTranMindEdx);	// #ly template <class InputIterator, class OutputIterator>  OutputIterator copy (InputIterator first, InputIterator last, OutputIterator result);
		std::copy(tmp1Tofbeta.begin(),tmp1Tofbeta.end(),TrkTranMinTofbeta);
		std::copy(tmp1Pt.begin(),tmp1Pt.end(),TrkTranMinPt);

		Spectrum_TranMaxPtvsLeadJetPt->Add(htmp2);
		Spectrum_TranMinPtvsLeadJetPt->Add(htmp1);
	}
	else {
		std::copy(tmp1dEdx.begin(),tmp1dEdx.end(),TrkTranMaxdEdx);	// #ly template <class InputIterator, class OutputIterator>  OutputIterator copy (InputIterator first, InputIterator last, OutputIterator result);
		std::copy(tmp1Tofbeta.begin(),tmp1Tofbeta.end(),TrkTranMaxTofbeta);
		std::copy(tmp1Pt.begin(),tmp1Pt.end(),TrkTranMaxPt);

		std::copy(tmp2dEdx.begin(),tmp2dEdx.end(),TrkTranMindEdx);	// #ly template <class InputIterator, class OutputIterator>  OutputIterator copy (InputIterator first, InputIterator last, OutputIterator result);
		std::copy(tmp2Tofbeta.begin(),tmp2Tofbeta.end(),TrkTranMinTofbeta);
		std::copy(tmp2Pt.begin(),tmp2Pt.end(),TrkTranMinPt);

		Spectrum_TranMaxPtvsLeadJetPt->Add(htmp1);
		Spectrum_TranMinPtvsLeadJetPt->Add(htmp2);
	}
	delete htmp1;		// #ly when 'new' is used for the variable, 'delete' is needed to free the memory 
	delete htmp2;

	//if(ntrktranmax>0) pttranmax/=ntrktranmax;
	//if(ntrktranmin>0) pttranmin/=ntrktranmin;

	//if(ntrklead>0) ptlead/=ntrklead;
	//if(ntrksublead>0) ptsublead/=ntrksublead;

	mLeadAreaPt=ptlead;
	mSubAreaPt=ptsublead;
	mTranMaxPt=pttranmax;
	mTranMinPt=pttranmin;
	mTranPt=pttran/2.;		// take the average of these two trans

	mLeadAreaNtrk=ntrklead;
	mSubAreaNtrk=ntrksublead;
	mTranMaxNtrk=ntrktranmax;
	mTranMinNtrk=ntrktranmin;
	mTranNtrk=ntrktran/2.;		// take the average of these two trans

	LeadJetNtrkvsLeadJetPt->Fill(Ptleadingjet,mLeadAreaNtrk);
	SubJetNtrkvsLeadJetPt->Fill(Ptleadingjet,mSubAreaNtrk);
	TranMaxNtrkvsLeadJetPt->Fill(Ptleadingjet,mTranMaxNtrk);
	TranMinNtrkvsLeadJetPt->Fill(Ptleadingjet,mTranMinNtrk);
	TranNtrkvsLeadJetPt->Fill(Ptleadingjet,mTranNtrk);

	LeadJetPtSumvsLeadJetPt->Fill(Ptleadingjet,mLeadAreaPt);
	SubJetPtSumvsLeadJetPt->Fill(Ptleadingjet,mSubAreaPt);
	TranMaxPtSumvsLeadJetPt->Fill(Ptleadingjet,mTranMaxPt);
	TranMinPtSumvsLeadJetPt->Fill(Ptleadingjet,mTranMinPt);
	TranPtSumvsLeadJetPt->Fill(Ptleadingjet,mTranPt);

	ResultTree->Fill();

	//std::cout<<"pJA = "<<pJA->GetBackgroundEstimator()<<std::endl;		// test

	nEventProcessed++;


	return 1;
}

int UnderlyingAna::Finish ( ) {

	fout->cd();
	ResultTree->Write();
	fout->Write();	// write histograms also

	std::cout << nEventProcessed << " events processed " << std::endl;

	std::cout << "Wrote to " << OutFileName << std::endl;

	fout->Close();

	return 1;

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

	evCuts->SetPVRankingCut ( 0 );		// Vertex ranking

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

	std::cout << "Using these tower cuts:" << std::endl;
	std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
	std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;

	// V0s: Turn off
	reader.SetProcessV0s(false);

	return reader;

}



