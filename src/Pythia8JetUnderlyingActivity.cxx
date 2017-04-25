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
#include <TVector2.h>

#include <utility>	// std::pair, std::make_pair
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
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


int stoi(string s) {
	int x = 0; 
	stringstream convert(s);
	convert >> x; 
	return x;
}

int main ( int argc, const char** argv ) {

	// Set up some convenient default
	// ------------------------------
	const char *defaults[] = {"PYTHIAJetUnderlyingActivity","/home/hep/caines/ly247/scratch60/pythia200/FullJet_TransCharged_pythiaSTARdefault.root","pythia","/home/hep/caines/ly247/scratch60/pythia2/pythia_STAR200default.root", "0", "0", "0", "1000000" };
	// {Code name, to be discard but needed since argv will use command name as the [0], output file name, triggername, intput file list, for variable IntTowScale to scale tower as systematics study, which effiencey file to use, start event, total event numbers }
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
	TString ChainName  = "tree";
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

	cout<<"Prepare to read pythia data"<<endl;

	// input
	int eventid;
	double weight;
	int nparticles;
	const int MaxSize=1000;
	int id[MaxSize], status[MaxSize], charge[MaxSize];
	float px[MaxSize], py[MaxSize], pz[MaxSize], energy[MaxSize];

	chain->SetBranchAddress("eventid",&eventid);
	chain->SetBranchAddress("eventweight",&weight);
	chain->SetBranchAddress("npart",&nparticles);
	chain->SetBranchAddress("id",id);
	chain->SetBranchAddress("status",status);
	chain->SetBranchAddress("charge",charge);
	chain->SetBranchAddress("px",px);
	chain->SetBranchAddress("py",py);
	chain->SetBranchAddress("pz",pz);
	chain->SetBranchAddress("energy",energy);

	// output
	vector<PseudoJet> Pparticles;		

	// Initialize analysis class
	// -------------------------
	//cout<<"initialize analysis"<<endl;		


	Long64_t nEvents=-1; // -1 for all
	nEvents = stoi(arguments.at(6));
	//nEvents=1000;	// test
	int startc = stoi(arguments.at(5));
	int count = 0; 
	if(startc>=0) {
		count = startc;
	}
	else {
		startc = 0;
	}
	TString NoFileTag = "";
	if(startc!=0&&nEvents>0)  { 
		NoFileTag+=Form("_from%d",startc);
		NoFileTag+=Form("_to%d",startc+nEvents);
	}
	else if(startc!=0){
		NoFileTag+=Form("_from%d",startc);
	}

	OutFileName.ReplaceAll(".root",NoFileTag+TString(".root"));

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

	ula->SetToMatchJetTrigger(false);

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
	PseudoJet pj;

	std::vector<EtaPhiPair> TrigLoc2Match;		// will be empty for our pythia8, no GEANT


	cout<<"Total "<<chain->GetEntries()<<" events"<<endl;
	if(nEvents>=0) {
		cout<<"Read from "<<count<<" to "<<startc+nEvents<<" events"<<endl;
	}
	else {
		cout<<"Read from "<<count<<endl;
	}
	
	try{
		while ( chain->GetEntry(count) ) {
			if(count%10000==0) cout<<"event "<<count<<endl;
			count++;
			if(nEvents>0 && count>startc+nEvents) break;


			TrigLoc2Match.clear();
			Pparticles.clear();

			// Load event particles
			// ----------
			for (int ip = 0; ip<nparticles ; ++ip ){
				float pt = sqrt(px[ip]*px[ip]+py[ip]*py[ip]);
				if(fabs(pt)<0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et. --> NOT in use anymore, moved cuts to UnderlyingAna class min_const_pt for all particles (tracks and towers)

				pj=PseudoJet( px[ip], py[ip], pz[ip], energy[ip] );
				pj.set_user_info ( new JetAnalysisUserInfo( charge[ip], 0, 0 ) );
				//cout<<"px = "<<pj.px()<<" py = "<<pj.py()<<" pz = "<<pj.pz()<<" eta = "<<pj.eta()<<" phi = "<<pj.phi()<<" pt = "<<pj.perp()<<endl;

				Pparticles.push_back ( pj );

				//}	      
			}    
			// Run analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	

			ula->AnalyzeAndFill( Pparticles, 
					eventid,
					0,
					nparticles,
					0,
					TrigLoc2Match,
					weight
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

