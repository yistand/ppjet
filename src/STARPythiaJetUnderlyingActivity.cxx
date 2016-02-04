//====================================================================================================
//
//
//	2016.02.01	Li Yi
//	Read Pythia
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

#include "STARPythia.h"
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
	const char *defaults[] = {"STARPythiaJetUnderlyingActivity","/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/FullJet_TransCharged_pythiaMB.root","pythia","/home/fas/caines/ly247/scratch/pythiadata/pythia_MB_pp200tree_151117.root", "0", "0","FullJet","Charge" };
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
	TString ChainName  = "genevents";
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


	cout<<"Prepare to read STAR pythia data"<<endl;
	STARPythia *mEv = new STARPythia(chain);

	//TString OutFileName = "test.root"; //test 
	//float R = 0.6;	// test
	UnderlyingAna *ula = new UnderlyingAna( R,
			AjParameters::max_track_rap,
			0.2,			// pt min for const.
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

	if(jetchargecode==2) ula->SetNetraulJetFracCut(true);			// whether apply neutral energy fraction in jet cut
	else ula->SetNetraulJetFracCut(false);

	ula->SetJetCharge(jetchargecode);			// Jet charge

	ula->SetUnderlyingParticleCharge(underlyingchargecode);			// underlying event charge: 0 for netural, 1 for charged, 2 for all

	ula->SetDiJetAngle(0);					// Use Dijet angle (1) or Monojet angle (0) 

	// Cycle through events
	// --------------------
	vector<PseudoJet> particles;		

	PseudoJet pj;

	Long64_t nEvents=-1; // -1 for all
	//nEvents=10000;	// test
	int count = 0;

	try{
		while ( mEv->GetEntry(count) ) {
			count++;

			if(count%10000==0) cout<<"event "<<count<<endl;


			// Make particle vector
			// --------------------
			particles.clear();

			for (int ip = 0; ip<mEv->mParticles_ ; ++ip ){

				if(mEv->mParticles_mStatus[ip]!=1) continue;		//	not final particle

				if(!(fabs(mEv->mParticles_mId[ip])==211||fabs(mEv->mParticles_mId[ip])==321 ||fabs(mEv->mParticles_mId[ip])==2212 || mEv->mParticles_mId[ip]==111||mEv->mParticles_mId[ip]==22)) continue;		// if not Pion, Kaon, Proton+-0, gamma


				if(fabs(pow(mEv->mParticles_mPx[ip],2)+pow(mEv->mParticles_mPy[ip],2))<0.2*0.2) continue;		// #ly CHECK!!!!!!!! minimum pT or Et

				pj=PseudoJet (mEv->mParticles_mPx[ip],mEv->mParticles_mPy[ip], mEv->mParticles_mPz[ip], mEv->mParticles_mEnergy[ip] );

				if( fabs(mEv->mParticles_mId[ip])==211||fabs(mEv->mParticles_mId[ip])==321 ||fabs(mEv->mParticles_mId[ip])==2212 ) {
					if(mEv->mParticles_mId[ip]>0) {
						pj.set_user_info ( new JetAnalysisUserInfo( 1 ) );
					}
					else {
						pj.set_user_info ( new JetAnalysisUserInfo( -1 ) );
					}
				}
				if (mEv->mParticles_mId[ip]==111||mEv->mParticles_mId[ip]==22) {
					pj.set_user_info ( new JetAnalysisUserInfo( 0 ) );
				}	

				particles.push_back ( pj );
			}    
			// Run analysis
			// ------------
			//cout<<"analyze and fill"<<endl; 	// test
			std::vector<std::pair<float,float> > ToMatch;
			ToMatch.clear();
			ula->AnalyzeAndFill( particles,
				mEv->mEventNumber,
				mEv->mRunNumber,
				mEv->mNumParticles,
				ToMatch
			 );

		} // while NextEvent
	} catch ( exception& e) {
		cerr << "Caught " << e.what() << endl;
		return -1;
	}
	cout << "##################################################################" << endl;


	// Close up shop
	// -------------
	ula->Finish();


	cout << "Bye." << endl;

	return 0;
}

