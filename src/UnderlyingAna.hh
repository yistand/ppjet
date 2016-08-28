//====================================================================================================
//
//	2015.09.22	Li Yi
//	modified from Kolja Kauder's Aj analysis to study underlying event activity dependence on jet
//	energy in pp 200 GeV
//
//====================================================================================================

/** 
    @author Kolja Kauder
    @version Revision 0.1
    @brief Class for A<SUB>J</SUB> analysis
    @details Uses JetAnalyzer objects to perform A<SUB>J</SUB> analysis.
    @date Mar 02, 2015
*/

#ifndef __UNDERLYINGANA_HH
#define __UNDERLYINGANA_HH

#include "AjParameters.hh"
#include "JetAnalyzer.hh"
//#include "ktTrackEff.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"


//// Not needed for analysis per se
//#include "TStarJetPicoReader.h"
//#include "TStarJetPicoEvent.h"
//#include "TStarJetPicoEventHeader.h"
//#include "TStarJetPicoEventCuts.h"
//
//#include "TStarJetPicoPrimaryTrack.h"
//#include "TStarJetPicoTrackCuts.h"
//#include "TStarJetPicoTowerCuts.h"
//
//#include "TStarJetVectorContainer.h"
//#include "TStarJetVector.h"
//#include "TStarJetPicoTriggerInfo.h"
//#include "TStarJetPicoUtils.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <string>

#define MAXARRAYLENGTH  5000

/**
   The main class
 */
class UnderlyingAna {
private :

  // count nEventProcessed;
  int nEventProcessed;

  // Match to trigger or not
  bool mNeedToMatchTrig;

  // Neutral/Total Pt of Jet fraction cut
  bool mNeutralJetFracCut;

  // Use Dijet angle (1) or Monojet angle (0)
  int mUseDijetAngle;

  // Underlying event particle: charged(1) or neutral(0) or all(2, or whatever val not equal 1 or 0)
  int mUnderlyingParticleCharge;

  // Jet: charged(1) or neutral(0) or all(2, or whatever val not equal 1 or 0)
  int mJetCharge;

  // Underlying event phi size: default 60 degree
  float mTranPhiSize;

  // For output 
  TFile *fout;
  TString OutFileName;
  // Tree and its variables
  TTree* ResultTree;

  // event header
  int eventid;
  int runid;
  double refmult; 
  double vz;


  // Jet-finding from fastjet
  TLorentzVector j1, jas, j2; 	// leading & away-side jet & subleading jet	(away-side jet passed dijet phi selection; subleading is the 2nd highest pt jet)
  float j1pt, jaspt, j2pt;
  float j1phi, jasphi, j2phi;
  float j1eta, jaseta, j2eta;
  float j1area,jasarea, j2area;
  float j1area_err, jasarea_err, j2area_err;
  float j1neutralfrac;
  float j1r1pt;			// leading jet with R = 1

  float rho, rhoerr;


  TLorentzVector j3, j4; 	// whether there are jets in transverse region
  float j3pt, j4pt;
  float j3phi, j4phi;
  float j3eta, j4eta;

  // underlying event info
  float mLeadAreaPt;
  float mSubAreaPt;
  float mTranMaxPt;
  float mTranMinPt;
  float mTranPt;

  int mLeadAreaNtrk;
  int mSubAreaNtrk;
  int mTranMaxNtrk;
  int mTranMinNtrk;
  int mTranNtrk;

  float TrkTranMaxdEdx[MAXARRAYLENGTH];		// for track in tranmax
  float TrkTranMaxTofbeta[MAXARRAYLENGTH];  	// for track in tranmax
  float TrkTranMaxPt[MAXARRAYLENGTH];		// for track in tranmax
  float TrkTranMaxPhi[MAXARRAYLENGTH];		// for track in tranmax
  float TrkTranMaxEta[MAXARRAYLENGTH];		// for track in tranmax
  float TrkTranMindEdx[MAXARRAYLENGTH];		// for track in tranmin
  float TrkTranMinTofbeta[MAXARRAYLENGTH];  	// for track in tranmin
  float TrkTranMinPt[MAXARRAYLENGTH];		// for track in tranmin
  float TrkTranMinPhi[MAXARRAYLENGTH];		// for track in tranmax
  float TrkTranMinEta[MAXARRAYLENGTH];		// for track in tranmax
  float TrkLeadAreadEdx[MAXARRAYLENGTH];		// for track in lead 
  float TrkLeadAreaTofbeta[MAXARRAYLENGTH];  		// for track in lead
  float TrkLeadAreaPt[MAXARRAYLENGTH];			// for track in lead
  float TrkLeadAreaPhi[MAXARRAYLENGTH];		// for track in lead
  float TrkLeadAreaEta[MAXARRAYLENGTH];		// for track in lead
  float TrkSubAreadEdx[MAXARRAYLENGTH];		// for track in sublead
  float TrkSubAreaTofbeta[MAXARRAYLENGTH];  	// for track in sublead
  float TrkSubAreaPt[MAXARRAYLENGTH];		// for track in sublead
  float TrkSubAreaPhi[MAXARRAYLENGTH];		// for track in sublead
  float TrkSubAreaEta[MAXARRAYLENGTH];		// for track in sublead
 
  // Histograms
  TH1D* LeadJetPt;
  TH1D* SubJetPt;
  
  TProfile* LeadJetNtrkvsLeadJetPt;
  TProfile* SubJetNtrkvsLeadJetPt;
  TProfile* TranMaxNtrkvsLeadJetPt;
  TProfile* TranMinNtrkvsLeadJetPt;
  TProfile* TranNtrkvsLeadJetPt;
  
  TProfile* LeadJetPtSumvsLeadJetPt;
  TProfile* SubJetPtSumvsLeadJetPt;
  TProfile* TranMaxPtSumvsLeadJetPt;
  TProfile* TranMinPtSumvsLeadJetPt;
  TProfile* TranPtSumvsLeadJetPt;
  
  TH2D* Spectrum_LeadJetPtvsLeadJetPt;
  TH2D* Spectrum_SubJetPtvsLeadJetPt;
  TH2D* Spectrum_TranMaxPtvsLeadJetPt;
  TH2D* Spectrum_TranMinPtvsLeadJetPt;
  TH2D* Spectrum_TranPtvsLeadJetPt;


  // efficiency  
  //ktTrackEff* tEff;


  // for jet finding
  double R;              ///< Resolution parameter ("jet radius")
  double max_rap;        ///< jet rapidity acceptance
  double ghost_maxrap;   ///< for ghosted area, should be >= max_rap + 2*R

  double max_const_rap;  ///< constituent rapidity cut
  double min_const_pt;   ///< constituent pt cut

  double dPhiCut;        ///< opening angle for dijet requirement. Accept only  |&phi;1 - &phi;2 - &pi;| < &Delta;&phi;.

  fastjet::JetDefinition jet_def;       ///< jet definition
  fastjet::JetDefinition other_jet_def; ///< jet definition with a different radius

  fastjet::Selector select_const_rap;   ///< constituent rapidity selector
  fastjet::Selector select_const_ptmin;   ///< constituent p<SUB>T</SUB>  selector

  fastjet::Selector sconst;                ///< compound selector for constituents for jet finding

  fastjet::Selector sUconst;                ///< compound selector for constituents for underlying event 

// Relevant jet candidates
  fastjet::Selector select_jet_rap;        ///< jet rapidity selector
  fastjet::Selector select_jet_pt_min;     ///< jet p<SUB>T</SUB> selector
  fastjet::Selector select_jet_pt_max;     ///< jet p<SUB>T</SUB> selector
  fastjet::Selector sjet;                  ///< compound jet selector

// Selector pseudojet not ghost
  fastjet::Selector NoGhosts;		///< non-ghost jet  

  fastjet::GhostedAreaSpec area_spec;      ///< ghosted area specification
  fastjet::AreaDefinition area_def;        ///< jet area definition

  JetAnalyzer* pJA;                      ///< JetAnalyzer object
  JetAnalyzer* pJA_bkgsub;                      ///< JetAnalyzer object with background subtraction
  
  std::vector<fastjet::PseudoJet> Jconstituents;     ///< constituents for jet finding

  std::vector<fastjet::PseudoJet> JAResult;  ///< Unaltered clustering result 
  
  std::vector<fastjet::PseudoJet> JAResult_bkgsub;  ///< Unaltered clustering result with background subtraction
  

  std::vector<fastjet::PseudoJet> DiJets;    ///< Dijet result 

  
  std::vector<fastjet::PseudoJet> Uconstituents;     ///< particles for underlying events

public:

  /** Standard constructor. Set up analysis parameters.
      \param R: jet resolution parameter (radius)
      \param jet_ptmin: minimum jet p<SUB>T</SUB>
      \param jet_ptmax: maximum jet p<SUB>T</SUB>
      \param LeadPtMin: leading jet minimum p<SUB>T</SUB>
      \param SubLeadPtMin: subleading jet minimum p<SUB>T</SUB>
      \param max_const_rap: constituent rapidity cut
      \param PtConsLo: constituent minimum p<SUB>T</SUB>
      \param PtConsHi: constituent maximum p<SUB>T</SUB>
      \param dPhiCut: opening angle for dijet requirement. Accept only  |&phi;1 - &phi;2 - &pi;| < dPhiCut.
   */
  UnderlyingAna ( double R = 0.4,
		//double jet_ptmin = 10.0, double jet_ptmax = 100.0,
		//double LeadPtMin = 20.0, double SubLeadPtMin = 10, 
		double max_const_rap = 1.0, //double PtConsLo=0.2, double PtConsHi=2.0,
		double min_const_pt = 0.2,
		double dPhiCut = 0.4,
		TString name = "underlyingoutput.root",
		std::string jetalgo = "antikt"	
	        );

  ~UnderlyingAna();

  int Init ();
  
  /** Main analysis routine.
      \param particles: Current event
      \param ToMatch: Optionally enforce matching of at least one of the dijets to a trigger
      \param EventClassifier: Used to separate between events, e.g. by RefmultBin
      \param UnmatchedAJ_hi: Dijet imbalance &Delta;p<SUB>T</SUB> / &Sigma;p<SUB>T</SUB> for all jets with high p<SUB>T</SUB> constituents.
      \param AJ_hi: Dijet imbalance &Delta;p<SUB>T</SUB> / &Sigma;p<SUB>T</SUB> for matched jets with high p<SUB>T</SUB> constituents.
      \param AJ_lo: Dijet imbalance &Delta;p<SUB>T</SUB> / &Sigma;p<SUB>T</SUB> for matched jets with low p<SUB>T</SUB> constituents.
      \param UnmatchedhPtHi: p<SUB>T</SUB><SUP>sub</SUP> vs. p<SUB>T</SUB><SUP>lead</SUP> spectrum for all jets with high p<SUB>T</SUB> constituents.
      \param hPtHi: p<SUB>T</SUB><SUP>sub</SUP> vs. p<SUB>T</SUB><SUP>lead</SUP> spectrum for matched jets with high p<SUB>T</SUB> constituents.
      \param hPtLo: p<SUB>T</SUB><SUP>sub</SUP> vs. p<SUB>T</SUB><SUP>lead</SUP> spectrum for matched jets with low p<SUB>T</SUB> constituents.
      \param UnmatchedhdPtHi: &Delta;p<SUB>T</SUB><SUP>sub</SUP> for all jets with high p<SUB>T</SUB> constituents.
      \param hdPtHi: &Delta;p<SUB>T</SUB><SUP>sub</SUP> for matched jets with high p<SUB>T</SUB> constituents.
      \param hdPtLo: &Delta;p<SUB>T</SUB><SUP>sub</SUP> for matched jets with low p<SUB>T</SUB> constituents.
      \param hdphiHi: Dijet angle for matched jets with high p<SUB>T</SUB> constituents.
      \param hdphiLo: Dijet angle for matched jets with low p<SUB>T</SUB> constituents.
      \param OtherAJ_lo: Dijet imbalance &Delta;p<SUB>T</SUB> / &Sigma;p<SUB>T</SUB> for matched jets with low p<SUB>T</SUB> constituents and different R
      \param OtherLeadPtLoss_lo: &Delta;p<SUB>T</SUB> between the two radii in the leading jet
      \param OtherSubLeadPtLoss_lo: &Delta;p<SUB>T</SUB> between the two radii in the sub-leading jet
      \param OtherR: Different radius to match to

      \param hdPtLead: Experimental
      \param hdPtSubLead: Experimental
      \param SpecialhdPtLead: Experimental
      \param SpecialhdPtSubLead: Experimental

   * Return value:
   *   - 0: No hard constituent dijet found, or not matched to ToMatch
   *   - 1: No soft constituent dijet found
   *   - 2: Soft constituent dijet found but not matched
   */

  //int AnalyzeAndFill ( std::vector<fastjet::PseudoJet>& particles, std::vector<fastjet::PseudoJet>& ToMatch,
  int AnalyzeAndFill (	const std::vector<fastjet::PseudoJet>& particles, 			
			//TStarJetVectorContainer<TStarJetVector>& container,
			int ineventid, int inrunid, double inrefmult, double invz,
			//TStarJetPicoReader &reader,
	 		//Int_t mEffUn,	
			const std::vector<std::pair<float,float> > &ToMatch
		      	//Double_t EventClassifier = 0

		       //TH1D* LeadJetPt=0, TH1D* SubJetPt=0, 
		       //TProfile* LeadJetNtrkvsLeadJetPt=0, TProfile* SubJetNtrkvsLeadJetPt=0, TProfile* TranMaxNtrkvsLeadJetPt=0, TProfile* TranMinNtrkvsLeadJetPt=0, TProfile* TranNtrkvsLeadJetPt=0, 
		       //TProfile* LeadJetPtvsLeadJetPt=0, TProfile* SubJetPtvsLeadJetPt=0, TProfile* TranMaxPtvsLeadJetPt=0, TProfile* TranMinPtvsLeadJetPt=0, TProfile* TranPtvsLeadJetPt =0,
		       //TProfile* TranPionNtrkvsLeadJetPt=0, TProfile* TranPionPtvsLeadJetPt=0, TProfile* TranProtonNtrkvsLeadJetPt=0, TProfile* TranProtonPtvsLeadJetPt=0,
		       //TH2D* Spectrum_LeadJetPtvsLeadJetPt=0, TH2D* Spectrum_SubJetPtvsLeadJetPt=0, TH2D* Spectrum_TranMaxPtvsLeadJetPt=0, TH2D* Spectrum_TranMinPtvsLeadJetPt=0, TH2D* Spectrum_TranPtvsLeadJetPt =0
	);

	int Finish();


  /** This little helper is true if there's at least one 10 GeV jet
   **/
  bool Has10Gev;
  /** This little helper is true if there's dijet
   **/
  bool HasDijet;

  // Getters and Setters
  // -------------------
  // Whether need to match jet found by fastjet with the location fired the trigger
  void SetToMatchJetTrigger(bool val) {mNeedToMatchTrig = val; }
  bool GetToMatchJetTrigger() {return mNeedToMatchTrig; }

  // Wether apply Neutral/Total Pt of Jet fraction cut
  void SetNetraulJetFracCut(bool val) {mNeutralJetFracCut = val; };

  // Use Dijet angle (1) or Monojet angle (0)
  void SetDiJetAngle(int val) {mUseDijetAngle = val; }

  // Underlying event particle: charged(1) or neutral(0) or all(2)
  void SetUnderlyingParticleCharge(int val);
  int GetUnderlyingParticleCharge() { return mUnderlyingParticleCharge; }

  // Jet: charged(1) or neutral(0) or all(2)
  void SetJetCharge(int val);
  int GetJetCharge() { return mJetCharge; }

  // Underlying phi size
  void SetTransversePhiSize(float val) { if(val<180&&val>0) {mTranPhiSize = val;} else {mTranPhiSize = 60;} }
  int GetTransversePhiSize() { return mTranPhiSize; }

  /// Get jet radius
  inline double GetR ( )                   { return R; }
  /// Set jet radius
  inline void   SetR ( const double newv ) { R=newv;   }


  /// Get jet rapidity acceptance
  inline double GetMax_rap ( )                   { return max_rap; }
  /// Set jet rapidity acceptance
  inline void   SetMax_rap ( const double newv ) { max_rap=newv;   }

  /// Get ghosted area rapidity cut, should be >= max_rap + 2*R
  inline double GetGhost_maxrap ( )                   { return ghost_maxrap; }
  /// Set ghosted area rapidity cut, should be >= max_rap + 2*R
  inline void   SetGhost_maxrap ( const double newv ) { ghost_maxrap=newv;   }

   /// Get dijet opening angle
  inline double GetDPhiCut ( )                   { return dPhiCut; }
  /// Set dijet opening angle
  inline void   SetDPhiCut ( const double newv ) { dPhiCut=newv;   }
  

  // Objects will be handed by _reference_! Obviates need for setter
  /// Handle to jet definition
  inline fastjet::JetDefinition& GetJet_def () { return jet_def; }
  /// Handle to selector for constituents
  inline fastjet::Selector& GetConsSelector () { return sconst; }
  
  /// Handle to selector for jet candidates
  inline fastjet::Selector& GetJetSelector () { return sjet; }

  /// Handle to ghosted area specification
  inline fastjet::GhostedAreaSpec& GetArea_spec () { return area_spec; }
  /// Handle to jet area definition
  inline fastjet::AreaDefinition& GetArea_def () { return area_def; }


  /// Handle to JetAnalyzer
  inline JetAnalyzer* GetJA() {return pJA; }

  /// Handle to unaltered clustering result with high pT constituents
  inline std::vector<fastjet::PseudoJet> GetJAResult() {return JAResult; }

  /// Handle to constituents
  inline std::vector<fastjet::PseudoJet> GetConstituents() {return Jconstituents; }

  /// Handle to transverse Ntrk, pt
  inline double GetLeadJetPt() { return mLeadAreaPt; } 
  inline double GetSubJetPt() { return mSubAreaPt; } 
  inline double GetTranMaxPt() { return mTranMaxPt; } 
  inline double GetTranMinPt() { return mTranMinPt; } 
  inline double GetTranPt() { return mTranPt; } 

  inline int GetLeadJetNtrk() { return mLeadAreaNtrk; } 
  inline int GetSubJetNtrk() { return mSubAreaNtrk; } 
  inline int GetTranMaxNtrk() { return mTranMaxNtrk; } 
  inline int GetTranMinNtrk() { return mTranMinNtrk; } 
  inline int GetTranNtrk() { return mTranNtrk; } 
  
  /// Handle to Dijet result 
  inline std::vector<fastjet::PseudoJet> GetDiJets() {return DiJets; };

};  

/** Helper to perform the TStarJetPicoReader initialization
 */
/*
TStarJetPicoReader GetReader ( TString ChainPattern="~putschke/Data/Pico_ppHT/*.root", 
			       TString TriggerString="ppHT",
			       TString ChainName="JetTree",
			       const double RefMultCut=0
			       );
*/
/** Slightly different, preferred version of GetReader
 */
/*
TStarJetPicoReader SetupReader ( TChain* chain, TString TriggerString, const double RefMultCut=0 );
*/


#endif // __UNDERLYINGANA_HH
